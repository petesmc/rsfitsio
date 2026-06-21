/*  This file, modkey.rs, contains routines that modify, insert, or update  */
/*  keywords in a FITS header.                                             */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use std::ffi::CStr;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::c_types::*;
use crate::fitscore::{
    ffc2s, ffcmrk_safe, ffgmsg_safe, ffiblk, ffmahd_safe, ffmkey, ffmkky_safe, ffpmrk_safe,
    ffpmsg_slice, ffpmsg_str, ffpsvc_safe, fftkey_safe, fits_strncasecmp,
};
use crate::getkey::{
    ffgcnt, ffgcrd_safe, ffghps_safe, ffgkcsl_safe, ffgkey_safe, ffgrec_safe, ffgskyc_safe,
    ffgstr_safe, ffmaky_safe,
};
use crate::nullable_slice_cstr;
use crate::putkey::{
    ffd2e, ffd2f, ffi2c, ffl2c, ffpkfc_safe, ffpkfm_safe, ffpkls_safe, ffpkyc_safe, ffpkyd_safe,
    ffpkye_safe, ffpkyf_safe, ffpkyg_safe, ffpkyj_safe, ffpkyl_safe, ffpkym_safe, ffpkys_safe,
    ffpkyu_safe, ffpkyuj_safe, ffprec_safe, ffr2e, ffr2f, ffs2c, ffu2c, fits_make_longstr_key_util,
};
use crate::wrappers::*;
use crate::{KeywordDatatype, fitsio2::*};
use crate::{bb, cs};
use crate::{buffers::*, raw_to_slice};
use crate::{fitsio::*, int_snprintf, slice_to_str};

/*--------------------------------------------------------------------------*/
/// Update the keyword, value and comment in the FITS header.
/// The datatype is specified by the 2nd argument.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffuky(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    datatype: c_int,        /* I - datatype of the value    */
    keyname: *const c_char, /* I - name of keyword to write */
    value: *const c_void,   /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        let datatype_with_value = KeywordDatatype::from_datatype(datatype, value);

        ffuky_safe(fptr, datatype_with_value, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Update the keyword, value and comment in the FITS header.
/// The datatype is specified by the 2nd argument.
///
/// Heavily modified for safety.
pub fn ffuky_safe(
    fptr: &mut fitsfile,       /* I - FITS file pointer        */
    datatype: KeywordDatatype, /* I - datatype of the value    */
    keyname: &[c_char],        /* I - name of keyword to write */
    comm: Option<&[c_char]>,   /* I - keyword comment          */
    status: &mut c_int,        /* IO - error status            */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    match datatype {
        KeywordDatatype::TSTRING(value) => {
            ffukys_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatype::TBYTE(value) => {
            ffukyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TSBYTE(value) => {
            ffukyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TUSHORT(value) => {
            ffukyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TSHORT(value) => {
            ffukyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TINT(value) => {
            ffukyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TUINT(value) => {
            ffukyg_safe(fptr, keyname, f64::from(*(value)), 0, comm, status);
        }
        KeywordDatatype::TLOGICAL(value) => {
            ffukyl_safe(fptr, keyname, *(value), comm, status);
        }
        KeywordDatatype::TULONG(value) => {
            ffukyg_safe(fptr, keyname, *value as f64, 0, comm, status);
        }
        KeywordDatatype::TLONG(value) => {
            ffukyj_safe(fptr, keyname, *(value) as LONGLONG, comm, status);
        }
        KeywordDatatype::TULONGLONG(value) => {
            ffukyuj_safe(fptr, keyname, *value, comm, status);
        }
        KeywordDatatype::TLONGLONG(value) => {
            ffukyj_safe(fptr, keyname, *(value), comm, status);
        }
        KeywordDatatype::TFLOAT(value) => {
            ffukye_safe(fptr, keyname, *(value), -7, comm, status);
        }
        KeywordDatatype::TDOUBLE(value) => {
            ffukyd_safe(fptr, keyname, *(value), -15, comm, status);
        }
        KeywordDatatype::TCOMPLEX(value) => {
            ffukyc_safe(fptr, keyname, value, -7, comm, status);
        }
        KeywordDatatype::TDBLCOMPLEX(value) => {
            ffukym_safe(fptr, keyname, value, -15, comm, status);
        }
        _ => {
            *status = BAD_DATATYPE;
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyu(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyu_safe(fptr, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyu_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyu_safe(fptr, keyname, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyu_safe(fptr, keyname, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukys(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const c_char,   /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(value);

        let comm: Option<&[c_char]> = match comm.is_null() {
            true => None,
            false => Some(cast_slice(CStr::from_ptr(comm).to_bytes_with_nul())),
        };

        ffukys_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukys_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[c_char],        /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkys_safe(fptr, keyname, value, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkys_safe(fptr, keyname, value, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukls(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const c_char,   /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(value);
        nullable_slice_cstr!(comm);

        ffukls_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukls_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[c_char],        /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    /* update a long string keyword */
    let mut junk: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkls_safe(fptr, keyname, value, comm, status) == KEY_NO_EXIST {
        /* since the ffmkls call failed, it wrote a bogus error message */

        ffgmsg_safe(&mut junk); /* clear the error message */

        *status = tstatus;
        ffpkls_safe(fptr, keyname, value, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyl(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: c_int,           /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyl_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyl_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: c_int,            /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyl_safe(fptr, keyname, value, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyl_safe(fptr, keyname, value, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyj(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: LONGLONG,        /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyj_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: LONGLONG,         /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyj_safe(fptr, keyname, value, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyj_safe(fptr, keyname, value, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyuj(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: ULONGLONG,       /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyuj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyuj_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: ULONGLONG,        /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyuj_safe(fptr, keyname, value, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyuj_safe(fptr, keyname, value, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyf(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f32,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyf_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyf_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f32,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyf_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyf_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukye(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f32,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukye_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukye_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f32,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkye_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkye_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyg(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f64,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyg_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyg_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f64,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyg_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyg_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyd(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f64,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyd_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyd_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f64,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyd_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyd_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukfc(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f32; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);
        ffukfc_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukfc_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f32; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkfc_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkfc_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukyc(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f32; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukyc_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukyc_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f32; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkyc_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkyc_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukfm(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f64; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukfm_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukfm_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f64; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkfm_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkfm_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffukym(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f64; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffukym_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffukym_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f64; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let tstatus = *status;

    if ffmkym_safe(fptr, keyname, value, decim, comm, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffpkym_safe(fptr, keyname, value, decim, comm, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffucrd(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    card: *const c_char,    /* I - card string value  */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        raw_to_slice!(keyname);
        raw_to_slice!(card);

        ffucrd_safe(fptr, keyname, card, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffucrd_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    keyname: &[c_char],  /* I - keyword name       */
    card: &[c_char],     /* I - card string value  */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut tstatus = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    tstatus = *status;

    if ffmcrd_safe(fptr, keyname, card, status) == KEY_NO_EXIST {
        *status = tstatus;
        ffprec_safe(fptr, card, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmrec(
    fptr: *mut fitsfile, /* I - FITS file pointer               */
    nkey: c_int,         /* I - number of the keyword to modify */
    card: *const c_char, /* I - card string value               */
    status: *mut c_int,  /* IO - error status                   */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(card);

        ffmrec_safe(fptr, nkey, card, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmrec_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer               */
    nkey: c_int,         /* I - number of the keyword to modify */
    card: &[c_char],     /* I - card string value               */
    status: &mut c_int,  /* IO - error status                   */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffmaky_safe(fptr, nkey + 1, status);
    ffmkey(fptr, card, status);
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmcrd(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    card: *const c_char,    /* I - card string value  */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(card);

        ffmcrd_safe(fptr, keyname, card, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmcrd_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    keyname: &[c_char],  /* I - keyword name       */
    card: &[c_char],     /* I - card string value  */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut tcard: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut valstring: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT]; // Original code as FLEN_CARD, likely an error
    let mut value: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut nextcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgcrd_safe(fptr, keyname, &mut tcard, status) > 0 {
        return *status;
    }

    ffmkey(fptr, card, status);

    /* calc position of keyword in header */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    let keypos = (((fptr.Fptr.nextkey) - (headstart[fptr.Fptr.curhdu as usize])) / 80) + 1;

    ffpsvc_safe(&tcard, &mut valstring, Some(&mut comm), status);

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check for string value which may be continued over multiple keywords */
    ffpmrk_safe(); /* put mark on message stack; erase any messages after this */
    ffc2s(&valstring, &mut value, status); /* remove quotes and trailing spaces */

    if *status == VALUE_UNDEFINED {
        ffcmrk_safe(); /* clear any spurious error messages, back to the mark */
        *status = 0;
    } else {
        let mut len = strlen_safe(&value);

        while len > 0 && value[len - 1] == bb(b'&') {
            /* ampersand used as continuation char */
            ffgcnt(fptr, &mut value, Some(&mut nextcomm), status);
            if value[0] != 0 {
                ffdrec_safe(fptr, keypos as c_int, status); /* delete the keyword */
                len = strlen_safe(&value);
            } else {
                /* a null valstring indicates no continuation */
                len = 0;
            }
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmnam(
    fptr: *mut fitsfile,    /* I - FITS file pointer     */
    oldname: *const c_char, /* I - existing keyword name */
    newname: *const c_char, /* I - new name for keyword  */
    status: *mut c_int,     /* IO - error status         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(oldname);
        raw_to_slice!(newname);

        ffmnam_safe(fptr, oldname, newname, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmnam_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer     */
    oldname: &[c_char],  /* I - existing keyword name */
    newname: &[c_char],  /* I - new name for keyword  */
    status: &mut c_int,  /* IO - error status         */
) -> c_int {
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, oldname, &mut value, Some(&mut comm), status) > 0 {
        return *status;
    }

    ffmkky_safe(newname, &value, Some(&comm), &mut card, status); /* construct the card */
    ffmkey(fptr, &card, status); /* rewrite with new name */

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmcom(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmcom_safe(fptr, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmcom_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut value, Some(&mut oldcomm), status) > 0 {
        return *status;
    }

    ffmkky_safe(keyname, &value, comm, &mut card, status); /* construct the card */
    ffmkey(fptr, &card, status); /* rewrite with new name */

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the units string into the comment field of the existing keyword.
///
/// This routine uses a  FITS convention  in which the units are enclosed in
/// square brackets following the '/' comment field delimiter, e.g.:
///
/// KEYWORD =                   12 / [kpc] comment string goes here
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpunt(
    fptr: *mut fitsfile,    /* I - FITS file pointer   */
    keyname: *const c_char, /* I - keyword name        */
    unit: *const c_char,    /* I - keyword unit string */
    status: *mut c_int,     /* IO - error status       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(unit);

        ffpunt_safe(fptr, keyname, unit, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the units string into the comment field of the existing keyword.
///
/// This routine uses a  FITS convention  in which the units are enclosed in
/// square brackets following the '/' comment field delimiter, e.g.:
///
/// KEYWORD =                   12 / [kpc] comment string goes here
pub fn ffpunt_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer   */
    keyname: &[c_char],  /* I - keyword name        */
    unit: &[c_char],     /* I - keyword unit string */
    status: &mut c_int,  /* IO - error status       */
) -> c_int {
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut newcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    let mut len = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut value, Some(&mut oldcomm), status) > 0 {
        return *status;
    }

    /* copy the units string to the new comment string if not null */
    if unit[0] != 0 {
        strcpy_safe(&mut newcomm, cs!(c"["));
        strncat_safe(&mut newcomm, unit, 45); /* max allowed length is about 45 chars */
        strcat_safe(&mut newcomm, cs!(c"] "));
        len = strlen_safe(&newcomm);
        len = FLEN_COMMENT - len - 1; /* amount of space left in the field */
    } else {
        newcomm[0] = 0;
        len = FLEN_COMMENT - 1;
    }

    if oldcomm[0] == bb(b'[')
    /* check for existing units field */
    {
        let loc = strchr_safe(&oldcomm, bb(b']')); /* look for the closing bracket */
        if let Some(mut loc_inner) = loc {
            loc_inner += 1;
            while oldcomm[loc_inner] == bb(b' ') {
                /* skip any blank spaces */
                loc_inner += 1;
            }

            strncat_safe(&mut newcomm, &oldcomm[loc_inner..], len); /* concat remainder of comment */
        } else {
            strncat_safe(&mut newcomm, &oldcomm, len); /* append old comment onto new */
        }
    } else {
        strncat_safe(&mut newcomm, &oldcomm, len);
    }

    ffmkky_safe(keyname, &value, Some(&newcomm), &mut card, status); /* construct the card */
    ffmkey(fptr, &card, status); /* rewrite with new units string */

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyu(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyu_safe(fptr, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyu_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }

    strcpy_safe(&mut valstring, cs!(c" ")); /* create a dummy value string */

    //if (!comm || comm[0] == '&') { /* preserve the current comment string */
    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }

    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkys(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const c_char,   /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(value);

        let comm: Option<&[c_char]> = match comm.is_null() {
            true => None,
            false => Some(cast_slice(CStr::from_ptr(comm).to_bytes_with_nul())),
        };

        ffmkys_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkys_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[c_char],        /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    /* NOTE: This routine does not support long continued strings */
    /*  It will correctly overwrite an existing long continued string, */
    /*  but it will not write a new long string.  */

    let mut oldval: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut nextcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut oldval, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }

    ffs2c(value, &mut valstring, status); /* convert value to a string */

    //if (!comm || comm[0] == '&')  /* preserve the current comment string */
    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }

    ffmkey(fptr, &card, status); /* overwrite the previous keyword */

    let headstart = fptr.Fptr.get_headstart_as_slice();
    let keypos =
        ((((fptr.Fptr.nextkey) - (headstart[fptr.Fptr.curhdu as usize])) / 80) + 1) as c_int;

    if *status > 0 {
        return *status;
    }

    /* check if old string value was continued over multiple keywords */
    ffpmrk_safe(); /* put mark on message stack; erase any messages after this */
    ffc2s(&oldval, &mut valstring, status); /* remove quotes and trailing spaces */

    if *status == VALUE_UNDEFINED {
        ffcmrk_safe(); /* clear any spurious error messages, back to the mark */
        *status = 0;
    } else {
        let mut len = strlen_safe(&valstring);

        while len > 0 && valstring[len - 1] == bb(b'&') {
            /* ampersand is continuation char */

            nextcomm[0] = 0;
            ffgcnt(fptr, &mut valstring, Some(&mut nextcomm), status);
            if valstring[0] != 0 || strlen_safe(&nextcomm) != 0 {
                ffdrec_safe(fptr, keypos, status); /* delete the continuation */
                len = strlen_safe(&valstring);
            } else {
                /* a null valstring indicates no continuation */
                len = 0;
            }
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Modify the value and optionally the comment of a long string keyword.
///
/// This routine supports the
/// HEASARC long string convention and can modify arbitrarily long string
/// keyword values.  The value is continued over multiple keywords that
/// have the name COMTINUE without an equal sign in column 9 of the card.
/// This routine also supports simple string keywords which are less than
/// 69 characters in length.
///
/// This routine is not very efficient, so it should be used sparingly.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkls(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    value: *const c_char,   /* I - keyword value            */
    incomm: *const c_char,  /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(value);

        nullable_slice_cstr!(incomm);

        ffmkls_safe(fptr, keyname, value, incomm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Modify the value and optionally the comment of a long string keyword.
///
/// This routine supports the
/// HEASARC long string convention and can modify arbitrarily long string
/// keyword values.  The value is continued over multiple keywords that
/// have the name COMTINUE without an equal sign in column 9 of the card.
/// This routine also supports simple string keywords which are less than
/// 69 characters in length.
///
/// This routine is not very efficient, so it should be used sparingly.
pub fn ffmkls_safe(
    fptr: &mut fitsfile,       /* I - FITS file pointer        */
    keyname: &[c_char],        /* I - name of keyword to write */
    value: &[c_char],          /* I - keyword value            */
    incomm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,        /* IO - error status            */
) -> c_int {
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut comm: Option<Vec<c_char>> = None;
    let mut nkeys: c_int = 0;
    let mut keypos: c_int = 0;
    let mut vlen: c_int = 0;
    let mut commlen: c_int = 0;
    let mut tmpvlen: c_int = 0;
    let mut tmpcommlen: c_int = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let incomm_empty_or_continue =
        incomm.is_none() || incomm.as_ref().is_some_and(|c| c[0] == bb(b'&'));

    /* preserve the old comment string */
    if incomm_empty_or_continue {
        ffghps_safe(fptr, Some(&mut nkeys), Some(&mut keypos), status); /* save current position */

        if ffgkcsl_safe(fptr, keyname, &mut vlen, &mut commlen, status) != 0 {
            return *status; /* keyword doesn't exist or is bad format */
        }

        let mut tmplongval: Vec<c_char> = vec![0; vlen as usize + 1];
        let mut c: Vec<c_char> = vec![0; commlen as usize + 1];

        ffgskyc_safe(
            fptr,
            keyname,
            1,
            vlen,
            commlen,
            Some(&mut tmplongval),
            &mut tmpvlen,
            Some(&mut c),
            &mut tmpcommlen,
            status,
        );

        comm = Some(c);

        /* move back to previous position to ensure that we delete */
        /* the right keyword in case there are more than one keyword */
        /* with this same name. */
        ffgrec_safe(fptr, keypos - 1, Some(&mut card), status);
    } else {
        /* copy the input comment string */
        let incomm_shadow = incomm.unwrap();
        commlen = strlen_safe(cast_slice(incomm_shadow)) as c_int;
        if commlen > 0 {
            let mut c = vec![0; commlen as usize + 1];
            strcpy_safe(&mut c, incomm_shadow);
            comm = Some(c);
        }
    }

    /* delete the old keyword */
    if ffdkey_safe(fptr, keyname, status) > 0 {
        return *status; /* keyword doesn't exist */
    }

    ffghps_safe(fptr, Some(&mut nkeys), Some(&mut keypos), status); /* save current position */

    fits_make_longstr_key_util(fptr, keyname, value, comm.as_deref(), keypos, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyl(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: c_int,           /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyl_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyl_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: c_int,            /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        /* get old comment */
        return *status;
    }

    ffl2c(value, &mut valstring, status); /* convert value to a string */

    // if (comm || comm[0] == bb(b'&')) {
    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }
    ffmkey(fptr, &card, status); /* rewrite with new name */

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyj(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: LONGLONG,        /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyj_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: LONGLONG,         /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }
    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }
    ffi2c(value, &mut valstring, status); /* convert value to a string */

    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }
    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyuj(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: ULONGLONG,       /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyuj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyuj_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: ULONGLONG,        /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }

    ffu2c(value, &mut valstring, status); /* convert value to a string */

    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }

    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyf(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f32,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyf_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyf_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f32,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        /* get old comment */
        return *status;
    }

    ffr2f(value, decim, &mut valstring, status); /* convert value to a string */

    // if (comm || comm[0] == bb(b'&')) {
    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }
    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkye(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f32,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkye_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkye_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f32,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        /* get old comment */
        return *status;
    }

    ffr2e(value, decim, &mut valstring, status); /* convert value to a string */

    // if (comm || comm[0] == bb(b'&')) {
    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }
    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyg(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f64,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyg_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyg_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f64,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        /* get old comment */
        return *status;
    }

    ffd2f(value, decim, &mut valstring, status); /* convert value to a string */

    // if (comm || comm[0] == bb(b'&')) {
    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }
    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyd(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f64,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyd_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyd_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f64,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        /* get old comment */
        return *status;
    }

    ffd2e(value, decim, &mut valstring, status); /* convert value to a string */

    // if (comm || comm[0] == bb(b'&')) {
    if comm.is_none() {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else if let Some(c) = comm
        && c[0] == bb(b'&')
    {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }
    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkfc(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f32; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkfc_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkfc_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f32; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffr2f(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkfc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffr2f(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkfc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    if comm.is_none() || (comm.is_some() && comm.unwrap()[0] == bb(b'&')) {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }

    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkyc(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f32; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkyc_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkyc_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f32; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffr2e(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkyc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffr2e(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkyc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    if comm.is_none() || (comm.is_some() && comm.unwrap()[0] == bb(b'&')) {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }

    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkfm(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f64; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkfm_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkfm_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f64; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffd2f(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkfm)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffd2f(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkfm)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    if comm.is_none() || (comm.is_some() && comm.unwrap()[0] == bb(b'&')) {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }

    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmkym(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f64; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffmkym_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffmkym_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f64; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut oldcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut oldcomm), status) > 0 {
        return *status; /* get old comment */
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffd2e(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkym)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffd2e(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffmkym)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    if comm.is_none() || (comm.is_some() && comm.unwrap()[0] == bb(b'&')) {
        /* preserve the current comment string */
        ffmkky_safe(keyname, &valstring, Some(&oldcomm), &mut card, status);
    } else {
        ffmkky_safe(keyname, &valstring, comm, &mut card, status);
    }

    ffmkey(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Insert a null-valued keyword and comment into the FITS header.  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikyu(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikyu_safer(fptr, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a null-valued keyword and comment into the FITS header.
pub fn ffikyu_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c" ")); /* create a dummy value string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikys(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const c_char,   /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(value);
        nullable_slice_cstr!(comm);

        ffikys_safer(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a string keyword into the FITS header at the current position.
pub fn ffikys_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[c_char],        /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffs2c(value, &mut valstring, status); /* put quotes around the string */

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Insert a long string keyword.  This routine supports the
/// HEASARC long string convention and can insert arbitrarily long string
/// keyword values.  The value is continued over multiple keywords that
/// have the name COMTINUE without an equal sign in column 9 of the card.
/// This routine also supports simple string keywords which are less than
/// 69 characters in length.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikls(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    value: *const c_char,   /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(value);
        nullable_slice_cstr!(comm);

        ffikls_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a long string keyword.  This routine supports the
/// HEASARC long string convention and can insert arbitrarily long string
/// keyword values.  The value is continued over multiple keywords that
/// have the name COMTINUE without an equal sign in column 9 of the card.
/// This routine also supports simple string keywords which are less than
/// 69 characters in length.
pub fn ffikls_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    value: &[c_char],        /* I - keyword value            */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut tmpkeyname: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut tstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tstatus = -1;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /*  construct the new keyword, and insert into header */
    let mut remain = strlen_safe(value) as c_int; /* number of characters to write out */
    let mut next = 0; /* pointer to next character to write */

    /* count the number of single quote characters in the string */
    let mut nquote = 0;
    for &ch in value.iter() {
        if ch == bb(b'\'') {
            nquote += 1;
        }
    }

    strncpy_safe(&mut tmpkeyname, keyname, 80);
    tmpkeyname[80] = 0;

    // Find first non-space character in keyname
    let mut cptr_offset = 0;
    while cptr_offset < tmpkeyname.len() && tmpkeyname[cptr_offset] == bb(b' ') {
        cptr_offset += 1;
    }

    /* determine the number of characters that will fit on the line */
    /* Note: each quote character is expanded to 2 quotes */

    let namelen = strlen_safe(&tmpkeyname[cptr_offset..]) as c_int;
    let mut nchar = if namelen <= 8 && (fftkey_safe(&tmpkeyname[cptr_offset..], &mut tstatus) <= 0)
    {
        /* This a normal 8-character FITS keyword */
        68 - nquote /*  max of 68 chars fit in a FITS string value */
    } else {
        80 - nquote - namelen - 5
    };

    let mut contin = 0;
    while remain > 0 {
        if nchar > FLEN_VALUE as c_int - 1 {
            ffpmsg_str("longstr keyword value is too long (ffikls)");
            *status = BAD_KEYCHAR;
            return *status;
        }

        // Copy substring to temporary buffer
        let copy_len = if nchar > remain { remain } else { nchar } as usize;
        for i in 0..copy_len {
            tstring[i] = value[(next as usize) + i];
        }
        tstring[copy_len] = 0;

        ffs2c(&tstring, &mut valstring, status); /* put quotes around the string */

        if remain > nchar
        /* if string is continued, put & as last char */
        {
            let vlen = strlen_safe(&valstring);
            nchar -= 1; /* outputting one less character now */

            if valstring[vlen - 2] != bb(b'\'') {
                valstring[vlen - 2] = bb(b'&'); /*  over write last char with &  */
            } else {
                /* last char was a pair of single quotes, so over write both */
                valstring[vlen - 3] = bb(b'&');
                valstring[vlen - 1] = 0;
            }
        }

        if contin != 0
        /* This is a CONTINUEd keyword */
        {
            ffmkky_safe(cs!(c"CONTINUE"), &valstring, comm, &mut card, status); /* make keyword */
            // Overwrite the '=' with spaces
            card[8] = bb(b' ');
            card[9] = bb(b' ');
        } else {
            ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* make keyword */
        }

        ffikey_safe(fptr, &card, status); /* insert the keyword */

        contin = 1;
        remain -= nchar;
        next += nchar;
        nchar = 68 - nquote;
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikyl(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: c_int,           /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikyl_safer(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a logical keyword into the FITS header at the current position.
pub fn ffikyl_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: c_int,            /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffl2c(value, &mut valstring, status); /* convert logical to 'T' or 'F' */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikyj(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: LONGLONG,        /* I - keyword value      */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikyj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffikyj_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: LONGLONG,         /* I - keyword value      */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffi2c(value, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
pub fn ffikyf_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f32,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    ffr2f(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikyf(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f32,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikyf_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikye(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f32,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikye_safer(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a float keyword into the FITS header at the current position.
pub fn ffikye_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f32,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffr2e(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikyg(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f64,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikyg_safer(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a double keyword into the FITS header at the current position.
pub fn ffikyg_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f64,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffd2f(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikyd(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: f64,             /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikyd_safer(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a double keyword into the FITS header at the current position.
pub fn ffikyd_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: f64,              /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffd2e(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikfc(
    fptr: *const fitsfile,  /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f32; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = (fptr as *mut fitsfile).as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikfc_safer(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a complex float keyword into the FITS header at the current position.
pub fn ffikfc_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f32; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffr2f(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikfc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffr2f(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikfc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikyc(
    fptr: *const fitsfile,  /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f32; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = (fptr as *mut fitsfile).as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikyc_safer(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a complex float keyword into the FITS header at the current position.
pub fn ffikyc_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f32; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffr2e(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikyc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffr2e(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikyc)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikfm(
    fptr: *const fitsfile,  /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f64; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = (fptr as *mut fitsfile).as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikfm_safer(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a complex double keyword into the FITS header at the current position.
pub fn ffikfm_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f64; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffd2f(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikfm)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffd2f(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikfm)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikym(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    value: *const [f64; 2], /* I - keyword value      */
    decim: c_int,           /* I - no of decimals     */
    comm: *const c_char,    /* I - keyword comment    */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffikym_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a complex double keyword into the FITS header at the current position.
pub fn ffikym_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer  */
    keyname: &[c_char],      /* I - keyword name       */
    value: &[f64; 2],        /* I - keyword value      */
    decim: c_int,            /* I - no of decimals     */
    comm: Option<&[c_char]>, /* I - keyword comment    */
    status: &mut c_int,      /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffd2e(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&tmpstring) + 3 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikym)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffd2e(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("complex key value too long (ffikym)");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffikey_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffirec(
    fptr: *mut fitsfile, /* I - FITS file pointer              */
    nkey: c_int,         /* I - position to insert new keyword */
    card: *const c_char, /* I - card string value              */
    status: *mut c_int,  /* IO - error status                  */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(card);

        ffirec_safe(fptr, nkey, card, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffirec_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer              */
    nkey: c_int,         /* I - position to insert new keyword */
    card: &[c_char],     /* I - card string value              */
    status: &mut c_int,  /* IO - error status                  */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffmaky_safe(fptr, nkey, status); /* move to insert position */
    ffikey_safe(fptr, card, status); /* insert the keyword card */

    *status
}

/*--------------------------------------------------------------------------*/
/// Insert a keyword at the position of (fptr->Fptr)->nextkey
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffikey(
    fptr: *mut fitsfile, /* I - FITS file pointer  */
    card: *const c_char, /* I - card string value  */
    status: *mut c_int,  /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(card);

        ffikey_safe(fptr, card, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a keyword at the position of (fptr->Fptr)->nextkey
pub fn ffikey_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    card: &[c_char],     /* I - card string value  */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut nblocks = 0;
    let buff1: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut buff2: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if (fptr.Fptr.datastart - fptr.Fptr.headend) == 80 {
        /* only room for END card */
        nblocks = 1;
        if ffiblk(fptr, nblocks, 0, status) > 0 {
            /* add new 2880-byte block*/
            return *status;
        }
    }

    /* no. keywords to shift */
    let nshift = ((fptr.Fptr.headend - fptr.Fptr.nextkey) / 80) as c_int;

    strncpy_safe(&mut buff2, card, 80); /* copy card to output buffer */
    buff2[80] = 0;

    let len = strlen_safe(&buff2);

    /* silently replace any illegal characters with a space */
    for ii in 0..len {
        if buff2[ii] < bb(b' ') || buff2[ii] > 126 {
            buff2[ii] = bb(b' ');
        }
    }

    for ii in len..80 {
        /* fill buffer with spaces if necessary */
        buff2[ii] = bb(b' ');
    }

    let mut keylength = strcspn_safe(&buff2, cs!(c"="));
    if keylength == 80 {
        keylength = 8;
    }

    /* test for the common commentary keywords which by definition have 8-char names */
    if fits_strncasecmp(cs!(c"COMMENT "), &buff2, 8) == 0
        || fits_strncasecmp(cs!(c"HISTORY "), &buff2, 8) == 0
        || fits_strncasecmp(cs!(c"        "), &buff2, 8) == 0
        || fits_strncasecmp(cs!(c"CONTINUE"), &buff2, 8) == 0
    {
        keylength = 8;
    }

    for ii in 0..(keylength as usize) {
        /* make sure keyword name is uppercase */
        buff2[ii] = toupper(buff2[ii]);
    }

    fftkey_safe(&buff2, status); /* test keyword name contains legal chars */

    /*  no need to do this any more, since any illegal characters have been removed
    fftrec(buff2, status);  */
    /* test rest of keyword for legal chars   */

    let mut inbuff = buff1;
    let mut outbuff = buff2;

    let mut bytepos = fptr.Fptr.nextkey; /* pointer to next keyword in header */
    ffmbyt_safe(fptr, bytepos, REPORT_EOF, status);

    for _ii in 0..(nshift as usize) {
        /* shift each keyword down one position */

        ffgbyt(fptr, 80, cast_slice_mut(&mut inbuff), status); /* read the current keyword */

        ffmbyt_safe(fptr, bytepos, REPORT_EOF, status); /* move back */
        ffpbyt(fptr, 80, cast_slice_mut(&mut outbuff), status); /* overwrite with other buffer */

        std::mem::swap(&mut inbuff, &mut outbuff);

        bytepos += 80;
    }

    ffpbyt(fptr, 80, cast_slice(&outbuff), status); /* write the final keyword */

    fptr.Fptr.headend += 80; /* increment the position of the END keyword */
    fptr.Fptr.nextkey += 80; /* increment the pointer to next keyword */

    *status
}

/*--------------------------------------------------------------------------*/
/// delete a specified header keyword
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdkey(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    keyname: *const c_char, /* I - keyword name       */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        ffdkey_safe(fptr, keyname, status)
    }
}

/*--------------------------------------------------------------------------*/
/// delete a specified header keyword
pub fn ffdkey_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    keyname: &[c_char],  /* I - keyword name       */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut nextcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut comm), status) > 0 {
        /* read keyword */

        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Could not find the {} keyword to delete (ffdkey)",
            slice_to_str!(keyname),
        );
        ffpmsg_slice(&message);
        return *status;
    }

    /* calc position of keyword in header */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    let keypos = (((fptr.Fptr.nextkey) - (headstart[fptr.Fptr.curhdu as usize])) / 80) as c_int;

    ffdrec_safe(fptr, keypos, status); /* delete the keyword */

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check for string value which may be continued over multiple keywords */
    ffpmrk_safe(); /* put mark on message stack; erase any messages after this */
    ffc2s(&valstring, &mut value, status); /* remove quotes and trailing spaces */

    if *status == VALUE_UNDEFINED {
        ffcmrk_safe(); /* clear any spurious error messages, back to the mark */
        *status = 0;
    } else {
        let mut len = strlen_safe(&value);

        while len != 0 && value[len - 1] == bb(b'&') {
            /* ampersand used as continuation char */

            nextcomm[0] = 0;
            ffgcnt(fptr, &mut value, Some(&mut nextcomm), status);
            if value[0] != 0 || strlen_safe(&nextcomm) != 0 {
                ffdrec_safe(fptr, keypos, status); /* delete the keyword */
                len = strlen_safe(&value);
            } else {
                /* a null valstring indicates no continuation */
                len = 0;
            }
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// delete a specified header keyword containing the input string
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdstr(
    fptr: *mut fitsfile,   /* I - FITS file pointer  */
    string: *const c_char, /* I - keyword name       */
    status: *mut c_int,    /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(string);

        ffdstr_safe(fptr, string, status)
    }
}

/*--------------------------------------------------------------------------*/
/// delete a specified header keyword containing the input string
pub fn ffdstr_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    string: &[c_char],   /* I - keyword string     */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut nextcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if ffgstr_safe(fptr, string, &mut card, status) > 0 {
        /* read keyword */
        int_snprintf!(
            message,
            FLEN_ERRMSG,
            "Could not find the {} keyword to delete (ffdkey)",
            CStr::from_bytes_until_nul(cast_slice(string))
                .unwrap()
                .to_str()
                .unwrap()
        );
        ffpmsg_slice(&message);
        return *status;
    }

    /* calc position of keyword in header */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    let keypos = (((fptr.Fptr.nextkey) - headstart[fptr.Fptr.curhdu as usize]) / 80) as c_int;

    ffdrec_safe(fptr, keypos, status); /* delete the keyword */

    /* check for string value which may be continued over multiple keywords */
    ffpsvc_safe(&card, &mut valstring, Some(&mut comm), status);

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check for string value which may be continued over multiple keywords */
    ffpmrk_safe(); /* put mark on message stack; erase any messages after this */
    ffc2s(&valstring, &mut value, status); /* remove quotes and trailing spaces */

    if *status == VALUE_UNDEFINED {
        ffcmrk_safe(); /* clear any spurious error messages, back to the mark */
        *status = 0;
    } else {
        let mut len = strlen_safe(&value) as c_int;

        while len > 0 && value[(len - 1) as usize] == bb(b'&')
        /* ampersand used as continuation char */
        {
            ffgcnt(fptr, &mut value, Some(&mut nextcomm), status);
            if value[0] != 0 {
                ffdrec_safe(fptr, keypos, status); /* delete the keyword */
                len = strlen_safe(&value) as c_int;
            } else {
                /* a null valstring indicates no continuation */
                len = 0;
            }
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Delete a header keyword at position keypos. The 1st keyword is at keypos=1.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdrec(
    fptr: *mut fitsfile, /* I - FITS file pointer  */
    keypos: c_int,       /* I - position in header of keyword to delete */
    status: *mut c_int,  /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffdrec_safe(fptr, keypos, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Delete a header keyword at position keypos. The 1st keyword is at keypos=1.
pub fn ffdrec_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    keypos: c_int,       /* I - position in header of keyword to delete */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut nshift = 0;
    let mut bytepos: LONGLONG = 0;

    //char *inbuff, *outbuff, *tmpbuff,
    let buff1: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut buff2: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let headstart = fptr.Fptr.get_headstart_as_slice();

    if keypos < 1
        || LONGLONG::from(keypos) > fptr.Fptr.headend - headstart[fptr.Fptr.curhdu as usize] / 80
    {
        *status = KEY_OUT_BOUNDS;
        return *status;
    }

    fptr.Fptr.nextkey = headstart[fptr.Fptr.curhdu as usize] + (LONGLONG::from(keypos) - 1) * 80;

    nshift = ((fptr.Fptr.headend - fptr.Fptr.nextkey) / 80) as c_int; /* no. keywords to shift */

    if nshift <= 0 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Cannot delete keyword number {}.  It does not exist.",
            keypos,
        );
        ffpmsg_slice(&message);
        *status = KEY_OUT_BOUNDS;
        return *status;
    }

    bytepos = fptr.Fptr.headend - 80; /* last keyword in header */

    /* construct a blank keyword */
    strcpy_safe(&mut buff2, cs!(c"                                        "));
    strcat_safe(&mut buff2, cs!(c"                                        "));
    let mut inbuff = buff1;
    let mut outbuff = buff2;
    for _ii in 0..(nshift as usize) {
        /* shift each keyword up one position */

        ffmbyt_safe(fptr, bytepos, REPORT_EOF, status);
        ffgbyt(fptr, 80, cast_slice_mut(&mut inbuff), status); /* read the current keyword */

        ffmbyt_safe(fptr, bytepos, REPORT_EOF, status);
        ffpbyt(fptr, 80, cast_slice_mut(&mut outbuff), status); /* overwrite with next keyword */

        std::mem::swap(&mut inbuff, &mut outbuff);

        bytepos -= 80;
    }

    fptr.Fptr.headend -= 80; /* decrement the position of the END keyword */
    *status
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CStr;
    use std::ptr;

    use crate::KeywordDatatype;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        BAD_DATATYPE, BAD_F2C, BYTE_IMG, FLEN_CARD, FLEN_VALUE, KEY_NO_EXIST, KEY_OUT_BOUNDS,
        LONGLONG, ULONGLONG, fitsfile,
    };
    use crate::helpers::testhelpers::{from_buf, to_buf, with_temp_file};
    use libc::{c_char, c_int, c_void};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Create a fresh file with an empty BYTE_IMG primary HDU and run `body`.
    fn with_byte_img<F>(body: F)
    where
        F: FnOnce(&mut fitsfile, &mut c_int),
    {
        // `with_temp_file` requires an `Fn` closure, but `body` is `FnOnce`,
        // so stash it in a `RefCell<Option<_>>` and take it out on the single
        // invocation.
        let body = std::cell::RefCell::new(Some(body));
        with_temp_file(|filename| {
            let body = body.borrow_mut().take().unwrap();
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut fptr: Option<Box<fitsfile>> = None;
            fits_create_file(&mut fptr, &name, &mut status);
            assert_eq!(status, 0, "ffinit failed");
            let f = fptr.as_deref_mut().unwrap();
            fits_write_imghdr(f, BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "ffphps failed");
            body(f, &mut status);
            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    /// Read a long string keyword and return an owned String, freeing the
    /// underlying allocation via fits_free_memory.
    fn read_longstr(f: &mut fitsfile, key: &[c_char], status: &mut c_int) -> String {
        let mut result: *mut c_char = ptr::null_mut();
        fits_read_key_longstr(f, key, &mut result, None, status);
        assert_eq!(*status, 0, "ffgkls failed");
        let s = unsafe { CStr::from_ptr(result) }
            .to_str()
            .unwrap()
            .to_string();
        fits_free_memory(result as *mut c_void, status);
        s
    }

    // ---------------- Update keyword tests ----------------

    #[test]
    fn test_ffuky_string() {
        with_byte_img(|f, status| {
            let mut value = [0 as c_char; FLEN_VALUE];
            // Update when keyword does not exist - should insert
            fits_update_key_str(
                f,
                &cc("TESTKEY"),
                &cc("value1"),
                Some(&cc("comment1")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_str(f, &cc("TESTKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "value1");
            // Update when keyword exists - should modify
            fits_update_key_str(
                f,
                &cc("TESTKEY"),
                &cc("value2"),
                Some(&cc("comment2")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_str(f, &cc("TESTKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "value2");
        });
    }

    #[test]
    fn test_ffuky_logical() {
        with_byte_img(|f, status| {
            let mut value: c_int = 0;
            fits_update_key_log(f, &cc("LOGKEY"), 1, Some(&cc("logical true")), status);
            assert_eq!(*status, 0);
            fits_read_key_log(f, &cc("LOGKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 1);
            fits_update_key_log(f, &cc("LOGKEY"), 0, Some(&cc("logical false")), status);
            assert_eq!(*status, 0);
            fits_read_key_log(f, &cc("LOGKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 0);
        });
    }

    #[test]
    fn test_ffuky_longlong() {
        with_byte_img(|f, status| {
            let mut value: LONGLONG = 0;
            fits_update_key_lng(
                f,
                &cc("LLKEY"),
                123456789012,
                Some(&cc("longlong value")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_lnglng(f, &cc("LLKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 123456789012);
            fits_update_key_lng(f, &cc("LLKEY"), -987654321012, Some(&cc("updated")), status);
            assert_eq!(*status, 0);
            fits_read_key_lnglng(f, &cc("LLKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, -987654321012);
        });
    }

    #[test]
    fn test_ffuky_ulonglong() {
        with_byte_img(|f, status| {
            let mut value: ULONGLONG = 0;
            fits_update_key_ulng(
                f,
                &cc("ULLKEY"),
                18446744073709551000,
                Some(&cc("ulonglong")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_ulnglng(f, &cc("ULLKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 18446744073709551000);
        });
    }

    #[test]
    fn test_ffuky_float_fixed() {
        with_byte_img(|f, status| {
            let mut value: f32 = 0.0;
            fits_update_key_fixflt(
                f,
                &cc("FLTKEY"),
                3.14159,
                5,
                Some(&cc("float fixed")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_flt(f, &cc("FLTKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 3.14159).abs() <= 0.0001);
        });
    }

    #[test]
    fn test_ffuky_float_exp() {
        with_byte_img(|f, status| {
            let mut value: f32 = 0.0;
            fits_update_key_flt(
                f,
                &cc("FLTEXP"),
                1.23e10,
                4,
                Some(&cc("float exponential")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_flt(f, &cc("FLTEXP"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 1.23e10).abs() / 1.23e10 <= 0.001);
        });
    }

    #[test]
    fn test_ffuky_double_fixed() {
        with_byte_img(|f, status| {
            let mut value: f64 = 0.0;
            fits_update_key_fixdbl(
                f,
                &cc("DBLKEY"),
                3.141592653589793,
                10,
                Some(&cc("double fixed")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dbl(f, &cc("DBLKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 3.141592653589793).abs() <= 1e-10);
        });
    }

    #[test]
    fn test_ffuky_double_exp() {
        with_byte_img(|f, status| {
            let mut value: f64 = 0.0;
            fits_update_key_dbl(
                f,
                &cc("DBLEXP"),
                1.23456789e20,
                8,
                Some(&cc("double exponential")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dbl(f, &cc("DBLEXP"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 1.23456789e20).abs() / 1.23456789e20 <= 1e-8);
        });
    }

    #[test]
    fn test_ffuky_undefined() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_update_key_null(f, &cc("UNDEF"), Some(&cc("undefined keyword")), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("UNDEF"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("UNDEF"));
        });
    }

    #[test]
    fn test_ffuky_complex_float() {
        with_byte_img(|f, status| {
            let cval: [f32; 2] = [1.5, 2.5];
            let mut result: [f32; 2] = [0.0; 2];
            fits_update_key_cmp(
                f,
                &cc("CMPLXF"),
                &cval,
                4,
                Some(&cc("complex float")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_cmp(f, &cc("CMPLXF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 1.5).abs() <= 0.01 && (result[1] - 2.5).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffuky_complex_double() {
        with_byte_img(|f, status| {
            let cval: [f64; 2] = [1.23456789, 9.87654321];
            let mut result: [f64; 2] = [0.0; 2];
            fits_update_key_dblcmp(
                f,
                &cc("CMPLXD"),
                &cval,
                8,
                Some(&cc("complex double")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dblcmp(f, &cc("CMPLXD"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 1.23456789).abs() <= 1e-6);
            assert!((result[1] - 9.87654321).abs() <= 1e-6);
        });
    }

    #[test]
    fn test_ffuky_complex_float_fixed() {
        with_byte_img(|f, status| {
            let cval: [f32; 2] = [10.5, 20.5];
            let mut result: [f32; 2] = [0.0; 2];
            fits_update_key_fixcmp(
                f,
                &cc("CMPLXFF"),
                &cval,
                2,
                Some(&cc("complex float fixed")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_cmp(f, &cc("CMPLXFF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 10.5).abs() <= 0.1 && (result[1] - 20.5).abs() <= 0.1);
        });
    }

    #[test]
    fn test_ffuky_complex_double_fixed() {
        with_byte_img(|f, status| {
            let cval: [f64; 2] = [100.125, 200.875];
            let mut result: [f64; 2] = [0.0; 2];
            fits_update_key_fixdblcmp(
                f,
                &cc("CMPLXDF"),
                &cval,
                3,
                Some(&cc("complex double fixed")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dblcmp(f, &cc("CMPLXDF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 100.125).abs() <= 0.01);
            assert!((result[1] - 200.875).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffuky_long_string() {
        with_byte_img(|f, status| {
            let longstr = "This is a very long string that exceeds the normal \
                FITS keyword value length limit of 68 characters and requires \
                CONTINUE keywords to store properly.";
            fits_update_key_longstr(
                f,
                &cc("LONGSTR"),
                &cc(longstr),
                Some(&cc("long string")),
                status,
            );
            assert_eq!(*status, 0);
            let result = read_longstr(f, &cc("LONGSTR"), status);
            assert_eq!(result, longstr);
        });
    }

    #[test]
    fn test_ffuky_update_card() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(
                f,
                &cc("OLDKEY"),
                &cc("oldvalue"),
                Some(&cc("old comment")),
                status,
            );
            fits_update_card(
                f,
                &cc("OLDKEY"),
                &cc("OLDKEY  = 'newvalue' / new comment"),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("OLDKEY"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("newvalue"));
        });
    }

    #[test]
    fn test_ffuky_generic() {
        with_byte_img(|f, status| {
            let sval = cc("teststring");
            let bval: u8 = 42;
            let sbval: c_char = -42;
            let usval: u16 = 1000;
            let shval: i16 = -1000;
            let ival: c_int = 123456;
            let uival: u32 = 4000000000;
            let lval: c_int = 1;
            let ulval: libc::c_ulong = 3000000000;
            let ullval: ULONGLONG = 18000000000000000000;
            let lngval: libc::c_long = 1234567890;
            let llval: LONGLONG = 9876543210;
            let fval: f32 = 3.14;
            let dval: f64 = 2.71828;
            let cval: [f32; 2] = [1.0, 2.0];
            let mval: [f64; 2] = [3.0, 4.0];

            let mut rstr = [0 as c_char; FLEN_VALUE];
            let mut rll: LONGLONG = 0;
            let mut rflt: f32 = 0.0;
            let mut rdbl: f64 = 0.0;
            let mut rint: c_int = 0;
            let mut rcmplx: [f32; 2] = [0.0; 2];
            let mut rdmplx: [f64; 2] = [0.0; 2];

            fits_update_key(
                f,
                KeywordDatatype::TSTRING(&sval),
                &cc("GENSTR"),
                Some(&cc("string")),
                status,
            );
            fits_read_key_str(f, &cc("GENSTR"), &mut rstr, None, status);
            assert_eq!(from_buf(&rstr), "teststring");

            fits_update_key(
                f,
                KeywordDatatype::TBYTE(&bval),
                &cc("GENBYTE"),
                Some(&cc("byte")),
                status,
            );
            fits_read_key_lnglng(f, &cc("GENBYTE"), &mut rll, None, status);
            assert_eq!(rll, 42);

            fits_update_key(
                f,
                KeywordDatatype::TSBYTE(&sbval),
                &cc("GENSBYTE"),
                Some(&cc("sbyte")),
                status,
            );
            fits_read_key_lnglng(f, &cc("GENSBYTE"), &mut rll, None, status);
            assert_eq!(rll, -42);

            fits_update_key(
                f,
                KeywordDatatype::TUSHORT(&usval),
                &cc("GENUSHRT"),
                Some(&cc("ushort")),
                status,
            );
            fits_read_key_lnglng(f, &cc("GENUSHRT"), &mut rll, None, status);
            assert_eq!(rll, 1000);

            fits_update_key(
                f,
                KeywordDatatype::TSHORT(&shval),
                &cc("GENSHORT"),
                Some(&cc("short")),
                status,
            );
            fits_read_key_lnglng(f, &cc("GENSHORT"), &mut rll, None, status);
            assert_eq!(rll, -1000);

            fits_update_key(
                f,
                KeywordDatatype::TINT(&ival),
                &cc("GENINT"),
                Some(&cc("int")),
                status,
            );
            fits_read_key_lnglng(f, &cc("GENINT"), &mut rll, None, status);
            assert_eq!(rll, 123456);

            fits_update_key(
                f,
                KeywordDatatype::TUINT(&uival),
                &cc("GENUINT"),
                Some(&cc("uint")),
                status,
            );
            fits_read_key_dbl(f, &cc("GENUINT"), &mut rdbl, None, status);
            assert!((rdbl - 4000000000.0).abs() <= 1.0);

            fits_update_key(
                f,
                KeywordDatatype::TLOGICAL(&lval),
                &cc("GENLOG"),
                Some(&cc("logical")),
                status,
            );
            fits_read_key_log(f, &cc("GENLOG"), &mut rint, None, status);
            assert_eq!(rint, 1);

            fits_update_key(
                f,
                KeywordDatatype::TULONG(&ulval),
                &cc("GENULONG"),
                Some(&cc("ulong")),
                status,
            );
            let mut rull: ULONGLONG = 0;
            fits_read_key_ulnglng(f, &cc("GENULONG"), &mut rull, None, status);
            assert_eq!(rull, 3000000000);

            fits_update_key(
                f,
                KeywordDatatype::TULONGLONG(&ullval),
                &cc("GENULL"),
                Some(&cc("ulonglong")),
                status,
            );
            fits_read_key_ulnglng(f, &cc("GENULL"), &mut rull, None, status);
            assert_eq!(rull, 18000000000000000000);

            fits_update_key(
                f,
                KeywordDatatype::TLONG(&lngval),
                &cc("GENLONG"),
                Some(&cc("long")),
                status,
            );
            fits_read_key_lnglng(f, &cc("GENLONG"), &mut rll, None, status);
            assert_eq!(rll, 1234567890);

            fits_update_key(
                f,
                KeywordDatatype::TLONGLONG(&llval),
                &cc("GENLL"),
                Some(&cc("longlong")),
                status,
            );
            fits_read_key_lnglng(f, &cc("GENLL"), &mut rll, None, status);
            assert_eq!(rll, 9876543210);

            fits_update_key(
                f,
                KeywordDatatype::TFLOAT(&fval),
                &cc("GENFLT"),
                Some(&cc("float")),
                status,
            );
            fits_read_key_flt(f, &cc("GENFLT"), &mut rflt, None, status);
            assert!((rflt - 3.14).abs() <= 0.01);

            fits_update_key(
                f,
                KeywordDatatype::TDOUBLE(&dval),
                &cc("GENDBL"),
                Some(&cc("double")),
                status,
            );
            fits_read_key_dbl(f, &cc("GENDBL"), &mut rdbl, None, status);
            assert!((rdbl - 2.71828).abs() <= 1e-5);

            fits_update_key(
                f,
                KeywordDatatype::TCOMPLEX(&cval),
                &cc("GENCMPLX"),
                Some(&cc("complex")),
                status,
            );
            fits_read_key_cmp(f, &cc("GENCMPLX"), &mut rcmplx, None, status);
            assert!((rcmplx[0] - 1.0).abs() <= 0.01);

            fits_update_key(
                f,
                KeywordDatatype::TDBLCOMPLEX(&mval),
                &cc("GENDCMPLX"),
                Some(&cc("dblcomplex")),
                status,
            );
            fits_read_key_dblcmp(f, &cc("GENDCMPLX"), &mut rdmplx, None, status);
            assert!((rdmplx[0] - 3.0).abs() <= 0.01);

            assert_eq!(*status, 0);
        });
    }

    #[test]
    fn test_ffuky_bad_datatype() {
        with_byte_img(|f, status| {
            let val: c_int = 1;
            fits_update_key(
                f,
                KeywordDatatype::INVALID(9999),
                &cc("BADTYPE"),
                Some(&cc("bad type")),
                status,
            );
            // C passes &val with datatype 9999; the Rust enum represents an
            // unknown datatype as INVALID(9999).
            let _ = val;
            assert_eq!(*status, BAD_DATATYPE);
            *status = 0;
        });
    }

    // ---------------- Modify keyword tests ----------------

    #[test]
    fn test_ffmrec() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(f, &cc("KEY1"), &cc("val1"), Some(&cc("comment1")), status);
            fits_modify_record(f, 7, &cc("KEY1    = 'modified' / new comment"), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("KEY1"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("modified"));
        });
    }

    #[test]
    fn test_ffmcrd() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(
                f,
                &cc("KEY1"),
                &cc("original"),
                Some(&cc("original comment")),
                status,
            );
            fits_modify_card(
                f,
                &cc("KEY1"),
                &cc("KEY1    = 'changed' / changed comment"),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("KEY1"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("changed"));
        });
    }

    #[test]
    fn test_ffmnam() {
        with_byte_img(|f, status| {
            let mut value = [0 as c_char; FLEN_VALUE];
            fits_write_key_str(
                f,
                &cc("OLDNAME"),
                &cc("testval"),
                Some(&cc("comment")),
                status,
            );
            fits_modify_name(f, &cc("OLDNAME"), &cc("NEWNAME"), status);
            assert_eq!(*status, 0);
            fits_read_key_str(f, &cc("NEWNAME"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "testval");
            // Old name should not exist
            fits_read_key_str(f, &cc("OLDNAME"), &mut value, None, status);
            assert_eq!(*status, KEY_NO_EXIST);
            *status = 0;
        });
    }

    #[test]
    fn test_ffmcom() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(
                f,
                &cc("TESTKEY"),
                &cc("value"),
                Some(&cc("old comment")),
                status,
            );
            fits_modify_comment(f, &cc("TESTKEY"), Some(&cc("new comment")), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("TESTKEY"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("new comment"));
        });
    }

    #[test]
    fn test_ffpunt() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_lng(f, &cc("DIST"), 100, Some(&cc("distance")), status);
            fits_write_key_unit(f, &cc("DIST"), &cc("kpc"), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("DIST"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("[kpc]"));
        });
    }

    #[test]
    fn test_ffpunt_replace_units() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_lng(f, &cc("DIST"), 100, Some(&cc("[pc] distance")), status);
            fits_write_key_unit(f, &cc("DIST"), &cc("kpc"), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("DIST"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("[kpc]"));
            assert!(!from_buf(&card).contains("[pc]"));
        });
    }

    #[test]
    fn test_ffpunt_empty() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_lng(f, &cc("DIST"), 100, Some(&cc("[pc] distance")), status);
            fits_write_key_unit(f, &cc("DIST"), &cc(""), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("DIST"), &mut card, status);
            assert_eq!(*status, 0);
            // Units should be removed
            assert!(!from_buf(&card).contains("["));
        });
    }

    #[test]
    fn test_ffpunt_no_existing_units() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_lng(f, &cc("DIST"), 100, Some(&cc("distance no units")), status);
            fits_write_key_unit(f, &cc("DIST"), &cc("kpc"), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("DIST"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("[kpc]"));
            assert!(from_buf(&card).contains("distance"));
        });
    }

    #[test]
    fn test_ffpunt_malformed_units() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            // Write a card with malformed units using ffprec
            fits_write_record(
                f,
                &cc("DIST    =                  100 / [malformed units no close bracket"),
                status,
            );
            fits_write_key_unit(f, &cc("DIST"), &cc("kpc"), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("DIST"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("[kpc]"));
        });
    }

    #[test]
    fn test_ffmkyu() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(f, &cc("KEY1"), &cc("value"), Some(&cc("comment")), status);
            fits_modify_key_null(f, &cc("KEY1"), Some(&cc("now undefined")), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("KEY1"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("now undefined"));
        });
    }

    #[test]
    fn test_ffmkyu_preserve_comment() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(
                f,
                &cc("KEY1"),
                &cc("value"),
                Some(&cc("original comment")),
                status,
            );
            fits_modify_key_null(f, &cc("KEY1"), None, status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("KEY1"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("original comment"));
        });
    }

    #[test]
    fn test_ffmkyu_ampersand_comment() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(f, &cc("KEY1"), &cc("value"), Some(&cc("keep this")), status);
            fits_modify_key_null(f, &cc("KEY1"), Some(&cc("&")), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("KEY1"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("keep this"));
        });
    }

    #[test]
    fn test_ffmkys() {
        with_byte_img(|f, status| {
            let mut value = [0 as c_char; FLEN_VALUE];
            fits_write_key_str(
                f,
                &cc("STRKEY"),
                &cc("oldval"),
                Some(&cc("old comment")),
                status,
            );
            fits_modify_key_str(
                f,
                &cc("STRKEY"),
                &cc("newval"),
                Some(&cc("new comment")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_str(f, &cc("STRKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "newval");
        });
    }

    #[test]
    fn test_ffmkys_preserve_comment() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(
                f,
                &cc("STRKEY"),
                &cc("oldval"),
                Some(&cc("preserve this")),
                status,
            );
            fits_modify_key_str(f, &cc("STRKEY"), &cc("newval"), Some(&cc("&")), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("STRKEY"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("preserve this"));
        });
    }

    #[test]
    fn test_ffmkys_delete_continuation() {
        with_byte_img(|f, status| {
            let longstr = "This is a very long string that exceeds the normal \
                FITS keyword value length and requires CONTINUE keywords to \
                store properly in the header.";
            let mut value = [0 as c_char; FLEN_VALUE];
            let mut nkeys_before: c_int = 0;
            let mut nkeys_after: c_int = 0;

            // Write long string with ffpkls - creates CONTINUE cards
            fits_write_key_longstr(
                f,
                &cc("CONTKEY"),
                &cc(longstr),
                Some(&cc("long string")),
                status,
            );
            fits_get_hdrspace(f, Some(&mut nkeys_before), None, status);
            // Use ffmkys to modify to short value - should delete CONTINUE
            fits_modify_key_str(
                f,
                &cc("CONTKEY"),
                &cc("short"),
                Some(&cc("now short")),
                status,
            );
            fits_get_hdrspace(f, Some(&mut nkeys_after), None, status);
            assert!(nkeys_after < nkeys_before);
            fits_read_key_str(f, &cc("CONTKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "short");
        });
    }

    #[test]
    fn test_ffmkyl() {
        with_byte_img(|f, status| {
            let mut value: c_int = 0;
            fits_write_key_log(f, &cc("LOGKEY"), 0, Some(&cc("was false")), status);
            fits_modify_key_log(f, &cc("LOGKEY"), 1, Some(&cc("now true")), status);
            assert_eq!(*status, 0);
            fits_read_key_log(f, &cc("LOGKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 1);
        });
    }

    #[test]
    fn test_ffmkyj() {
        with_byte_img(|f, status| {
            let mut value: LONGLONG = 0;
            fits_write_key_lng(f, &cc("INTKEY"), 100, Some(&cc("original")), status);
            fits_modify_key_lng(f, &cc("INTKEY"), 200, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_lnglng(f, &cc("INTKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 200);
        });
    }

    #[test]
    fn test_ffmkyuj() {
        with_byte_img(|f, status| {
            let mut value: ULONGLONG = 0;
            fits_write_key_ulng(f, &cc("ULLKEY"), 1000, Some(&cc("original")), status);
            fits_modify_key_ulng(f, &cc("ULLKEY"), 2000, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_ulnglng(f, &cc("ULLKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 2000);
        });
    }

    #[test]
    fn test_ffmkyuj_preserve_comment() {
        with_byte_img(|f, status| {
            let mut value: ULONGLONG = 0;
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_ulng(
                f,
                &cc("ULLKEY"),
                1000,
                Some(&cc("keep this comment")),
                status,
            );
            fits_modify_key_ulng(f, &cc("ULLKEY"), 2000, Some(&cc("&")), status);
            assert_eq!(*status, 0);
            fits_read_key_ulnglng(f, &cc("ULLKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 2000);
            fits_read_card(f, &cc("ULLKEY"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("keep this"));
        });
    }

    #[test]
    fn test_ffmkyf() {
        with_byte_img(|f, status| {
            let mut value: f32 = 0.0;
            fits_write_key_fixflt(f, &cc("FLTKEY"), 1.5, 2, Some(&cc("original")), status);
            fits_modify_key_fixflt(f, &cc("FLTKEY"), 2.5, 2, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_flt(f, &cc("FLTKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 2.5).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffmkye() {
        with_byte_img(|f, status| {
            let mut value: f32 = 0.0;
            fits_write_key_flt(f, &cc("FLTEXP"), 1.5e5, 3, Some(&cc("original")), status);
            fits_modify_key_flt(f, &cc("FLTEXP"), 2.5e5, 3, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_flt(f, &cc("FLTEXP"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 2.5e5).abs() / 2.5e5 <= 0.01);
        });
    }

    #[test]
    fn test_ffmkyg() {
        with_byte_img(|f, status| {
            let mut value: f64 = 0.0;
            fits_write_key_fixdbl(f, &cc("DBLKEY"), 1.234, 4, Some(&cc("original")), status);
            fits_modify_key_fixdbl(f, &cc("DBLKEY"), 5.678, 4, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_dbl(f, &cc("DBLKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 5.678).abs() <= 0.001);
        });
    }

    #[test]
    fn test_ffmkyd() {
        with_byte_img(|f, status| {
            let mut value: f64 = 0.0;
            fits_write_key_dbl(f, &cc("DBLEXP"), 1.234e10, 6, Some(&cc("original")), status);
            fits_modify_key_dbl(f, &cc("DBLEXP"), 5.678e10, 6, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_dbl(f, &cc("DBLEXP"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 5.678e10).abs() / 5.678e10 <= 1e-6);
        });
    }

    #[test]
    fn test_ffmkfc() {
        with_byte_img(|f, status| {
            let cval: [f32; 2] = [1.0, 2.0];
            let newval: [f32; 2] = [3.0, 4.0];
            let mut result: [f32; 2] = [0.0; 2];
            fits_write_key_fixcmp(f, &cc("CMPLXF"), &cval, 2, Some(&cc("original")), status);
            fits_modify_key_fixcmp(f, &cc("CMPLXF"), &newval, 2, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_cmp(f, &cc("CMPLXF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 3.0).abs() <= 0.01);
            assert!((result[1] - 4.0).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffmkfc_preserve_comment() {
        with_byte_img(|f, status| {
            let cval: [f32; 2] = [1.0, 2.0];
            let newval: [f32; 2] = [3.0, 4.0];
            let mut result: [f32; 2] = [0.0; 2];
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_fixcmp(f, &cc("CMPLXF"), &cval, 2, Some(&cc("keep this")), status);
            fits_modify_key_fixcmp(f, &cc("CMPLXF"), &newval, 2, Some(&cc("&")), status);
            assert_eq!(*status, 0);
            fits_read_key_cmp(f, &cc("CMPLXF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 3.0).abs() <= 0.01);
            fits_read_card(f, &cc("CMPLXF"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("keep this"));
        });
    }

    #[test]
    fn test_ffmkyc() {
        with_byte_img(|f, status| {
            let cval: [f32; 2] = [1.0, 2.0];
            let newval: [f32; 2] = [5.0, 6.0];
            let mut result: [f32; 2] = [0.0; 2];
            fits_write_key_cmp(f, &cc("CMPLXE"), &cval, 3, Some(&cc("original")), status);
            fits_modify_key_cmp(f, &cc("CMPLXE"), &newval, 3, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_cmp(f, &cc("CMPLXE"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 5.0).abs() <= 0.01);
            assert!((result[1] - 6.0).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffmkfm() {
        with_byte_img(|f, status| {
            let cval: [f64; 2] = [1.0, 2.0];
            let newval: [f64; 2] = [7.0, 8.0];
            let mut result: [f64; 2] = [0.0; 2];
            fits_write_key_fixdblcmp(f, &cc("CMPLXDF"), &cval, 3, Some(&cc("original")), status);
            fits_modify_key_fixdblcmp(f, &cc("CMPLXDF"), &newval, 3, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_dblcmp(f, &cc("CMPLXDF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 7.0).abs() <= 0.01);
            assert!((result[1] - 8.0).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffmkym() {
        with_byte_img(|f, status| {
            let cval: [f64; 2] = [1.0, 2.0];
            let newval: [f64; 2] = [9.0, 10.0];
            let mut result: [f64; 2] = [0.0; 2];
            fits_write_key_dblcmp(f, &cc("CMPLXDE"), &cval, 4, Some(&cc("original")), status);
            fits_modify_key_dblcmp(f, &cc("CMPLXDE"), &newval, 4, Some(&cc("modified")), status);
            assert_eq!(*status, 0);
            fits_read_key_dblcmp(f, &cc("CMPLXDE"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 9.0).abs() <= 0.01);
            assert!((result[1] - 10.0).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffmkls() {
        with_byte_img(|f, status| {
            let longstr = "This is a very long string that exceeds the normal \
                FITS keyword value length limit and requires CONTINUE keywords.";
            let newstr = "This is a completely different long string value that \
                also requires CONTINUE keywords for proper storage in the header.";
            fits_write_key_longstr(
                f,
                &cc("LONGSTR"),
                &cc(longstr),
                Some(&cc("original")),
                status,
            );
            fits_modify_key_longstr(
                f,
                &cc("LONGSTR"),
                &cc(newstr),
                Some(&cc("modified")),
                status,
            );
            assert_eq!(*status, 0);
            let result = read_longstr(f, &cc("LONGSTR"), status);
            assert_eq!(result, newstr);
        });
    }

    // ---------------- Insert keyword tests ----------------

    #[test]
    fn test_ffirec() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            fits_write_key_str(f, &cc("KEY1"), &cc("val1"), Some(&cc("first key")), status);
            fits_write_key_str(f, &cc("KEY3"), &cc("val3"), Some(&cc("third key")), status);
            // Insert at position 5 (after KEY1, before KEY3)
            fits_insert_record(f, 5, &cc("KEY2    = 'val2' / second key"), status);
            assert_eq!(*status, 0);
            // Read card at position 5 to verify order
            fits_read_record(f, 5, Some(&mut card), status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("KEY2"));
        });
    }

    #[test]
    fn test_ffikey() {
        with_byte_img(|f, status| {
            let mut card = [0 as c_char; FLEN_CARD];
            let mut nkeys: c_int = 0;
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            // Move to position after existing keywords
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_card(f, &cc("NEWKEY  = 'inserted' / inserted card"), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("NEWKEY"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("inserted"));
        });
    }

    #[test]
    fn test_ffikyu() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut card = [0 as c_char; FLEN_CARD];
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_null(f, &cc("UNDEF"), Some(&cc("undefined value")), status);
            assert_eq!(*status, 0);
            fits_read_card(f, &cc("UNDEF"), &mut card, status);
            assert_eq!(*status, 0);
            assert!(from_buf(&card).contains("UNDEF"));
        });
    }

    #[test]
    fn test_ffikys() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut value = [0 as c_char; FLEN_VALUE];
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_str(
                f,
                &cc("STRINS"),
                &cc("inserted"),
                Some(&cc("string insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_str(f, &cc("STRINS"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "inserted");
        });
    }

    #[test]
    fn test_ffikls() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let longstr = "This is a very long string value that needs to be \
                continued over multiple CONTINUE keywords in the FITS header.";
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_longstr(
                f,
                &cc("LONGINS"),
                &cc(longstr),
                Some(&cc("long insert")),
                status,
            );
            assert_eq!(*status, 0);
            let result = read_longstr(f, &cc("LONGINS"), status);
            assert_eq!(result, longstr);
        });
    }

    #[test]
    fn test_ffikyl() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut value: c_int = 0;
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_log(f, &cc("LOGINS"), 1, Some(&cc("logical insert")), status);
            assert_eq!(*status, 0);
            fits_read_key_log(f, &cc("LOGINS"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 1);
        });
    }

    #[test]
    fn test_ffikyj() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut value: LONGLONG = 0;
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_lng(
                f,
                &cc("INTINS"),
                12345678901,
                Some(&cc("integer insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_lnglng(f, &cc("INTINS"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(value, 12345678901);
        });
    }

    #[test]
    fn test_ffikyf() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut value: f32 = 0.0;
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_fixflt(
                f,
                &cc("FLTINS"),
                1.234,
                3,
                Some(&cc("float insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_flt(f, &cc("FLTINS"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 1.234).abs() <= 0.001);
        });
    }

    #[test]
    fn test_ffikye() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut value: f32 = 0.0;
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_flt(
                f,
                &cc("FLTEXP"),
                1.234e5,
                3,
                Some(&cc("float exp insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_flt(f, &cc("FLTEXP"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 1.234e5).abs() / 1.234e5 <= 0.01);
        });
    }

    #[test]
    fn test_ffikyg() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut value: f64 = 0.0;
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_fixdbl(
                f,
                &cc("DBLINS"),
                1.23456789,
                8,
                Some(&cc("double insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dbl(f, &cc("DBLINS"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 1.23456789).abs() <= 1e-8);
        });
    }

    #[test]
    fn test_ffikyd() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let mut value: f64 = 0.0;
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_dbl(
                f,
                &cc("DBLEXP"),
                1.234567e15,
                6,
                Some(&cc("double exp insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dbl(f, &cc("DBLEXP"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert!((value - 1.234567e15).abs() / 1.234567e15 <= 1e-6);
        });
    }

    #[test]
    fn test_ffikfc() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let cval: [f32; 2] = [1.5, 2.5];
            let mut result: [f32; 2] = [0.0; 2];
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_fixcmp(
                f,
                &cc("CMPLXF"),
                &cval,
                2,
                Some(&cc("complex fixed insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_cmp(f, &cc("CMPLXF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 1.5).abs() <= 0.01);
            assert!((result[1] - 2.5).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffikyc() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let cval: [f32; 2] = [1.5e3, 2.5e3];
            let mut result: [f32; 2] = [0.0; 2];
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_cmp(
                f,
                &cc("CMPLXE"),
                &cval,
                3,
                Some(&cc("complex exp insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_cmp(f, &cc("CMPLXE"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 1.5e3).abs() / 1.5e3 <= 0.01);
        });
    }

    #[test]
    fn test_ffikfm() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let cval: [f64; 2] = [1.5, 2.5];
            let mut result: [f64; 2] = [0.0; 2];
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_fixdblcmp(
                f,
                &cc("DCMPLXF"),
                &cval,
                3,
                Some(&cc("dcomplex fixed insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dblcmp(f, &cc("DCMPLXF"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 1.5).abs() <= 0.01);
            assert!((result[1] - 2.5).abs() <= 0.01);
        });
    }

    #[test]
    fn test_ffikym() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let cval: [f64; 2] = [1.5e10, 2.5e10];
            let mut result: [f64; 2] = [0.0; 2];
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_dblcmp(
                f,
                &cc("DCMPLXE"),
                &cval,
                4,
                Some(&cc("dcomplex exp insert")),
                status,
            );
            assert_eq!(*status, 0);
            fits_read_key_dblcmp(f, &cc("DCMPLXE"), &mut result, None, status);
            assert_eq!(*status, 0);
            assert!((result[0] - 1.5e10).abs() / 1.5e10 <= 0.001);
        });
    }

    // ---------------- Delete keyword tests ----------------

    #[test]
    fn test_ffdkey() {
        with_byte_img(|f, status| {
            let mut value = [0 as c_char; FLEN_VALUE];
            fits_write_key_str(
                f,
                &cc("DELKEY"),
                &cc("todelete"),
                Some(&cc("to be deleted")),
                status,
            );
            fits_read_key_str(f, &cc("DELKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "todelete");
            fits_delete_key(f, &cc("DELKEY"), status);
            assert_eq!(*status, 0);
            fits_read_key_str(f, &cc("DELKEY"), &mut value, None, status);
            assert_eq!(*status, KEY_NO_EXIST);
            *status = 0;
        });
    }

    #[test]
    fn test_ffdkey_not_exist() {
        with_byte_img(|f, status| {
            fits_delete_key(f, &cc("NOEXIST"), status);
            assert_eq!(*status, KEY_NO_EXIST);
            *status = 0;
        });
    }

    #[test]
    fn test_ffdstr() {
        with_byte_img(|f, status| {
            fits_write_record(f, &cc("COMMENT This is a test comment to delete"), status);
            fits_delete_str(f, &cc("test comment"), status);
            assert_eq!(*status, 0);
            fits_delete_str(f, &cc("test comment"), status);
            assert_eq!(*status, KEY_NO_EXIST);
            *status = 0;
        });
    }

    #[test]
    fn test_ffdrec() {
        with_byte_img(|f, status| {
            let mut value = [0 as c_char; FLEN_VALUE];
            let mut nkeys: c_int = 0;
            fits_write_key_str(f, &cc("KEY1"), &cc("val1"), Some(&cc("first")), status);
            fits_write_key_str(f, &cc("KEY2"), &cc("val2"), Some(&cc("second")), status);
            fits_write_key_str(f, &cc("KEY3"), &cc("val3"), Some(&cc("third")), status);
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            // Delete KEY2 (at position 8, after SIMPLE, BITPIX, NAXIS, EXTEND, 2 COMMENT, KEY1)
            fits_delete_record(f, 8, status);
            assert_eq!(*status, 0);
            // KEY3 should still exist
            fits_read_key_str(f, &cc("KEY3"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "val3");
            // KEY2 should not exist
            fits_read_key_str(f, &cc("KEY2"), &mut value, None, status);
            assert_eq!(*status, KEY_NO_EXIST);
            *status = 0;
        });
    }

    #[test]
    fn test_ffdrec_out_of_bounds() {
        with_byte_img(|f, status| {
            fits_delete_record(f, 0, status);
            assert_eq!(*status, KEY_OUT_BOUNDS);
            *status = 0;
            fits_delete_record(f, 1000, status);
            assert_eq!(*status, KEY_OUT_BOUNDS);
            *status = 0;
        });
    }

    #[test]
    fn test_ffdkey_continued_string() {
        with_byte_img(|f, status| {
            let longstr = "This is a very long string that requires CONTINUE \
                keywords to store properly in the FITS header format.";
            fits_write_key_longstr(
                f,
                &cc("LONGKEY"),
                &cc(longstr),
                Some(&cc("long string")),
                status,
            );
            fits_delete_key(f, &cc("LONGKEY"), status);
            assert_eq!(*status, 0);
            // Verify the key and all CONTINUE cards are gone
            fits_delete_key(f, &cc("LONGKEY"), status);
            assert_eq!(*status, KEY_NO_EXIST);
            *status = 0;
        });
    }

    // ---------------- Additional coverage tests ----------------

    #[test]
    fn test_ffmcrd_undefined_value() {
        with_byte_img(|f, status| {
            let mut value = [0 as c_char; FLEN_VALUE];
            // Write a card with no value (just keyword and comment)
            fits_write_record(
                f,
                &cc("TESTKEY =                      / comment only"),
                status,
            );
            // Modify the card
            fits_modify_card(
                f,
                &cc("TESTKEY"),
                &cc("TESTKEY = 'newvalue' / new comment"),
                status,
            );
            assert_eq!(*status, 0);
            // Verify the modification
            fits_read_key_str(f, &cc("TESTKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "newvalue");
        });
    }

    #[test]
    fn test_ffmcrd_continued_string() {
        with_byte_img(|f, status| {
            let longstr = "This is a very long string that is more than \
                sixty-eight characters and requires CONTINUE keywords to \
                store in the FITS header properly.";
            let mut value = [0 as c_char; FLEN_VALUE];
            let mut nkeys_before: c_int = 0;
            let mut nkeys_after: c_int = 0;

            // Write a long string that uses continuation
            fits_write_key_longstr(
                f,
                &cc("LONGKEY"),
                &cc(longstr),
                Some(&cc("long string")),
                status,
            );
            // Count keywords - should include CONTINUE cards
            fits_get_hdrspace(f, Some(&mut nkeys_before), None, status);
            // Verify long string was written correctly
            let result = read_longstr(f, &cc("LONGKEY"), status);
            assert_eq!(result, longstr);
            // Modify the keyword to a short value - this should delete CONTINUE
            fits_modify_card(
                f,
                &cc("LONGKEY"),
                &cc("LONGKEY = 'short' / now short"),
                status,
            );
            // Count keywords again - should be fewer now
            fits_get_hdrspace(f, Some(&mut nkeys_after), None, status);
            assert!(nkeys_after < nkeys_before);
            fits_read_key_str(f, &cc("LONGKEY"), &mut value, None, status);
            assert_eq!(*status, 0);
            assert_eq!(from_buf(&value), "short");
        });
    }

    #[test]
    fn test_error_status_returns() {
        with_byte_img(|f, status| {
            let card = cc("TEST    = 'value'");
            let fcval: [f32; 2] = [1.0, 2.0];
            let dcval: [f64; 2] = [1.0, 2.0];

            // Test ffuky* functions with error status
            *status = 1;
            fits_update_key_str(f, &cc("KEY"), &cc("val"), None, status);
            assert_eq!(*status, 1);
            fits_update_key_log(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_update_key_lng(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_update_key_ulng(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_update_key_fixflt(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_update_key_flt(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_update_key_fixdbl(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_update_key_dbl(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_update_key_null(f, &cc("KEY"), None, status);
            assert_eq!(*status, 1);
            fits_update_key_cmp(f, &cc("KEY"), &fcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_update_key_dblcmp(f, &cc("KEY"), &dcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_update_key_fixcmp(f, &cc("KEY"), &fcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_update_key_fixdblcmp(f, &cc("KEY"), &dcval, 2, None, status);
            assert_eq!(*status, 1);

            // Test ffmky* functions with error status
            fits_modify_record(f, 1, &card, status);
            assert_eq!(*status, 1);
            fits_modify_card(f, &cc("KEY"), &card, status);
            assert_eq!(*status, 1);
            fits_modify_name(f, &cc("KEY"), &cc("NEW"), status);
            assert_eq!(*status, 1);
            fits_modify_comment(f, &cc("KEY"), Some(&cc("comment")), status);
            assert_eq!(*status, 1);
            fits_modify_key_null(f, &cc("KEY"), None, status);
            assert_eq!(*status, 1);
            fits_modify_key_str(f, &cc("KEY"), &cc("val"), None, status);
            assert_eq!(*status, 1);
            fits_modify_key_log(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_lng(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_ulng(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_fixflt(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_flt(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_fixdbl(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_dbl(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_fixcmp(f, &cc("KEY"), &fcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_cmp(f, &cc("KEY"), &fcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_fixdblcmp(f, &cc("KEY"), &dcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_modify_key_dblcmp(f, &cc("KEY"), &dcval, 2, None, status);
            assert_eq!(*status, 1);

            // Test ffiky* functions with error status
            fits_insert_record(f, 1, &card, status);
            assert_eq!(*status, 1);
            fits_insert_card(f, &card, status);
            assert_eq!(*status, 1);
            fits_insert_key_null(f, &cc("KEY"), None, status);
            assert_eq!(*status, 1);
            fits_insert_key_str(f, &cc("KEY"), &cc("val"), None, status);
            assert_eq!(*status, 1);
            fits_insert_key_log(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_lng(f, &cc("KEY"), 1, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_fixflt(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_flt(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_fixdbl(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_dbl(f, &cc("KEY"), 1.0, 2, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_fixcmp(f, &cc("KEY"), &fcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_cmp(f, &cc("KEY"), &fcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_fixdblcmp(f, &cc("KEY"), &dcval, 2, None, status);
            assert_eq!(*status, 1);
            fits_insert_key_dblcmp(f, &cc("KEY"), &dcval, 2, None, status);
            assert_eq!(*status, 1);

            // Test ffd* functions with error status
            fits_delete_key(f, &cc("KEY"), status);
            assert_eq!(*status, 1);
            fits_delete_record(f, 1, status);
            assert_eq!(*status, 1);

            *status = 0;
        });
    }

    #[test]
    fn test_complex_value_overflow() {
        with_byte_img(|f, status| {
            let fcval: [f32; 2] = [1.0, 2.0];
            let dcval: [f64; 2] = [1.0, 2.0];

            // Create a complex key first
            fits_write_key_cmp(
                f,
                &cc("CPLXKEY"),
                &fcval,
                2,
                Some(&cc("test complex")),
                status,
            );
            assert_eq!(*status, 0);

            // Try to modify with very high precision - should cause overflow
            fits_modify_key_fixcmp(f, &cc("CPLXKEY"), &fcval, 50, None, status);
            assert_eq!(*status, BAD_F2C);
            *status = 0;

            // Similar for ffmkyc
            fits_modify_key_cmp(f, &cc("CPLXKEY"), &fcval, 50, None, status);
            assert_eq!(*status, BAD_F2C);
            *status = 0;

            // Create a double complex key
            fits_write_key_dblcmp(
                f,
                &cc("DCPLXKEY"),
                &dcval,
                2,
                Some(&cc("test complex double")),
                status,
            );
            assert_eq!(*status, 0);

            // Try to modify with very high precision
            fits_modify_key_fixdblcmp(f, &cc("DCPLXKEY"), &dcval, 50, None, status);
            assert_eq!(*status, BAD_F2C);
            *status = 0;

            fits_modify_key_dblcmp(f, &cc("DCPLXKEY"), &dcval, 50, None, status);
            assert_eq!(*status, BAD_F2C);
            *status = 0;
        });
    }

    #[test]
    fn test_insert_long_string_with_quotes() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            // String with embedded quotes - tests quote counting loop
            let value = "This string has 'single quotes' and more 'quotes' here.";
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_longstr(f, &cc("LONGSTR"), &cc(value), None, status);
            assert_eq!(*status, 0);
            let result = read_longstr(f, &cc("LONGSTR"), status);
            assert_eq!(result, value);
        });
    }

    #[test]
    fn test_insert_long_string_ending_with_quote() {
        with_byte_img(|f, status| {
            let mut nkeys: c_int = 0;
            let value = "This is a very long string that will need to \
                continue onto the next card and has 'quote' near the end";
            fits_get_hdrspace(f, Some(&mut nkeys), None, status);
            fits_movabs_key(f, nkeys + 1, status);
            fits_insert_key_longstr(f, &cc("QTEST"), &cc(value), None, status);
            assert_eq!(*status, 0);
            let result = read_longstr(f, &cc("QTEST"), status);
            assert_eq!(result, value);
        });
    }
}
