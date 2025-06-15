/*  This file, modkey.c, contains routines that modify, insert, or update  */
/*  keywords in a FITS header.                                             */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use std::cmp;
use std::ffi::CStr;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::c_types::*;
use crate::fitscore::{
    ffc2s, ffcmrk_safe, ffgmsg_safe, ffiblk, ffmahd_safe, ffmkey, ffmkky_safe, ffpmrk_safe,
    ffpmsg_slice, ffpsvc_safe, fftkey_safe, fits_strncasecmp,
};
use crate::getkey::{ffgcnt, ffgcrd_safe, ffgkey_safe, ffmaky_safe};
use crate::nullable_slice_cstr;
use crate::putkey::{
    ffd2e, ffd2f, ffi2c, ffl2c, ffpkfc_safe, ffpkfm_safe, ffpkls_safe, ffpkyc_safe, ffpkyd_safe,
    ffpkye_safe, ffpkyf_safe, ffpkyg_safe, ffpkyj_safe, ffpkyl_safe, ffpkym_safe, ffpkys_safe,
    ffpkyu_safe, ffpkyuj_safe, ffprec_safe, ffr2e, ffr2f, ffs2c, ffu2c,
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
            ffukyj_safe(fptr, keyname, *(value) as LONGLONG, comm, status);
        }
        KeywordDatatype::TSBYTE(value) => {
            ffukyj_safe(fptr, keyname, *(value) as LONGLONG, comm, status);
        }
        KeywordDatatype::TUSHORT(value) => {
            ffukyj_safe(fptr, keyname, *(value) as LONGLONG, comm, status);
        }
        KeywordDatatype::TSHORT(value) => {
            ffukyj_safe(fptr, keyname, *(value) as LONGLONG, comm, status);
        }
        KeywordDatatype::TINT(value) => {
            ffukyj_safe(fptr, keyname, *(value) as LONGLONG, comm, status);
        }
        KeywordDatatype::TUINT(value) => {
            ffukyg_safe(fptr, keyname, *(value) as f64, 0, comm, status);
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
            if valstring[0] == 0 || strlen_safe(&nextcomm) != 0 {
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
    todo!();
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
    todo!();
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
    todo!();
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
    todo!();
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
    todo!();
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
    todo!();
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
    todo!();
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
    todo!();
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
    todo!();
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
    let valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    todo!();

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
    let keylength = 0;
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

    for ii in 0..(nshift as usize) {
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
    todo!()
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
        || (keypos as LONGLONG) > fptr.Fptr.headend - headstart[fptr.Fptr.curhdu as usize] / 80
    {
        *status = KEY_OUT_BOUNDS;
        return *status;
    }

    fptr.Fptr.nextkey = headstart[fptr.Fptr.curhdu as usize] + (keypos as LONGLONG - 1) * 80;

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
    for ii in 0..(nshift as usize) {
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
