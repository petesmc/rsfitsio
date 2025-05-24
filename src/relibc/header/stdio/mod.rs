use crate::c_types::{c_char, c_int, size_t};
use bytemuck::cast_slice;
use std::ffi::CStr;

use crate::relibc::platform;

mod lookaheadreader;
mod printf;
mod scanf;

pub(crate) unsafe extern "C" fn snprintf(
    s: *mut c_char,
    n: size_t,
    format: *const c_char,
    mut __valist: ...
) -> c_int {
    unsafe {
        printf::printf(
            &mut platform::StringWriter(s as *mut u8, n),
            CStr::from_ptr(format),
            __valist.as_va_list(),
        )
    }
}
#[allow(improper_ctypes_definitions)] // Needs to be extern to get access to variadics
pub(crate) unsafe extern "C" fn snprintf_safer(
    s: &mut [c_char],
    n: size_t,
    format: &[c_char],
    mut __valist: ...
) -> c_int {
    unsafe {
        printf::printf(
            // TODO: Why does SliceWriter not work here?
            //&mut platform::SliceWriter(cast_slice_mut(s), n),
            &mut platform::UnsafeStringWriter(s.as_mut_ptr() as *mut u8),
            CStr::from_bytes_until_nul(cast_slice(format)).unwrap(),
            __valist.as_va_list(),
        )
    }
}

pub(crate) unsafe extern "C" fn sprintf(
    s: *mut c_char,
    format: *const c_char,
    mut __valist: ...
) -> c_int {
    unsafe {
        printf::printf(
            &mut platform::UnsafeStringWriter(s as *mut u8),
            CStr::from_ptr(format),
            __valist.as_va_list(),
        )
    }
}

#[allow(improper_ctypes_definitions)] // Needs to be extern to get access to variadics
pub(crate) unsafe extern "C" fn sprintf_safer(
    s: &mut [c_char],
    format: &[c_char],
    mut __valist: ...
) -> c_int {
    unsafe {
        let l = s.len();
        printf::printf(
            // TODO: Why does SliceWriter not work here?
            // &mut platform::SliceWriter(cast_slice_mut(s), l),
            &mut platform::UnsafeStringWriter(s.as_mut_ptr() as *mut u8),
            CStr::from_bytes_until_nul(cast_slice(format)).unwrap(),
            __valist.as_va_list(),
        )
    }
}

pub(crate) fn sprintf_f64(s: &mut [c_char], format: &[c_char], val: f64) -> c_int {
    unsafe { sprintf_safer(s, format, val) }
}

pub(crate) fn sprintf_string_width(
    s: &mut [c_char],
    format: &[c_char],
    width: c_int,
    val: &[c_char],
) -> c_int {
    unsafe { sprintf_safer(s, format, width, val.as_ptr()) }
}

pub(crate) fn sprintf_string(s: &mut [c_char], format: &[c_char], val: &[c_char]) -> c_int {
    unsafe { sprintf_safer(s, format, val.as_ptr()) }
}

pub(crate) fn snprintf_f64(s: &mut [c_char], n: size_t, format: &[c_char], val: f64) -> c_int {
    unsafe { snprintf_safer(s, n, format, val) }
}

pub(crate) fn snprintf_cint(s: &mut [c_char], n: size_t, format: &[c_char], val: c_int) -> c_int {
    unsafe { snprintf_safer(s, n, format, val) }
}

pub(crate) unsafe extern "C" fn sscanf(
    s: *const c_char,
    format: *const c_char,
    mut __valist: ...
) -> c_int {
    unsafe {
        let reader = (s as *const u8).into();
        scanf::scanf(reader, format, __valist.as_va_list())
    }
}
