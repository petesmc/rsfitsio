use crate::c_types::{c_char, c_int, size_t};
use bytemuck::cast_slice;
use printf::{CustomVaList, VaArg};
use std::ffi::{c_void, CStr};

use crate::relibc::platform;

mod lookaheadreader;
mod printf;
mod scanf;


pub(crate) fn sprintf_f64(s: &mut [c_char], format: &[c_char], val: f64) -> c_int {

    unsafe {
        let n = s.len();
        let mut valist = CustomVaList::new();
        valist.push(VaArg::c_double(val));

        printf::printf(
            &mut platform::StringWriter(s.as_mut_ptr() as *mut u8, n),
            CStr::from_bytes_until_nul(cast_slice(format)).unwrap(),
            valist,
        )
    }
}

pub(crate) fn sprintf_string_width(
    s: &mut [c_char],
    format: &[c_char],
    width: c_int,
    val: &[c_char],
) -> c_int {
    unsafe {
        let n = s.len();
        let mut valist = CustomVaList::new();
        valist.push(VaArg::c_int(width));
        valist.push(VaArg::pointer(val.as_ptr() as *const c_void));

        printf::printf(
            &mut platform::StringWriter(s.as_mut_ptr() as *mut u8, n),
            CStr::from_bytes_until_nul(cast_slice(format)).unwrap(),
            valist,
        )
    }
}

pub(crate) fn snprintf_f64(s: &mut [c_char], n: size_t, format: &[c_char], val: f64) -> c_int {
    unsafe {
        let mut valist = CustomVaList::new();
        valist.push(VaArg::c_double(val));

        printf::printf(
            &mut platform::StringWriter(s.as_mut_ptr() as *mut u8, n),
            CStr::from_bytes_until_nul(cast_slice(format)).unwrap(),
            valist,
        )
    }
}

pub(crate) fn snprintf_cint(s: &mut [c_char], n: size_t, format: &[c_char], val: c_int) -> c_int {
    unsafe {
        let mut valist = CustomVaList::new();
        valist.push(VaArg::c_int(val));

        printf::printf(
            &mut platform::StringWriter(s.as_mut_ptr() as *mut u8, n),
            CStr::from_bytes_until_nul(cast_slice(format)).unwrap(),
            valist,
        )
    }
}

pub(crate) fn snprintf_f64_decim(
    s: &mut [c_char],
    n: size_t,
    format: &[c_char],
    decim: c_int,
    val: f64,
) -> c_int {
    unsafe {
        let mut valist = CustomVaList::new();
        valist.push(VaArg::c_int(decim));
        valist.push(VaArg::c_double(val));

        printf::printf(
            &mut platform::StringWriter(s.as_mut_ptr() as *mut u8, n),
            CStr::from_bytes_until_nul(cast_slice(format)).unwrap(),
            valist,
        )
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_printf() {
        let format = c"%3d";
        let mut buffer: [c_char; 100] = [0; 100];
        let result = snprintf_cint(&mut buffer, 100, cast_slice(format.to_bytes_with_nul()), 42);
        
        assert_eq!(&buffer[..(result + 1) as usize], cast_slice(c" 42".to_bytes_with_nul()));
    }
}
