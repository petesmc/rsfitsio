use std::{cmp, ffi::CStr, ptr, str::FromStr};

use bytemuck::cast_slice;
use cbitset::BitSet256;

use crate::c_types::{c_char, c_int, c_uchar, size_t};
use libc::atoi;

use crate::bb;

// Several wrappers in this file were taken from the Redoc Relibc
// project. The original source code can be found at:
// https://gitlab.redox-os.org/redox-os/relibc

pub(crate) unsafe fn strcpy(dst: *mut c_char, src: *const c_char) -> *mut c_char {
    unsafe {
        let mut i = 0;

        loop {
            let byte = *src.offset(i);
            *dst.offset(i) = byte;

            if byte == 0 {
                break;
            }

            i += 1;
        }

        dst
    }
}

pub(crate) unsafe fn strncpy(dst: *mut c_char, src: *const c_char, n: size_t) -> *mut c_char {
    unsafe {
        let mut i = 0;

        while *src.add(i) != 0 && i < n {
            *dst.add(i) = *src.add(i);
            i += 1;
        }

        for i in i..n {
            *dst.add(i) = 0;
        }

        dst
    }
}

pub fn strncpy_safe(dst: &mut [c_char], src: &[c_char], n: usize) {
    let mut i = 0;

    assert!(n <= dst.len());

    while src[i] != 0 && i < n {
        dst[i] = src[i];
        i += 1;
    }

    for i in i..n {
        dst[i] = 0;
    }
}

pub(crate) fn strcspn_safe(cs: &[c_char], ct: &[c_char]) -> usize {
    strcspn_inner(cs, ct, false)
}

fn strcspn_inner(cs: &[c_char], ct: &[c_char], cmp: bool) -> usize {
    let mut it = 0;
    while ct[it] != 0 {
        it += 1;
    }

    let mut is = 0;
    while cs[is] != 0 {
        if ct[..it].contains(&cs[is]) != cmp {
            return is;
        }
        is += 1;
    }

    is
}

pub(crate) fn strspn_safe(cs: &[c_char], ct: &[c_char]) -> usize {
    strcspn_inner(cs, ct, true)
}

pub fn strcpy_safe(dst: &mut [c_char], src: &[c_char]) {
    let mut i = 0;
    loop {
        dst[i] = src[i];
        if src[i] == 0 {
            break;
        }
        i += 1;
    }
}

pub(crate) fn strcmp_safe(cs: &[c_char], ct: &[c_char]) -> c_int {
    strncmp_safe(cs, ct, cmp::max(cs.len(), ct.len()))
}

pub fn strncmp_safe(cs: &[c_char], ct: &[c_char], n: usize) -> c_int {
    let min_len = cmp::min(cs.len(), ct.len());
    let min_n = cmp::min(n, min_len);
    for i in 0..min_n {
        let a = cs[i];
        let b = ct[i];

        if a != b || a == 0 {
            return (a as c_int) - (b as c_int);
        }
    }

    // If we reached here, all bytes so far are equal but we haven't reached
    // the end of either string. So we need to check if one of them is
    // shorter than the other.
    if min_len < n {
        // We reached the end of one string, but not the other
        return if cs.len() < ct.len() {
            -(ct[min_len] as c_int)
        } else {
            cs[min_len] as c_int
        };
    }

    // If we reached here, all bytes are equal and we reached the end of both

    0
}

pub(crate) unsafe fn strcmp(s1: *const c_char, s2: *const c_char) -> c_int {
    unsafe { strncmp(s1, s2, 1024) }
}

pub(crate) unsafe fn strncmp(s1: *const c_char, s2: *const c_char, n: size_t) -> c_int {
    unsafe {
        let s1 = core::slice::from_raw_parts(s1 as *const c_uchar, n);
        let s2 = core::slice::from_raw_parts(s2 as *const c_uchar, n);

        for (&a, &b) in s1.iter().zip(s2.iter()) {
            let val = (a as c_int) - (b as c_int);
            if a != b || a == 0 {
                return val;
            }
        }

        0
    }
}

pub(crate) fn strlen_safe_cstr(str: &CStr) -> usize {
    str.to_bytes().len()
}

pub(crate) fn strlen_safe(cs: &[c_char]) -> usize {
    let mut ii = 0;
    let len = cs.len();
    while ii < len {
        if cs[ii] == 0 {
            return ii;
        }
        ii += 1;
    }

    panic!("Invalid C Style String");
    len
}

pub(crate) unsafe fn strlen(s: *const c_char) -> size_t {
    unsafe { strnlen(s, usize::MAX) }
}

pub(crate) unsafe fn strnlen(s: *const c_char, size: size_t) -> size_t {
    unsafe {
        let mut i = 0;
        while i < size {
            if *s.add(i) == 0 {
                break;
            }
            i += 1;
        }
        i as size_t
    }
}

pub(crate) fn strtol_safe<F: FromStr>(input: &[c_char]) -> Result<(F, usize), <F as FromStr>::Err> {
    let strlen = input.len() - 1;
    let input: &[u8] = cast_slice(input);

    // Find the first non-numeric character or the end of the string
    let mut start = 0;
    let mut end = 0;

    while start < strlen && (input[start] as char).is_whitespace() {
        start += 1;
    }

    while start < strlen
        && !(input[start] as char).is_numeric()
        && (input[start] as char) != '+'
        && (input[start] as char) != '-'
    {
        start += 1;
    }

    end = start + 1;

    while end < strlen && (input[end] as char).is_numeric() {
        end += 1;
    }

    let str = str::from_utf8(&input[start..end]).unwrap();

    let res = str.parse::<F>()?;

    Ok((res, end))
}

pub(crate) fn strtod_safe(s: &[c_char], endp: &mut usize) -> f64 {
    strto_float_impl(s, endp)
}

pub(crate) unsafe fn strchr(mut s: *const c_char, c: c_int) -> *mut c_char {
    unsafe {
        let c = c as c_char;
        while *s != 0 {
            if *s == c {
                return s as *mut c_char;
            }
            s = s.offset(1);
        }
        ptr::null_mut()
    }
}

pub fn strchr_safe(cs: &[c_char], c: c_char) -> Option<usize> {
    cs.iter().position(|&x| x == c)
}

unsafe fn inner_strstr(
    mut haystack: *const c_char,
    needle: *const c_char,
    mask: c_char,
) -> *mut c_char {
    unsafe {
        while *haystack != 0 {
            let mut i = 0;
            loop {
                if *needle.offset(i) == 0 {
                    // We reached the end of the needle, everything matches this far
                    return haystack as *mut c_char;
                }
                if *haystack.offset(i) & mask != *needle.offset(i) & mask {
                    break;
                }

                i += 1;
            }

            haystack = haystack.offset(1);
        }
        ptr::null_mut()
    }
}

pub(crate) unsafe fn strstr(haystack: *const c_char, needle: *const c_char) -> *mut c_char {
    unsafe { inner_strstr(haystack, needle, !0) }
}

pub(crate) fn strstr_safe(cs: &[c_char], ct: &[c_char]) -> Option<usize> {
    unsafe {
        let t = strstr(cs.as_ptr(), ct.as_ptr());

        if t.is_null() {
            None
        } else {
            Some(t.offset_from(cs.as_ptr()) as usize)
        }
    }
}

pub(crate) unsafe fn strcat(s1: *mut c_char, s2: *const c_char) -> *mut c_char {
    unsafe { strncat(s1, s2, usize::MAX) }
}

pub(crate) unsafe fn strncat(s1: *mut c_char, s2: *const c_char, n: size_t) -> *mut c_char {
    unsafe {
        let len = strlen(s1 as *const c_char);
        let mut i = 0;
        while i < n {
            let b = *s2.add(i);
            if b == 0 {
                break;
            }

            *s1.add(len + i) = b;
            i += 1;
        }
        *s1.add(len + i) = 0;

        s1
    }
}

pub fn strcat_safe(s: &mut [c_char], ct: &[c_char]) {
    let ct_len = strlen_safe(ct);

    strncat_safe(s, ct, ct_len);
}

pub(crate) fn strncat_safe(s: &mut [c_char], ct: &[c_char], n: usize) {
    let s_len = strlen_safe(s);
    let ct_len = strlen_safe(ct);

    let n = cmp::min(n, ct_len);
    let mut i = 0;
    while i < n {
        let b = ct[i];
        if b == 0 {
            break;
        }

        s[s_len + i] = b;
        i += 1;
    }
    s[s_len + i] = 0;
}

pub(crate) fn toupper(c: c_char) -> c_char {
    if islower(c) {
        return c & 0x5f;
    }
    c
}

pub(crate) fn islower(c: c_char) -> bool {
    (c >= 97) && (c <= 122)
}

pub(crate) fn isupper(c: c_char) -> bool {
    (c >= 65) && (c <= 90)
}

pub(crate) fn isspace(c: c_char) -> bool {
    let c = c as c_char;
    c == bb(b' ') || c == bb(b'\t') || c == bb(b'\n') || c == bb(b'\r') || c == 0x0b || c == 0x0c
}

pub(crate) fn atoi_safe(cs: &[c_char]) -> c_int {
    unsafe { atoi(cs.as_ptr()) }
}

pub(crate) fn atof_safe(cs: &[c_char]) -> f64 {
    let mut dummy = 0;
    strtod_safe(cs, &mut dummy)
}

pub(crate) fn isdigit_safe(c: c_char) -> bool {
    (b'0' as c_char) <= c && c <= (b'9' as c_char)
}

pub(crate) fn ffstrtok(
    str: *mut c_char,
    delim: *const c_char,
    saveptr: *mut *mut c_char,
) -> *mut c_char {
    unsafe { strtok_r(str, delim, saveptr) }
}

// Copied from Redox's Relibc: https://gitlab.redox-os.org/redox-os/relibc/-/blob/master/src/header/string/mod.rs
pub(crate) unsafe fn strtok_r(
    s: *mut c_char,
    delimiter: *const c_char,
    lasts: *mut *mut c_char,
) -> *mut c_char {
    unsafe {
        // Loosely based on GLIBC implementation
        let mut haystack = s;
        if haystack.is_null() {
            if (*lasts).is_null() {
                return ptr::null_mut();
            }
            haystack = *lasts;
        }

        // Skip past any extra delimiter left over from previous call
        haystack = haystack.add(strspn(haystack, delimiter));
        if *haystack == 0 {
            *lasts = ptr::null_mut();
            return ptr::null_mut();
        }

        // Build token by injecting null byte into delimiter
        let token = haystack;
        haystack = strpbrk(token, delimiter);
        if !haystack.is_null() {
            haystack.write(0);
            haystack = haystack.add(1);
            *lasts = haystack;
        } else {
            *lasts = ptr::null_mut();
        }

        token
    }
}

// Copied from Redox's Relibc: https://gitlab.redox-os.org/redox-os/relibc/-/blob/master/src/header/string/mod.rs
unsafe fn strpbrk(s1: *const c_char, s2: *const c_char) -> *mut c_char {
    unsafe {
        let p = s1.add(strcspn(s1, s2));
        if *p != 0 {
            p as *mut c_char
        } else {
            ptr::null_mut()
        }
    }
}

// Copied from Redox's Relibc: https://gitlab.redox-os.org/redox-os/relibc/-/blob/master/src/header/string/mod.rs
unsafe fn strspn(s1: *const c_char, s2: *const c_char) -> size_t {
    unsafe { inner_strspn(s1, s2, true) }
}

// Copied from Redox's Relibc: https://gitlab.redox-os.org/redox-os/relibc/-/blob/master/src/header/string/mod.rs
unsafe fn strcspn(s1: *const c_char, s2: *const c_char) -> size_t {
    unsafe { inner_strspn(s1, s2, false) }
}

// Copied from Redox's Relibc: https://gitlab.redox-os.org/redox-os/relibc/-/blob/master/src/header/string/mod.rs
unsafe fn inner_strspn(s1: *const c_char, s2: *const c_char, cmp: bool) -> size_t {
    unsafe {
        let mut s1 = s1 as *const u8;
        let mut s2 = s2 as *const u8;

        // The below logic is effectively ripped from the musl implementation. It
        // works by placing each byte as it's own bit in an array of numbers. Each
        // number can hold up to 8 * mem::size_of::<usize>() bits. We need 256 bits
        // in total, to fit one byte.

        let mut set = BitSet256::new();

        while *s2 != 0 {
            set.insert(*s2 as usize);
            s2 = s2.offset(1);
        }

        let mut i = 0;
        while *s1 != 0 {
            if set.contains(*s1 as usize) != cmp {
                break;
            }
            i += 1;
            s1 = s1.offset(1);
        }
        i
    }
}

// Copied and modified from Redox's Relibc: https://gitlab.redox-os.org/redox-os/relibc/-/blob/master/src/macros.rs
fn strto_float_impl(s: &[c_char], endptr: &mut usize) -> f64 {
    let mut s = s;

    let mut si = 0;
    while isspace(s[si]) {
        si += 1;
    }

    let mut result: f64 = 0.0;
    let mut exponent: Option<f64> = None;
    let mut radix = 10;

    let result_sign = match s[si] as u8 {
        b'-' => {
            si += 1;
            -1.0
        }
        b'+' => {
            si += 1;
            1.0
        }
        _ => 1.0,
    };

    let rust_s = unsafe { CStr::from_ptr(s.as_ptr()).to_string_lossy() };

    // detect NaN, Inf
    if rust_s.to_lowercase().starts_with("inf") {
        result = f64::INFINITY;
        si += 3;
    } else if rust_s.to_lowercase().starts_with("nan") {
        // we cannot signal negative NaN in LLVM backed languages
        // https://github.com/rust-lang/rust/issues/73328 , https://github.com/rust-lang/rust/issues/81261
        result = f64::NAN;
        si += 3;
    } else {
        if s[si] as u8 == b'0' && s[si + 1] as u8 == b'x' {
            si += 2;
            radix = 16;
        }

        while let Some(digit) = (s[si] as u8 as char).to_digit(radix) {
            result *= radix as f64;
            result += digit as f64;
            si += 1;
        }

        if s[si] as u8 == b'.' {
            si += 1;

            let mut i = 1.0;
            while let Some(digit) = (s[si] as u8 as char).to_digit(radix) {
                i *= radix as f64;
                result += digit as f64 / i;
                si += 1;
            }
        }

        let s_before_exponent = s;

        exponent = match (s[si] as u8, radix) {
            (b'e' | b'E', 10) | (b'p' | b'P', 16) => {
                si += 1;

                let is_exponent_positive = match s[si] as u8 {
                    b'-' => {
                        si += 1;
                        false
                    }
                    b'+' => {
                        si += 1;
                        true
                    }
                    _ => true,
                };

                // Exponent digits are always in base 10.
                if (s[si] as u8 as char).is_ascii_digit() {
                    let mut exponent_value = 0;

                    while let Some(digit) = (s[si] as u8 as char).to_digit(10) {
                        exponent_value *= 10;
                        exponent_value += digit;
                        si += 1;
                    }

                    let exponent_base = match radix {
                        10 => 10u128,
                        16 => 2u128,
                        _ => unreachable!(),
                    };

                    if is_exponent_positive {
                        Some(exponent_base.pow(exponent_value) as f64)
                    } else {
                        Some(1.0 / (exponent_base.pow(exponent_value) as f64))
                    }
                } else {
                    // Exponent had no valid digits after 'e'/'p' and '+'/'-', rollback
                    s = s_before_exponent;
                    None
                }
            }
            _ => None,
        };
    }

    *endptr = si;

    if let Some(exponent) = exponent {
        result_sign * result * exponent
    } else {
        result_sign * result
    }
}

#[cfg(test)]
mod tests {
    use libc::c_longlong;

    use crate::wrappers::*;

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_strtol() {
        let s: &[libc::c_char] = bytemuck::cast_slice(b"12345 XC");
        let mut endpu: *mut libc::c_char = std::ptr::null_mut();

        let ru = unsafe { libc::strtol(s.as_ptr(), &mut endpu, 10) };
        let (rs, endps): (c_longlong, usize) = strtol_safe(s).unwrap();

        assert_eq!(ru, 12345);
        unsafe {
            assert_eq!(endpu.offset_from(s.as_ptr()), 5);
        }
        assert_eq!(rs, 12345);
        assert_eq!(endps, 5)
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_strtol_safer_vs_strtol() {
        // Test with leadng whitespace
        let s: &[c_char] = bytemuck::cast_slice(b"   12345 XC\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let (ru, endp) = strtol_safe(s).unwrap();
        assert_eq!(ru, 12345);
        assert_eq!(endp, 8); // 3 spaces + 5 digits
        assert_eq!(rs, ru);
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with leading zeros
        let s: &[c_char] = bytemuck::cast_slice(b"00012345 XC\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let (ru, endp) = strtol_safe(s).unwrap();
        assert_eq!(ru, 12345);
        assert_eq!(endp, 8); // 3 zeros + 5 digits
        assert_eq!(rs, ru);
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with negative number
        let s: &[c_char] = bytemuck::cast_slice(b"-12345 XC\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let (ru, endp) = strtol_safe(s).unwrap();
        assert_eq!(ru, -12345);
        assert_eq!(endp, 6); // 1 minus + 5 digits
        assert_eq!(rs, ru);
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with invalid characters
        let s: &[c_char] = bytemuck::cast_slice(b"12345a XC\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let (ru, endp) = strtol_safe(s).unwrap();
        assert_eq!(ru, 12345);
        assert_eq!(endp, 5); // 5 digits
        assert_eq!(rs, ru);
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with empty string
        let s: &[c_char] = bytemuck::cast_slice(b"\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let r = strtol_safe::<c_longlong>(s);
        let endp = 0;
        assert!(r.is_err());
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with only whitespace
        let s: &[c_char] = bytemuck::cast_slice(b"   \0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let r = strtol_safe::<c_longlong>(s);
        let endp = 0;
        assert!(r.is_err());
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with only invalid characters
        let s: &[c_char] = bytemuck::cast_slice(b"abcde\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let r = strtol_safe::<c_longlong>(s);
        let endp = 0;
        assert!(r.is_err());
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with leading zeros and invalid characters
        let s: &[c_char] = bytemuck::cast_slice(b"00012345a XC\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let (ru, endp) = strtol_safe(s).unwrap();
        assert_eq!(ru, 12345);
        assert_eq!(endp, 8); // 3 zeros + 5 digits
        assert_eq!(rs, ru);
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with negative number and invalid characters
        let s: &[c_char] = bytemuck::cast_slice(b"-12345a XC\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let (ru, endp) = strtol_safe(s).unwrap();
        assert_eq!(ru, -12345);
        assert_eq!(endp, 6); // 1 minus + 5 digits
        assert_eq!(rs, ru);
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);

        // Test with leading whitespace and invalid characters
        let s: &[c_char] = bytemuck::cast_slice(b"   12345a XC\0");
        let mut endps: *mut libc::c_char = std::ptr::null_mut();
        let rs = unsafe { libc::strtol(s.as_ptr(), &mut endps, 10) };
        let (ru, endp) = strtol_safe(s).unwrap();
        assert_eq!(ru, 12345);
        assert_eq!(endp, 8); // 3 spaces + 5 digits
        assert_eq!(rs, ru);
        assert_eq!(endp, unsafe { endps.offset_from(s.as_ptr()) } as usize);
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_strlen() {
        let s: &[libc::c_char] = bytemuck::cast_slice(b"12345 XC\0");

        let ru = unsafe { libc::strlen(s.as_ptr()) };
        let rs = strlen_safe(s);
        assert_eq!(ru, rs);
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_toupper() {
        let ru = unsafe { libc::toupper(b'a' as libc::c_int) as libc::c_char };
        let rs = toupper(b'a' as libc::c_char);

        assert_eq!(ru, rs);
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_islower() {
        let ru = unsafe { libc::islower(b'a' as libc::c_int) > 0 };
        let rs = islower(b'a' as libc::c_char);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::islower(b'A' as libc::c_int) > 0 };
        let rs = islower(b'A' as libc::c_char);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::islower(b'.' as libc::c_int) > 0 };
        let rs = islower(b'.' as libc::c_char);
        assert_eq!(ru, rs);
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_strncmp_safe() {
        let s1: &[libc::c_char] = bytemuck::cast_slice(b"hello\0");
        let s2: &[libc::c_char] = bytemuck::cast_slice(b"hello\0");
        let s3: &[libc::c_char] = bytemuck::cast_slice(b"world\0");
        let s4: &[libc::c_char] = bytemuck::cast_slice(b"worl\0");

        let rs = strncmp_safe(s1, s2, 5);
        assert_eq!(rs, 0);

        let rs = strncmp_safe(s1, s3, 5);
        assert!(rs < 0);

        let rs = strncmp_safe(s3, s1, 5);
        assert!(rs > 0);

        let rs = strncmp_safe(s1, s2, 3);
        assert_eq!(rs, 0);

        let rs = strncmp_safe(s1, s3, 3);
        assert!(rs < 0);

        let rs = strncmp_safe(s3, s4, 6);
        assert!(rs > 0);
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_compare_strncmp_strncmp_safe() {
        let s1: &[libc::c_char] = bytemuck::cast_slice(b"hello\0");
        let s2: &[libc::c_char] = bytemuck::cast_slice(b"hello\0");
        let s3: &[libc::c_char] = bytemuck::cast_slice(b"world\0");
        let s4: &[libc::c_char] = bytemuck::cast_slice(b"worl\0");

        let ru = unsafe { libc::strncmp(s1.as_ptr(), s2.as_ptr(), 5) };
        let rs = strncmp_safe(s1, s2, 5);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strncmp(s1.as_ptr(), s3.as_ptr(), 5) };
        let rs = strncmp_safe(s1, s3, 5);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strncmp(s3.as_ptr(), s1.as_ptr(), 5) };
        let rs = strncmp_safe(s3, s1, 5);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strncmp(s1.as_ptr(), s2.as_ptr(), 3) };
        let rs = strncmp_safe(s1, s2, 3);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strncmp(s1.as_ptr(), s3.as_ptr(), 3) };
        let rs = strncmp_safe(s1, s3, 3);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strncmp(s3.as_ptr(), s4.as_ptr(), 6) };
        let rs = strncmp_safe(s3, s4, 6);
        assert_eq!(ru, rs);
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_compare_strcmp_strcmp_safe() {
        let s1: &[libc::c_char] = bytemuck::cast_slice(b"hello\0");
        let s2: &[libc::c_char] = bytemuck::cast_slice(b"hello\0");
        let s3: &[libc::c_char] = bytemuck::cast_slice(b"world\0");
        let s4: &[libc::c_char] = bytemuck::cast_slice(b"worl\0");

        let ru = unsafe { libc::strcmp(s1.as_ptr(), s2.as_ptr()) };
        let rs = strcmp_safe(s1, s2);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strcmp(s1.as_ptr(), s3.as_ptr()) };
        let rs = strcmp_safe(s1, s3);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strcmp(s3.as_ptr(), s1.as_ptr()) };
        let rs = strcmp_safe(s3, s1);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strcmp(s1.as_ptr(), s2.as_ptr()) };
        let rs = strcmp_safe(s1, s2);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strcmp(s1.as_ptr(), s3.as_ptr()) };
        let rs = strcmp_safe(s1, s3);
        assert_eq!(ru, rs);

        let ru = unsafe { libc::strcmp(s3.as_ptr(), s4.as_ptr()) };
        let rs = strcmp_safe(s3, s4);
        assert_eq!(ru, rs);
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_compare_strcmp_strcmp_safe_shorter() {
        let s1: &[libc::c_char] = bytemuck::cast_slice(b"..\0");
        let s2: &[libc::c_char] = bytemuck::cast_slice(b".\0");

        // Without nulls
        let s3: &[libc::c_char] = bytemuck::cast_slice(b"..");
        let s4: &[libc::c_char] = bytemuck::cast_slice(b".");

        let ru = unsafe { libc::strcmp(s1.as_ptr(), s2.as_ptr()) };
        let rs = strcmp_safe(s3, s4);
        assert_eq!(ru, rs);
    }

    // Write test cases for strto_float_impl
    #[test]
    fn test_strto_float_impl() {
        let mut endptr: usize = 0;

        // Test with valid float
        let s: &[libc::c_char] = bytemuck::cast_slice(b"123.456\0");
        let result = strto_float_impl(s, &mut endptr);
        assert_eq!(result, 123.456);
        assert_eq!(endptr, 7); // "123.456" + null terminator

        // Test with negative float
        let s: &[libc::c_char] = bytemuck::cast_slice(b"-123.456\0");
        let result = strto_float_impl(s, &mut endptr);
        assert_eq!(result, -123.456);
        assert_eq!(endptr, 8); // "-123.456" + null terminator

        // Test with scientific notation
        let s: &[libc::c_char] = bytemuck::cast_slice(b"1.23e4\0");
        let result = strto_float_impl(s, &mut endptr);
        assert_eq!(result, 12300.0);
        assert_eq!(endptr, 6); // "1.23e4" + null terminator

        // Test with invalid float
        let s: &[libc::c_char] = bytemuck::cast_slice(b"abc\0");
        let result = strto_float_impl(s, &mut endptr);
        assert_eq!(result, 0.0);
        assert_eq!(endptr, 0); // No valid float parsed

        // Test with empty string
        let s: &[libc::c_char] = bytemuck::cast_slice(b"\0");
        let result = strto_float_impl(s, &mut endptr);
        assert_eq!(result, 0.0);
        assert_eq!(endptr, 0); // No valid float parsed

        // Test with Nan
        let s: &[libc::c_char] = bytemuck::cast_slice(b"nan\0");
        let result = strto_float_impl(s, &mut endptr);
        assert!(result.is_nan());
        assert_eq!(endptr, 3); // NaN parsed
    }
}
