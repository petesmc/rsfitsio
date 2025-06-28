use crate::c_types::{c_char, c_double, c_int, c_long, c_uint, size_t};
use bytemuck::cast_slice;
use printf::{CustomVaList, VaArg};
use std::ffi::{CStr, c_void};

use crate::relibc::platform;

mod lookaheadreader;
mod printf;
mod scanf;

#[cfg(test)]
mod printf_tests;
#[cfg(test)]
mod sscanf_tests;

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

fn sscanf_internal(s: *const c_char, format: *const c_char, valist: CustomVaList) -> c_int {
    unsafe {
        let reader = (s as *const u8).into();
        scanf::scanf_custom(reader, format, valist)
    }
}

pub(crate) fn sscanf_ld(s: &[c_char], format: &[c_char], val: *mut c_long) -> c_int {
    let mut valist = CustomVaList::new();
    valist.push(VaArg::pointer(val as *const c_void));

    sscanf_internal(s.as_ptr(), format.as_ptr(), valist)
}

pub(crate) fn sscanf_d(s: &[c_char], format: &[c_char], val: *mut c_int) -> c_int {
    let mut valist = CustomVaList::new();
    valist.push(VaArg::pointer(val as *const c_void));

    sscanf_internal(s.as_ptr(), format.as_ptr(), valist)
}

pub(crate) fn sscanf_u(s: &[c_char], format: &[c_char], val: *mut c_uint) -> c_int {
    let mut valist = CustomVaList::new();
    valist.push(VaArg::pointer(val as *const c_void));

    sscanf_internal(s.as_ptr(), format.as_ptr(), valist)
}

pub(crate) fn sscanf_lf(s: &[c_char], format: &[c_char], val: *mut c_double) -> c_int {
    let mut valist = CustomVaList::new();
    valist.push(VaArg::pointer(val as *const c_void));

    sscanf_internal(s.as_ptr(), format.as_ptr(), valist)
}

#[cfg(test)]
mod tests {
    use super::*;
    use libc::{self, c_float};

    #[test]
    fn test_printf() {
        let format = c"%3d";
        let mut buffer: [c_char; 100] = [0; 100];
        let result = snprintf_cint(&mut buffer, 100, cast_slice(format.to_bytes_with_nul()), 42);

        assert_eq!(
            &buffer[..(result + 1) as usize],
            cast_slice(c" 42".to_bytes_with_nul())
        );
    }

    #[test]
    fn test_sscanf_helper_functions() {
        // Test sscanf_d for %d format
        let input = c"42";
        let format = c"%d";
        let mut value: c_int = 0;
        let result = sscanf_d(
            cast_slice(input.to_bytes_with_nul()),
            cast_slice(format.to_bytes_with_nul()),
            &mut value,
        );
        assert_eq!(result, 1);
        assert_eq!(value, 42);

        // Test sscanf_ld for %ld format
        let input = c"123456789";
        let format = c"%ld";
        let mut value: c_long = 0;
        let result = sscanf_ld(
            cast_slice(input.to_bytes_with_nul()),
            cast_slice(format.to_bytes_with_nul()),
            &mut value,
        );
        assert_eq!(result, 1);
        assert_eq!(value, 123456789);

        // Test sscanf_u for %u format
        let input = c"4294967295";
        let format = c"%u";
        let mut value: c_uint = 0;
        let result = sscanf_u(
            cast_slice(input.to_bytes_with_nul()),
            cast_slice(format.to_bytes_with_nul()),
            &mut value,
        );
        assert_eq!(result, 1);
        assert_eq!(value, 4294967295);

        // Test sscanf_lf for %lf format
        let input = c"3.24159";
        let format = c"%lf";
        let mut value: c_double = 0.0;
        let result = sscanf_lf(
            cast_slice(input.to_bytes_with_nul()),
            cast_slice(format.to_bytes_with_nul()),
            &mut value,
        );
        assert_eq!(result, 1);
        assert!((value - 3.24159).abs() < 0.00001);
    }

    #[test]
    fn test_sscanf_basic_debug() {
        // First test: Check if our original sscanf works
        #[cfg(not(target_family = "windows"))]
        unsafe {
            let input = c"42";
            let format = c"%d";
            let mut value: c_int = 0;
            let result = libc::sscanf(input.as_ptr(), format.as_ptr(), &mut value as *mut c_int);
            println!("Original sscanf result: {result}, value: {value}");
        }

        // Second test: Check sscanf_internal directly

        let input = c"42";
        let format = c"%d";
        let mut value: c_int = 0;
        let mut valist = CustomVaList::new();
        valist.push(VaArg::pointer(&mut value as *mut c_int as *const c_void));
        let result = sscanf_internal(input.as_ptr(), format.as_ptr(), valist);
        println!("sscanf_internal result: {result}, value: {value}");

        // Verify expected values
        assert_eq!(result, 1);
        assert_eq!(value, 42);
    }

    #[test]
    fn test_sscanf_internal_vs_libc() {
        // Test 1: Basic integer parsing with %d
        unsafe {
            let input = c"42";
            let format = c"%d";

            // Test with our implementation
            let mut rust_value: c_int = 0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_int as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert_eq!(rust_value, 42, "Expected value to be 42");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_int = 0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_int,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for basic integer: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_value, libc_value,
                    "Value mismatch for basic integer: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 2: Long integer parsing with %ld
        unsafe {
            let input = c"1234567890";
            let format = c"%ld";

            // Test with our implementation
            let mut rust_value: c_long = 0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_long as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert_eq!(rust_value, 1234567890, "Expected value to be 1234567890");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_long = 0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_long,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for long integer: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_value, libc_value,
                    "Value mismatch for long integer: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 3: Unsigned integer parsing with %u
        unsafe {
            let input = c"4294967295";
            let format = c"%u";

            // Test with our implementation
            let mut rust_value: c_uint = 0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_uint as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert_eq!(rust_value, 4294967295, "Expected value to be 4294967295");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_uint = 0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_uint,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for unsigned integer: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_value, libc_value,
                    "Value mismatch for unsigned integer: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 4: Double parsing with %lf
        unsafe {
            let input = c"3.14159265359";
            let format = c"%lf";

            // Test with our implementation
            let mut rust_value: c_double = 0.0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_double as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert!(
                (rust_value - 3.14159265359).abs() < 1e-10,
                "Expected value to be close to 3.14159265359, got {rust_value}"
            );

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_double = 0.0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_double,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for double: rust={rust_result}, libc={libc_result}"
                );
                assert!(
                    (rust_value - libc_value).abs() < 1e-15,
                    "Value mismatch for double: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 5: Float parsing with %f
        unsafe {
            let input = c"2.71828";
            let format = c"%f";

            // Test with our implementation
            let mut rust_value: c_float = 0.0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_float as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert!(
                (rust_value - 2.71828).abs() < 1e-5,
                "Expected value to be close to 2.71828, got {rust_value}"
            );

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_float = 0.0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_float,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for float: rust={rust_result}, libc={libc_result}"
                );
                assert!(
                    (rust_value - libc_value).abs() < 1e-6,
                    "Value mismatch for float: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 6: Negative integer
        unsafe {
            let input = c"-12345";
            let format = c"%d";

            // Test with our implementation
            let mut rust_value: c_int = 0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_int as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert_eq!(rust_value, -12345, "Expected value to be -12345");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_int = 0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_int,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for negative integer: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_value, libc_value,
                    "Value mismatch for negative integer: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 7: Scientific notation
        unsafe {
            let input = c"1.23e-4";
            let format = c"%lf";

            // Test with our implementation
            let mut rust_value: c_double = 0.0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_double as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert!(
                (rust_value - 1.23e-4).abs() < 1e-15,
                "Expected value to be close to 1.23e-4, got {rust_value}"
            );

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_double = 0.0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_double,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for scientific notation: rust={rust_result}, libc={libc_result}"
                );
                assert!(
                    (rust_value - libc_value).abs() < 1e-15,
                    "Value mismatch for scientific notation: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 8: Multiple values
        unsafe {
            let input = c"10 20 30";
            let format = c"%d %d %d";

            // Test with our implementation
            let mut rust_a: c_int = 0;
            let mut rust_b: c_int = 0;
            let mut rust_c: c_int = 0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(&mut rust_a as *mut c_int as *const c_void));
            rust_valist.push(VaArg::pointer(&mut rust_b as *mut c_int as *const c_void));
            rust_valist.push(VaArg::pointer(&mut rust_c as *mut c_int as *const c_void));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 3, "Expected 3 successful conversions");
            assert_eq!(rust_a, 10, "First value should be 10");
            assert_eq!(rust_b, 20, "Second value should be 20");
            assert_eq!(rust_c, 30, "Third value should be 30");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_a: c_int = 0;
                let mut libc_b: c_int = 0;
                let mut libc_c: c_int = 0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_a as *mut c_int,
                    &mut libc_b as *mut c_int,
                    &mut libc_c as *mut c_int,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for multiple values: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_a, libc_a,
                    "First value mismatch: rust={rust_a}, libc={libc_a}"
                );
                assert_eq!(
                    rust_b, libc_b,
                    "Second value mismatch: rust={rust_b}, libc={libc_b}"
                );
                assert_eq!(
                    rust_c, libc_c,
                    "Third value mismatch: rust={rust_c}, libc={libc_c}"
                );
            }
        }

        // Test 9: Empty input
        unsafe {
            let input = c"";
            let format = c"%d";

            // Test with our implementation
            let mut rust_value: c_int = 999;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_int as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, -1, "Expected -1 (EOF) for empty input");
            assert_eq!(rust_value, 999, "Value should remain unchanged (999)");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_int = 999;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_int,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for empty input: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_value, libc_value,
                    "Value mismatch for empty input: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 10: Invalid format
        unsafe {
            let input = c"abc";
            let format = c"%d";

            // Test with our implementation
            let mut rust_value: c_int = 999;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_int as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(
                rust_result, 0,
                "Expected 0 conversions for non-numeric input"
            );
            assert_eq!(rust_value, 999, "Value should remain unchanged (999)");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_int = 999;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_int,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for invalid format: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_value, libc_value,
                    "Value mismatch for invalid format: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 11: Whitespace handling
        unsafe {
            let input = c"   42   ";
            let format = c"%d";

            // Test with our implementation
            let mut rust_value: c_int = 0;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(
                &mut rust_value as *mut c_int as *const c_void,
            ));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(rust_result, 1, "Expected 1 successful conversion");
            assert_eq!(rust_value, 42, "Expected value to be 42");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_value: c_int = 0;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_value as *mut c_int,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for whitespace: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_value, libc_value,
                    "Value mismatch for whitespace: rust={rust_value}, libc={libc_value}"
                );
            }
        }

        // Test 12: Partial match
        unsafe {
            let input = c"123 abc";
            let format = c"%d %d";

            // Test with our implementation
            let mut rust_a: c_int = 0;
            let mut rust_b: c_int = 999;
            let mut rust_valist = CustomVaList::new();
            rust_valist.push(VaArg::pointer(&mut rust_a as *mut c_int as *const c_void));
            rust_valist.push(VaArg::pointer(&mut rust_b as *mut c_int as *const c_void));
            let rust_result = sscanf_internal(input.as_ptr(), format.as_ptr(), rust_valist);

            // Expected values
            assert_eq!(
                rust_result, 1,
                "Expected 1 successful conversion (partial match)"
            );
            assert_eq!(rust_a, 123, "First value should be 123");
            assert_eq!(rust_b, 999, "Second value should remain unchanged (999)");

            // Compare with libc on non-Windows platforms
            #[cfg(not(target_family = "windows"))]
            {
                let mut libc_a: c_int = 0;
                let mut libc_b: c_int = 999;
                let libc_result = libc::sscanf(
                    input.as_ptr(),
                    format.as_ptr(),
                    &mut libc_a as *mut c_int,
                    &mut libc_b as *mut c_int,
                );

                assert_eq!(
                    rust_result, libc_result,
                    "Result mismatch for partial match: rust={rust_result}, libc={libc_result}"
                );
                assert_eq!(
                    rust_a, libc_a,
                    "First value mismatch for partial match: rust={rust_a}, libc={libc_a}"
                );
                assert_eq!(
                    rust_b, libc_b,
                    "Second value mismatch for partial match: rust={rust_b}, libc={libc_b}"
                );
            }
        }
    }
}
