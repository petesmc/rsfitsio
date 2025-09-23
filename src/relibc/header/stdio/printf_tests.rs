#[cfg(test)]
mod tests {

    use crate::{
        c_types::{c_char, c_int},
        relibc::header::stdio::{
            snprintf_cint, snprintf_f64, snprintf_f64_decim, sprintf_f64, sprintf_string_width,
        },
    };
    use bytemuck::cast_slice;
    use libc;
    use std::ffi::CStr;

    // Helper to convert string literals to c_char arrays
    macro_rules! c_str {
        ($s:literal) => {{
            const S: &str = concat!($s, "\0");
            std::mem::transmute::<*const u8, *const c_char>(S.as_ptr())
        }};
    }

    // Test the existing snprintf_cint function against libc
    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_snprintf_cint() {
        unsafe {
            let value: c_int = 42;
            let mut rust_buffer: [c_char; 100] = [0; 100];

            let rust_result = snprintf_cint(&mut rust_buffer, 100, cast_slice(b"%d\0"), value);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            assert_eq!(rust_result, 2, "Expected return value of 2 for '42'");

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%d");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(rust_str, libc_str, "Integer formatting mismatch for %d");
                assert_eq!(rust_result, libc_result, "Return value mismatch for %d");
            }

            assert_eq!(rust_str, "42");
        }
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_snprintf_f64() {
        unsafe {
            let value: f64 = 3.24159;
            let mut rust_buffer: [c_char; 100] = [0; 100];

            let rust_result = snprintf_f64(&mut rust_buffer, 100, cast_slice(b"%f\0"), value);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            assert_eq!(
                rust_str, "3.241590",
                "Expected float formatting for 3.24159"
            );

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%f");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(rust_str, libc_str, "Float formatting mismatch for %f");
                assert_eq!(rust_result, libc_result, "Return value mismatch for %f");
            }
        }
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_sprintf_f64() {
        unsafe {
            let value: f64 = 2.81828;
            let mut rust_buffer: [c_char; 100] = [0; 100];

            let rust_result = sprintf_f64(&mut rust_buffer, cast_slice(b"%f\0"), value);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            assert_eq!(
                rust_str, "2.818280",
                "Expected float formatting for 2.81828"
            );

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%f");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::sprintf(libc_buffer.as_mut_ptr(), format, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(
                    rust_str, libc_str,
                    "Float formatting mismatch for sprintf_f64"
                );
                assert_eq!(
                    rust_result, libc_result,
                    "Return value mismatch for sprintf_f64"
                );
            }
        }
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_snprintf_f64_decim() {
        unsafe {
            let value: f64 = 3.24159265;
            let decim: c_int = 2;
            let mut rust_buffer: [c_char; 100] = [0; 100];

            let rust_result =
                snprintf_f64_decim(&mut rust_buffer, 100, cast_slice(b"%.*f\0"), decim, value);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            assert_eq!(rust_result, 4, "Expected return value of 4 for '3.24'");

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%.*f");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result =
                    libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, decim, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(rust_str, libc_str, "Precision float formatting mismatch");
                assert_eq!(
                    rust_result, libc_result,
                    "Return value mismatch for precision float"
                );
            }

            assert_eq!(rust_str, "3.24");
        }
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_sprintf_string_width() {
        unsafe {
            let width: c_int = 10;
            let mut rust_buffer: [c_char; 100] = [0; 100];

            let rust_result = sprintf_string_width(
                &mut rust_buffer,
                cast_slice(b"%*s\0"),
                width,
                cast_slice(b"test\0"),
            );
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            assert_eq!(
                rust_result, 10,
                "Expected return value of 10 for '      test'"
            );

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%*s");
                let value = c_str!("test");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::sprintf(libc_buffer.as_mut_ptr(), format, width, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(rust_str, libc_str, "String width formatting mismatch");
                assert_eq!(
                    rust_result, libc_result,
                    "Return value mismatch for string width"
                );
            }

            assert_eq!(rust_str, "      test");
        }
    }

    // Test edge cases
    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_sprintf_edge_cases() {
        unsafe {
            // Test zero value
            let value: c_int = 0;
            let mut rust_buffer: [c_char; 100] = [0; 100];

            let rust_result = snprintf_cint(&mut rust_buffer, 100, cast_slice(b"%d\0"), value);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            assert_eq!(rust_result, 1, "Expected return value of 1 for '0'");

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%d");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(rust_str, libc_str, "Zero value formatting mismatch");
                assert_eq!(rust_result, libc_result, "Return value mismatch for zero");
            }

            assert_eq!(rust_str, "0");

            // Test negative value
            let value: c_int = -42;
            rust_buffer = [0; 100];

            let rust_result = snprintf_cint(&mut rust_buffer, 100, cast_slice(b"%d\0"), value);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            assert_eq!(rust_result, 3, "Expected return value of 3 for '-42'");

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%d");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(rust_str, libc_str, "Negative value formatting mismatch");
                assert_eq!(
                    rust_result, libc_result,
                    "Return value mismatch for negative"
                );
            }

            assert_eq!(rust_str, "-42");
        }
    }

    // Comprehensive integer tests
    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_snprintf_integers_comprehensive() {
        unsafe {
            let test_values = [
                (1, "1"),
                (42, "42"),
                (123, "123"),
                (-1, "-1"),
                (-42, "-42"),
                (-123, "-123"),
                (2147483647, "2147483647"),   // INT_MAX
                (-2147483648, "-2147483648"), // INT_MIN
                (1000, "1000"),
                (-1000, "-1000"),
                (12345, "12345"),
                (-12345, "-12345"),
                (999999, "999999"),
                (-999999, "-999999"),
            ];

            for (value, expected) in test_values.iter() {
                let mut rust_buffer: [c_char; 100] = [0; 100];

                let rust_result = snprintf_cint(&mut rust_buffer, 100, cast_slice(b"%d\0"), *value);
                let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(
                    rust_result as usize,
                    expected.len(),
                    "Expected return value to match string length for {value}"
                );

                #[cfg(not(target_family = "windows"))]
                {
                    let format = c_str!("%d");
                    let mut libc_buffer: [c_char; 100] = [0; 100];
                    let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, *value);
                    let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                    assert_eq!(rust_str, libc_str, "Integer {value} formatting mismatch");
                    assert_eq!(
                        rust_result, libc_result,
                        "Return value mismatch for {value}"
                    );
                }

                assert_eq!(rust_str, *expected, "Expected output mismatch for {value}");
            }
        }
    }

    // Comprehensive float tests
    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_snprintf_floats_comprehensive() {
        unsafe {
            let test_values = [
                (0.0, "0.000000"),
                (1.0, "1.000000"),
                (-1.0, "-1.000000"),
                (3.24159, "3.241590"),
                (-3.24159, "-3.241590"),
                (2.81828, "2.818280"),
                (123.456, "123.456000"),
                (-123.456, "-123.456000"),
                (0.000001, "0.000001"),
                (-0.000001, "-0.000001"),
                (1000000.0, "1000000.000000"),
                (-1000000.0, "-1000000.000000"),
            ];

            for (value, expected) in test_values.iter() {
                let mut rust_buffer: [c_char; 100] = [0; 100];

                let rust_result = snprintf_f64(&mut rust_buffer, 100, cast_slice(b"%f\0"), *value);
                let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(
                    rust_result as usize,
                    expected.len(),
                    "Expected return value to match string length for {value}"
                );

                #[cfg(not(target_family = "windows"))]
                {
                    let format = c_str!("%f");
                    let mut libc_buffer: [c_char; 100] = [0; 100];
                    let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, *value);
                    let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                    assert_eq!(rust_str, libc_str, "Float {value} formatting mismatch");
                    assert_eq!(
                        rust_result, libc_result,
                        "Return value mismatch for {value}"
                    );
                }

                assert_eq!(rust_str, *expected, "Expected output mismatch for {value}");
            }
        }
    }

    // Test precision specifiers
    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_snprintf_precision_comprehensive() {
        unsafe {
            let test_cases = [
                (3.24159, 0, "3"),
                (3.24159, 1, "3.2"),
                (3.24159, 2, "3.24"),
                (3.24159, 3, "3.242"),
                (3.24159, 4, "3.2416"),
                (3.24159, 5, "3.24159"),
                (123.456, 0, "123"),
                (123.456, 1, "123.5"),
                (123.456, 2, "123.46"),
                (0.0, 0, "0"),
                (0.0, 3, "0.000"),
                (-2.5, 1, "-2.5"),
            ];

            for (value, precision, expected) in test_cases.iter() {
                let mut rust_buffer: [c_char; 100] = [0; 100];

                let rust_result = snprintf_f64_decim(
                    &mut rust_buffer,
                    100,
                    cast_slice(b"%.*f\0"),
                    *precision,
                    *value,
                );
                let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

                #[cfg(not(target_family = "windows"))]
                {
                    let format = c_str!("%.*f");
                    let mut libc_buffer: [c_char; 100] = [0; 100];
                    let libc_result =
                        libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, *precision, *value);
                    let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                    assert_eq!(
                        rust_str, libc_str,
                        "Precision {precision} for {value} formatting mismatch"
                    );
                    assert_eq!(
                        rust_result, libc_result,
                        "Return value mismatch for precision {precision} value {value}"
                    );
                }

                assert_eq!(
                    rust_result as usize,
                    expected.len(),
                    "Expected return value to match string length for precision {precision} value {value}"
                );

                assert_eq!(
                    rust_str, *expected,
                    "Expected output mismatch for precision {precision} value {value}"
                );
            }
        }
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_many_test_cases() {
        unsafe {
            // Test case: 233 with %d should produce "233"
            let value: c_int = 233;
            let mut rust_buffer: [c_char; 100] = [0; 100];

            let rust_result = snprintf_cint(&mut rust_buffer, 100, cast_slice(b"%d\0"), value);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%d");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, value);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(rust_str, libc_str, "test case mismatch for %d with 233");
                assert_eq!(
                    rust_result, libc_result,
                    "Return value mismatch for %d test"
                );
            }

            assert_eq!(rust_result, 3, "Expected return value of 3 for '233'");
            assert_eq!(rust_str, "233");

            // Test case: 233000L with %ld should produce "233000"
            // Note: Testing with the available integer function
            let value_long: c_int = 233000; // Using c_int for compatibility with snprintf_cint
            rust_buffer = [0; 100];

            let rust_result = snprintf_cint(&mut rust_buffer, 100, cast_slice(b"%d\0"), value_long);
            let rust_str = CStr::from_ptr(rust_buffer.as_ptr()).to_str().unwrap();

            #[cfg(not(target_family = "windows"))]
            {
                let format = c_str!("%d");
                let mut libc_buffer: [c_char; 100] = [0; 100];
                let libc_result = libc::snprintf(libc_buffer.as_mut_ptr(), 100, format, value_long);
                let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();

                assert_eq!(
                    rust_str, libc_str,
                    "Z88dk test case mismatch for large integer"
                );
                assert_eq!(
                    rust_result, libc_result,
                    "Return value mismatch for large integer"
                );
            }

            assert_eq!(rust_result, 6, "Expected return value of 6 for '233000'");
            assert_eq!(rust_str, "233000");
        }
    }
}
