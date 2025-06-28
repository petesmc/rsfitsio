#[cfg(test)]
mod tests {

    use crate::{
        c_types::{c_char, c_double, c_int, c_long, c_uint},
        relibc::header::stdio::{
            printf::{CustomVaList, VaArg},
            sscanf_internal,
        },
    };

    use libc;
    use std::ffi::{CStr, CString, c_void};

    // Helper to create null-terminated string from &str
    fn to_c_string(s: &str) -> Vec<c_char> {
        let mut result: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        result.push(0);
        result
    }

    // Test helper for sscanf_internal with comparison to libc
    unsafe fn test_sscanf_internal(
        input: &str,
        format: &str,
        valist: CustomVaList,
        expected_matches: c_int,
    ) -> c_int {
        let input_cstr = to_c_string(input);
        let format_cstr = to_c_string(format);

        let result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);
        if result != expected_matches {
            println!(
                "FAILED: input='{input}', format='{format}', result={result}, expected={expected_matches}"
            );
        }
        assert_eq!(
            result, expected_matches,
            "sscanf_internal failed for input '{input}' with format '{format}'"
        );
        result
    }

    // Helper to compare our sscanf implementation with libc for c_int
    unsafe fn test_sscanf_int_with_libc(input: &str, format: &str) -> (c_int, c_int) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_int = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_int,
                );

                // Test with our implementation
                let mut our_val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result > 0 {
                    assert_eq!(
                        our_val, libc_val,
                        "Parsed values should match for input='{input}', format='{format}': our={our_val}, libc={libc_val}"
                    );
                }

                (our_result, our_val)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for c_uint
    unsafe fn test_sscanf_uint_with_libc(input: &str, format: &str) -> (c_int, c_uint) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_uint = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_uint,
                );

                // Test with our implementation
                let mut our_val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_uint as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result > 0 {
                    assert_eq!(
                        our_val, libc_val,
                        "Parsed values should match for input='{input}', format='{format}': our={our_val}, libc={libc_val}"
                    );
                }

                (our_result, our_val)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_uint as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for c_long
    unsafe fn test_sscanf_long_with_libc(input: &str, format: &str) -> (c_int, c_long) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_long = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_long,
                );

                // Test with our implementation
                let mut our_val: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_long as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result > 0 {
                    assert_eq!(
                        our_val, libc_val,
                        "Parsed values should match for input='{input}', format='{format}': our={our_val}, libc={libc_val}"
                    );
                }

                (our_result, our_val)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_long as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for c_double
    unsafe fn test_sscanf_double_with_libc(input: &str, format: &str) -> (c_int, c_double) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_double = 0.0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_double,
                );

                // Test with our implementation
                let mut our_val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val as *mut c_double as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result > 0 {
                    assert!(
                        (our_val - libc_val).abs() < 1e-10,
                        "Parsed values should match for input='{input}', format='{format}': our={our_val}, libc={libc_val}"
                    );
                }

                (our_result, our_val)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val as *mut c_double as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for c_char
    unsafe fn test_sscanf_char_with_libc(input: &str, format: &str) -> (c_int, c_char) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_char = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_char,
                );

                // Test with our implementation
                let mut our_val: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_char as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result > 0 {
                    assert_eq!(
                        our_val, libc_val,
                        "Parsed values should match for input='{input}', format='{format}': our={our_val}, libc={libc_val}"
                    );
                }

                (our_result, our_val)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_char as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for string (%s)
    unsafe fn test_sscanf_string_with_libc(
        input: &str,
        format: &str,
        buffer_size: usize,
    ) -> (c_int, String) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_buffer: Vec<c_char> = vec![0; buffer_size];
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    libc_buffer.as_mut_ptr(),
                );

                // Test with our implementation
                let mut our_buffer: Vec<c_char> = vec![0; buffer_size];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(our_buffer.as_mut_ptr() as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );

                if libc_result > 0 {
                    // Convert buffers to strings for comparison
                    let libc_str = CStr::from_ptr(libc_buffer.as_ptr()).to_str().unwrap();
                    let our_str = CStr::from_ptr(our_buffer.as_ptr()).to_str().unwrap();
                    assert_eq!(
                        our_str, libc_str,
                        "Parsed strings should match for input='{input}', format='{format}': our='{our_str}', libc='{libc_str}'"
                    );
                    (our_result, our_str.to_string())
                } else {
                    (our_result, String::new())
                }
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_buffer: Vec<c_char> = vec![0; buffer_size];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(our_buffer.as_mut_ptr() as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                if our_result > 0 {
                    let our_str = CStr::from_ptr(our_buffer.as_ptr()).to_str().unwrap();
                    (our_result, our_str.to_string())
                } else {
                    (our_result, String::new())
                }
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for literal tests (no arguments)
    unsafe fn test_sscanf_literal_with_libc(input: &str, format: &str) -> c_int {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let libc_result = libc::sscanf(input_cstring.as_ptr(), format_cstring.as_ptr());

                // Test with our implementation
                let valist = CustomVaList::new();
                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );

                our_result
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let valist = CustomVaList::new();
                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                our_result
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for two integers
    unsafe fn test_sscanf_two_ints_with_libc(input: &str, format: &str) -> (c_int, c_int, c_int) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val1: c_int = 0;
                let mut libc_val2: c_int = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val1 as *mut c_int,
                    &mut libc_val2 as *mut c_int,
                );

                // Test with our implementation
                let mut our_val1: c_int = 0;
                let mut our_val2: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut our_val2 as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result >= 1 {
                    assert_eq!(
                        our_val1, libc_val1,
                        "First parsed value should match for input='{input}', format='{format}': our={our_val1}, libc={libc_val1}"
                    );
                }
                if libc_result >= 2 {
                    assert_eq!(
                        our_val2, libc_val2,
                        "Second parsed value should match for input='{input}', format='{format}': our={our_val2}, libc={libc_val2}"
                    );
                }

                (our_result, our_val1, our_val2)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val1: c_int = 0;
                let mut our_val2: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut our_val2 as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val1, our_val2)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for two unsigned integers
    unsafe fn test_sscanf_two_uints_with_libc(
        input: &str,
        format: &str,
    ) -> (c_int, c_uint, c_uint) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val1: c_uint = 0;
                let mut libc_val2: c_uint = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val1 as *mut c_uint,
                    &mut libc_val2 as *mut c_uint,
                );

                // Test with our implementation
                let mut our_val1: c_uint = 0;
                let mut our_val2: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_uint as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_uint as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result >= 1 {
                    assert_eq!(
                        our_val1, libc_val1,
                        "First parsed value should match for input='{input}', format='{format}': our={our_val1}, libc={libc_val1}"
                    );
                }
                if libc_result >= 2 {
                    assert_eq!(
                        our_val2, libc_val2,
                        "Second parsed value should match for input='{input}', format='{format}': our={our_val2}, libc={libc_val2}"
                    );
                }

                (our_result, our_val1, our_val2)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val1: c_uint = 0;
                let mut our_val2: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_uint as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_uint as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val1, our_val2)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for two longs
    unsafe fn test_sscanf_two_longs_with_libc(
        input: &str,
        format: &str,
    ) -> (c_int, c_long, c_long) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val1: c_long = 0;
                let mut libc_val2: c_long = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val1 as *mut c_long,
                    &mut libc_val2 as *mut c_long,
                );

                // Test with our implementation
                let mut our_val1: c_long = 0;
                let mut our_val2: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_long as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_long as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result >= 1 {
                    assert_eq!(
                        our_val1, libc_val1,
                        "First parsed value should match for input='{input}', format='{format}': our={our_val1}, libc={libc_val1}"
                    );
                }
                if libc_result >= 2 {
                    assert_eq!(
                        our_val2, libc_val2,
                        "Second parsed value should match for input='{input}', format='{format}': our={our_val2}, libc={libc_val2}"
                    );
                }

                (our_result, our_val1, our_val2)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val1: c_long = 0;
                let mut our_val2: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_long as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_long as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val1, our_val2)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for hex unsigned integer
    unsafe fn test_sscanf_hex_with_libc(input: &str, format: &str) -> (c_int, c_uint) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_uint = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_uint,
                );

                // Test with our implementation
                let mut our_val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_uint as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result > 0 {
                    assert_eq!(
                        our_val, libc_val,
                        "Parsed values should match for input='{input}', format='{format}': our={our_val}, libc={libc_val}"
                    );
                }

                (our_result, our_val)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_uint as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for two hex unsigned integers
    unsafe fn test_sscanf_two_hex_with_libc(input: &str, format: &str) -> (c_int, c_uint, c_uint) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val1: c_uint = 0;
                let mut libc_val2: c_uint = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val1 as *mut c_uint,
                    &mut libc_val2 as *mut c_uint,
                );

                // Test with our implementation
                let mut our_val1: c_uint = 0;
                let mut our_val2: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_uint as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_uint as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result >= 1 {
                    assert_eq!(
                        our_val1, libc_val1,
                        "First parsed value should match for input='{input}', format='{format}': our={our_val1}, libc={libc_val1}"
                    );
                }
                if libc_result >= 2 {
                    assert_eq!(
                        our_val2, libc_val2,
                        "Second parsed value should match for input='{input}', format='{format}': our={our_val2}, libc={libc_val2}"
                    );
                }

                (our_result, our_val1, our_val2)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val1: c_uint = 0;
                let mut our_val2: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_uint as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_uint as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val1, our_val2)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for octal unsigned integer
    unsafe fn test_sscanf_octal_with_libc(input: &str, format: &str) -> (c_int, c_uint) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_uint = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_uint,
                );

                // Test with our implementation
                let mut our_val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_uint as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result > 0 {
                    assert_eq!(
                        our_val, libc_val,
                        "Parsed values should match for input='{input}', format='{format}': our={our_val}, libc={libc_val}"
                    );
                }

                (our_result, our_val)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_uint as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val)
            }
        }
    }

    // Helper to compare our sscanf implementation with libc for two octal unsigned integers
    unsafe fn test_sscanf_two_octal_with_libc(
        input: &str,
        format: &str,
    ) -> (c_int, c_uint, c_uint) {
        unsafe {
            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val1: c_uint = 0;
                let mut libc_val2: c_uint = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val1 as *mut c_uint,
                    &mut libc_val2 as *mut c_uint,
                );

                // Test with our implementation
                let mut our_val1: c_uint = 0;
                let mut our_val2: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_uint as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_uint as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(
                    our_result, libc_result,
                    "Return values should match for input='{input}', format='{format}': our={our_result}, libc={libc_result}"
                );
                if libc_result >= 1 {
                    assert_eq!(
                        our_val1, libc_val1,
                        "First parsed value should match for input='{input}', format='{format}': our={our_val1}, libc={libc_val1}"
                    );
                }
                if libc_result >= 2 {
                    assert_eq!(
                        our_val2, libc_val2,
                        "Second parsed value should match for input='{input}', format='{format}': our={our_val2}, libc={libc_val2}"
                    );
                }

                (our_result, our_val1, our_val2)
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val1: c_uint = 0;
                let mut our_val2: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(
                    &mut our_val1 as *mut c_uint as *const c_void,
                ));
                valist.push(VaArg::pointer(
                    &mut our_val2 as *mut c_uint as *const c_void,
                ));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                (our_result, our_val1, our_val2)
            }
        }
    }

    // COMPREHENSIVE SCANF TESTS - INSPIRED

    #[test]
    fn test_scanf_d_comprehensive() {
        unsafe {
            // Test 1: Simple decimal integer
            let (result, val) = test_sscanf_int_with_libc("42", "%d");
            assert_eq!(val, 42);
            assert_eq!(result, 1);

            // Test 2: Negative integer
            let (result, val) = test_sscanf_int_with_libc("-123", "%d");
            assert_eq!(val, -123);
            assert_eq!(result, 1);

            // Test 3: Zero
            let (result, val) = test_sscanf_int_with_libc("0", "%d");
            assert_eq!(val, 0);
            assert_eq!(result, 1);

            // Test 4: Large positive number
            let (result, val) = test_sscanf_int_with_libc("2147483647", "%d");
            assert_eq!(val, 2147483647);
            assert_eq!(result, 1);

            // Test 5: Large negative number
            let (result, val) = test_sscanf_int_with_libc("-2147483648", "%d");
            assert_eq!(val, -2147483648);
            assert_eq!(result, 1);

            // Test 6: Leading whitespace
            let (result, val) = test_sscanf_int_with_libc("   123", "%d");
            assert_eq!(val, 123);
            assert_eq!(result, 1);

            // Test 7: Leading zeros
            let (result, val) = test_sscanf_int_with_libc("000123", "%d");
            assert_eq!(val, 123);
            assert_eq!(result, 1);

            // Test 8: Plus sign
            let (result, val) = test_sscanf_int_with_libc("+456", "%d");
            assert_eq!(val, 456);
            assert_eq!(result, 1);

            // Test 9: Multiple integers with spaces
            let (result, val1, val2) = test_sscanf_two_ints_with_libc("12 34", "%d %d");
            assert_eq!(val1, 12);
            assert_eq!(val2, 34);
            assert_eq!(result, 2);

            // Test 10: Multiple integers with tabs and newlines
            let (result, val1, val2) = test_sscanf_two_ints_with_libc("12\t\n34", "%d %d");
            assert_eq!(val1, 12);
            assert_eq!(val2, 34);
            assert_eq!(result, 2);
        }
    }

    #[test]
    fn test_scanf_u_comprehensive() {
        unsafe {
            // Test 1: Simple unsigned integer
            let (result, val) = test_sscanf_uint_with_libc("42", "%u");
            assert_eq!(val, 42);
            assert_eq!(result, 1);

            // Test 2: Zero
            let (result, val) = test_sscanf_uint_with_libc("0", "%u");
            assert_eq!(val, 0);
            assert_eq!(result, 1);

            // Test 3: Large unsigned number
            let (result, val) = test_sscanf_uint_with_libc("4294967295", "%u");
            assert_eq!(val, 4294967295);
            assert_eq!(result, 1);

            // Test 4: Leading whitespace
            let (result, val) = test_sscanf_uint_with_libc("   123", "%u");
            assert_eq!(val, 123);
            assert_eq!(result, 1);

            // Test 5: Leading zeros
            let (result, val) = test_sscanf_uint_with_libc("000456", "%u");
            assert_eq!(val, 456);
            assert_eq!(result, 1);

            // Test 6: Multiple unsigned integers
            let (result, val1, val2) = test_sscanf_two_uints_with_libc("123 456", "%u %u");
            assert_eq!(val1, 123);
            assert_eq!(val2, 456);
            assert_eq!(result, 2);
        }
    }

    #[test]
    fn test_scanf_x_comprehensive() {
        unsafe {
            // Test 1: Simple hex number
            let (result, val) = test_sscanf_hex_with_libc("FF", "%x");
            assert_eq!(val, 255);
            assert_eq!(result, 1);

            // Test 2: Lowercase hex
            let (result, val) = test_sscanf_hex_with_libc("ff", "%x");
            assert_eq!(val, 255);
            assert_eq!(result, 1);

            // Test 3: Mixed case hex
            let (result, val) = test_sscanf_hex_with_libc("aBcD", "%x");
            assert_eq!(val, 0xABCD);
            assert_eq!(result, 1);

            // Test 4: Hex with 0x prefix
            let (result, val) = test_sscanf_hex_with_libc("0xABCD", "%x");
            assert_eq!(val, 0xABCD);
            assert_eq!(result, 1);

            // Test 5: Hex with 0X prefix
            let (result, val) = test_sscanf_hex_with_libc("0X1234", "%x");
            assert_eq!(val, 0x1234);
            assert_eq!(result, 1);

            // Test 6: Zero hex
            let (result, val) = test_sscanf_hex_with_libc("0", "%x");
            assert_eq!(val, 0);
            assert_eq!(result, 1);

            // Test 7: Single digit hex
            let (result, val) = test_sscanf_hex_with_libc("A", "%x");
            assert_eq!(val, 10);
            assert_eq!(result, 1);

            // Test 8: Maximum hex value
            let (result, val) = test_sscanf_hex_with_libc("FFFFFFFF", "%x");
            assert_eq!(val, 0xFFFFFFFF);
            assert_eq!(result, 1);

            // Test 9: Multiple hex values
            let (result, val1, val2) = test_sscanf_two_hex_with_libc("1A 2B", "%x %x");
            assert_eq!(val1, 0x1A);
            assert_eq!(val2, 0x2B);
            assert_eq!(result, 2);

            // Test 10: Uppercase X format
            let (result, val) = test_sscanf_hex_with_libc("BEEF", "%X");
            assert_eq!(val, 0xBEEF);
            assert_eq!(result, 1);
        }
    }

    #[test]
    fn test_scanf_o_comprehensive() {
        unsafe {
            // Test 1: Simple octal number
            let (result, val) = test_sscanf_octal_with_libc("777", "%o");
            assert_eq!(val, 0o777);
            assert_eq!(result, 1);

            // Test 2: Zero octal
            let (result, val) = test_sscanf_octal_with_libc("0", "%o");
            assert_eq!(val, 0);
            assert_eq!(result, 1);

            // Test 3: Single digit octal
            let (result, val) = test_sscanf_octal_with_libc("7", "%o");
            assert_eq!(val, 7);
            assert_eq!(result, 1);

            // Test 4: Leading zeros in octal
            let (result, val) = test_sscanf_octal_with_libc("0123", "%o");
            assert_eq!(val, 0o123);
            assert_eq!(result, 1);

            // Test 5: Large octal number
            let (result, val) = test_sscanf_octal_with_libc("37777777777", "%o");
            assert_eq!(val, 0o37777777777);
            assert_eq!(result, 1);

            // Test 6: Multiple octal values
            let (result, val1, val2) = test_sscanf_two_octal_with_libc("123 456", "%o %o");
            assert_eq!(val1, 0o123);
            assert_eq!(val2, 0o456);
            assert_eq!(result, 2);

            // Test 7: Octal with whitespace
            let (result, val) = test_sscanf_octal_with_libc("   755", "%o");
            assert_eq!(val, 0o755);
            assert_eq!(result, 1);
        }
    }

    #[test]
    fn test_scanf_f_comprehensive() {
        unsafe {
            // Test 1: Simple float
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("3.24", "%lf", valist, 1);
                assert!((val - 3.24).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 2: Integer as float
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("42", "%lf", valist, 1);
                assert!((val - 42.0).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 3: Negative float
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("-2.5", "%lf", valist, 1);
                assert!((val - (-2.5)).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 4: Scientific notation
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("1.23e-4", "%lf", valist, 1);
                assert!((val - 1.23e-4).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 5: Scientific notation with E
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("2.5E+3", "%lf", valist, 1);
                assert!((val - 2500.0).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 6: Zero float
            {
                let mut val: c_double = 999.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("0.0", "%lf", valist, 1);
                assert!((val - 0.0).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 7: Very small float
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("0.000001", "%lf", valist, 1);
                assert!((val - 0.000001).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 8: Large float
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("123456.789", "%lf", valist, 1);
                assert!((val - 123456.789).abs() < 1e-6);
                assert_eq!(result, 1);
            }

            // Test 9: Float with plus sign
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("+3.24", "%lf", valist, 1);
                assert!((val - 3.24).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 10: Multiple floats
            {
                let mut val1: c_double = 0.0;
                let mut val2: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_double as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_double as *const c_void));

                let result = test_sscanf_internal("1.5 2.5", "%lf %lf", valist, 2);
                assert!((val1 - 1.5).abs() < 1e-10);
                assert!((val2 - 2.5).abs() < 1e-10);
                assert_eq!(result, 2);
            }
        }
    }

    #[test]
    fn test_scanf_ld_comprehensive() {
        unsafe {
            // Test 1: Simple long integer
            let (result, val) = test_sscanf_long_with_libc("123456", "%ld");
            assert_eq!(val, 123456);
            assert_eq!(result, 1);

            // Test 2: Negative long
            let (result, val) = test_sscanf_long_with_libc("-987654", "%ld");
            assert_eq!(val, -987654);
            assert_eq!(result, 1);

            // Test 3: Zero long
            let (result, val) = test_sscanf_long_with_libc("0", "%ld");
            assert_eq!(val, 0);
            assert_eq!(result, 1);

            // Test 4: Very large long - 32bit
            let (result, val) = test_sscanf_long_with_libc("2147483647", "%ld");
            assert_eq!(val, 2147483647);
            assert_eq!(result, 1);

            // Test 5: Multiple long integers
            {
                let mut val1: c_long = 0;
                let mut val2: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_long as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_long as *const c_void));

                let result = test_sscanf_internal("1000000 2000000", "%ld %ld", valist, 2);
                assert_eq!(val1, 1000000);
                assert_eq!(val2, 2000000);
                assert_eq!(result, 2);
            }

            // Test 6: Long hex
            {
                let mut val: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_long as *const c_void));

                let result = test_sscanf_internal("DEAD", "%lx", valist, 1);
                assert_eq!(val, 0xDEAD);
                assert_eq!(result, 1);
            }

            // Test 7: Long octal
            {
                let mut val: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_long as *const c_void));

                let result = test_sscanf_internal("177777", "%lo", valist, 1);
                assert_eq!(val, 0o177777);
                assert_eq!(result, 1);
            }

            // Test 8: Long with plus sign
            {
                let mut val: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_long as *const c_void));

                let result = test_sscanf_internal("+456789", "%ld", valist, 1);
                assert_eq!(val, 456789);
                assert_eq!(result, 1);
            }
        }
    }

    #[test]
    fn test_scanf_c_comprehensive() {
        unsafe {
            // Test 1: Single character
            {
                let mut val: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_char as *const c_void));

                let result = test_sscanf_internal("A", "%c", valist, 1);
                assert_eq!(val as u8, b'A');
                assert_eq!(result, 1);
            }

            // Test 2: Digit character
            {
                let mut val: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_char as *const c_void));

                let result = test_sscanf_internal("5", "%c", valist, 1);
                assert_eq!(val as u8, b'5');
                assert_eq!(result, 1);
            }

            // Test 3: Space character
            {
                let mut val: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_char as *const c_void));

                let result = test_sscanf_internal(" ", "%c", valist, 1);
                assert_eq!(val as u8, b' ');
                assert_eq!(result, 1);
            }

            // Test 4: Multiple characters
            {
                let mut val1: c_char = 0;
                let mut val2: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_char as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_char as *const c_void));

                let result = test_sscanf_internal("AB", "%c%c", valist, 2);
                assert_eq!(val1 as u8, b'A');
                assert_eq!(val2 as u8, b'B');
                assert_eq!(result, 2);
            }

            // Test 5: Characters with space format
            {
                let mut val1: c_char = 0;
                let mut val2: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_char as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_char as *const c_void));

                let result = test_sscanf_internal("A B", "%c %c", valist, 2);
                assert_eq!(val1 as u8, b'A');
                assert_eq!(val2 as u8, b'B');
                assert_eq!(result, 2);
            }

            // Test 6: Special characters
            {
                let mut val: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_char as *const c_void));

                let result = test_sscanf_internal("@", "%c", valist, 1);
                assert_eq!(val as u8, b'@');
                assert_eq!(result, 1);
            }
        }
    }

    #[test]
    fn test_scanf_s_comprehensive() {
        unsafe {
            // Test 1: Simple string
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("hello", "%s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "hello");
                assert_eq!(result, 1);
            }

            // Test 2: String with trailing text
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("hello world", "%s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "hello");
                assert_eq!(result, 1);
            }

            // Test 3: Multiple strings
            {
                let mut buffer1: [c_char; 100] = [0; 100];
                let mut buffer2: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer1.as_mut_ptr() as *const c_void));
                valist.push(VaArg::pointer(buffer2.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("hello world", "%s %s", valist, 2);
                let s1 = CStr::from_ptr(buffer1.as_ptr()).to_str().unwrap();
                let s2 = CStr::from_ptr(buffer2.as_ptr()).to_str().unwrap();
                assert_eq!(s1, "hello");
                assert_eq!(s2, "world");
                assert_eq!(result, 2);
            }

            // Test 4: String with numbers
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("test123", "%s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "test123");
                assert_eq!(result, 1);
            }

            // Test 5: Single character as string
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("A", "%s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "A");
                assert_eq!(result, 1);
            }

            // Test 6: String with special characters
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("@#$%", "%s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "@#$%");
                assert_eq!(result, 1);
            }
        }
    }

    #[test]
    fn test_scanf_mixed_comprehensive() {
        unsafe {
            // Test 1: Integer and string
            {
                let mut val: c_int = 0;
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("42 hello", "%d %s", valist, 2);
                assert_eq!(val, 42);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "hello");
                assert_eq!(result, 2);
            }

            // Test 2: Float and integer
            {
                let mut fval: c_double = 0.0;
                let mut ival: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut fval as *mut c_double as *const c_void));
                valist.push(VaArg::pointer(&mut ival as *mut c_int as *const c_void));

                let result = test_sscanf_internal("3.24 42", "%lf %d", valist, 2);
                assert!((fval - 3.24).abs() < 1e-10);
                assert_eq!(ival, 42);
                assert_eq!(result, 2);
            }

            // Test 3: Hex, decimal, and string
            {
                let mut hex_val: c_uint = 0;
                let mut dec_val: c_int = 0;
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut hex_val as *mut c_uint as *const c_void));
                valist.push(VaArg::pointer(&mut dec_val as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("FF 255 test", "%x %d %s", valist, 3);
                assert_eq!(hex_val, 255);
                assert_eq!(dec_val, 255);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "test");
                assert_eq!(result, 3);
            }

            // Test 4: Character and integer
            {
                let mut cval: c_char = 0;
                let mut ival: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut cval as *mut c_char as *const c_void));
                valist.push(VaArg::pointer(&mut ival as *mut c_int as *const c_void));

                let result = test_sscanf_internal("A 65", "%c %d", valist, 2);
                assert_eq!(cval as u8, b'A');
                assert_eq!(ival, 65);
                assert_eq!(result, 2);
            }

            // Test 5: Multiple formats with complex whitespace
            {
                let mut val1: c_int = 0;
                let mut val2: c_double = 0.0;
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_double as *const c_void));
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result =
                    test_sscanf_internal("  123  \t\n  2.5   word  ", "%d %lf %s", valist, 3);
                assert_eq!(val1, 123);
                assert!((val2 - 2.5).abs() < 1e-10);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "word");
                assert_eq!(result, 3);
            }
        }
    }

    #[test]
    fn test_scanf_edge_cases() {
        unsafe {
            // Test 1: Empty input
            {
                let mut val: c_int = 999;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                let result = test_sscanf_internal("", "%d", valist, -1);
                assert_eq!(result, -1); // EOF
            }

            // Test 2: Whitespace only input
            {
                let mut val: c_int = 999;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                let result = test_sscanf_internal("   ", "%d", valist, -1);
                assert_eq!(result, -1); // EOF
            }

            // Test 3: Invalid number format
            {
                let mut val: c_int = 999;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                let result = test_sscanf_internal("abc", "%d", valist, 0);
                assert_eq!(result, 0); // No matches
                assert_eq!(val, 999); // Value unchanged
            }

            // Test 4: Partial match
            {
                let mut val1: c_int = 999;
                let mut val2: c_int = 888;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_int as *const c_void));

                let result = test_sscanf_internal("42 abc", "%d %d", valist, 1);
                assert_eq!(result, 1); // Only first match
                assert_eq!(val1, 42);
                assert_eq!(val2, 888); // Second value unchanged
            }

            // Test 5: Overflow behavior (large number)
            {
                let mut val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                // This should handle overflow gracefully
                let result = test_sscanf_internal("999999999999999999999", "%d", valist, 1);
                assert_eq!(result, 1);
                // Value might wrap around due to overflow
            }

            // Test 6: Hex without valid hex digits
            {
                let mut val: c_uint = 999;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_uint as *const c_void));

                let result = test_sscanf_internal("xyz", "%x", valist, 0);
                assert_eq!(result, 0);
                assert_eq!(val, 999); // Value unchanged
            }

            // Test 7: Float with invalid format
            {
                let mut val: c_double = 999.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("not.a.number", "%lf", valist, 0);
                assert_eq!(result, 0);
                assert_eq!(val, 999.0); // Value unchanged
            }

            // Test 8: Missing second operand
            {
                let mut val1: c_int = 0;
                let mut val2: c_int = 999;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_int as *const c_void));

                let result = test_sscanf_internal("42", "%d %d", valist, 1);
                assert_eq!(result, 1);
                assert_eq!(val1, 42);
                assert_eq!(val2, 999); // Second value unchanged
            }
        }
    }

    #[test]
    fn test_scanf_width_specifiers() {
        unsafe {
            // Test 1: Width specifier for integer
            {
                let mut val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                let result = test_sscanf_internal("12345", "%3d", valist, 1);
                assert_eq!(val, 123);
                assert_eq!(result, 1);
            }

            // Test 2: Width specifier for string
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("hellothere", "%5s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "hello");
                assert_eq!(result, 1);
            }

            // Test 3: Width specifier for hex
            {
                let mut val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_uint as *const c_void));

                let result = test_sscanf_internal("ABCDEF", "%4x", valist, 1);
                assert_eq!(val, 0xABCD);
                assert_eq!(result, 1);
            }

            // Test 4: Width specifier for character
            {
                let mut buffer: [c_char; 10] = [0; 10];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("ABCDE", "%3c", valist, 1);
                assert_eq!(buffer[0] as u8, b'A');
                assert_eq!(buffer[1] as u8, b'B');
                assert_eq!(buffer[2] as u8, b'C');
                assert_eq!(result, 1);
            }

            // Test 5: Multiple width specifiers
            {
                let mut val1: c_int = 0;
                let mut val2: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_int as *const c_void));

                let result = test_sscanf_internal("123456789", "%3d%4d", valist, 2);
                assert_eq!(val1, 123);
                assert_eq!(val2, 4567);
                assert_eq!(result, 2);
            }
        }
    }

    #[test]
    fn test_scanf_percent_literal() {
        unsafe {
            // Test 1: Single percent literal
            {
                let mut val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                let result = test_sscanf_internal("%42", "%%%d", valist, 1);
                assert_eq!(val, 42);
                assert_eq!(result, 1);
            }

            // Test 2: Percent literal with multiple values
            {
                let mut val1: c_int = 0;
                let mut val2: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_int as *const c_void));

                let result = test_sscanf_internal("42%100", "%d%%%d", valist, 2);
                assert_eq!(val1, 42);
                assert_eq!(val2, 100);
                assert_eq!(result, 2);
            }

            // Test 3: Mismatched percent
            {
                let mut val: c_int = 999;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                let result = test_sscanf_internal("42", "%%%d", valist, 0);
                assert_eq!(result, 0); // Should fail to match
                assert_eq!(val, 999); // Value unchanged
            }
        }
    }

    #[test]
    fn test_scanf_many() {
        unsafe {
            // Test 1: Range of integers -32767 to 32767 (sampling)
            for i in [-32767, -1000, -1, 0, 1, 1000, 32767].iter() {
                let mut val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_int as *const c_void));

                let input = format!("{i}");
                let result = test_sscanf_internal(&input, "%d", valist, 1);
                assert_eq!(val, *i);
                assert_eq!(result, 1);
            }

            // Test 2: Whitespace handling
            {
                let mut val1: c_int = 0;
                let mut val2: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_int as *const c_void));

                let result = test_sscanf_internal("12 \t\n32", "%d %d", valist, 2);
                assert_eq!(val1, 12);
                assert_eq!(val2, 32);
                assert_eq!(result, 2);
            }

            // Test 3: Character parsing
            {
                let mut val1: c_char = 0;
                let mut val2: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_char as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_char as *const c_void));

                let result = test_sscanf_internal("a b", "%c %c", valist, 2);
                assert_eq!(val1 as u8, b'a');
                assert_eq!(val2 as u8, b'b');
                assert_eq!(result, 2);
            }

            // Test 4: String parsing
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("hellothere", "%s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "hellothere");
                assert_eq!(result, 1);
            }

            // Test 5: String with width specifier
            {
                let mut buffer: [c_char; 100] = [0; 100];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal("hellothere", "%5s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "hello");
                assert_eq!(result, 1);
            }

            // Test 6: Mixed numeric formats
            {
                let mut dec_val: c_int = 0;
                let mut hex_val: c_uint = 0;
                let mut oct_val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut dec_val as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut hex_val as *mut c_uint as *const c_void));
                valist.push(VaArg::pointer(&mut oct_val as *mut c_uint as *const c_void));

                let result = test_sscanf_internal("255 FF 377", "%d %x %o", valist, 3);
                assert_eq!(dec_val, 255);
                assert_eq!(hex_val, 255);
                assert_eq!(oct_val, 255);
                assert_eq!(result, 3);
            }

            // Test 7: Long integers
            {
                let mut val: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_long as *const c_void));

                let result = test_sscanf_internal("1234567890", "%ld", valist, 1);
                assert_eq!(val, 1234567890);
                assert_eq!(result, 1);
            }

            // Test 8: Floating point
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("3.24159", "%lf", valist, 1);
                assert!((val - 3.24159).abs() < 1e-10);
                assert_eq!(result, 1);
            }

            // Test 9: Scientific notation
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("1.5e10", "%lf", valist, 1);
                assert!((val - 1.5e10).abs() < 1e5);
                assert_eq!(result, 1);
            }

            // Test 10: Complex mixed format
            {
                let mut ival: c_int = 0;
                let mut buffer: [c_char; 100] = [0; 100];
                let mut fval: c_double = 0.0;
                let mut cval: c_char = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut ival as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));
                valist.push(VaArg::pointer(&mut fval as *mut c_double as *const c_void));
                valist.push(VaArg::pointer(&mut cval as *mut c_char as *const c_void));

                let result = test_sscanf_internal("42 hello 3.24 X", "%d %s %lf %c", valist, 4);
                assert_eq!(ival, 42);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, "hello");
                assert!((fval - 3.24).abs() < 1e-10);
                assert_eq!(cval as u8, b'X');
                assert_eq!(result, 4);
            }
        }
    }

    #[test]
    fn test_scanf_extreme_edge_cases() {
        unsafe {
            // Test 1: Maximum string length
            {
                let long_string = "a".repeat(500);
                let mut buffer: [c_char; 1000] = [0; 1000];
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(buffer.as_mut_ptr() as *const c_void));

                let result = test_sscanf_internal(&long_string, "%s", valist, 1);
                let s = CStr::from_ptr(buffer.as_ptr()).to_str().unwrap();
                assert_eq!(s, long_string);
                assert_eq!(result, 1);
            }

            // Test 2: Very long number
            {
                let mut val: c_long = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_long as *const c_void));

                let result = test_sscanf_internal("2147483647", "%ld", valist, 1);
                assert_eq!(val, 2147483647);
                assert_eq!(result, 1);
            }

            // Test 3: Many consecutive format specifiers
            {
                let mut vals: [c_int; 10] = [0; 10];
                let mut valist = CustomVaList::new();
                for i in 0..10 {
                    valist.push(VaArg::pointer(&mut vals[i] as *mut c_int as *const c_void));
                }

                let result = test_sscanf_internal(
                    "1 2 3 4 5 6 7 8 9 10",
                    "%d %d %d %d %d %d %d %d %d %d",
                    valist,
                    10,
                );
                for i in 0..10 {
                    assert_eq!(vals[i], (i + 1) as c_int);
                }
                assert_eq!(result, 10);
            }

            // Test 4: Extremely small float
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("1e-100", "%lf", valist, 1);
                assert!((val - 1e-100).abs() < 1e-105);
                assert_eq!(result, 1);
            }

            // Test 5: Extremely large float
            {
                let mut val: c_double = 0.0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_double as *const c_void));

                let result = test_sscanf_internal("1e50", "%lf", valist, 1);
                assert!((val - 1e50).abs() < 1e45);
                assert_eq!(result, 1);
            }

            // Test 6: All hex digits
            {
                let mut val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_uint as *const c_void));

                let result = test_sscanf_internal("01234567", "%x", valist, 1);
                assert_eq!(val, 0x01234567);
                assert_eq!(result, 1);
            }

            // Test 7: All octal digits
            {
                let mut val: c_uint = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val as *mut c_uint as *const c_void));

                let result = test_sscanf_internal("01234567", "%o", valist, 1);
                assert_eq!(val, 0o1234567);
                assert_eq!(result, 1);
            }

            // Test 8: Mixed with many whitespace variations
            {
                let mut val1: c_int = 0;
                let mut val2: c_int = 0;
                let mut val3: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut val1 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val2 as *mut c_int as *const c_void));
                valist.push(VaArg::pointer(&mut val3 as *mut c_int as *const c_void));

                let result = test_sscanf_internal(
                    "  \t\n  1   \t\n\t   2  \n\n  3  \t",
                    "%d %d %d",
                    valist,
                    3,
                );
                assert_eq!(val1, 1);
                assert_eq!(val2, 2);
                assert_eq!(val3, 3);
                assert_eq!(result, 3);
            }
        }
    }

    // Comparison tests with libc::sscanf for failing cases
    #[test]
    fn test_libc_comparison_overflow() {
        unsafe {
            // Test the overflow case that's failing
            let input = "999999999999999999999";
            let format = "%d";

            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_int = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_int,
                );

                // Test with our implementation
                let mut our_val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(our_result, libc_result, "Return values should match");
                if libc_result > 0 {
                    assert_eq!(our_val, libc_val, "Parsed values should match");
                }
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // On Windows, just verify our implementation runs without panicking
                assert!(
                    our_result >= 0,
                    "Our implementation should handle overflow gracefully"
                );
            }
        }
    }

    #[test]
    fn test_libc_comparison_percent_literal() {
        unsafe {
            // Test the percent literal case that's failing
            let input = "%42";
            let format = "%%%d";

            #[cfg(not(target_family = "windows"))]
            {
                // Test with libc::sscanf
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_int = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_int,
                );

                // Test with our implementation
                let mut our_val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // Compare results
                assert_eq!(our_result, libc_result, "Return values should match");
                if libc_result > 0 {
                    assert_eq!(our_val, libc_val, "Parsed values should match");
                }
            }
            #[cfg(target_family = "windows")]
            {
                // On Windows, just test our implementation
                let mut our_val: c_int = 0;
                let mut valist = CustomVaList::new();
                valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                let input_cstr = to_c_string(input);
                let format_cstr = to_c_string(format);
                let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                // On Windows, just verify our implementation runs
                assert!(
                    our_result >= 0,
                    "Our implementation should handle percent literals"
                );
            }
        }
    }

    #[test]
    fn test_debug_minimal_percent() {
        unsafe {
            // Very simple test: just %% followed by %d
            let input = "%1";
            let format = "%%%d";

            let mut our_val: c_int = 0;
            let mut valist = CustomVaList::new();
            valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

            let input_cstr = to_c_string(input);
            let format_cstr = to_c_string(format);
            let our_result = sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

            println!(
                "MINIMAL TEST: input='{input}', format='{format}', result={our_result}, value={our_val}"
            );

            #[cfg(not(target_family = "windows"))]
            {
                // Test libc too
                let input_cstring = CString::new(input).unwrap();
                let format_cstring = CString::new(format).unwrap();
                let mut libc_val: c_int = 0;
                let libc_result = libc::sscanf(
                    input_cstring.as_ptr(),
                    format_cstring.as_ptr(),
                    &mut libc_val as *mut c_int,
                );
                println!("LIBC: result={libc_result}, value={libc_val}");
            }
            #[cfg(target_family = "windows")]
            {
                println!("WINDOWS: Only testing our implementation");
            }
        }
    }

    #[test]
    fn test_debug_percent_step_by_step() {
        unsafe {
            // Test step 1: Just matching a single %
            {
                let input = "%";
                let format = "%%";

                #[cfg(not(target_family = "windows"))]
                {
                    let input_cstring = CString::new(input).unwrap();
                    let format_cstring = CString::new(format).unwrap();
                    let libc_result = libc::sscanf(input_cstring.as_ptr(), format_cstring.as_ptr());

                    let input_cstr = to_c_string(input);
                    let format_cstr = to_c_string(format);
                    let valist = CustomVaList::new(); // no arguments for %%
                    let our_result =
                        sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                    println!(
                        "Step 1 - input: '{input}', format: '{format}' - libc: {libc_result}, ours: {our_result}"
                    );
                }
                #[cfg(target_family = "windows")]
                {
                    let input_cstr = to_c_string(input);
                    let format_cstr = to_c_string(format);
                    let valist = CustomVaList::new(); // no arguments for %%
                    let our_result =
                        sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                    println!("Step 1 - input: '{input}', format: '{format}' - ours: {our_result}");
                }
            }

            // Test step 2: % followed by a number
            {
                let input = "%42";
                let format = "%%";

                #[cfg(not(target_family = "windows"))]
                {
                    let input_cstring = CString::new(input).unwrap();
                    let format_cstring = CString::new(format).unwrap();
                    let libc_result = libc::sscanf(input_cstring.as_ptr(), format_cstring.as_ptr());

                    let input_cstr = to_c_string(input);
                    let format_cstr = to_c_string(format);
                    let valist = CustomVaList::new(); // no arguments for %%
                    let our_result =
                        sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                    println!(
                        "Step 2 - input: '{input}', format: '{format}' - libc: {libc_result}, ours: {our_result}"
                    );
                }
                #[cfg(target_family = "windows")]
                {
                    let input_cstr = to_c_string(input);
                    let format_cstr = to_c_string(format);
                    let valist = CustomVaList::new(); // no arguments for %%
                    let our_result =
                        sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                    println!("Step 2 - input: '{input}', format: '{format}' - ours: {our_result}");
                }
            }

            // Test step 3: Just parse the number after %
            {
                let input = "42";
                let format = "%d";

                #[cfg(not(target_family = "windows"))]
                {
                    let input_cstring = CString::new(input).unwrap();
                    let format_cstring = CString::new(format).unwrap();
                    let mut libc_val: c_int = 0;
                    let libc_result = libc::sscanf(
                        input_cstring.as_ptr(),
                        format_cstring.as_ptr(),
                        &mut libc_val as *mut c_int,
                    );

                    let mut our_val: c_int = 0;
                    let mut valist = CustomVaList::new();
                    valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                    let input_cstr = to_c_string(input);
                    let format_cstr = to_c_string(format);
                    let our_result =
                        sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                    println!(
                        "Step 3 - input: '{input}', format: '{format}' - libc: {libc_result} (val: {libc_val}), ours: {our_result} (val: {our_val})"
                    );
                }
                #[cfg(target_family = "windows")]
                {
                    let mut our_val: c_int = 0;
                    let mut valist = CustomVaList::new();
                    valist.push(VaArg::pointer(&mut our_val as *mut c_int as *const c_void));

                    let input_cstr = to_c_string(input);
                    let format_cstr = to_c_string(format);
                    let our_result =
                        sscanf_internal(input_cstr.as_ptr(), format_cstr.as_ptr(), valist);

                    println!(
                        "Step 3 - input: '{input}', format: '{format}' - ours: {our_result} (val: {our_val})"
                    );
                }
            }
        }
    }
}
