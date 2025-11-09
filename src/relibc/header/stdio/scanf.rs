use super::lookaheadreader::LookAheadReader;
use super::printf::{CustomVaList, VaArg};
use crate::c_types::{
    c_char, c_double, c_float, c_int, c_long, c_longlong, c_short, c_uchar, c_uint, c_ulong,
    c_ulonglong, c_ushort, intmax_t, ptrdiff_t, size_t, ssize_t, uintmax_t,
};

use std::string::String;
use std::vec::Vec;

#[derive(PartialEq, Eq)]
enum IntKind {
    Byte,
    Short,
    Int,
    Long,
    LongLong,
    IntMax,
    PtrDiff,
    Size,
}

/// Helper function for progressing a C string
unsafe fn next_byte(string: &mut *const c_char) -> Result<u8, c_int> {
    unsafe {
        let c = **string as u8;
        *string = string.offset(1);
        if c == 0 { Err(-1) } else { Ok(c) }
    }
}

pub unsafe fn scanf_custom(
    r: LookAheadReader,
    format: *const c_char,
    mut ap: CustomVaList,
) -> c_int {
    unsafe {
        match inner_scanf_custom(r, format, &mut ap) {
            Ok(n) => n,
            Err(n) => n,
        }
    }
}

unsafe fn inner_scanf_custom(
    mut r: LookAheadReader,
    mut format: *const c_char,
    ap: &mut CustomVaList,
) -> Result<c_int, c_int> {
    unsafe {
        let mut matched = 0;
        let mut byte = 0;
        let mut skip_read = false;
        let mut first_read = true;

        macro_rules! read {
            () => {{
                match r.lookahead1() {
                    Ok(None) => false,
                    Ok(Some(b)) => {
                        byte = b;
                        true
                    }
                    Err(x) => return Err(x),
                }
            }};
        }

        macro_rules! maybe_read {
            () => {
                maybe_read!(inner false);
            };
            (noreset) => {
                maybe_read!(inner);
            };
            (inner $($placeholder:expr)*) => {
                if !skip_read && !read!() {
                    // If this is the first read attempt and we have empty input,
                    // return EOF (-1) instead of 0 to match libc behavior
                    if first_read && matched == 0 {
                        return Err(-1);
                    }
                    return Ok(matched);
                }
                $(else {
                    // Hacky way of having this optional
                    skip_read = $placeholder;
                })*
                first_read = false;
            }
        }

        // Follow the same structure as the original inner_scanf
        while *format != 0 {
            let mut c = next_byte(&mut format)?;

            if c == b' ' {
                maybe_read!(noreset);

                while (byte as char).is_whitespace() {
                    if !read!() {
                        return Ok(matched);
                    }
                }

                skip_read = true;
            } else if c != b'%' {
                maybe_read!();
                if c != byte {
                    return Ok(matched);
                }
                r.commit();
            } else {
                c = next_byte(&mut format)?;

                let mut ignore = false;
                if c == b'*' {
                    ignore = true;
                    c = next_byte(&mut format)?;
                }

                let mut width = String::new();
                while c >= b'0' && c <= b'9' {
                    width.push(c as char);
                    c = next_byte(&mut format)?;
                }
                let mut width = if width.is_empty() {
                    None
                } else {
                    match width.parse::<usize>() {
                        Ok(n) => Some(n),
                        Err(_) => return Err(-1),
                    }
                };

                // When an EOF occurs, eof is set, stuff is marked matched
                // as usual, and finally it is returned
                let mut eof = false;

                let mut kind = IntKind::Int;

                // Handle length modifiers (mostly following original)
                match c {
                    b'h' => {
                        kind = IntKind::Short;
                        c = next_byte(&mut format)?;
                        if c == b'h' {
                            kind = IntKind::Byte;
                            c = next_byte(&mut format)?;
                        }
                    }
                    b'l' => {
                        kind = IntKind::Long;
                        c = next_byte(&mut format)?;
                        if c == b'l' {
                            kind = IntKind::LongLong;
                            c = next_byte(&mut format)?;
                        }
                    }
                    b'j' => {
                        kind = IntKind::IntMax;
                        c = next_byte(&mut format)?;
                    }
                    b'z' => {
                        kind = IntKind::Size;
                        c = next_byte(&mut format)?;
                    }
                    b't' => {
                        kind = IntKind::PtrDiff;
                        c = next_byte(&mut format)?;
                    }
                    _ => {}
                }

                if c != b'n' {
                    maybe_read!(noreset);
                }
                match c {
                    b'%' => {
                        while (byte as char).is_whitespace() {
                            if !read!() {
                                // If no items have been matched yet, return EOF
                                if matched == 0 {
                                    return Err(-1);
                                }
                                return Ok(matched);
                            }
                        }

                        if byte != b'%' {
                            return Ok(matched);
                        }
                        r.commit();
                        skip_read = false; // Reset skip_read so next format specifier reads fresh input
                    }
                    b'd' | b'i' | b'o' | b'u' | b'x' | b'X' | b'f' | b'e' | b'g' | b'E' | b'a'
                    | b'p' => {
                        while (byte as char).is_whitespace() {
                            if !read!() {
                                // If no items have been matched yet, return EOF
                                if matched == 0 {
                                    return Err(-1);
                                }
                                return Ok(matched);
                            }
                        }

                        let pointer = c == b'p';
                        // Pointers aren't automatic, but we do want to parse "0x"
                        let auto = c == b'i' || pointer;
                        let float = c == b'f' || c == b'e' || c == b'g' || c == b'E' || c == b'a';

                        let mut radix = match c {
                            b'o' => 8,
                            b'x' | b'X' | b'p' => 16,
                            _ => 10,
                        };

                        let mut n = String::new();
                        let mut dot = false;
                        let mut negative = false;
                        let mut has_sign = false;

                        // Handle leading sign for numeric conversions
                        if !float && (byte == b'-' || byte == b'+') {
                            negative = byte == b'-';
                            has_sign = true;
                            r.commit();
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                return Ok(matched);
                            }
                        } else if float && (byte == b'-' || byte == b'+') {
                            n.push(byte as char);
                            r.commit();
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                return Ok(matched);
                            }
                        }

                        // Skip 0x prefix for %x format
                        if (c == b'x' || c == b'X')
                            && byte == b'0'
                            && width.map(|w| w > 0).unwrap_or(true)
                        {
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && read!() {
                                if byte == b'x' || byte == b'X' {
                                    width = width.map(|w| w - 1);
                                    if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                        return Ok(matched);
                                    }
                                } else {
                                    // It was just a 0, put it back
                                    n.push('0');
                                    r.commit();
                                }
                            } else {
                                // Just a 0 at EOF
                                n.push('0');
                                r.commit();
                                eof = true;
                            }
                        }

                        while width.map(|w| w > 0).unwrap_or(true)
                            && ((byte >= b'0' && byte <= b'7')
                                || (radix >= 10 && (byte >= b'8' && byte <= b'9'))
                                || (float && !dot && byte == b'.')
                                || (radix == 16
                                    && ((byte >= b'a' && byte <= b'f')
                                        || (byte >= b'A' && byte <= b'F'))))
                        {
                            if auto
                                && n.is_empty()
                                && byte == b'0'
                                && width.map(|w| w > 0).unwrap_or(true)
                            {
                                if !pointer {
                                    radix = 8;
                                }
                                width = width.map(|w| w - 1);
                                r.commit();
                                if !read!() {
                                    n.push('0');
                                    break;
                                }
                                if width.map(|w| w > 0).unwrap_or(true)
                                    && (byte == b'x' || byte == b'X')
                                {
                                    radix = 16;
                                    width = width.map(|w| w - 1);
                                    r.commit();
                                    if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                        return Ok(matched);
                                    }
                                } else {
                                    skip_read = true;
                                }
                                continue;
                            }
                            if byte == b'.' {
                                // Don't allow another dot
                                dot = true;
                            }
                            n.push(byte as char);
                            r.commit();
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                eof = true;
                                break;
                            }
                        }

                        // Handle scientific notation for floats
                        if float
                            && width.map(|w| w > 0).unwrap_or(true)
                            && !n.is_empty()
                            && (byte == b'e' || byte == b'E')
                        {
                            n.push(byte as char);
                            r.commit();
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && read!() {
                                if byte == b'+' || byte == b'-' {
                                    n.push(byte as char);
                                    r.commit();
                                    width = width.map(|w| w - 1);
                                    if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                        break;
                                    }
                                }
                                // Read exponent digits
                                while width.map(|w| w > 0).unwrap_or(true)
                                    && byte >= b'0'
                                    && byte <= b'9'
                                {
                                    n.push(byte as char);
                                    r.commit();
                                    width = width.map(|w| w - 1);
                                    if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                        break;
                                    }
                                }
                            }
                        }

                        // Custom parse_type macro for CustomVaList
                        macro_rules! parse_type_custom {
                            (noformat $type:ident) => {{
                                if n.is_empty() {
                                    return Ok(matched);
                                }
                                let n = n.parse::<$type>().map_err(|_| matched)?;
                                if !ignore {
                                    let ptr = match ap.arg() {
                                        VaArg::pointer(p) => p as *mut $type,
                                        _ => panic!("Expected pointer for conversion"),
                                    };
                                    *ptr = n;
                                    matched += 1;
                                }
                            }};
                            ($type:ident) => {
                                parse_type_custom!($type, $type)
                            };
                            ($type:ident, $final:ty) => {{
                                if n.is_empty() {
                                    return Ok(matched);
                                }
                                // Handle special case for minimum signed integer values and overflow
                                let val = if negative && has_sign && radix == 10 {
                                    // For negative numbers, parse as positive then negate
                                    match $type::from_str_radix(&n, radix) {
                                        Ok(v) => v.wrapping_neg(),
                                        Err(_) => {
                                            // Try parsing with a minus sign for edge cases like "-2147483648"
                                            let neg_n = format!("-{}", n);
                                            match $type::from_str_radix(&neg_n, radix) {
                                                Ok(v) => v,
                                                Err(_) => {
                                                    // Overflow case - return minimum value like libc
                                                    $type::MIN
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    match $type::from_str_radix(&n, radix) {
                                        Ok(v) => v,
                                        Err(_) => {
                                            // Overflow case - mimic libc behavior
                                            if $type::MIN == 0 {
                                                // unsigned type
                                                $type::MAX // wrap around for unsigned
                                            } else {
                                                // signed type
                                                (-1i32) as $type // -1 for signed overflow
                                            }
                                        }
                                    }
                                };
                                if !ignore {
                                    let ptr = match ap.arg() {
                                        VaArg::pointer(p) => p as *mut $final,
                                        _ => panic!("Expected pointer for conversion"),
                                    };
                                    *ptr = val as $final;
                                    matched += 1;
                                }
                            }};
                        }

                        if float {
                            if kind == IntKind::Long || kind == IntKind::LongLong {
                                parse_type_custom!(noformat c_double);
                            } else {
                                parse_type_custom!(noformat c_float);
                            }
                        } else {
                            let unsigned = c == b'o' || c == b'u' || c == b'x' || c == b'X';

                            match kind {
                                IntKind::Byte => {
                                    if unsigned {
                                        parse_type_custom!(c_uchar);
                                    } else {
                                        parse_type_custom!(c_char);
                                    }
                                }
                                IntKind::Short => {
                                    if unsigned {
                                        parse_type_custom!(c_ushort)
                                    } else {
                                        parse_type_custom!(c_short)
                                    }
                                }
                                IntKind::Int => {
                                    if unsigned {
                                        parse_type_custom!(c_uint)
                                    } else {
                                        parse_type_custom!(c_int)
                                    }
                                }
                                IntKind::Long => {
                                    if unsigned {
                                        parse_type_custom!(c_ulong)
                                    } else {
                                        parse_type_custom!(c_long)
                                    }
                                }
                                IntKind::LongLong => {
                                    if unsigned {
                                        parse_type_custom!(c_ulonglong)
                                    } else {
                                        parse_type_custom!(c_longlong)
                                    }
                                }
                                IntKind::IntMax => {
                                    if unsigned {
                                        parse_type_custom!(uintmax_t)
                                    } else {
                                        parse_type_custom!(intmax_t)
                                    }
                                }
                                IntKind::Size => {
                                    if unsigned {
                                        parse_type_custom!(size_t)
                                    } else {
                                        parse_type_custom!(ssize_t)
                                    }
                                }
                                IntKind::PtrDiff => {
                                    parse_type_custom!(ptrdiff_t)
                                }
                            }
                        }

                        // Don't return early due to eof here - we might still need to process %n or other
                        // format specifiers that don't consume input

                        if width != Some(0) && c != b'n' {
                            // It didn't hit the width, so an extra character was read and matched.
                            // But this character did not match so let's reuse it.
                            skip_read = true;
                        }
                    }
                    b's' => {
                        while (byte as char).is_whitespace() {
                            if !read!() {
                                // If no items have been matched yet, return EOF
                                if matched == 0 {
                                    return Err(-1);
                                }
                                return Ok(matched);
                            }
                        }

                        let mut ptr: Option<*mut c_char> = if ignore {
                            None
                        } else {
                            match ap.arg() {
                                VaArg::pointer(p) => Some(p as *mut c_char),
                                _ => panic!("Expected pointer for %s conversion"),
                            }
                        };

                        while width.map(|w| w > 0).unwrap_or(true)
                            && !(byte as char).is_whitespace()
                        {
                            if let Some(ref mut ptr) = ptr {
                                **ptr = byte as c_char;
                                *ptr = ptr.offset(1);
                            }
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                eof = true;
                                break;
                            }
                        }

                        if let Some(ptr) = ptr {
                            *ptr = 0;
                            matched += 1;
                            r.commit();
                        }
                    }
                    b'c' => {
                        let ptr: Option<*mut c_char> = if ignore {
                            None
                        } else {
                            match ap.arg() {
                                VaArg::pointer(p) => Some(p as *mut c_char),
                                _ => panic!("Expected pointer for %c conversion"),
                            }
                        };

                        for i in 0..width.unwrap_or(1) {
                            if let Some(ptr) = ptr {
                                *ptr.add(i) = byte as c_char;
                            }
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                eof = true;
                                break;
                            }
                        }

                        if ptr.is_some() {
                            matched += 1;
                            r.commit();
                        }
                    }
                    b'[' => {
                        c = next_byte(&mut format)?;

                        let mut matches = Vec::new();
                        let invert = if c == b'^' {
                            c = next_byte(&mut format)?;
                            true
                        } else {
                            false
                        };

                        let mut prev;
                        loop {
                            matches.push(c);
                            prev = c;
                            c = next_byte(&mut format)?;
                            if c == b'-' {
                                if prev == b']' {
                                    continue;
                                }
                                c = next_byte(&mut format)?;
                                if c == b']' {
                                    matches.push(b'-');
                                    break;
                                }
                                prev += 1;
                                while prev < c {
                                    matches.push(prev);
                                    prev += 1;
                                }
                            } else if c == b']' {
                                break;
                            }
                        }

                        let mut ptr: Option<*mut c_char> = if ignore {
                            None
                        } else {
                            match ap.arg() {
                                VaArg::pointer(p) => Some(p as *mut c_char),
                                _ => panic!("Expected pointer for %[ conversion"),
                            }
                        };

                        // While we haven't used up all the width, and it matches
                        let mut data_stored = false;
                        while width.map(|w| w > 0).unwrap_or(true)
                            && invert != matches.contains(&byte)
                        {
                            if let Some(ref mut ptr) = ptr {
                                **ptr = byte as c_char;
                                *ptr = ptr.offset(1);
                                data_stored = true;
                            }
                            r.commit();
                            // Decrease the width, and read a new character unless the width is 0
                            width = width.map(|w| w - 1);
                            if width.map(|w| w > 0).unwrap_or(true) && !read!() {
                                // Reading a new character has failed, return after
                                // actually marking this as matched
                                break;
                            }
                        }

                        if data_stored {
                            *ptr.unwrap() = 0;
                            matched += 1;
                        }
                    }
                    b'n' => {
                        if !ignore {
                            let ptr = match ap.arg() {
                                VaArg::pointer(p) => p as *mut c_int,
                                _ => panic!("Expected pointer for %n conversion"),
                            };
                            *ptr = r.position() as c_int;
                        }
                    }
                    _ => return Err(-1),
                }

                // If we hit EOF and we've successfully parsed at least one item, return.
                // But allow %n to be processed even at EOF since it doesn't consume input.
                if eof {
                    // Check if the next format character is %n
                    let next_is_n = unsafe {
                        *format == b'%' as c_char && {
                            let next_char = *format.offset(1);
                            // Skip any width specifiers and modifiers to find the actual format char
                            let mut temp_format = format.offset(1);
                            let mut temp_c = *temp_format as u8;
                            // Skip '*' (ignore flag)
                            if temp_c == b'*' {
                                temp_format = temp_format.offset(1);
                                temp_c = *temp_format as u8;
                            }
                            // Skip width digits
                            while temp_c >= b'0' && temp_c <= b'9' {
                                temp_format = temp_format.offset(1);
                                temp_c = *temp_format as u8;
                            }
                            // Skip length modifiers (h, l, j, z, t)
                            while temp_c == b'h' || temp_c == b'l' || temp_c == b'j' || temp_c == b'z' || temp_c == b't' {
                                temp_format = temp_format.offset(1);
                                temp_c = *temp_format as u8;
                            }
                            temp_c == b'n'
                        }
                    };
                    if !next_is_n {
                        return Ok(matched);
                    }
                }

                if width != Some(0) && c != b'n' && c != b'%' {
                    // It didn't hit the width, so an extra character was read and matched.
                    // But this character did not match so let's reuse it.
                    // Note: Exclude % literal case as it should not affect skip_read
                    skip_read = true;
                }
            }
        }

        Ok(matched)
    }
}
