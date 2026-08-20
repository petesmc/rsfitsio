#[cfg(not(windows))]
use libc::{c_char, c_int};
#[cfg(not(windows))]
use rsfitsio::drvrsmem::{
    shared_getaddr, shared_init, shared_list, shared_recover, shared_uncond_delete,
};
use std::ptr;
use std::{ffi::CStr, process::ExitCode};

#[cfg(windows)]
pub fn main() -> ExitCode {
    println!("smem not supported on windows");
    ExitCode::from(10)
}

/// `sscanf(s, "%d", &out)`: skip leading whitespace, take an optional sign and
/// one or more digits, ignore any trailing garbage. `None` when nothing was
/// converted, which is the `!= 1` return the C checks for.
#[cfg(not(windows))]
fn scan_int(s: &str) -> Option<c_int> {
    let s = s.trim_start();
    let (sign, rest) = match s.strip_prefix('-') {
        Some(r) => (-1, r),
        None => (1, s.strip_prefix('+').unwrap_or(s)),
    };

    let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
    if digits.is_empty() {
        return None;
    }

    digits.parse::<c_int>().ok().map(|v| sign * v)
}

#[cfg(not(windows))]
pub fn main() -> ExitCode {
    let mut cmdok: bool = true;
    let mut listmode: bool = false;
    let _longlistmode: bool;
    let mut recovermode: bool = false;
    let mut deletemode: bool = false;
    let mut id: c_int = -1;
    let mut status: c_int;

    let args = std::env::args();
    let args: Vec<String> = args.collect();

    let argc = args.len();

    unsafe {
        match argc {
            1 => {
                listmode = true;
            }

            2 => match args[1].as_str() {
                "-l" => _longlistmode = true,
                "-r" => recovermode = true,
                "-d" => deletemode = true,
                _ => cmdok = false,
            },
            3 => {
                let mut c = true;

                match args[1].as_str() {
                    "-r" => recovermode = true,
                    "-d" => deletemode = true,
                    _ => {
                        cmdok = false; /* signal invalid cmd line syntax */
                        c = false;
                    }
                }

                /* C: sscanf(argv[2], "%d", &id) != 1 */
                match c.then(|| scan_int(&args[2])).flatten() {
                    Some(v) => id = v,
                    None => {
                        if c {
                            cmdok = false;
                        }
                    }
                }
            }
            _ => {
                cmdok = false;
            }
        }

        if !cmdok {
            print!("usage :\n\n");
            println!("smem            - list all shared memory segments");
            println!("\t!\tcouldn't obtain RDONLY lock - info unreliable");
            println!("\tIdx\thandle of shared memory segment (visible by application)");
            println!("\tKey\tcurrent system key of shared memory segment. Key");
            println!("\t\tchanges whenever shmem segment is reallocated. Use");
            println!("\t\tipcs (or ipcs -a) to view all shmem segments");
            println!("\tNproc\tnumber of processes attached to segment");
            println!("\tSize\tsize of shmem segment in bytes");
            println!("\tFlags\tRESIZABLE - realloc allowed, PERSIST - segment is not");
            println!("\t\tdeleted after shared_free called by last process attached");
            println!("\t\tto it.");
            println!("smem -d         - delete all shared memory segments (may block)");
            println!("smem -d id      - delete specific shared memory segment (may block)");
            println!("smem -r         - unconditionally reset all shared memory segments");
            println!("\t\t(does not block, recovers zombie handles left by kill -9)");
            println!("smem -r id      - unconditionally reset specific shared memory segment");
        }

        if shared_init(0) != 0 {
            println!("couldn't initialize shared memory, aborting ...");
            return ExitCode::from(10);
        }

        if listmode {
            shared_list(id);
        } else if recovermode {
            shared_recover(id);
        } else if deletemode {
            shared_uncond_delete(id);
        }

        for id in 0..16 {
            let mut address: *mut c_char = ptr::null_mut();

            status = shared_getaddr(id, &mut address);

            if status == 0 {
                let address_str = CStr::from_ptr(address).to_str().unwrap();
                println!("id, status, address {id} {status} {address:p} {address_str:30}");
            }
        }
        ExitCode::from(0)
    }
}

#[cfg(all(test, not(windows)))]
mod tests {
    use super::scan_int;

    /* Expected values are C sscanf("%d") semantics, hand-derived: leading
    whitespace skipped, optional sign, digits, trailing garbage ignored, and
    a 0 return (here None) when no digits were converted. */
    #[test]
    fn scan_int_matches_sscanf_percent_d() {
        assert_eq!(scan_int("0"), Some(0));
        assert_eq!(scan_int("7"), Some(7));
        assert_eq!(scan_int("15"), Some(15));
        assert_eq!(scan_int("+15"), Some(15));
        assert_eq!(scan_int("-3"), Some(-3));
        assert_eq!(scan_int("  42"), Some(42));
        assert_eq!(scan_int("\t\n5"), Some(5));
        assert_eq!(scan_int("007"), Some(7));
        /* trailing garbage is left in the stream, the conversion still counts */
        assert_eq!(scan_int("12abc"), Some(12));
        assert_eq!(scan_int("3 4"), Some(3));

        /* no digits converted -> sscanf returns 0, not 1 */
        assert_eq!(scan_int(""), None);
        assert_eq!(scan_int("   "), None);
        assert_eq!(scan_int("abc"), None);
        assert_eq!(scan_int("-"), None);
        assert_eq!(scan_int("+"), None);
        assert_eq!(scan_int("-abc"), None);
        /* whitespace between sign and digits is not accepted by %d */
        assert_eq!(scan_int("- 5"), None);
    }
}
