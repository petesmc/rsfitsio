#![allow(deprecated)]

#[cfg(not(windows))]
use libc::{c_char, c_int, sscanf, strcmp};
#[cfg(not(windows))]
use rsfitsio::drvrsmem::{
    shared_getaddr, shared_init, shared_list, shared_recover, shared_uncond_delete,
};
use std::ffi::CString;
use std::ptr;
use std::{ffi::CStr, process::ExitCode};

#[cfg(windows)]
pub fn main() -> ExitCode {
    println!("smem not supported on windows");
    ExitCode::from(10)
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

            2 => {
                let a1 = CString::new(args[1].as_str()).unwrap().as_c_str().as_ptr();

                if 0 == strcmp(c"-l".as_ptr(), a1) {
                    _longlistmode = true;
                } else if 0 == strcmp(c"-r".as_ptr(), a1) {
                    recovermode = true;
                } else if 0 == strcmp(c"-d".as_ptr(), a1) {
                    deletemode = true;
                } else {
                    cmdok = false;
                }
            }
            3 => {
                let mut c = true;

                let a1 = CString::new(args[1].as_str()).unwrap().as_c_str().as_ptr();
                let a2 = CString::new(args[2].as_str()).unwrap().as_c_str().as_ptr();

                if 0 == strcmp(c"-r".as_ptr(), a1) {
                    recovermode = true;
                } else if 0 == strcmp(c"-d".as_ptr(), a1) {
                    deletemode = true;
                } else {
                    cmdok = false; /* signal invalid cmd line syntax */
                    c = false;
                }

                if c && (1 != sscanf(a2, c"%d".as_ptr(), &mut id)) {
                    cmdok = false;
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

            // NOTE: shared_getaddr leaves each successfully-probed segment
            // attached and locked (it deliberately omits smem_close) — see its
            // lock/lifetime contract in src/drvrsmem.rs. This loop never frees,
            // so it leaks one lock+attach per live segment. That matches the C
            // `smem` tool's behaviour; this short-lived CLI exits immediately
            // afterwards, so the OS reclaims everything.
            status = shared_getaddr(id, &mut address);

            if status == 0 {
                let address_str = CStr::from_ptr(address).to_str().unwrap();
                println!("id, status, address {id} {status} {address:p} {address_str:30}");
            }
        }
        ExitCode::from(0)
    }
}
