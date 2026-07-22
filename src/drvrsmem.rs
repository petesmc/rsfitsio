#![allow(static_mut_refs)]

/* configuration parameters */

use core::slice;
use std::{
    cmp,
    ffi::{CStr, c_void},
    mem::MaybeUninit,
    ptr,
};

use bytemuck::{cast_slice, cast_slice_mut};
use libc::{
    EACCES, EAGAIN, F_RDLCK, F_SETLK, F_SETLKW, F_UNLCK, F_WRLCK, GETVAL, IPC_CREAT, IPC_EXCL,
    IPC_RMID, IPC_STAT, O_CREAT, O_EXCL, O_RDWR, O_TRUNC, SEM_UNDO, atexit, atoi, close, fcntl,
    flock, getenv, memcpy, mode_t, open, sembuf, semctl, semget, semid_ds, semop, shmat, shmctl,
    shmdt, shmget, shmid_ds, umask,
};

use crate::{
    BL,
    fitsio::{
        FLEN_FILENAME, LONGLONG, READWRITE, SHARED_AGAIN, SHARED_BADARG, SHARED_ERRBASE,
        SHARED_IPCERR, SHARED_NOFILE, SHARED_NOMEM, SHARED_NORESIZE, SHARED_NOTINIT, SHARED_NULPTR,
    },
};
use crate::{
    c_types::{c_char, c_int, c_long, c_short, c_ulong, c_ushort},
    cs, int_snprintf,
    relibc::header::stdio::sscanf_d,
    wrappers::strcpy_safe,
};

const SHARED_MAXSEG_VAL: c_int = 16; /* maximum number of shared memory blocks */

const SHARED_KEYBASE: c_int = 14011963; /* base for shared memory keys, may be overriden by getenv */
const SHARED_FDNAME: &CStr = c"/tmp/.shmem-lockfile"; /* template for lock file name */

const SHARED_ENV_KEYBASE: &CStr = c"SHMEM_LIB_KEYBASE"; /* name of environment variable */
const SHARED_ENV_MAXSEG: &CStr = c"SHMEM_LIB_MAXSEG"; /* name of environment variable */

/* useful constants */

const SHARED_RDONLY: c_int = 0; /* flag for shared_(un)lock, lock for read */
const SHARED_RDWRITE: c_int = 1; /* flag for shared_(un)lock, lock for write */
const SHARED_WAIT: c_int = 0; /* flag for shared_lock, block if cannot lock immediate */
const SHARED_NOWAIT: c_int = 2; /* flag for shared_lock, fail if cannot lock immediate */
const SHARED_NOLOCK: c_int = 0x100; /* flag for shared_validate function */

const SHARED_RESIZE: c_int = 4; /* flag for shared_malloc, object is resizeable */
const SHARED_PERSIST: c_int = 8; /* flag for shared_malloc, object is not deleted after last proc detaches */

const SHARED_INVALID: c_int = -1; /* invalid handle for semaphore/shared memory */

const SHARED_EMPTY: c_int = 0; /* entries for shared_used table */
const SHARED_USED: c_int = 1;

const SHARED_GRANUL: c_int = 16384; /* granularity of shared_malloc allocation = phys page size, system dependent */

/* checkpoints in shared memory segments - might be omitted */

const SHARED_ID_0: c_char = b'J' as c_char; /* first byte of identifier in BLKHEAD */
const SHARED_ID_1: c_char = b'B' as c_char; /* second byte of identifier in BLKHEAD */

const BLOCK_REG: c_char = 0; /* value for tflag member of BLKHEAD */
const BLOCK_SHARED: c_char = 1; /* value for tflag member of BLKHEAD */

/* generic error codes */

const SHARED_OK: c_int = 0;

const SHARED_ERR_MIN_IDX: c_int = SHARED_BADARG;
const SHARED_ERR_MAX_IDX: c_int = SHARED_NORESIZE;

const DAL_SHM_FREE: c_int = 0;
const DAL_SHM_USED: c_int = 1;

const DAL_SHM_ID0: u8 = b'D';
const DAL_SHM_ID1: u8 = b'S';
const DAL_SHM_ID2: u8 = b'M';

const DAL_SHM_SEGHEAD_ID: c_int = 0x19630114;

/* data types */

/* BLKHEAD object is placed at the beginning of every memory segment (both
shared and regular) to allow automatic recognition of segments type */

#[repr(C, align(8))]
#[derive(Debug, Copy, Clone)]
struct BLKHEAD {
    ID: [c_char; 2], /* ID = 'JB', just as a checkpoint */
    tflag: c_char,   /* is it shared memory or regular one ? */
    handle: c_int,   /* this is not necessary, used only for non-resizeable objects via ptr */
}

type SHARED_P = *const c_void; /* generic type of shared memory pointer */

#[repr(C)]
struct SHARED_GTAB /* data type used in global table */ {
    sem: c_int,        /* access semaphore (1 field): process count */
    semkey: c_int,     /* key value used to generate semaphore handle */
    key: c_int,        /* key value used to generate shared memory handle (realloc changes it) */
    handle: c_int,     /* handle of shared memory segment */
    size: c_int,       /* size of shared memory segment */
    nprocdebug: c_int, /* attached proc counter, helps remove zombie segments */
    attr: c_char,      /* attributes of shared memory object */
}

#[derive(Default, Clone)]
#[repr(C)]
struct SHARED_LTAB /* data type used in local table */ {
    p: Option<*const BLKHEAD>, /* pointer to segment (may be null) */
    tcnt: c_int,               /* number of threads in this process attached to segment */
    lkcnt: c_int,              /* >=0 <- number of read locks, -1 - write lock */
    seekpos: c_long,           /* current pointer position, read/write/seek operations change it */
}

unsafe impl Send for SHARED_LTAB {}

/* system dependent definitions */

type flock_t = flock;

#[derive(Clone, Copy)]
#[repr(C)]
union semun {
    val: c_int,
    buf: *mut semid_ds,
    array: *mut c_ushort,
}

type DAL_SHM_SEGHEAD = DAL_SHM_SEGHEAD_STRUCT;

#[repr(C)]
struct DAL_SHM_SEGHEAD_STRUCT {
    ID: c_int,      /* ID for debugging */
    h: c_int,       /* handle of sh. mem */
    size: c_int,    /* size of data area */
    nodeidx: c_int, /* offset of root object (node struct typically) */
}

/*              S H A R E D   M E M O R Y   D R I V E R
                =======================================

                  by Jerzy.Borkowski@obs.unige.ch

09-Mar-98 : initial version 1.0 released
23-Mar-98 : shared_malloc now accepts new handle as an argument
23-Mar-98 : shmem://0, shmem://1, etc changed to shmem://h0, etc due to bug
            in url parser.
10-Apr-98 : code cleanup
13-May-99 : delayed initialization added, global table deleted on exit when
            no shmem segments remain, and last process terminates
*/

static mut SHARED_MAXSEG: c_int = 0; /* max number of shared memory blocks */
static mut SHARED_RANGE: c_int = 0; /* max number of tried entries */
static mut SHARED_FD: c_int = SHARED_INVALID; /* handle of global access lock file */
static mut SHARED_GT_H: c_int = SHARED_INVALID; /* handle of global table segment */

static mut SHARED_KBASE: c_int = 0; /* base for shared memory handles */
static mut SHARED_DEBUG: bool = false; /* simple debugging tool, set to 0 to disable messages */
static mut SHARED_CREATE_MODE: c_int = 0o666; /* permission flags for created objects */
static mut SHARED_INIT_CALLED: bool = false; /* flag whether shared_init() has been called, used for delayed init */

static mut SHARED_GT: &mut [SHARED_GTAB] = &mut []; /* global table pointer */
static mut SHARED_GT_PTR: *mut SHARED_GTAB = ptr::null_mut();
static mut SHARED_LT: Vec<SHARED_LTAB> = Vec::new(); /* local table pointer */

unsafe fn shared_clear_entry(idx: usize) -> c_int /* unconditionally clear entry */ {
    unsafe {
        if (idx < 0) || (idx >= SHARED_MAXSEG.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        SHARED_GT[idx].key = SHARED_INVALID; /* clear entries in global table */
        SHARED_GT[idx].handle = SHARED_INVALID;
        SHARED_GT[idx].sem = SHARED_INVALID;
        SHARED_GT[idx].semkey = SHARED_INVALID;
        SHARED_GT[idx].nprocdebug = 0;
        SHARED_GT[idx].size = 0;
        SHARED_GT[idx].attr = 0;

        SHARED_OK
    }
}

unsafe fn shared_destroy_entry(idx: usize) -> c_int /* unconditionally destroy sema & shseg and clear entry */
{
    unsafe {
        let mut r: c_int = 0;
        let mut r2: c_int = 0;
        let filler: semun = semun { val: 0 };

        if (idx < 0) || (idx >= SHARED_MAXSEG.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        r2 = SHARED_OK;
        r = SHARED_OK;

        if SHARED_INVALID != SHARED_GT[idx].sem {
            r = semctl(SHARED_GT[idx].sem, 0, IPC_RMID, filler); /* destroy semaphore */
        }
        if SHARED_INVALID != SHARED_GT[idx].handle {
            r2 = shmctl(SHARED_GT[idx].handle, IPC_RMID, std::ptr::null_mut()); /* destroy shared memory segment */
        }
        if SHARED_OK == r {
            r = r2; /* accumulate error code in r, free r2 */
        }
        r2 = shared_clear_entry(idx);
        if SHARED_OK == r { r2 } else { r }
    }
}

/// This must (should) be called during exit/abort
#[cfg_attr(not(test), unsafe(no_mangle))]
pub extern "C" fn shared_cleanup() {
    unsafe {
        let i: c_int = 0;
        let j: c_int = 0;
        let mut r: c_int = 0;
        let mut oktodelete: bool = false;
        let mut filelocked: bool = false;
        let mut segmentspresent: bool = false;

        let mut flk: flock_t = flock {
            l_type: 0,
            l_whence: 0,
            l_start: 0,
            l_len: 0,
            l_pid: 0,
        };

        let mut ds: MaybeUninit<shmid_ds> = MaybeUninit::uninit();

        if SHARED_DEBUG {
            print!("shared_cleanup:");
        }

        if !SHARED_LT.is_empty() {
            if SHARED_DEBUG {
                print!(" deleting segments:");
            }

            for i in 0..SHARED_MAXSEG {
                if 0 == SHARED_LT[i as usize].tcnt {
                    continue; /* we're not using this segment, skip this ... */
                }
                if -1 != SHARED_LT[i as usize].lkcnt {
                    continue; /* seg not R/W locked by us, skip this ... */
                }

                r = shared_destroy_entry(i.try_into().unwrap()); /* destroy unconditionally sema & segment */
                if SHARED_DEBUG {
                    if SHARED_OK == r {
                        print!(" [{i}]");
                    } else {
                        print!(" [error on {i} !!!!]");
                    }
                }
            }

            SHARED_LT = Vec::new(); /* free local table */
        }

        /* detach global index table */
        if !SHARED_GT.len() == 0 {
            oktodelete = false;
            filelocked = false;
            if SHARED_DEBUG {
                print!(" detaching globalsharedtable");
            }

            if SHARED_INVALID != SHARED_FD {
                flk.l_type = F_WRLCK as c_short; /* lock whole lock file */

                flk.l_whence = 0;
                flk.l_start = 0;
                flk.l_len = SHARED_MAXSEG as c_long;

                if -1 != fcntl(SHARED_FD, F_SETLK, &mut flk) {
                    filelocked = true; /* success, scan global table, to see if there are any segs */
                    segmentspresent = false; /* assume, there are no segs in the system */
                    for j in 0..SHARED_MAXSEG {
                        if SHARED_INVALID != SHARED_GT[j as usize].key {
                            segmentspresent = true; /* yes, there is at least one */
                            break;
                        }
                    }

                    if !segmentspresent {
                        /* if there are no segs ... */

                        if 0 == shmctl(SHARED_GT_H, IPC_STAT, ds.as_mut_ptr()) {
                            /* get number of processes attached to table */

                            let ds: shmid_ds = ds.assume_init();

                            if ds.shm_nattch <= 1 {
                                oktodelete = true; /* if only one (we), then it is safe (but see text 4 lines later) to unlink */
                            }
                        }
                    }
                }
            }

            let _ = shmdt(SHARED_GT.as_ptr() as *const _); /* detach global table */

            /* delete global table from system, if no shm seg present */
            if oktodelete {
                shmctl(SHARED_GT_H, IPC_RMID, ptr::null_mut()); /* there is a race condition here - time window between shmdt and shmctl */
                SHARED_GT_H = SHARED_INVALID;
            }

            SHARED_GT = &mut [];

            /* if we locked, we need to unlock */
            if filelocked {
                flk.l_type = F_UNLCK.try_into().unwrap();
                flk.l_whence = 0;
                flk.l_start = 0;
                flk.l_len = SHARED_MAXSEG.into();
                fcntl(SHARED_FD, F_SETLK, &flk);
            }
        }
        SHARED_GT_H = SHARED_INVALID;

        /* close lock file */
        if SHARED_INVALID != SHARED_FD {
            if SHARED_DEBUG {
                print!(" closing lockfile");
            }
            close(SHARED_FD);
            SHARED_FD = SHARED_INVALID;
        }

        SHARED_KBASE = 0;
        SHARED_MAXSEG = 0;
        SHARED_RANGE = 0;
        SHARED_INIT_CALLED = false;

        if SHARED_DEBUG {
            println!(" <<done>>");
        }
    }
}

/// Initialize shared memory stuff, you have to call this routine once
#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_init(debug_msgs: c_int) -> c_int {
    unsafe { shared_init_safer(debug_msgs) }
}

/// Initialize shared memory stuff, you have to call this routine once
pub unsafe fn shared_init_safer(debug_msgs: c_int) -> c_int {
    unsafe {
        let i: c_int = 0;
        let mut buf: [c_char; 1000] = [0; 1000];
        let mut p: *mut c_char = ptr::null_mut();
        let mut oldumask: mode_t = 0;

        SHARED_INIT_CALLED = true; /* tell everybody no need to call us for the 2nd time */
        SHARED_DEBUG = debug_msgs != 0; /* set required debug mode */

        if SHARED_DEBUG {
            print!("shared_init:");
        }

        SHARED_KBASE = 0; /* adapt to current env. settings */
        p = getenv(SHARED_ENV_KEYBASE.as_ptr());

        if !p.is_null() {
            SHARED_KBASE = atoi(p);
        }

        if 0 == SHARED_KBASE {
            SHARED_KBASE = SHARED_KEYBASE;
        }

        if SHARED_DEBUG {
            let skb = SHARED_KBASE;
            print!(" keybase={}", skb);
        }

        SHARED_MAXSEG = 0;

        (p = getenv(SHARED_ENV_MAXSEG.as_ptr()));

        if !p.is_null() {
            SHARED_MAXSEG = atoi(p);
        }

        if 0 == SHARED_MAXSEG {
            SHARED_MAXSEG = SHARED_MAXSEG_VAL;
        }

        if SHARED_DEBUG {
            let sms = SHARED_MAXSEG;
            print!(" maxseg={}", sms);
        }

        SHARED_RANGE = 3 * SHARED_MAXSEG;

        /* create rw locking file (this file is never deleted) */
        if SHARED_INVALID == SHARED_FD {
            if SHARED_DEBUG {
                print!(" lockfileinit=");
            }

            let skb = SHARED_KBASE;
            let sms = SHARED_MAXSEG;
            int_snprintf!(
                buf,
                1000,
                "{}.{}.{}",
                SHARED_FDNAME.to_str().unwrap(),
                skb,
                sms,
            );

            oldumask = umask(0);

            SHARED_FD = open(
                buf.as_ptr(),
                O_TRUNC | O_EXCL | O_CREAT | O_RDWR,
                SHARED_CREATE_MODE,
            );

            umask(oldumask);

            /* or just open rw locking file, in case it already exists */
            if SHARED_INVALID == SHARED_FD {
                SHARED_FD = open(buf.as_mut_ptr(), O_TRUNC | O_RDWR, SHARED_CREATE_MODE);

                if SHARED_INVALID == SHARED_FD {
                    return SHARED_NOFILE;
                }

                if SHARED_DEBUG {
                    print!("slave");
                }
            } else if SHARED_DEBUG {
                print!("master");
            }
        }

        /* global table not attached, try to create it in shared memory */
        if SHARED_INVALID == SHARED_GT_H {
            if SHARED_DEBUG {
                print!(" globalsharedtableinit=");
            }

            SHARED_GT_H = shmget(
                SHARED_KBASE,
                (SHARED_MAXSEG * size_of::<SHARED_GTAB>() as c_int)
                    .try_into()
                    .unwrap(),
                IPC_CREAT | IPC_EXCL | SHARED_CREATE_MODE,
            ); /* try open as a master */

            /* if failed, try to open as a slave */
            if SHARED_INVALID == SHARED_GT_H {
                let segment_size_bytes = (SHARED_MAXSEG * size_of::<SHARED_GTAB>() as c_int)
                    .try_into()
                    .unwrap();

                SHARED_GT_H = shmget(SHARED_KBASE, segment_size_bytes, SHARED_CREATE_MODE);

                if SHARED_INVALID == SHARED_GT_H {
                    return SHARED_IPCERR; /* means deleted ID residing in system, shared mem unusable ... */
                }

                SHARED_GT_PTR = shmat(SHARED_GT_H, ptr::null(), 0) as *mut SHARED_GTAB; /* attach segment */

                if (SHARED_INVALID as *mut SHARED_GTAB) == SHARED_GT_PTR {
                    return SHARED_IPCERR;
                }

                SHARED_GT = slice::from_raw_parts_mut(SHARED_GT_PTR, SHARED_MAXSEG as usize);

                if SHARED_DEBUG {
                    print!("slave");
                }
            } else {
                SHARED_GT_PTR = shmat(SHARED_GT_H, ptr::null(), 0) as *mut SHARED_GTAB; /* attach segment */

                if (SHARED_INVALID as *mut SHARED_GTAB) == SHARED_GT_PTR {
                    return SHARED_IPCERR;
                }

                SHARED_GT = slice::from_raw_parts_mut(SHARED_GT_PTR, SHARED_MAXSEG as usize);

                for i in 0..SHARED_MAXSEG {
                    shared_clear_entry(i.try_into().unwrap()); /* since we are master, init data */
                }

                if SHARED_DEBUG {
                    print!("master");
                }
            }
        }

        /* initialize local table */
        if SHARED_LT.is_empty() {
            if SHARED_DEBUG {
                print!(" localtableinit=");
            }

            let mut v = Vec::new();

            if v.try_reserve_exact(SHARED_MAXSEG as usize).is_err() {
                return SHARED_NOMEM; /* not enough memory to allocate local table */
            } else {
                v.resize(SHARED_MAXSEG as usize, SHARED_LTAB::default());
            }

            for i in 0..SHARED_MAXSEG as usize {
                v[i].p = None; /* not mapped */
                v[i].tcnt = 0; /* unused (or zero threads using this seg) */
                v[i].lkcnt = 0; /* segment is unlocked */
                v[i].seekpos = 0; /* r/w pointer at the beginning of file */
            }

            SHARED_LT = v;

            if SHARED_DEBUG {
                print!("ok");
            }
        }

        atexit(shared_cleanup); /* we want shared_cleanup to be called at exit or abort */

        if SHARED_DEBUG {
            println!(" <<done>>");
        }

        SHARED_OK
    }
}

/// try to recover dormant segments after applic crash
#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_recover(id: c_int) -> c_int {
    unsafe {
        let i: c_int = 0;
        let mut r: c_int = 0;
        let mut r2: c_int = 0;

        if SHARED_GT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        r = SHARED_OK;

        // let shared_lt = shared_lt.as_mut().unwrap();

        for i in 0..SHARED_MAXSEG {
            if -1 != id && i != id {
                continue;
            }

            if SHARED_LT[i as usize].tcnt != 0 {
                continue; /* somebody (we) is using it */
            }

            if SHARED_INVALID == SHARED_GT[i as usize].key {
                continue; /* unused slot */
            }

            if shared_mux(i.try_into().unwrap(), SHARED_NOWAIT | SHARED_RDWRITE) != 0 {
                continue; /* acquire exclusive access to segment, but do not wait */
            }

            r2 = shared_process_count(SHARED_GT[i as usize].sem);
            if (SHARED_GT[i as usize].nprocdebug > r2) || (0 == r2) {
                if SHARED_DEBUG {
                    print!(
                        "Bogus handle={} nproc={} sema={}:",
                        i, SHARED_GT[i as usize].nprocdebug, r2,
                    );
                }
                r = shared_destroy_entry(i.try_into().unwrap());
                if SHARED_DEBUG {
                    print!(
                        "{}",
                        if r != 0 {
                            "error couldn't clear handle"
                        } else {
                            "handle cleared"
                        },
                    );
                }
            }
            shared_demux(i.try_into().unwrap(), SHARED_RDWRITE);
        }
        r /* table full */
    }
}

/* API routines - mutexes and locking */

/// obtain exclusive access to specified segment
unsafe fn shared_mux(idx: usize, mode: c_int) -> c_int {
    unsafe {
        let mut flk: flock_t = flock {
            l_type: 0,
            l_whence: 0,
            l_start: 0,
            l_len: 0,
            l_pid: 0,
        };

        let mut r: c_int = 0;

        /* delayed initialization */
        if !SHARED_INIT_CALLED {
            r = shared_init(0);
            if SHARED_OK != r {
                return r;
            }
        }

        if SHARED_INVALID == SHARED_FD {
            return SHARED_NOTINIT;
        }

        if (idx < 0) || (idx >= SHARED_MAXSEG.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        flk.l_type = if (mode & SHARED_RDWRITE) != 0 {
            F_WRLCK.try_into().unwrap()
        } else {
            F_RDLCK.try_into().unwrap()
        };
        flk.l_whence = 0;
        flk.l_start = idx as c_long;
        flk.l_len = 1;

        if SHARED_DEBUG {
            print!(" [mux ({idx}): ");
        }

        if -1
            == fcntl(
                SHARED_FD,
                if (mode & SHARED_NOWAIT) != 0 {
                    F_SETLK
                } else {
                    F_SETLKW
                },
                &flk,
            )
        {
            let errno = std::io::Error::last_os_error().raw_os_error().unwrap_or(0);

            match errno {
                EAGAIN | EACCES => {
                    if SHARED_DEBUG {
                        print!("again]");
                    }
                    return SHARED_AGAIN;
                }
                _ => {
                    if SHARED_DEBUG {
                        print!("err]");
                    }
                    return SHARED_IPCERR;
                }
            }
        }

        if SHARED_DEBUG {
            print!("ok]");
        }

        SHARED_OK
    }
}

unsafe fn shared_demux(idx: usize, mode: c_int) -> c_int /* free exclusive access to specified segment */
{
    unsafe {
        let mut flk: flock_t = flock {
            l_type: 0,
            l_whence: 0,
            l_start: 0,
            l_len: 0,
            l_pid: 0,
        };

        if SHARED_INVALID == SHARED_FD {
            return SHARED_NOTINIT;
        }

        if (idx < 0) || (idx >= SHARED_MAXSEG.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        flk.l_type = F_UNLCK.try_into().unwrap();
        flk.l_whence = 0;
        flk.l_start = idx as c_long;
        flk.l_len = 1;

        if SHARED_DEBUG {
            print!(" [demux ({idx}): ");
        }

        if -1 == fcntl(SHARED_FD, F_SETLKW, &flk) {
            let errno = std::io::Error::last_os_error().raw_os_error().unwrap_or(0);

            match errno {
                EAGAIN | EACCES => {
                    if SHARED_DEBUG {
                        print!("again]");
                    }
                    return SHARED_AGAIN;
                }
                _ => {
                    if SHARED_DEBUG {
                        print!("err]");
                    }
                    return SHARED_IPCERR;
                }
            }
        }

        if SHARED_DEBUG {
            print!("mode={mode} ok]");
        }

        SHARED_OK
    }
}

unsafe fn shared_process_count(sem: c_int) -> c_int /* valid only for time of invocation */ {
    let su: semun = semun { val: 0 };

    unsafe { semctl(sem, 0, GETVAL, su) } /* su is unused here */
}

unsafe fn shared_delta_process(sem: c_int, delta: c_int) -> c_int /* change number of processes hanging on segment */
{
    unsafe {
        let mut sb = sembuf {
            sem_num: 0,
            sem_op: 0,
            sem_flg: 0,
        };

        if SHARED_INVALID == sem {
            return SHARED_BADARG; /* semaphore not attached */
        }

        sb.sem_num = 0;
        sb.sem_op = delta as c_short;
        sb.sem_flg = SEM_UNDO as c_short;

        if -1 == semop(sem, &mut sb, 1) {
            SHARED_IPCERR
        } else {
            SHARED_OK
        }
    }
}

unsafe fn shared_attach_process(sem: c_int) -> c_int {
    unsafe {
        if SHARED_DEBUG {
            print!(" [attach process]");
        }
        shared_delta_process(sem, 1)
    }
}

unsafe fn shared_detach_process(sem: c_int) -> c_int {
    unsafe {
        if SHARED_DEBUG {
            print!(" [detach process]");
        }
        shared_delta_process(sem, -1)
    }
}

/* API routines - hashing and searching */

unsafe fn shared_get_free_entry(newhandle: c_int) -> c_int /* get newhandle, or -1, entry is set rw locked */
{
    unsafe {
        if SHARED_GT.is_empty() {
            return -1; /* not initialized */
        }

        if SHARED_LT.is_empty() {
            return -1; /* not initialized */
        }

        if newhandle < 0 {
            return -1;
        }

        if newhandle >= SHARED_MAXSEG {
            return -1;
        }

        if !SHARED_LT.is_empty() && SHARED_LT[newhandle as usize].tcnt != 0 {
            return -1; /* somebody (we) is using it */
        }

        if shared_mux(
            newhandle.try_into().unwrap(),
            SHARED_NOWAIT | SHARED_RDWRITE,
        ) != 0
        {
            return -1; /* used by others */
        }

        if SHARED_INVALID == SHARED_GT[newhandle as usize].key {
            return newhandle; /* we have found free slot, lock it and return index */
        }

        shared_demux(newhandle.try_into().unwrap(), SHARED_RDWRITE);

        if SHARED_DEBUG {
            print!("[free_entry - ERROR - entry unusable]");
        }
        -1 /* table full */
    }
}

/// return hash value for malloc
unsafe fn shared_get_hash(size: c_int, idx: usize) -> c_int {
    unsafe {
        static mut COUNTER: c_int = 0;
        let mut hash: c_int = 0;

        hash = (COUNTER + size * idx as c_int) % SHARED_RANGE;
        COUNTER = (COUNTER + 1) % SHARED_RANGE;
        hash
    }
}

fn shared_adjust_size(size: usize) -> c_int {
    ((size + size_of::<BLKHEAD>()).div_ceil(SHARED_GRANUL as usize) * SHARED_GRANUL as usize)
        as c_int
}

/* API routines - core : malloc/realloc/free/attach/detach/lock/unlock */

/// return idx or SHARED_INVALID
#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_malloc(size: c_long, mode: c_int, newhandle: c_int) -> c_int {
    let mut h: c_int = 0;
    let mut i: c_int = 0;
    let mut r: c_int = 0;
    let mut idx: c_int = 0;
    let mut key: c_int = 0;
    let filler: semun = semun { val: 0 };

    unsafe {
        /* delayed initialization */
        if !SHARED_INIT_CALLED {
            r = shared_init(0);
            if SHARED_OK != r {
                return r;
            }
        }

        if SHARED_DEBUG {
            print!("malloc (size = {size}, mode = {mode}):");
        }

        if size < 0 {
            return SHARED_INVALID;
        }

        (idx = shared_get_free_entry(newhandle));
        if -1 == idx {
            return SHARED_INVALID;
        }

        if SHARED_DEBUG {
            print!(" idx={idx}");
        }

        i = 0;
        loop {
            /* table full, signal error & exit */
            if i >= SHARED_RANGE {
                shared_demux(idx.try_into().unwrap(), SHARED_RDWRITE);
                return SHARED_INVALID;
            }

            key = SHARED_KBASE
                + ((i + shared_get_hash(size.try_into().unwrap(), idx.try_into().unwrap()))
                    % SHARED_RANGE);

            if SHARED_DEBUG {
                print!(" key={key}");
            }

            h = shmget(
                key,
                shared_adjust_size(size.try_into().unwrap())
                    .try_into()
                    .unwrap(),
                IPC_CREAT | IPC_EXCL | SHARED_CREATE_MODE,
            );

            if SHARED_DEBUG {
                print!(" handle={h}");
            }

            if SHARED_INVALID == h {
                i += 1;
                continue; /* segment already accupied */
            }

            let bp = shmat(h, std::ptr::null(), 0) as *mut BLKHEAD; /* try attach */

            if SHARED_DEBUG {
                print!(" p={bp:p}");
            }

            /* cannot attach, delete segment, try with another key */
            if (SHARED_INVALID as *mut BLKHEAD) == bp {
                shmctl(h, IPC_RMID, std::ptr::null_mut());
                i += 1;
                continue;
            } /* now create semaphor counting number of processes attached */

            SHARED_GT[idx as usize].sem = semget(key, 1, IPC_CREAT | IPC_EXCL | SHARED_CREATE_MODE);
            if SHARED_INVALID == SHARED_GT[idx as usize].sem {
                shmdt(bp as *const c_void); /* cannot create segment, delete everything */
                shmctl(h, IPC_RMID, ptr::null_mut());
                i += 1;
                continue; /* try with another key */
            }

            if SHARED_DEBUG {
                print!(" sem={}", SHARED_GT[idx as usize].sem);
            }

            /* try attach process */
            if shared_attach_process(SHARED_GT[idx as usize].sem) != 0 {
                semctl(SHARED_GT[idx as usize].sem, 0, IPC_RMID, filler); /* destroy semaphore */
                shmdt(bp as *const c_void); /* detach shared mem segment */
                shmctl(h, IPC_RMID, std::ptr::null_mut()); /* destroy shared mem segment */
                i += 1;
                continue; /* try with another key */
            }

            (*bp).tflag = BLOCK_SHARED; /* fill in data in segment's header (this is really not necessary) */
            (*bp).ID[0] = SHARED_ID_0;
            (*bp).ID[1] = SHARED_ID_1;
            (*bp).handle = idx; /* used in yorick */

            if !SHARED_LT.is_empty() {
                if mode & SHARED_RESIZE != 0 {
                    if shmdt(bp as *const c_void) != 0 {
                        r = SHARED_IPCERR; /* if segment is resizable, then detach segment */
                    }
                    SHARED_LT[idx as usize].p = None;
                } else {
                    SHARED_LT[idx as usize].p = Some(bp);
                }

                SHARED_LT[idx as usize].tcnt = 1; /* one thread using segment */
                SHARED_LT[idx as usize].lkcnt = 0; /* no locks at the moment */
                SHARED_LT[idx as usize].seekpos = 0; /* r/w pointer positioned at beg of block */
                SHARED_GT[idx as usize].handle = h; /* fill in data in global table */
                SHARED_GT[idx as usize].size = size as c_int;
                SHARED_GT[idx as usize].attr = mode as c_char;
                SHARED_GT[idx as usize].semkey = key;
                SHARED_GT[idx as usize].key = key;
                SHARED_GT[idx as usize].nprocdebug = 0;
            }

            break;
        }

        shared_demux(idx.try_into().unwrap(), SHARED_RDWRITE); /* hope this will not fail */

        idx.try_into().unwrap()
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_attach(idx: usize) -> c_int {
    unsafe {
        let mut r: c_int = 0;
        let mut r2: c_int = 0;

        r = shared_mux(idx, SHARED_RDWRITE | SHARED_WAIT);
        if SHARED_OK != r {
            return r;
        }

        r = shared_map(idx);
        if SHARED_OK != r {
            shared_demux(idx, SHARED_RDWRITE);
            return r;
        }

        // Assume shared_lt is initialized and has the same length as shared_gt
        // let shared_lt = shared_lt.as_mut().unwrap();

        /* try attach process */
        if shared_attach_process(SHARED_GT[idx].sem) != 0 {
            shmdt(SHARED_LT[idx].p.unwrap() as *const c_void); /* cannot attach process, detach everything */
            SHARED_LT[idx].p = None;
            shared_demux(idx, SHARED_RDWRITE);
            return SHARED_BADARG;
        }

        SHARED_LT[idx].tcnt += 1; /* one more thread is using segment */

        /* if resizeable, detach and return special pointer */
        if SHARED_GT[idx].attr as c_int & SHARED_RESIZE != 0 {
            if shmdt(SHARED_LT[idx].p.unwrap() as *const c_void) != 0 {
                r = SHARED_IPCERR; /* if segment is resizable, then detach segment */
            }
            SHARED_LT[idx].p = None;
        }

        SHARED_LT[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
        r2 = shared_demux(idx, SHARED_RDWRITE);
        if r != 0 { r } else { r2 }
    }
}

unsafe fn shared_check_locked_index(idx: usize) -> c_int /* verify that given idx is valid */ {
    unsafe {
        let mut r: c_int = 0;

        /* delayed initialization */
        if !SHARED_INIT_CALLED {
            r = shared_init(0);
            if SHARED_OK != r {
                return r;
            }
        }

        if (idx < 0) || (idx >= SHARED_MAXSEG.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_deref().unwrap();

        if SHARED_LT[idx].p.is_none() {
            return SHARED_BADARG; /* NULL pointer, not attached ?? */
        }

        if 0 == SHARED_LT[idx].lkcnt {
            return SHARED_BADARG; /* not locked ?? */
        }

        if (SHARED_ID_0 != (*SHARED_LT[idx].p.unwrap()).ID[0])
            || (SHARED_ID_1 != (*SHARED_LT[idx].p.unwrap()).ID[1])
            || (BLOCK_SHARED != (*SHARED_LT[idx].p.unwrap()).tflag)
        {
            /* invalid data in segment */
            return SHARED_BADARG;
        }
        SHARED_OK
    }
}

unsafe fn shared_map(idx: usize) -> c_int /* map all tables for given idx, check for validity */ {
    unsafe {
        /* have to obtain excl. access before calling shared_map */

        if (idx < 0) || (idx >= SHARED_MAXSEG.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        if SHARED_INVALID == SHARED_GT[idx].key {
            return SHARED_BADARG;
        }

        let h: c_int = shmget(SHARED_GT[idx].key, 1, SHARED_CREATE_MODE);

        if SHARED_INVALID == h {
            return SHARED_BADARG;
        }

        let bp = shmat(h, std::ptr::null(), 0) as *mut BLKHEAD;

        if (SHARED_INVALID as *const BLKHEAD) == bp {
            return SHARED_BADARG;
        }
        if (SHARED_ID_0 != (*bp).ID[0])
            || (SHARED_ID_1 != (*bp).ID[1])
            || (BLOCK_SHARED != (*bp).tflag)
            || (h != SHARED_GT[idx].handle)
        {
            shmdt(bp as *const c_void); /* invalid segment, detach everything */
            return SHARED_BADARG;
        }

        /* check if sema is still there */
        if SHARED_GT[idx].sem != semget(SHARED_GT[idx].semkey, 1, SHARED_CREATE_MODE) {
            shmdt(bp as *const c_void); /* cannot attach semaphore, detach everything */
            return SHARED_BADARG;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        SHARED_LT[idx].p = Some(bp); /* store pointer to shmem data */
        SHARED_OK
    }
}

unsafe fn shared_validate(idx: usize, mode: c_int) -> c_int /* use intrnally inside crit.sect !!! */
{
    unsafe {
        let mut r: c_int = shared_mux(idx, mode);

        if SHARED_OK != r {
            return r; /* idx checked by shared_mux */
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if SHARED_LT[idx].p.is_none() {
            r = shared_map(idx);
            if SHARED_OK != r {
                shared_demux(idx, mode);
                return r;
            }
        }

        if (SHARED_ID_0 != (*SHARED_LT[idx].p.unwrap()).ID[0])
            || (SHARED_ID_1 != (*SHARED_LT[idx].p.unwrap()).ID[1])
            || (BLOCK_SHARED != (*SHARED_LT[idx].p.unwrap()).tflag)
        {
            shared_demux(idx, mode);
            return r;
        }

        SHARED_OK
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_realloc(idx: usize, newsize: c_int) -> SHARED_P /* realloc shared memory segment */
{
    unsafe {
        let mut h: c_int = 0;
        let mut key: c_int = 0;
        let mut i: c_int = 0;
        let mut r: c_int = 0;
        let mut transfersize: c_long = 0;

        let mut bp: *mut BLKHEAD = ptr::null_mut();

        r = SHARED_OK;
        if newsize < 0 {
            return ptr::null();
        }

        if shared_check_locked_index(idx) != 0 {
            return ptr::null();
        }

        if 0 == (SHARED_GT[idx].attr as c_int & SHARED_RESIZE) {
            return ptr::null();
        }

        if SHARED_LT.is_empty() {
            return ptr::null(); /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if -1 != SHARED_LT[idx].lkcnt {
            return ptr::null(); /* check for RW lock */
        }

        if shared_adjust_size(SHARED_GT[idx].size.try_into().unwrap())
            == shared_adjust_size(newsize.try_into().unwrap())
        {
            SHARED_GT[idx].size = newsize;

            return ((SHARED_LT[idx].p.unwrap()).add(1)) as SHARED_P;
        }

        loop {
            if i >= SHARED_RANGE {
                return ptr::null(); /* table full, signal error & exit */
            }
            key = SHARED_KBASE + ((i + shared_get_hash(newsize, idx)) % SHARED_RANGE);
            h = shmget(
                key,
                shared_adjust_size(newsize.try_into().unwrap())
                    .try_into()
                    .unwrap(),
                IPC_CREAT | IPC_EXCL | SHARED_CREATE_MODE,
            );
            if SHARED_INVALID == h {
                i += 1;
                continue; /* segment already accupied */
            }

            bp = shmat(h, std::ptr::null(), 0) as *mut BLKHEAD; /* try attach */

            /* cannot attach, delete segment, try with another key */
            if (SHARED_INVALID as *mut BLKHEAD) == bp {
                shmctl(h, IPC_RMID, std::ptr::null_mut());
                i += 1;
                continue;
            }

            *bp = *(SHARED_LT[idx].p.unwrap()); /* copy header, then data */

            transfersize = if newsize < SHARED_GT[idx].size {
                newsize.into()
            } else {
                SHARED_GT[idx].size.into()
            };

            if transfersize > 0 {
                memcpy(
                    bp.add(1) as *mut c_void,
                    (SHARED_LT[idx].p.unwrap()).add(1) as *const c_void,
                    transfersize.try_into().unwrap(),
                );
            }

            if shmdt(SHARED_LT[idx].p.unwrap() as *const c_void) != 0 {
                r = SHARED_IPCERR; /* try to detach old segment */
            }

            if shmctl(SHARED_GT[idx].handle, IPC_RMID, std::ptr::null_mut()) != 0 && SHARED_OK == r
            {
                r = SHARED_IPCERR; /* destroy old shared memory segment */
            }

            SHARED_GT[idx].size = newsize; /* signal new size */
            SHARED_GT[idx].handle = h; /* signal new handle */
            SHARED_GT[idx].key = key; /* signal new key */
            SHARED_LT[idx].p = Some(bp);
            break;
        }

        (bp.add(1)) as SHARED_P
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_free(idx: usize) -> c_int /* detach segment, if last process & !PERSIST, destroy segment */
{
    unsafe {
        let mut cnt: c_int = 0;
        let r: c_int = 0;
        let mut r2: c_int = 0;

        let mut r = shared_validate(idx, SHARED_RDWRITE | SHARED_WAIT);
        if SHARED_OK != r {
            return r;
        }

        r = shared_detach_process(SHARED_GT[idx].sem);

        /* update number of processes using segment */
        if SHARED_OK != r {
            shared_demux(idx, SHARED_RDWRITE);
            return r;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        SHARED_LT[idx].tcnt -= 1; /* update number of threads using segment */

        if SHARED_LT[idx].tcnt > 0 {
            return shared_demux(idx, SHARED_RDWRITE); /* if more threads are using segment we are done */
        }

        /* if, we are the last thread, try to detach segment */
        if shmdt(SHARED_LT[idx].p.unwrap() as *const c_void) != 0 {
            shared_demux(idx, SHARED_RDWRITE);
            return SHARED_IPCERR;
        }

        SHARED_LT[idx].p = None; /* clear entry in local table */
        SHARED_LT[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
        (cnt = shared_process_count(SHARED_GT[idx].sem));

        /* get number of processes hanging on segment */
        if -1 == cnt {
            shared_demux(idx, SHARED_RDWRITE);
            return SHARED_IPCERR;
        }

        if (0 == cnt) && (0 == (SHARED_GT[idx].attr as c_int & SHARED_PERSIST)) {
            r = shared_destroy_entry(idx); /* no procs on seg, destroy it */
        }

        r2 = shared_demux(idx, SHARED_RDWRITE);
        if r != 0 { r } else { r2 }
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_lock(idx: usize, mode: c_int) -> SHARED_P /* lock given segment for exclusive access */
{
    unsafe {
        let mut r: c_int = 0;

        if shared_mux(idx, mode) != 0 {
            return ptr::null(); /* idx checked by shared_mux */
        }

        if SHARED_LT.is_empty() {
            return ptr::null(); /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if 0 != SHARED_LT[idx].lkcnt {
            /* are we already locked ?? */
            r = shared_map(idx);
            if SHARED_OK != r {
                shared_demux(idx, mode);
                return ptr::null();
            }
        }
        if SHARED_LT[idx].p.is_none() {
            /* stupid pointer ?? */
            r = shared_map(idx);
            if SHARED_OK != r {
                shared_demux(idx, mode);
                return ptr::null();
            }
        }

        if (SHARED_ID_0 != (*SHARED_LT[idx].p.unwrap()).ID[0])
            || (SHARED_ID_1 != (*SHARED_LT[idx].p.unwrap()).ID[1])
            || (BLOCK_SHARED != (*SHARED_LT[idx].p.unwrap()).tflag)
        {
            shared_demux(idx, mode);
            return ptr::null();
        }

        if (mode & SHARED_RDWRITE) != 0 {
            SHARED_LT[idx].lkcnt = -1;

            SHARED_GT[idx].nprocdebug += 1;
        } else {
            SHARED_LT[idx].lkcnt += 1;
        }

        SHARED_LT[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
        (SHARED_LT[idx].p.unwrap().add(1)) as SHARED_P
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_unlock(idx: usize) -> c_int /* unlock given segment, assumes seg is locked !! */
{
    unsafe {
        let mut r: c_int = 0;
        let mut r2: c_int = 0;
        let mut mode: c_int = 0;
        r = shared_check_locked_index(idx);

        if SHARED_OK != r {
            return r;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if SHARED_LT[idx].lkcnt > 0 {
            SHARED_LT[idx].lkcnt -= 1; /* unlock read lock */
            mode = SHARED_RDONLY;
        } else {
            SHARED_LT[idx].lkcnt = 0; /* unlock write lock */
            SHARED_GT[idx].nprocdebug -= 1;
            mode = SHARED_RDWRITE;
        }

        if 0 == SHARED_LT[idx].lkcnt && SHARED_GT[idx].attr as c_int & SHARED_RESIZE != 0 {
            if shmdt(SHARED_LT[idx].p.unwrap() as *const c_void) != 0 {
                r = SHARED_IPCERR; /* segment is resizable, then detach segment */
            }
            SHARED_LT[idx].p = None; /* signal detachment in local table */
        }
        r2 = shared_demux(idx, mode); /* unlock segment, rest is only parameter checking */
        if r != 0 { r } else { r2 }
    }
}

/* API routines - support and info routines */

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_attr(idx: usize) -> c_int /* get the attributes of the shared memory segment */
{
    unsafe {
        if shared_check_locked_index(idx) != 0 {
            return SHARED_INVALID;
        }

        let r: c_int = SHARED_GT[idx].attr as c_int;

        r
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_set_attr(idx: usize, newattr: c_int) -> c_int /* get the attributes of the shared memory segment */
{
    unsafe {
        if shared_check_locked_index(idx) != 0 {
            return SHARED_INVALID;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if -1 != SHARED_LT[idx].lkcnt {
            return SHARED_INVALID; /* ADDED - check for RW lock */
        }

        let r: c_int = SHARED_GT[idx].attr as c_int;

        SHARED_GT[idx].attr = newattr as c_char;
        r
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_set_debug(mode: c_int) -> c_int /* set/reset debug mode */ {
    unsafe {
        let r: c_int = SHARED_DEBUG as c_int;

        SHARED_DEBUG = mode != 0;
        r
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_set_createmode(mode: c_int) -> c_int /* set/reset debug mode */ {
    unsafe {
        let r: c_int = SHARED_CREATE_MODE;

        SHARED_CREATE_MODE = mode;
        r
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_list(id: c_int) -> c_int {
    unsafe {
        let i: c_int = 0;
        let mut r: c_int = 0;

        if SHARED_GT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if SHARED_DEBUG {
            print!("shared_list:");
        }

        r = SHARED_OK;
        println!(" Idx    Key   Nproc   Size   Flags");
        println!("==============================================");
        for i in 0..SHARED_MAXSEG {
            if -1 != id && i != id {
                continue;
            }

            if SHARED_INVALID == SHARED_GT[i as usize].key {
                continue; /* unused slot */
            }

            /* acquire exclusive access to segment, but do not wait */
            match shared_mux(i.try_into().unwrap(), SHARED_NOWAIT | SHARED_RDONLY) {
                SHARED_AGAIN => {
                    print!(
                        "!{:3} {:08x} {:4}  {:8}",
                        i,
                        SHARED_GT[i as usize].key as c_ulong,
                        SHARED_GT[i as usize].nprocdebug,
                        SHARED_GT[i as usize].size
                    );
                    if SHARED_RESIZE & SHARED_GT[i as usize].attr as c_int != 0 {
                        print!(" RESIZABLE");
                    }
                    if SHARED_PERSIST & SHARED_GT[i as usize].attr as c_int != 0 {
                        print!(" PERSIST");
                    }
                    println!();
                }
                SHARED_OK => {
                    print!(
                        " {:3} {:08x} {:4}  {:8}",
                        i,
                        SHARED_GT[i as usize].key as c_ulong,
                        SHARED_GT[i as usize].nprocdebug,
                        SHARED_GT[i as usize].size
                    );
                    if SHARED_RESIZE & SHARED_GT[i as usize].attr as c_int != 0 {
                        print!(" RESIZABLE");
                    }
                    if SHARED_PERSIST & SHARED_GT[i as usize].attr as c_int != 0 {
                        print!(" PERSIST");
                    }
                    println!();
                    shared_demux(i.try_into().unwrap(), SHARED_RDONLY);
                }
                _ => continue,
            }
        }
        if SHARED_DEBUG {
            println!(" done");
        }
        r /* table full */
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_getaddr(id: c_int, address: &mut *mut c_char) -> c_int {
    unsafe {
        let mut i: c_int = 0;
        let mut segname: [c_char; 10] = [0; 10];

        if SHARED_GT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        strcpy_safe(&mut segname, cs!(c"h"));

        int_snprintf!(&mut segname[1..], 9, "{}", id);

        if smem_open(&mut segname, 0, &mut i) != 0 {
            return SHARED_BADARG;
        }

        *address = (((SHARED_LT[i as usize].p.unwrap().add(1)) as *const DAL_SHM_SEGHEAD).add(1))
            as *mut c_char;
        /*  smem_close(i); */
        SHARED_OK
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_uncond_delete(id: c_int) -> c_int {
    unsafe {
        let i: c_int = 0;
        let mut r: c_int = 0;

        if SHARED_GT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if SHARED_DEBUG {
            print!("shared_uncond_delete:");
        }

        r = SHARED_OK;
        for i in 0..SHARED_MAXSEG {
            if -1 != id && i != id {
                continue;
            }

            if shared_attach(i.try_into().unwrap()) != 0 {
                if -1 != id {
                    println!("no such handle");
                }
                continue;
            }

            print!("handle {i}:");

            if shared_lock(i.try_into().unwrap(), SHARED_RDWRITE | SHARED_NOWAIT).is_null() {
                println!(" cannot lock in RW mode, not deleted");
                continue;
            }

            if shared_set_attr(i.try_into().unwrap(), SHARED_RESIZE) >= SHARED_ERRBASE {
                print!(" cannot clear PERSIST attribute");
            }

            if shared_free(i.try_into().unwrap()) != 0 {
                println!(" delete failed");
            } else {
                println!(" deleted");
            }
        }

        if SHARED_DEBUG {
            println!(" done");
        }
        r /* table full */
    }
}

/************************* CFITSIO DRIVER FUNCTIONS ***************************/

pub(crate) fn smem_init() -> c_int {
    0
}

pub(crate) fn smem_shutdown() -> c_int {
    unsafe {
        if SHARED_INIT_CALLED {
            shared_cleanup();
        }
        0
    }
}

pub(crate) fn smem_setoptions(mut option: c_int) -> c_int {
    option = 0;
    0
}

pub(crate) fn smem_getoptions(options: &mut c_int) -> c_int {
    *options = 0;
    0
}

pub(crate) fn smem_getversion(version: &mut c_int) -> c_int {
    *version = 10;
    0
}

pub(crate) fn smem_open(filename: &mut [c_char], rwmode: c_int, driverhandle: &mut c_int) -> c_int {
    unsafe {
        let mut h: c_int = 0;
        let mut nitems: c_int = 0;

        if filename.is_empty() {
            return SHARED_NULPTR;
        }

        nitems = sscanf_d(filename, cs!(c"h%d"), &mut h);

        if 1 != nitems {
            return SHARED_BADARG;
        }

        let r: c_int = shared_attach(h.try_into().unwrap());

        if SHARED_OK != r {
            return r;
        }

        let sp: *mut DAL_SHM_SEGHEAD = shared_lock(
            h.try_into().unwrap(),
            if READWRITE == rwmode {
                SHARED_RDWRITE
            } else {
                SHARED_RDONLY
            },
        ) as *mut DAL_SHM_SEGHEAD;

        if sp.is_null() {
            shared_free(h.try_into().unwrap());
            return SHARED_BADARG;
        }

        if (h != (*sp).h) || (DAL_SHM_SEGHEAD_ID != (*sp).ID) {
            shared_unlock(h.try_into().unwrap());
            shared_free(h.try_into().unwrap());

            return SHARED_BADARG;
        }

        *driverhandle = h;
        0
    }
}

pub(crate) fn smem_create(
    filename: &mut [c_char; FLEN_FILENAME],
    driverhandle: &mut c_int,
) -> c_int {
    unsafe {
        let mut sz: c_int = 0;
        let mut nitems: c_int = 0;
        let mut h: c_int = 0;

        if filename.is_empty() {
            return SHARED_NULPTR; /* currently ignored */
        }

        nitems = sscanf_d(filename, cs!(c"h%d"), &mut h);
        if 1 != nitems {
            return SHARED_BADARG;
        }

        sz = BL!() + size_of::<DAL_SHM_SEGHEAD>() as c_int;

        h = shared_malloc(sz.into(), SHARED_RESIZE | SHARED_PERSIST, h);

        if SHARED_INVALID == h {
            return SHARED_NOMEM;
        }

        let sp: *mut DAL_SHM_SEGHEAD =
            shared_lock(h.try_into().unwrap(), SHARED_RDWRITE) as *mut DAL_SHM_SEGHEAD;
        if sp.is_null() {
            shared_free(h.try_into().unwrap());
            return SHARED_BADARG;
        }

        (*sp).ID = DAL_SHM_SEGHEAD_ID;
        (*sp).h = h;
        (*sp).size = sz;
        (*sp).nodeidx = -1;

        *driverhandle = h;

        0
    }
}

pub(crate) fn smem_close(driverhandle: c_int) -> c_int {
    let mut r: c_int = 0;

    r = unsafe { shared_unlock(driverhandle.try_into().unwrap()) };
    if SHARED_OK != r {
        return r;
    }
    unsafe { shared_free(driverhandle.try_into().unwrap()) }
}

pub(crate) fn smem_remove(filename: &[c_char]) -> c_int {
    unsafe {
        let mut h: c_int = 0;
        let mut nitems: c_int = 0;
        let mut r: c_int = 0;

        if filename.is_empty() {
            return SHARED_NULPTR;
        }

        nitems = sscanf_d(filename, cs!(c"h%d"), &mut h);

        if 1 != nitems {
            return SHARED_BADARG;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        /* are we locked ? */
        if 0 == shared_check_locked_index(h.try_into().unwrap()) {
            if -1 != SHARED_LT[h as usize].lkcnt
            /* are we locked RO ? */
            {
                r = shared_unlock(h.try_into().unwrap());
                if SHARED_OK != r {
                    return r; /* yes, so relock in RW */
                }
                if shared_lock(h.try_into().unwrap(), SHARED_RDWRITE).is_null() {
                    return SHARED_BADARG;
                }
            }
        } else {
            /* not locked */

            // WARNING: This is bad! We are just converting a immutable slice to a mutable
            // SAFETY: Absolutely none.
            let f = slice::from_raw_parts_mut(filename.as_ptr() as *mut c_char, filename.len());

            r = smem_open(f, READWRITE, &mut h);
            if SHARED_OK != r {
                return r; /* so open in RW mode */
            }
        }

        shared_set_attr(h.try_into().unwrap(), SHARED_RESIZE); /* delete PERSIST attribute */
        smem_close(h.try_into().unwrap()) /* detach segment (this will delete it) */
    }
}

pub(crate) fn smem_size(driverhandle: c_int, size: &mut usize) -> c_int {
    // Hack
    let size = Some(size);

    unsafe {
        match size {
            Some(sz) => {
                if shared_check_locked_index(driverhandle.try_into().unwrap()) != 0 {
                    return SHARED_INVALID;
                }

                *sz = (SHARED_GT[driverhandle as usize].size
                    - size_of::<DAL_SHM_SEGHEAD>() as c_int) as usize;
            }

            None => {
                return SHARED_NULPTR;
            }
        }

        0
    }
}

pub(crate) fn smem_flush(driverhandle: c_int) -> c_int {
    if unsafe { shared_check_locked_index(driverhandle as usize) } != 0 {
        return SHARED_INVALID;
    }
    0
}

pub(crate) fn smem_seek(driverhandle: c_int, offset: LONGLONG) -> c_int {
    unsafe {
        if offset < 0 {
            return SHARED_BADARG;
        }

        if shared_check_locked_index(driverhandle.try_into().unwrap()) != 0 {
            return SHARED_INVALID;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        SHARED_LT[driverhandle as usize].seekpos = offset;
        0
    }
}

pub(crate) fn smem_read(driverhandle: c_int, buffer: &mut [u8], nbytes: usize) -> c_int {
    unsafe {
        if buffer.is_empty() {
            return SHARED_NULPTR;
        }

        if shared_check_locked_index(driverhandle.try_into().unwrap()) != 0 {
            return SHARED_INVALID;
        }

        if nbytes < 0 {
            return SHARED_BADARG;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if (SHARED_LT[driverhandle as usize].seekpos + nbytes as c_long)
            > SHARED_GT[driverhandle as usize].size.into()
        {
            return SHARED_BADARG; /* read beyond EOF */
        }

        memcpy(
            buffer.as_mut_ptr() as *mut c_void,
            ((((SHARED_LT[driverhandle as usize].p.unwrap().add(1)) as *const DAL_SHM_SEGHEAD)
                .add(1)) as *const c_char)
                .add(SHARED_LT[driverhandle as usize].seekpos as usize)
                as *const c_void,
            nbytes,
        );

        SHARED_LT[driverhandle as usize].seekpos += nbytes as c_long;
        0
    }
}

pub(crate) fn smem_write(driverhandle: c_int, buffer: &[u8], nbytes: usize) -> c_int {
    unsafe {
        if buffer.is_empty() {
            return SHARED_NULPTR;
        }

        if shared_check_locked_index(driverhandle.try_into().unwrap()) != 0 {
            return SHARED_INVALID;
        }

        if SHARED_LT.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if -1 != SHARED_LT[driverhandle as usize].lkcnt {
            return SHARED_INVALID; /* are we locked RW ? */
        }

        if nbytes < 0 {
            return SHARED_BADARG;
        }

        if (SHARED_LT[driverhandle as usize].seekpos + nbytes as c_long) as c_ulong
            > (SHARED_GT[driverhandle as usize].size - size_of::<DAL_SHM_SEGHEAD>() as c_int)
                as c_ulong
        {
            /* need to realloc shmem */
            if shared_realloc(
                driverhandle.try_into().unwrap(),
                (SHARED_LT[driverhandle as usize].seekpos
                    + nbytes as c_long
                    + size_of::<DAL_SHM_SEGHEAD>() as c_long)
                    .try_into()
                    .unwrap(),
            )
            .is_null()
            {
                return SHARED_NOMEM;
            }
        }

        memcpy(
            ((((SHARED_LT[driverhandle as usize].p.unwrap().add(1)) as *const DAL_SHM_SEGHEAD)
                .add(1)) as *const c_char)
                .add(SHARED_LT[driverhandle as usize].seekpos as usize) as *mut c_void,
            buffer.as_ptr() as *mut c_void,
            nbytes,
        );

        SHARED_LT[driverhandle as usize].seekpos += nbytes as c_long;
        0
    }
}

// The whole C file is inside `#ifdef HAVE_SHMEM_SERVICES`.  In Rust the module
// `src/drvrsmem.rs` is gated `#[cfg(all(feature = "shared_mem", not(target_os =
// "windows")))]`, so this entire test module is gated the same way.  Without the
// feature this file compiles to nothing and is harmless.
//
// IMPORTANT: shared-memory segments are PROCESS-GLOBAL and persist for the life
// of the process.  The tests use fixed segment names h0..h14, so this test MUST
// be run single-threaded to avoid cross-test interference:
//
//     cargo test --features shared_mem --test test_drvrsmem -- --test-threads=1
//
// Each test deletes its own segment(s) at the end, mirroring the C exactly.
#[cfg(test)]
#[cfg(all(feature = "shared_mem", not(target_os = "windows")))]
mod tests {
    use super::*;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        BYTE_IMG, FLEN_COMMENT, FLEN_FILENAME, FLEN_VALUE, LONGLONG, READONLY, READWRITE,
        SHARED_BADARG, SHORT_IMG, TBYTE, TLONG, TSHORT, fitsfile,
    };
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{c_char, c_int, c_long};

    // Serializes the shmem tests w.r.t. EACH OTHER only; the rest of the suite
    // still runs in parallel. The shmem driver is process-global (fixed segment
    // names h0..h14, global static state), so two shmem tests running at once
    // collide. Recover on poison so one failing test doesn't cascade into the
    // rest. Equivalent to `serial_test`'s `#[serial]`, without adding a dep.
    static SHMEM_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    fn shmem_guard() -> std::sync::MutexGuard<'static, ()> {
        SHMEM_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Build a NUL-terminated fixed-size `[c_char; FLEN_FILENAME]` from a `&str`.
    fn cc_fixed(s: &str) -> [c_char; FLEN_FILENAME] {
        let mut buf = [0 as c_char; FLEN_FILENAME];
        for (i, b) in s.bytes().enumerate() {
            buf[i] = b as c_char;
        }
        buf
    }

    /// Read a NUL-terminated `[c_char]` into a Rust `String`.
    fn from_buf(buf: &[c_char]) -> String {
        let bytes: Vec<u8> = buf
            .iter()
            .take_while(|&&c| c != 0)
            .map(|&c| c as u8)
            .collect();
        String::from_utf8_lossy(&bytes).into_owned()
    }

    /// Create and populate the `shmem://h0` segment exactly as
    /// `test_create_shmem_file` does in the C source.
    ///
    /// In C these tests run sequentially from `main()` and share the single h0
    /// segment (create -> read -> keywords -> delete).  Rust runs each `#[test]`
    /// independently (and, even single-threaded, in alphabetical rather than
    /// source order), so the create/read/keywords/delete tests cannot rely on a
    /// previous test having left h0 in place.  Each h0-using test therefore sets
    /// up its own copy with this helper and tears it down at the end.
    fn setup_h0() {
        let mut status: c_int = 0;
        let naxes: [c_long; 2] = [10, 10];
        let mut data = [0i16; 100];
        for (i, d) in data.iter_mut().enumerate() {
            *d = i as i16;
        }

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h0"), &mut status);
        assert_eq!(status, 0, "setup ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 2, &naxes, &mut status);
        assert_eq!(status, 0, "setup ffphps failed");
        fits_write_img(
            f.as_deref_mut().unwrap(),
            TSHORT,
            1,
            100,
            cast_slice(&data),
            &mut status,
        );
        assert_eq!(status, 0, "setup ffppr failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "setup ffclos failed");
    }

    /// Delete the `shmem://h0` segment (open RW then ffdelt).
    fn cleanup_h0() {
        let mut status: c_int = 0;
        let mut f: Option<Box<fitsfile>> = None;
        fits_open_file(&mut f, &cc("shmem://h0"), READWRITE, &mut status);
        assert_eq!(status, 0, "cleanup ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "cleanup ffdelt failed");
    }

    #[test]
    fn test_create_shmem_file() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 2] = [10, 10];
        let mut data = [0i16; 100];
        for (i, d) in data.iter_mut().enumerate() {
            *d = i as i16;
        }

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h0"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 2, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_write_img(
            f.as_deref_mut().unwrap(),
            TSHORT,
            1,
            100,
            cast_slice(&data),
            &mut status,
        );
        assert_eq!(status, 0, "ffppr failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        // Tear down h0 so this test is self-contained (see setup_h0 comment).
        cleanup_h0();
    }

    #[test]
    fn test_open_and_read_shmem() {
        let _g = shmem_guard();
        setup_h0();

        let mut status: c_int = 0;
        let mut data = [0i16; 100];
        let mut anynull: c_int = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_open_file(&mut f, &cc("shmem://h0"), READONLY, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_read_img(
            f.as_deref_mut().unwrap(),
            TSHORT,
            1,
            100,
            None,
            cast_slice_mut(&mut data),
            Some(&mut anynull),
            &mut status,
        );
        assert_eq!(status, 0, "ffgpv failed");
        assert!(!(data[0] != 0));
        assert!(!(data[50] != 50));
        assert!(!(data[99] != 99));
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        cleanup_h0();
    }

    #[test]
    fn test_shmem_keywords() {
        let _g = shmem_guard();
        setup_h0();

        let mut status: c_int = 0;
        let mut strval = [0 as c_char; FLEN_VALUE];
        let mut comment = [0 as c_char; FLEN_COMMENT];
        let mut longval: c_long = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_open_file(&mut f, &cc("shmem://h0"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_write_key_str(
            f.as_deref_mut().unwrap(),
            &cc("TESTKEY"),
            &cc("testval"),
            Some(&cc("test comment")),
            &mut status,
        );
        assert_eq!(status, 0, "ffpkys failed");
        fits_write_key_lng(
            f.as_deref_mut().unwrap(),
            &cc("INTKEY"),
            12345,
            Some(&cc("integer key")),
            &mut status,
        );
        assert_eq!(status, 0, "ffpkyj failed");
        fits_read_key_str(
            f.as_deref_mut().unwrap(),
            &cc("TESTKEY"),
            &mut strval,
            Some(&mut comment),
            &mut status,
        );
        assert_eq!(status, 0, "ffgkys failed");
        assert!(!(from_buf(&strval) != "testval"));
        fits_read_key(
            f.as_deref_mut().unwrap(),
            crate::KeywordDatatypeMut::TLONG(&mut longval),
            &cc("INTKEY"),
            Some(&mut comment),
            &mut status,
        );
        assert_eq!(status, 0, "ffgky failed");
        assert!(!(longval != 12345));
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        // ensure TLONG datatype constant is referenced (mirrors C's TLONG usage)
        let _ = TLONG;

        cleanup_h0();
    }

    #[test]
    fn test_delete_shmem() {
        let _g = shmem_guard();
        setup_h0();

        let mut status: c_int = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_open_file(&mut f, &cc("shmem://h0"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_create_second_segment() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [5];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h1"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_smem_option_functions() {
        let _g = shmem_guard();
        let mut options: c_int = 0;
        let mut version: c_int = 0;

        assert!(!(smem_setoptions(42) != 0));
        assert!(!(smem_getoptions(&mut options) != 0));
        assert!(!(options != 0));
        // C also asserts smem_getoptions(NULL) == SHARED_NULPTR, but the Rust
        // signature takes `&mut c_int` (no NULL possible), so that case cannot
        // be exercised through this API.
        assert!(!(smem_getversion(&mut version) != 0));
        assert!(!(version != 10));
        // Likewise smem_getversion(NULL) == SHARED_NULPTR is not testable here.
    }

    #[test]
    fn test_shared_utility_functions() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [5];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h6"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        unsafe {
            assert!(!(shared_list(-1) != 0));
            assert!(!(shared_list(6) != 0));
            let mut addr: *mut c_char = std::ptr::null_mut();
            assert!(!(shared_getaddr(6, &mut addr) != 0));
            assert!(!(addr.is_null()));
            assert!(!(shared_recover(-1) != 0));
            assert!(!(shared_recover(6) != 0));
            assert!(!(shared_uncond_delete(6) != 0));
        }
    }

    #[test]
    fn test_read_beyond_eof() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [2];
        let mut data = [0u8; 1000];
        let mut anynull: c_int = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h3"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        fits_open_file(&mut f, &cc("shmem://h3"), READONLY, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        let rc = fits_read_img(
            f.as_deref_mut().unwrap(),
            TBYTE,
            1,
            1000,
            None,
            &mut data,
            Some(&mut anynull),
            &mut status,
        );
        assert!(!(rc == 0));

        status = 0;
        fits_close_file(f.take().unwrap(), &mut status);

        status = 0;
        fits_open_file(&mut f, &cc("shmem://h3"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_smem_read_beyond_eof() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [2];
        let mut data = [0u8; 1000];
        let mut handle: c_int = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h4"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        let mut name = cc("h4");
        assert!(!(smem_open(&mut name, READONLY, &mut handle) != 0));
        assert!(!(smem_seek(handle, 99999) != 0));
        assert!(!(smem_read(handle, &mut data, 1000) != SHARED_BADARG));
        assert!(!(smem_close(handle) != 0));

        fits_open_file(&mut f, &cc("shmem://h4"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_cleanup_locked_segment() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [3];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h2"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        smem_shutdown();
        // C deliberately leaks the open fitsfile here (cleanup deletes the
        // segment underneath it); we do the same by forgetting the handle so we
        // don't attempt to close a segment that has been torn down.
        std::mem::forget(f.take());
    }

    #[test]
    fn test_cleanup_with_debug() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [3];

        unsafe {
            shared_set_debug(1);
        }
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h5"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        assert!(!(smem_shutdown() != 0));
        // restore debug mode and avoid closing a torn-down segment
        unsafe {
            shared_set_debug(0);
        }
        std::mem::forget(f.take());
    }

    #[test]
    fn test_smem_size() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 2] = [10, 10];
        let mut handle: c_int = 0;
        let mut size: usize = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h7"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 2, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        let mut name = cc("h7");
        assert!(!(smem_open(&mut name, READONLY, &mut handle) != 0));
        assert!(!(smem_size(handle, &mut size) != 0));
        assert!(!((size as LONGLONG) < 2880));
        assert!(!(smem_close(handle) != 0));

        fits_open_file(&mut f, &cc("shmem://h7"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_smem_flush() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [5];
        let mut handle: c_int = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h8"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        let mut name = cc("h8");
        assert!(!(smem_open(&mut name, READWRITE, &mut handle) != 0));
        assert!(!(smem_flush(handle) != 0));
        assert!(!(smem_close(handle) != 0));

        fits_open_file(&mut f, &cc("shmem://h8"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_smem_write() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [10];
        let mut handle: c_int = 0;
        let data: [u8; 4] = [0xFF, 0xFE, 0xFD, 0xFC];
        let mut readbuf = [0u8; 4];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h9"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        let mut name = cc("h9");
        assert!(!(smem_open(&mut name, READWRITE, &mut handle) != 0));
        assert!(!(smem_seek(handle, 2880) != 0));
        assert!(!(smem_write(handle, &data, 4) != 0));
        assert!(!(smem_seek(handle, 2880) != 0));
        assert!(!(smem_read(handle, &mut readbuf, 4) != 0));
        assert!(!(readbuf[0] != 0xFF));
        assert!(!(readbuf[1] != 0xFE));
        assert!(!(smem_close(handle) != 0));

        fits_open_file(&mut f, &cc("shmem://h9"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_smem_remove() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [3];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h10"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        // smem_remove should fail on open segment or work on closed
        assert!(!(smem_remove(&cc("h10")) != 0));
    }

    #[test]
    fn test_shared_set_createmode() {
        let _g = shmem_guard();
        unsafe {
            let oldmode = shared_set_createmode(0o666);
            shared_set_createmode(oldmode);
        }
    }

    #[test]
    fn test_smem_open_nonexistent() {
        let _g = shmem_guard();
        let mut handle: c_int = 0;
        // Opening non-existent segment should fail
        let mut name = cc("nonexistent_segment_xyz");
        assert!(!(smem_open(&mut name, READONLY, &mut handle) == 0));
    }

    #[test]
    fn test_smem_create_and_delete() {
        let _g = shmem_guard();
        let mut handle: c_int = 0;
        let mut name = cc_fixed("h11");
        assert!(!(smem_create(&mut name, &mut handle) != 0));
        assert!(!(smem_close(handle) != 0));
        let mut openname = cc("h11");
        assert!(!(smem_open(&mut openname, READWRITE, &mut handle) != 0));
        assert!(!(smem_close(handle) != 0));
        assert!(!(smem_remove(&cc("h11")) != 0));
    }

    #[test]
    fn test_shared_attr() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [10];
        let mut handle: c_int = 0;

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h13"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        fits_close_file(f.take().unwrap(), &mut status);
        assert_eq!(status, 0, "ffclos failed");

        let mut name = cc("h13");
        assert!(!(smem_open(&mut name, READONLY, &mut handle) != 0));
        let _attr = unsafe { shared_attr(handle as usize) };
        assert!(!(smem_close(handle) != 0));

        fits_open_file(&mut f, &cc("shmem://h13"), READWRITE, &mut status);
        assert_eq!(status, 0, "ffopen failed");
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }

    #[test]
    fn test_list_with_segments() {
        let _g = shmem_guard();
        let mut status: c_int = 0;
        let naxes: [c_long; 1] = [5];

        // Create a segment and list it
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("shmem://h14"), &mut status);
        assert_eq!(status, 0, "ffinit failed");
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
        assert_eq!(status, 0, "ffphps failed");
        unsafe {
            assert!(!(shared_list(-1) != 0));
        }
        fits_delete_file(&mut f, &mut status);
        assert_eq!(status, 0, "ffdelt failed");
    }
}
