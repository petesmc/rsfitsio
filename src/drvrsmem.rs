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

const SHARED_MAXSEG: c_int = 16; /* maximum number of shared memory blocks */

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

const SHARED_ID_0: u8 = b'J'; /* first byte of identifier in BLKHEAD */
const SHARED_ID_1: u8 = b'B'; /* second byte of identifier in BLKHEAD */

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

static mut shared_maxseg: c_int = 0; /* max number of shared memory blocks */
static mut shared_range: c_int = 0; /* max number of tried entries */
static mut shared_fd: c_int = SHARED_INVALID; /* handle of global access lock file */
static mut shared_gt_h: c_int = SHARED_INVALID; /* handle of global table segment */

static mut shared_kbase: c_int = 0; /* base for shared memory handles */
static mut shared_debug: bool = false; /* simple debugging tool, set to 0 to disable messages */
static mut shared_create_mode: c_int = 0o666; /* permission flags for created objects */
static mut shared_init_called: bool = false; /* flag whether shared_init() has been called, used for delayed init */

static mut shared_gt: &mut [SHARED_GTAB] = &mut []; /* global table pointer */
static mut shared_gt_ptr: *mut SHARED_GTAB = ptr::null_mut();
static mut shared_lt: Vec<SHARED_LTAB> = Vec::new(); /* local table pointer */

unsafe fn shared_clear_entry(idx: usize) -> c_int /* unconditionally clear entry */ {
    unsafe {
        if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        shared_gt[idx].key = SHARED_INVALID; /* clear entries in global table */
        shared_gt[idx].handle = SHARED_INVALID;
        shared_gt[idx].sem = SHARED_INVALID;
        shared_gt[idx].semkey = SHARED_INVALID;
        shared_gt[idx].nprocdebug = 0;
        shared_gt[idx].size = 0;
        shared_gt[idx].attr = 0;

        SHARED_OK
    }
}

unsafe fn shared_destroy_entry(idx: usize) -> c_int /* unconditionally destroy sema & shseg and clear entry */
{
    unsafe {
        let mut r: c_int = 0;
        let mut r2: c_int = 0;
        let filler: semun = semun { val: 0 };

        if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        r2 = SHARED_OK;
        r = SHARED_OK;

        if SHARED_INVALID != shared_gt[idx].sem {
            r = semctl(shared_gt[idx].sem, 0, IPC_RMID, filler); /* destroy semaphore */
        }
        if SHARED_INVALID != shared_gt[idx].handle {
            r2 = shmctl(shared_gt[idx].handle, IPC_RMID, std::ptr::null_mut()); /* destroy shared memory segment */
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

        if shared_debug {
            print!("shared_cleanup:");
        }

        if !shared_lt.is_empty() {
            if shared_debug {
                print!(" deleting segments:");
            }

            for i in 0..shared_maxseg {
                if 0 == shared_lt[i as usize].tcnt {
                    continue; /* we're not using this segment, skip this ... */
                }
                if -1 != shared_lt[i as usize].lkcnt {
                    continue; /* seg not R/W locked by us, skip this ... */
                }

                r = shared_destroy_entry(i.try_into().unwrap()); /* destroy unconditionally sema & segment */
                if shared_debug {
                    if SHARED_OK == r {
                        print!(" [{i}]");
                    } else {
                        print!(" [error on {i} !!!!]");
                    }
                }
            }

            shared_lt = Vec::new(); /* free local table */
        }

        /* detach global index table */
        if !shared_gt.len() == 0 {
            oktodelete = false;
            filelocked = false;
            if shared_debug {
                print!(" detaching globalsharedtable");
            }

            if SHARED_INVALID != shared_fd {
                flk.l_type = F_WRLCK as c_short; /* lock whole lock file */

                flk.l_whence = 0;
                flk.l_start = 0;
                flk.l_len = shared_maxseg as c_long;

                if -1 != fcntl(shared_fd, F_SETLK, &mut flk) {
                    filelocked = true; /* success, scan global table, to see if there are any segs */
                    segmentspresent = false; /* assume, there are no segs in the system */
                    for j in 0..shared_maxseg {
                        if SHARED_INVALID != shared_gt[j as usize].key {
                            segmentspresent = true; /* yes, there is at least one */
                            break;
                        }
                    }

                    if !segmentspresent {
                        /* if there are no segs ... */

                        if 0 == shmctl(shared_gt_h, IPC_STAT, ds.as_mut_ptr()) {
                            /* get number of processes attached to table */

                            let ds: shmid_ds = ds.assume_init();

                            if ds.shm_nattch <= 1 {
                                oktodelete = true; /* if only one (we), then it is safe (but see text 4 lines later) to unlink */
                            }
                        }
                    }
                }
            }

            unsafe { shmdt(shared_gt.as_ptr() as *const _) }; /* detach global table */

            /* delete global table from system, if no shm seg present */
            if oktodelete {
                shmctl(shared_gt_h, IPC_RMID, ptr::null_mut()); /* there is a race condition here - time window between shmdt and shmctl */
                shared_gt_h = SHARED_INVALID;
            }

            shared_gt = &mut [];

            /* if we locked, we need to unlock */
            if filelocked {
                flk.l_type = F_UNLCK.try_into().unwrap();
                flk.l_whence = 0;
                flk.l_start = 0;
                flk.l_len = shared_maxseg.into();
                fcntl(shared_fd, F_SETLK, &flk);
            }
        }
        shared_gt_h = SHARED_INVALID;

        /* close lock file */
        if SHARED_INVALID != shared_fd {
            if shared_debug {
                print!(" closing lockfile");
            }
            close(shared_fd);
            shared_fd = SHARED_INVALID;
        }

        shared_kbase = 0;
        shared_maxseg = 0;
        shared_range = 0;
        shared_init_called = false;

        if shared_debug {
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

        shared_init_called = true; /* tell everybody no need to call us for the 2nd time */
        shared_debug = debug_msgs != 0; /* set required debug mode */

        if shared_debug {
            print!("shared_init:");
        }

        shared_kbase = 0; /* adapt to current env. settings */
        p = getenv(SHARED_ENV_KEYBASE.as_ptr());

        if !p.is_null() {
            shared_kbase = atoi(p);
        }

        if 0 == shared_kbase {
            shared_kbase = SHARED_KEYBASE;
        }

        if shared_debug {
            let skb = shared_kbase;
            print!(" keybase={}", skb);
        }

        shared_maxseg = 0;

        (p = getenv(SHARED_ENV_MAXSEG.as_ptr()));

        if !p.is_null() {
            shared_maxseg = atoi(p);
        }

        if 0 == shared_maxseg {
            shared_maxseg = SHARED_MAXSEG;
        }

        if shared_debug {
            let sms = shared_maxseg;
            print!(" maxseg={}", sms);
        }

        shared_range = 3 * shared_maxseg;

        /* create rw locking file (this file is never deleted) */
        if SHARED_INVALID == shared_fd {
            if shared_debug {
                print!(" lockfileinit=");
            }

            let skb = shared_kbase;
            let sms = shared_maxseg;
            int_snprintf!(
                buf,
                1000,
                "{}.{}.{}",
                SHARED_FDNAME.to_str().unwrap(),
                skb,
                sms,
            );

            oldumask = umask(0);

            shared_fd = open(
                buf.as_ptr(),
                O_TRUNC | O_EXCL | O_CREAT | O_RDWR,
                shared_create_mode,
            );

            umask(oldumask);

            /* or just open rw locking file, in case it already exists */
            if SHARED_INVALID == shared_fd {
                shared_fd = open(buf.as_mut_ptr(), O_TRUNC | O_RDWR, shared_create_mode);

                if SHARED_INVALID == shared_fd {
                    return SHARED_NOFILE;
                }

                if shared_debug {
                    print!("slave");
                }
            } else if shared_debug {
                print!("master");
            }
        }

        /* global table not attached, try to create it in shared memory */
        if SHARED_INVALID == shared_gt_h {
            if shared_debug {
                print!(" globalsharedtableinit=");
            }

            shared_gt_h = shmget(
                shared_kbase,
                (shared_maxseg * size_of::<SHARED_GTAB>() as c_int)
                    .try_into()
                    .unwrap(),
                IPC_CREAT | IPC_EXCL | shared_create_mode,
            ); /* try open as a master */

            /* if failed, try to open as a slave */
            if SHARED_INVALID == shared_gt_h {
                let segment_size_bytes = (shared_maxseg * size_of::<SHARED_GTAB>() as c_int)
                    .try_into()
                    .unwrap();

                shared_gt_h = shmget(shared_kbase, segment_size_bytes, shared_create_mode);

                if SHARED_INVALID == shared_gt_h {
                    return SHARED_IPCERR; /* means deleted ID residing in system, shared mem unusable ... */
                }

                shared_gt_ptr = unsafe { shmat(shared_gt_h, ptr::null(), 0) } as *mut SHARED_GTAB; /* attach segment */

                if (SHARED_INVALID as *mut SHARED_GTAB) == shared_gt_ptr {
                    return SHARED_IPCERR;
                }

                shared_gt = slice::from_raw_parts_mut(shared_gt_ptr, shared_maxseg as usize);

                if shared_debug {
                    print!("slave");
                }
            } else {
                shared_gt_ptr = unsafe { shmat(shared_gt_h, ptr::null(), 0) } as *mut SHARED_GTAB; /* attach segment */

                if (SHARED_INVALID as *mut SHARED_GTAB) == shared_gt_ptr {
                    return SHARED_IPCERR;
                }

                shared_gt = slice::from_raw_parts_mut(shared_gt_ptr, shared_maxseg as usize);

                for i in 0..shared_maxseg {
                    shared_clear_entry(i.try_into().unwrap()); /* since we are master, init data */
                }

                if shared_debug {
                    print!("master");
                }
            }
        }

        /* initialize local table */
        if shared_lt.is_empty() {
            if shared_debug {
                print!(" localtableinit=");
            }

            let mut v = Vec::new();

            if v.try_reserve_exact(shared_maxseg as usize).is_err() {
                return SHARED_NOMEM; /* not enough memory to allocate local table */
            } else {
                v.resize(shared_maxseg as usize, SHARED_LTAB::default());
            }

            for i in 0..shared_maxseg as usize {
                v[i].p = None; /* not mapped */
                v[i].tcnt = 0; /* unused (or zero threads using this seg) */
                v[i].lkcnt = 0; /* segment is unlocked */
                v[i].seekpos = 0; /* r/w pointer at the beginning of file */
            }

            shared_lt = v;

            if shared_debug {
                print!("ok");
            }
        }

        atexit(shared_cleanup); /* we want shared_cleanup to be called at exit or abort */

        if shared_debug {
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

        if shared_gt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        r = SHARED_OK;

        // let shared_lt = shared_lt.as_mut().unwrap();

        for i in 0..shared_maxseg {
            if -1 != id && i != id {
                continue;
            }

            if shared_lt[i as usize].tcnt != 0 {
                continue; /* somebody (we) is using it */
            }

            if SHARED_INVALID == shared_gt[i as usize].key {
                continue; /* unused slot */
            }

            if shared_mux(i.try_into().unwrap(), SHARED_NOWAIT | SHARED_RDWRITE) != 0 {
                continue; /* acquire exclusive access to segment, but do not wait */
            }

            r2 = shared_process_count(shared_gt[i as usize].sem);
            if (shared_gt[i as usize].nprocdebug > r2) || (0 == r2) {
                if shared_debug {
                    print!(
                        "Bogus handle={} nproc={} sema={}:",
                        i, shared_gt[i as usize].nprocdebug, r2,
                    );
                }
                r = shared_destroy_entry(i.try_into().unwrap());
                if shared_debug {
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
        if !shared_init_called {
            r = shared_init(0);
            if SHARED_OK != r {
                return r;
            }
        }

        if SHARED_INVALID == shared_fd {
            return SHARED_NOTINIT;
        }

        if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
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

        if shared_debug {
            print!(" [mux ({idx}): ");
        }

        if -1
            == fcntl(
                shared_fd,
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
                    if shared_debug {
                        print!("again]");
                    }
                    return SHARED_AGAIN;
                }
                _ => {
                    if shared_debug {
                        print!("err]");
                    }
                    return SHARED_IPCERR;
                }
            }
        }

        if shared_debug {
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

        if SHARED_INVALID == shared_fd {
            return SHARED_NOTINIT;
        }

        if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        flk.l_type = F_UNLCK.try_into().unwrap();
        flk.l_whence = 0;
        flk.l_start = idx as c_long;
        flk.l_len = 1;

        if shared_debug {
            print!(" [demux ({idx}): ");
        }

        if -1 == fcntl(shared_fd, F_SETLKW, &flk) {
            let errno = std::io::Error::last_os_error().raw_os_error().unwrap_or(0);

            match errno {
                EAGAIN | EACCES => {
                    if shared_debug {
                        print!("again]");
                    }
                    return SHARED_AGAIN;
                }
                _ => {
                    if shared_debug {
                        print!("err]");
                    }
                    return SHARED_IPCERR;
                }
            }
        }

        if shared_debug {
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
        if shared_debug {
            print!(" [attach process]");
        }
        shared_delta_process(sem, 1)
    }
}

unsafe fn shared_detach_process(sem: c_int) -> c_int {
    unsafe {
        if shared_debug {
            print!(" [detach process]");
        }
        shared_delta_process(sem, -1)
    }
}

/* API routines - hashing and searching */

unsafe fn shared_get_free_entry(newhandle: c_int) -> c_int /* get newhandle, or -1, entry is set rw locked */
{
    unsafe {
        if shared_gt.is_empty() {
            return -1; /* not initialized */
        }

        if shared_lt.is_empty() {
            return -1; /* not initialized */
        }

        if newhandle < 0 {
            return -1;
        }

        if newhandle >= shared_maxseg {
            return -1;
        }

        if !shared_lt.is_empty() && shared_lt[newhandle as usize].tcnt != 0 {
            return -1; /* somebody (we) is using it */
        }

        if shared_mux(
            newhandle.try_into().unwrap(),
            SHARED_NOWAIT | SHARED_RDWRITE,
        ) != 0
        {
            return -1; /* used by others */
        }

        if SHARED_INVALID == shared_gt[newhandle as usize].key {
            return newhandle; /* we have found free slot, lock it and return index */
        }

        shared_demux(newhandle.try_into().unwrap(), SHARED_RDWRITE);

        if shared_debug {
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

        hash = (COUNTER + size * idx as c_int) % shared_range;
        COUNTER = (COUNTER + 1) % shared_range;
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
        if !shared_init_called {
            r = shared_init(0);
            if SHARED_OK != r {
                return r;
            }
        }

        if shared_debug {
            print!("malloc (size = {size}, mode = {mode}):");
        }

        if size < 0 {
            return SHARED_INVALID;
        }

        (idx = shared_get_free_entry(newhandle));
        if -1 == idx {
            return SHARED_INVALID;
        }

        if shared_debug {
            print!(" idx={idx}");
        }

        i = 0;
        loop {
            /* table full, signal error & exit */
            if i >= shared_range {
                shared_demux(idx.try_into().unwrap(), SHARED_RDWRITE);
                return SHARED_INVALID;
            }

            key = shared_kbase
                + ((i + shared_get_hash(size.try_into().unwrap(), idx.try_into().unwrap()))
                    % shared_range);

            if shared_debug {
                print!(" key={key}");
            }

            h = shmget(
                key,
                shared_adjust_size(size.try_into().unwrap())
                    .try_into()
                    .unwrap(),
                IPC_CREAT | IPC_EXCL | shared_create_mode,
            );

            if shared_debug {
                print!(" handle={h}");
            }

            if SHARED_INVALID == h {
                i += 1;
                continue; /* segment already accupied */
            }

            let bp = shmat(h, std::ptr::null(), 0) as *mut BLKHEAD; /* try attach */

            if shared_debug {
                print!(" p={bp:p}");
            }

            /* cannot attach, delete segment, try with another key */
            if (SHARED_INVALID as *mut BLKHEAD) == bp {
                shmctl(h, IPC_RMID, std::ptr::null_mut());
                i += 1;
                continue;
            } /* now create semaphor counting number of processes attached */

            shared_gt[idx as usize].sem = semget(key, 1, IPC_CREAT | IPC_EXCL | shared_create_mode);
            if SHARED_INVALID == shared_gt[idx as usize].sem {
                shmdt(bp as *const c_void); /* cannot create segment, delete everything */
                shmctl(h, IPC_RMID, ptr::null_mut());
                i += 1;
                continue; /* try with another key */
            }

            if shared_debug {
                print!(" sem={}", shared_gt[idx as usize].sem);
            }

            /* try attach process */
            if shared_attach_process(shared_gt[idx as usize].sem) != 0 {
                semctl(shared_gt[idx as usize].sem, 0, IPC_RMID, filler); /* destroy semaphore */
                shmdt(bp as *const c_void); /* detach shared mem segment */
                shmctl(h, IPC_RMID, std::ptr::null_mut()); /* destroy shared mem segment */
                i += 1;
                continue; /* try with another key */
            }

            (*bp).tflag = BLOCK_SHARED; /* fill in data in segment's header (this is really not necessary) */
            (*bp).ID[0] = SHARED_ID_0;
            (*bp).ID[1] = SHARED_ID_1;
            (*bp).handle = idx; /* used in yorick */

            if !shared_lt.is_empty() {
                if mode & SHARED_RESIZE != 0 {
                    if unsafe { shmdt(bp as *const c_void) != 0 } {
                        r = SHARED_IPCERR; /* if segment is resizable, then detach segment */
                    }
                    shared_lt[idx as usize].p = None;
                } else {
                    shared_lt[idx as usize].p = Some(bp);
                }

                shared_lt[idx as usize].tcnt = 1; /* one thread using segment */
                shared_lt[idx as usize].lkcnt = 0; /* no locks at the moment */
                shared_lt[idx as usize].seekpos = 0; /* r/w pointer positioned at beg of block */
                shared_gt[idx as usize].handle = h; /* fill in data in global table */
                shared_gt[idx as usize].size = size as c_int;
                shared_gt[idx as usize].attr = mode as c_char;
                shared_gt[idx as usize].semkey = key;
                shared_gt[idx as usize].key = key;
                shared_gt[idx as usize].nprocdebug = 0;
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
        if shared_attach_process(shared_gt[idx].sem) != 0 {
            shmdt(shared_lt[idx].p.unwrap() as *const c_void); /* cannot attach process, detach everything */
            shared_lt[idx].p = None;
            shared_demux(idx, SHARED_RDWRITE);
            return SHARED_BADARG;
        }

        shared_lt[idx].tcnt += 1; /* one more thread is using segment */

        /* if resizeable, detach and return special pointer */
        if shared_gt[idx].attr as c_int & SHARED_RESIZE != 0 {
            if shmdt(shared_lt[idx].p.unwrap() as *const c_void) != 0 {
                r = SHARED_IPCERR; /* if segment is resizable, then detach segment */
            }
            shared_lt[idx].p = None;
        }

        shared_lt[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
        r2 = shared_demux(idx, SHARED_RDWRITE);
        if r != 0 { r } else { r2 }
    }
}

unsafe fn shared_check_locked_index(idx: usize) -> c_int /* verify that given idx is valid */ {
    unsafe {
        let mut r: c_int = 0;

        /* delayed initialization */
        if !shared_init_called {
            r = shared_init(0);
            if SHARED_OK != r {
                return r;
            }
        }

        if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_deref().unwrap();

        if shared_lt[idx].p.is_none() {
            return SHARED_BADARG; /* NULL pointer, not attached ?? */
        }

        if 0 == shared_lt[idx].lkcnt {
            return SHARED_BADARG; /* not locked ?? */
        }

        if (SHARED_ID_0 != (*shared_lt[idx].p.unwrap()).ID[0])
            || (SHARED_ID_1 != (*shared_lt[idx].p.unwrap()).ID[1])
            || (BLOCK_SHARED != (*shared_lt[idx].p.unwrap()).tflag)
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

        if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
            return SHARED_BADARG;
        }

        if SHARED_INVALID == shared_gt[idx].key {
            return SHARED_BADARG;
        }

        let h: c_int = shmget(shared_gt[idx].key, 1, shared_create_mode);

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
            || (h != shared_gt[idx].handle)
        {
            shmdt(bp as *const c_void); /* invalid segment, detach everything */
            return SHARED_BADARG;
        }

        /* check if sema is still there */
        if shared_gt[idx].sem != semget(shared_gt[idx].semkey, 1, shared_create_mode) {
            shmdt(bp as *const c_void); /* cannot attach semaphore, detach everything */
            return SHARED_BADARG;
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        shared_lt[idx].p = Some(bp); /* store pointer to shmem data */
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

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if shared_lt[idx].p.is_none() {
            r = shared_map(idx);
            if SHARED_OK != r {
                shared_demux(idx, mode);
                return r;
            }
        }

        if (SHARED_ID_0 != (*shared_lt[idx].p.unwrap()).ID[0])
            || (SHARED_ID_1 != (*shared_lt[idx].p.unwrap()).ID[1])
            || (BLOCK_SHARED != (*shared_lt[idx].p.unwrap()).tflag)
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

        if 0 == (shared_gt[idx].attr as c_int & SHARED_RESIZE) {
            return ptr::null();
        }

        if shared_lt.is_empty() {
            return ptr::null(); /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if -1 != shared_lt[idx].lkcnt {
            return ptr::null(); /* check for RW lock */
        }

        if shared_adjust_size(shared_gt[idx].size.try_into().unwrap())
            == shared_adjust_size(newsize.try_into().unwrap())
        {
            shared_gt[idx].size = newsize;

            return ((shared_lt[idx].p.unwrap()).add(1)) as SHARED_P;
        }

        loop {
            if i >= shared_range {
                return ptr::null(); /* table full, signal error & exit */
            }
            key = shared_kbase + ((i + shared_get_hash(newsize, idx)) % shared_range);
            h = shmget(
                key,
                shared_adjust_size(newsize.try_into().unwrap())
                    .try_into()
                    .unwrap(),
                IPC_CREAT | IPC_EXCL | shared_create_mode,
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

            *bp = *(shared_lt[idx].p.unwrap()); /* copy header, then data */

            transfersize = if newsize < shared_gt[idx].size {
                newsize.into()
            } else {
                shared_gt[idx].size.into()
            };

            if transfersize > 0 {
                memcpy(
                    bp.add(1) as *mut c_void,
                    (shared_lt[idx].p.unwrap()).add(1) as *const c_void,
                    transfersize.try_into().unwrap(),
                );
            }

            if shmdt(shared_lt[idx].p.unwrap() as *const c_void) != 0 {
                r = SHARED_IPCERR; /* try to detach old segment */
            }

            if shmctl(shared_gt[idx].handle, IPC_RMID, std::ptr::null_mut()) != 0 && SHARED_OK == r
            {
                r = SHARED_IPCERR; /* destroy old shared memory segment */
            }

            shared_gt[idx].size = newsize; /* signal new size */
            shared_gt[idx].handle = h; /* signal new handle */
            shared_gt[idx].key = key; /* signal new key */
            shared_lt[idx].p = Some(bp);
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

        r = shared_detach_process(shared_gt[idx].sem);

        /* update number of processes using segment */
        if SHARED_OK != r {
            shared_demux(idx, SHARED_RDWRITE);
            return r;
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        shared_lt[idx].tcnt -= 1; /* update number of threads using segment */

        if shared_lt[idx].tcnt > 0 {
            return shared_demux(idx, SHARED_RDWRITE); /* if more threads are using segment we are done */
        }

        /* if, we are the last thread, try to detach segment */
        if shmdt(shared_lt[idx].p.unwrap() as *const c_void) != 0 {
            shared_demux(idx, SHARED_RDWRITE);
            return SHARED_IPCERR;
        }

        shared_lt[idx].p = None; /* clear entry in local table */
        shared_lt[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
        (cnt = shared_process_count(shared_gt[idx].sem));

        /* get number of processes hanging on segment */
        if -1 == cnt {
            shared_demux(idx, SHARED_RDWRITE);
            return SHARED_IPCERR;
        }

        if (0 == cnt) && (0 == (shared_gt[idx].attr as c_int & SHARED_PERSIST)) {
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

        if shared_lt.is_empty() {
            return ptr::null(); /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if 0 != shared_lt[idx].lkcnt {
            /* are we already locked ?? */
            r = shared_map(idx);
            if SHARED_OK != r {
                shared_demux(idx, mode);
                return ptr::null();
            }
        }
        if shared_lt[idx].p.is_none() {
            /* stupid pointer ?? */
            r = shared_map(idx);
            if SHARED_OK != r {
                shared_demux(idx, mode);
                return ptr::null();
            }
        }

        if (SHARED_ID_0 != (*shared_lt[idx].p.unwrap()).ID[0])
            || (SHARED_ID_1 != (*shared_lt[idx].p.unwrap()).ID[1])
            || (BLOCK_SHARED != (*shared_lt[idx].p.unwrap()).tflag)
        {
            shared_demux(idx, mode);
            return ptr::null();
        }

        if (mode & SHARED_RDWRITE) != 0 {
            shared_lt[idx].lkcnt = -1;

            shared_gt[idx].nprocdebug += 1;
        } else {
            shared_lt[idx].lkcnt += 1;
        }

        shared_lt[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
        (shared_lt[idx].p.unwrap().add(1)) as SHARED_P
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

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if shared_lt[idx].lkcnt > 0 {
            shared_lt[idx].lkcnt -= 1; /* unlock read lock */
            mode = SHARED_RDONLY;
        } else {
            shared_lt[idx].lkcnt = 0; /* unlock write lock */
            shared_gt[idx].nprocdebug -= 1;
            mode = SHARED_RDWRITE;
        }

        if 0 == shared_lt[idx].lkcnt && shared_gt[idx].attr as c_int & SHARED_RESIZE != 0 {
            if shmdt(shared_lt[idx].p.unwrap() as *const c_void) != 0 {
                r = SHARED_IPCERR; /* segment is resizable, then detach segment */
            }
            shared_lt[idx].p = None; /* signal detachment in local table */
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

        let r: c_int = shared_gt[idx].attr as c_int;

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

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if -1 != shared_lt[idx].lkcnt {
            return SHARED_INVALID; /* ADDED - check for RW lock */
        }

        let r: c_int = shared_gt[idx].attr as c_int;

        shared_gt[idx].attr = newattr as u8;
        r
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_set_debug(mode: c_int) -> c_int /* set/reset debug mode */ {
    unsafe {
        let r: c_int = shared_debug as c_int;

        shared_debug = mode != 0;
        r
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_set_createmode(mode: c_int) -> c_int /* set/reset debug mode */ {
    unsafe {
        let r: c_int = shared_create_mode;

        shared_create_mode = mode;
        r
    }
}

#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_list(id: c_int) -> c_int {
    unsafe {
        let i: c_int = 0;
        let mut r: c_int = 0;

        if shared_gt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if shared_debug {
            print!("shared_list:");
        }

        r = SHARED_OK;
        println!(" Idx    Key   Nproc   Size   Flags");
        println!("==============================================");
        for i in 0..shared_maxseg {
            if -1 != id && i != id {
                continue;
            }

            if SHARED_INVALID == shared_gt[i as usize].key {
                continue; /* unused slot */
            }

            /* acquire exclusive access to segment, but do not wait */
            match shared_mux(i.try_into().unwrap(), SHARED_NOWAIT | SHARED_RDONLY) {
                SHARED_AGAIN => {
                    print!(
                        "!{:3} {:08x} {:4}  {:8}",
                        i,
                        shared_gt[i as usize].key as c_ulong,
                        shared_gt[i as usize].nprocdebug,
                        shared_gt[i as usize].size
                    );
                    if SHARED_RESIZE & shared_gt[i as usize].attr as c_int != 0 {
                        print!(" RESIZABLE");
                    }
                    if SHARED_PERSIST & shared_gt[i as usize].attr as c_int != 0 {
                        print!(" PERSIST");
                    }
                    println!();
                }
                SHARED_OK => {
                    print!(
                        " {:3} {:08x} {:4}  {:8}",
                        i,
                        shared_gt[i as usize].key as c_ulong,
                        shared_gt[i as usize].nprocdebug,
                        shared_gt[i as usize].size
                    );
                    if SHARED_RESIZE & shared_gt[i as usize].attr as c_int != 0 {
                        print!(" RESIZABLE");
                    }
                    if SHARED_PERSIST & shared_gt[i as usize].attr as c_int != 0 {
                        print!(" PERSIST");
                    }
                    println!();
                    shared_demux(i.try_into().unwrap(), SHARED_RDONLY);
                }
                _ => continue,
            }
        }
        if shared_debug {
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

        if shared_gt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        strcpy_safe(&mut segname, cs!(c"h"));

        int_snprintf!(&mut segname[1..], 9, "{}", id);

        if smem_open(&mut segname, 0, &mut i) != 0 {
            return SHARED_BADARG;
        }

        *address = (((shared_lt[i as usize].p.unwrap().add(1)) as *const DAL_SHM_SEGHEAD).add(1))
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

        if shared_gt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        if shared_debug {
            print!("shared_uncond_delete:");
        }

        r = SHARED_OK;
        for i in 0..shared_maxseg {
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

        if shared_debug {
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
        if shared_init_called {
            unsafe { shared_cleanup() };
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

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        /* are we locked ? */
        if 0 == shared_check_locked_index(h.try_into().unwrap()) {
            if -1 != shared_lt[h as usize].lkcnt
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

                *sz = (shared_gt[driverhandle as usize].size
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

        if unsafe { shared_check_locked_index(driverhandle.try_into().unwrap()) } != 0 {
            return SHARED_INVALID;
        }

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        shared_lt[driverhandle as usize].seekpos = offset;
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

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if (shared_lt[driverhandle as usize].seekpos + nbytes as c_long)
            > shared_gt[driverhandle as usize].size.into()
        {
            return SHARED_BADARG; /* read beyond EOF */
        }

        memcpy(
            buffer.as_mut_ptr() as *mut c_void,
            ((((shared_lt[driverhandle as usize].p.unwrap().add(1)) as *const DAL_SHM_SEGHEAD)
                .add(1)) as *const c_char)
                .add(shared_lt[driverhandle as usize].seekpos as usize)
                as *const c_void,
            nbytes,
        );

        shared_lt[driverhandle as usize].seekpos += nbytes as c_long;
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

        if shared_lt.is_empty() {
            return SHARED_NOTINIT; /* not initialized */
        }

        // let shared_lt = shared_lt.as_mut().unwrap();

        if -1 != shared_lt[driverhandle as usize].lkcnt {
            return SHARED_INVALID; /* are we locked RW ? */
        }

        if nbytes < 0 {
            return SHARED_BADARG;
        }

        if (shared_lt[driverhandle as usize].seekpos + nbytes as c_long) as c_ulong
            > (shared_gt[driverhandle as usize].size - size_of::<DAL_SHM_SEGHEAD>() as c_int)
                as c_ulong
        {
            /* need to realloc shmem */
            if shared_realloc(
                driverhandle.try_into().unwrap(),
                (shared_lt[driverhandle as usize].seekpos
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
            ((((shared_lt[driverhandle as usize].p.unwrap().add(1)) as *const DAL_SHM_SEGHEAD)
                .add(1)) as *const c_char)
                .add(shared_lt[driverhandle as usize].seekpos as usize) as *mut c_void,
            buffer.as_ptr() as *mut c_void,
            nbytes,
        );

        shared_lt[driverhandle as usize].seekpos += nbytes as c_long;
        0
    }
}
