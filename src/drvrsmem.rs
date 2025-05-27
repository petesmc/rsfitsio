/* configuration parameters */

use std::{
    ffi::{CStr, c_void},
    ptr,
};

use libc::{
    EACCES, EAGAIN, F_RDLCK, F_SETLK, F_SETLKW, F_UNLCK, F_WRLCK, GETVAL, IPC_CREAT, IPC_EXCL,
    IPC_RMID, IPC_STAT, O_CREAT, O_EXCL, O_RDWR, O_TRUNC, SEM_UNDO, atexit, atoi, close, fcntl,
    flock, free, getenv, malloc, memcpy, mode_t, open, sembuf, semctl, semget, semid_ds, semop,
    shmat, shmctl, shmdt, shmget, umask,
};

use crate::c_types::{c_char, c_int, c_long, c_ushort};
use crate::{
    BL,
    fitsio::{
        LONGLONG, READWRITE, SHARED_AGAIN, SHARED_BADARG, SHARED_ERRBASE, SHARED_IPCERR,
        SHARED_NOFILE, SHARED_NOMEM, SHARED_NORESIZE, SHARED_NOTINIT, SHARED_NULPTR,
    },
    relibc::header::stdio::{snprintf, sscanf},
    wrappers::strcpy,
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

#[repr(C)]
#[derive(Debug, Copy, Clone)]
struct BLKHEADstruct {
    ID: [c_char; 2], /* ID = 'JB', just as a checkpoint */
    tflag: c_char,   /* is it shared memory or regular one ? */
    handle: c_int,   /* this is not necessary, used only for non-resizeable objects via ptr */
}

#[repr(C)]
union BLKHEAD {
    s: BLKHEADstruct,
    d: f64, /* for proper alignment on every machine */
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

#[repr(C)]
struct SHARED_LTAB<'a> /* data type used in local table */ {
    p: Option<&'a mut BLKHEAD>, /* pointer to segment (may be null) */
    tcnt: c_int,                /* number of threads in this process attached to segment */
    lkcnt: c_int,               /* >=0 <- number of read locks, -1 - write lock */
    seekpos: c_long,            /* current pointer position, read/write/seek operations change it */
}

/* system dependent definitions */

type flock_t = flock;

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
    nodeidx: usize, /* offset of root object (node struct typically) */
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

static mut shared_kbase: c_int = 0; /* base for shared memory handles */
static shared_maxseg: c_int = 0; /* max number of shared memory blocks */
static shared_range: c_int = 0; /* max number of tried entries */
static shared_fd: c_int = SHARED_INVALID; /* handle of global access lock file */
static shared_gt_h: c_int = SHARED_INVALID; /* handle of global table segment */
static mut shared_lt: Vec<SHARED_LTAB> = Vec::new(); /* local table pointer */
static mut shared_gt: Vec<SHARED_GTAB> = Vec::new(); /* global table pointer */
static shared_create_mode: c_int = 0666; /* permission flags for created objects */
static shared_debug: bool = true; /* simple debugging tool, set to 0 to disable messages */
static shared_init_called: bool = false; /* flag whether shared_init() has been called, used for delayed init */

unsafe fn shared_clear_entry(idx: usize) -> c_int /* unconditionally clear entry */ {
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

    return SHARED_OK;
}

unsafe fn shared_destroy_entry(idx: usize) -> c_int /* unconditionally destroy sema & shseg and clear entry */
{
    let mut r: c_int = 0;
    let mut r2: c_int = 0;
    let mut filler: semun = semun { val: 0 };

    if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
        return SHARED_BADARG;
    }
    r2 = SHARED_OK;
    r = SHARED_OK;
    filler.val = 0; /* this is to make cc happy (warning otherwise) */
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
    return if SHARED_OK == r { r2 } else { r };
}

unsafe fn shared_cleanup() /* this must (should) be called during exit/abort */
{
    let mut i: c_int = 0;
    let mut j: c_int = 0;
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
    let mut ds: semid_ds = semid_ds::default();

    if shared_debug {
        print!("shared_cleanup:");
    }
    if !shared_lt.is_null() {
        if shared_debug {
            print!(" deleting segments:");
        }
        for i in 0..shared_maxseg {
            if 0 == shared_lt[i].tcnt {
                continue; /* we're not using this segment, skip this ... */
            }
            if -1 != shared_lt[i].lkcnt {
                continue; /* seg not R/W locked by us, skip this ... */
            }

            r = shared_destroy_entry(i.try_into().unwrap()); /* destroy unconditionally sema & segment */
            if shared_debug {
                if SHARED_OK == r {
                    print!(" [{}]", i);
                } else {
                    print!(" [error on {} !!!!]", i);
                }
            }
        }
        free(shared_lt); /* free local table */
        shared_lt = ptr::null_mut();
    }
    if !shared_gt.is_null()
    /* detach global index table */
    {
        oktodelete = false;
        filelocked = false;
        if shared_debug {
            print!(" detaching globalsharedtable");
        }

        // WARNING: Original C code looks odd here
        if SHARED_INVALID != shared_fd {
            flk.l_type = F_WRLCK; /* lock whole lock file */
        }
        flk.l_whence = 0;
        flk.l_start = 0;
        flk.l_len = shared_maxseg;
        if -1 != fcntl(shared_fd, F_SETLK, &flk) {
            filelocked = true; /* success, scan global table, to see if there are any segs */
            segmentspresent = false; /* assume, there are no segs in the system */
            for j in 0..shared_maxseg {
                if SHARED_INVALID != shared_gt[j].key {
                    segmentspresent = true; /* yes, there is at least one */
                    break;
                }
            }
            if !segmentspresent {
                /* if there are no segs ... */
                if 0 == shmctl(shared_gt_h, IPC_STAT, &ds)
                /* get number of processes attached to table */
                {
                    if ds.shm_nattch <= 1 {
                        oktodelete = true; /* if only one (we), then it is safe (but see text 4 lines later) to unlink */
                    }
                }
            }
        }
        shmdt(shared_gt as *const _); /* detach global table */
        if oktodelete
        /* delete global table from system, if no shm seg present */
        {
            shmctl(shared_gt_h, IPC_RMID, ptr::null_mut()); /* there is a race condition here - time window between shmdt and shmctl */
            shared_gt_h = SHARED_INVALID;
        }
        shared_gt = ptr::null_mut();
        if filelocked
        /* if we locked, we need to unlock */
        {
            flk.l_type = F_UNLCK.try_into().unwrap();
            flk.l_whence = 0;
            flk.l_start = 0;
            flk.l_len = shared_maxseg.try_into().unwrap();
            fcntl(shared_fd, F_SETLK, &flk);
        }
    }
    shared_gt_h = SHARED_INVALID;

    if SHARED_INVALID != shared_fd
    /* close lock file */
    {
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
        print!(" <<done>>\n");
    }
    return;
}

unsafe fn shared_init(debug_msgs: c_int) -> c_int /* initialize shared memory stuff, you have to call this routine once */
{
    let mut i: c_int = 0;
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
        print!(" keybase={}", shared_kbase);
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
        print!(" maxseg={}", shared_maxseg);
    }

    shared_range = 3 * shared_maxseg;

    if SHARED_INVALID == shared_fd
    /* create rw locking file (this file is never deleted) */
    {
        if shared_debug {
            print!(" lockfileinit=");
        }
        snprintf(
            buf.as_mut_ptr(),
            1000,
            c"%s.%d.%d".as_ptr(),
            SHARED_FDNAME,
            shared_kbase,
            shared_maxseg,
        );
        oldumask = umask(0);

        shared_fd = open(buf, O_TRUNC | O_EXCL | O_CREAT | O_RDWR, shared_create_mode);
        umask(oldumask);
        if SHARED_INVALID == shared_fd
        /* or just open rw locking file, in case it already exists */
        {
            shared_fd = open(buf.as_mut_ptr(), O_TRUNC | O_RDWR, shared_create_mode);
            if SHARED_INVALID == shared_fd {
                return SHARED_NOFILE;
            }
            if shared_debug {
                print!("slave");
            }
        } else {
            if shared_debug {
                print!("master");
            }
        }
    }

    if SHARED_INVALID == shared_gt_h
    /* global table not attached, try to create it in shared memory */
    {
        if shared_debug {
            print!(" globalsharedtableinit=");
        }
        shared_gt_h = shmget(
            shared_kbase,
            (shared_maxseg * size_of::<SHARED_GTAB>())
                .try_into()
                .unwrap(),
            IPC_CREAT | IPC_EXCL | shared_create_mode,
        ); /* try open as a master */
        if SHARED_INVALID == shared_gt_h
        /* if failed, try to open as a slave */
        {
            shared_gt_h = shmget(
                shared_kbase,
                (shared_maxseg * size_of::<SHARED_GTAB>())
                    .try_into()
                    .unwrap(),
                shared_create_mode,
            );
            if SHARED_INVALID == shared_gt_h {
                return SHARED_IPCERR; /* means deleted ID residing in system, shared mem unusable ... */
            }
            shared_gt = shmat(shared_gt_h, ptr::null(), 0) as *mut SHARED_GTAB; /* attach segment */
            if (SHARED_INVALID as *mut SHARED_GTAB) == shared_gt {
                return SHARED_IPCERR;
            }
            if shared_debug {
                print!("slave");
            }
        } else {
            shared_gt = shmat(shared_gt_h, ptr::null(), 0) as *mut SHARED_GTAB; /* attach segment */
            if (SHARED_INVALID as *mut SHARED_GTAB) == shared_gt {
                return SHARED_IPCERR;
            }
            for i in 0..shared_maxseg {
                shared_clear_entry(i.try_into().unwrap()); /* since we are master, init data */
            }
            if shared_debug {
                print!("master");
            }
        }
    }

    /* initialize local table */
    if shared_lt.is_null() {
        if shared_debug {
            print!(" localtableinit=");
        }
        if (shared_lt = malloc(
            (shared_maxseg * size_of::<SHARED_LTAB>())
                .try_into()
                .unwrap(),
        ) as *const SHARED_LTAB)
            .is_null()
        {
            return SHARED_NOMEM;
        }
        for i in 0..shared_maxseg as usize {
            shared_lt[i].p = ptr::null_mut(); /* not mapped */
            shared_lt[i].tcnt = 0; /* unused (or zero threads using this seg) */
            shared_lt[i].lkcnt = 0; /* segment is unlocked */
            shared_lt[i].seekpos = 0; /* r/w pointer at the beginning of file */
        }
        if shared_debug {
            print!("ok");
        }
    }

    atexit(shared_cleanup); /* we want shared_cleanup to be called at exit or abort */

    if shared_debug {
        print!(" <<done>>\n");
    }
    return SHARED_OK;
}

unsafe fn shared_recover(id: c_int) -> c_int /* try to recover dormant segments after applic crash */
{
    let mut i: c_int = 0;
    let mut r: c_int = 0;
    let mut r2: c_int = 0;

    if shared_gt.is_null() {
        return SHARED_NOTINIT; /* not initialized */
    }
    if shared_lt.is_null() {
        return SHARED_NOTINIT; /* not initialized */
    }
    r = SHARED_OK;
    for i in 0..shared_maxseg {
        if -1 != id {
            if i != id {
                continue;
            }
        }
        if shared_lt[i].tcnt {
            continue; /* somebody (we) is using it */
        }
        if SHARED_INVALID == shared_gt[i].key {
            continue; /* unused slot */
        }
        if shared_mux(i.try_into().unwrap(), SHARED_NOWAIT | SHARED_RDWRITE) != 0 {
            continue; /* acquire exclusive access to segment, but do not wait */
        }
        r2 = shared_process_count(shared_gt[i].sem);
        if (shared_gt[i].nprocdebug > r2) || (0 == r2) {
            if shared_debug {
                print!(
                    "Bogus handle={} nproc={} sema={}:",
                    i, shared_gt[i].nprocdebug, r2,
                );
            }
            r = shared_destroy_entry(i.try_into().unwrap());
            if shared_debug {
                print!(
                    "{}",
                    if r {
                        "error couldn't clear handle"
                    } else {
                        "handle cleared"
                    },
                );
            }
        }
        shared_demux(i.try_into().unwrap(), SHARED_RDWRITE);
    }
    return r; /* table full */
}

/* API routines - mutexes and locking */

unsafe fn shared_mux(idx: usize, mode: c_int) -> c_int /* obtain exclusive access to specified segment */
{
    let mut flk: flock_t = flock {
        l_type: 0,
        l_whence: 0,
        l_start: 0,
        l_len: 0,
        l_pid: 0,
    };

    let mut r: c_int = 0;

    if !shared_init_called
    /* delayed initialization */
    {
        (r = shared_init(0));
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
    flk.l_start = idx;
    flk.l_len = 1;
    if shared_debug {
        print!(" [mux ({}): ", idx);
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
    return SHARED_OK;
}

unsafe fn shared_demux(idx: usize, mode: c_int) -> c_int /* free exclusive access to specified segment */
{
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
    flk.l_start = idx.into();
    flk.l_len = 1;
    if shared_debug {
        print!(" [demux ({}): ", idx);
    }
    if -1 == fcntl(shared_fd, F_SETLKW, &flk) {
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
        print!("mode={} ok]", mode);
    }
    return SHARED_OK;
}

unsafe fn shared_process_count(sem: c_int) -> c_int /* valid only for time of invocation */ {
    let su: semun = semun { val: 0 };

    return unsafe { semctl(sem, 0, GETVAL, su) }; /* su is unused here */
}

unsafe fn shared_delta_process(sem: c_int, delta: c_int) -> c_int /* change number of processes hanging on segment */
{
    let mut sb = sembuf {
        sem_num: 0,
        sem_op: 0,
        sem_flg: 0,
    };

    if SHARED_INVALID == sem {
        return SHARED_BADARG; /* semaphore not attached */
    }
    sb.sem_num = 0;
    sb.sem_op = delta;
    sb.sem_flg = SEM_UNDO;
    return if -1 == semop(sem, &sb, 1) {
        SHARED_IPCERR
    } else {
        SHARED_OK
    };
}

unsafe fn shared_attach_process(sem: c_int) -> c_int {
    if shared_debug {
        print!(" [attach process]");
    }
    return shared_delta_process(sem, 1);
}

unsafe fn shared_detach_process(sem: c_int) -> c_int {
    if shared_debug {
        print!(" [detach process]");
    }
    return shared_delta_process(sem, -1);
}

/* API routines - hashing and searching */

unsafe fn shared_get_free_entry(newhandle: c_int) -> c_int /* get newhandle, or -1, entry is set rw locked */
{
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
    if shared_lt[newhandle].tcnt {
        return -1; /* somebody (we) is using it */
    }
    if shared_mux(
        newhandle.try_into().unwrap(),
        SHARED_NOWAIT | SHARED_RDWRITE,
    ) != 0
    {
        return -1; /* used by others */
    }
    if SHARED_INVALID == shared_gt[newhandle].key {
        return newhandle; /* we have found free slot, lock it and return index */
    }
    shared_demux(newhandle.try_into().unwrap(), SHARED_RDWRITE);
    if shared_debug {
        print!("[free_entry - ERROR - entry unusable]");
    }
    return -1; /* table full */
}

unsafe fn shared_get_hash(size: c_int, idx: usize) -> c_int /* return hash value for malloc */ {
    static mut counter: c_int = 0;
    let mut hash: c_int = 0;

    hash = (counter + size * idx as c_int) % shared_range;
    counter = (counter + 1) % shared_range;
    return hash;
}

fn shared_adjust_size(size: c_int) -> c_int /* size must be >= 0 !!! */ {
    return ((size + size_of::<BLKHEAD>() as c_int + SHARED_GRANUL - 1) / SHARED_GRANUL)
        * SHARED_GRANUL;
}

/* API routines - core : malloc/realloc/free/attach/detach/lock/unlock */

unsafe fn shared_malloc(size: c_long, mode: c_int, newhandle: c_int) -> c_int /* return idx or SHARED_INVALID */
{
    let mut h: c_int = 0;
    let mut i: c_int = 0;
    let mut r: c_int = 0;
    let mut idx: usize = 0;
    let mut key: c_int = 0;
    let mut filler: semun = semun { val: 0 };
    let mut bp: BLKHEAD = BLKHEAD {
        s: BLKHEADstruct {
            ID: [0; 2],
            tflag: 0,
            handle: 0,
        },
    };

    if !shared_init_called
    /* delayed initialization */
    {
        (r = shared_init(0));
        if SHARED_OK != r {
            return r;
        }
    }
    if shared_debug {
        print!("malloc (size = {}, mode = {}):", size, mode);
    }
    if size < 0 {
        return SHARED_INVALID;
    }

    (idx = shared_get_free_entry(newhandle));
    if -1 == idx {
        return SHARED_INVALID;
    }
    if shared_debug {
        print!(" idx={}", idx);
    }
    i = 0;
    loop {
        if i >= shared_range
        /* table full, signal error & exit */
        {
            shared_demux(idx, SHARED_RDWRITE);
            return SHARED_INVALID;
        }
        key = shared_kbase + ((i + shared_get_hash(size.try_into().unwrap(), idx)) % shared_range);
        if shared_debug {
            print!(" key={}", key);
        }
        h = shmget(
            key,
            shared_adjust_size(size.try_into().unwrap())
                .try_into()
                .unwrap(),
            IPC_CREAT | IPC_EXCL | shared_create_mode,
        );
        if shared_debug {
            print!(" handle={}", h);
        }
        if SHARED_INVALID == h {
            i += 1;
            continue; /* segment already accupied */
        }
        bp = shmat(h, std::ptr::null(), 0) as *const BLKHEAD; /* try attach */
        if shared_debug {
            print!(" p={:p}", bp);
        }
        if (SHARED_INVALID as *const BLKHEAD) == bp
        /* cannot attach, delete segment, try with another key */
        {
            shmctl(h, IPC_RMID, std::ptr::null_mut());
            i += 1;
            continue;
        } /* now create semaphor counting number of processes attached */
        if SHARED_INVALID
            == (shared_gt[idx].sem = semget(key, 1, IPC_CREAT | IPC_EXCL | shared_create_mode))
        {
            shmdt(bp); /* cannot create segment, delete everything */
            shmctl(h, IPC_RMID, ptr::null_mut());
            i += 1;
            continue; /* try with another key */
        }
        if shared_debug {
            print!(" sem={}", shared_gt[idx].sem);
        }
        if shared_attach_process(shared_gt[idx].sem)
        /* try attach process */
        {
            semctl(shared_gt[idx].sem, 0, IPC_RMID, filler); /* destroy semaphore */
            shmdt(bp); /* detach shared mem segment */
            shmctl(h, IPC_RMID, std::ptr::null_mut()); /* destroy shared mem segment */
            i += 1;
            continue; /* try with another key */
        }
        bp.s.tflag = BLOCK_SHARED; /* fill in data in segment's header (this is really not necessary) */
        bp.s.ID[0] = SHARED_ID_0;
        bp.s.ID[1] = SHARED_ID_1;
        bp.s.handle = idx; /* used in yorick */
        if mode & SHARED_RESIZE {
            if shmdt(bp) {
                r = SHARED_IPCERR; /* if segment is resizable, then detach segment */
            }
            shared_lt[idx].p = ptr::null_mut();
        } else {
            shared_lt[idx].p = bp;
        }
        shared_lt[idx].tcnt = 1; /* one thread using segment */
        shared_lt[idx].lkcnt = 0; /* no locks at the moment */
        shared_lt[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
        shared_gt[idx].handle = h; /* fill in data in global table */
        shared_gt[idx].size = size;
        shared_gt[idx].attr = mode;
        shared_gt[idx].semkey = key;
        shared_gt[idx].key = key;
        shared_gt[idx].nprocdebug = 0;

        break;
    }
    shared_demux(idx, SHARED_RDWRITE); /* hope this will not fail */
    return idx.try_into().unwrap();
}

unsafe fn shared_attach(idx: usize) -> c_int {
    let mut r: c_int = 0;
    let mut r2: c_int = 0;

    (r = shared_mux(idx, SHARED_RDWRITE | SHARED_WAIT));
    if SHARED_OK != r {
        return r;
    }

    (r = shared_map(idx));
    if SHARED_OK != r {
        shared_demux(idx, SHARED_RDWRITE);
        return r;
    }
    if shared_attach_process(shared_gt[idx].sem) != 0
    /* try attach process */
    {
        shmdt(shared_lt[idx].p); /* cannot attach process, detach everything */
        shared_lt[idx].p = ptr::null_mut();
        shared_demux(idx, SHARED_RDWRITE);
        return SHARED_BADARG;
    }
    shared_lt[idx].tcnt += 1; /* one more thread is using segment */
    if shared_gt[idx].attr & SHARED_RESIZE
    /* if resizeable, detach and return special pointer */
    {
        if shmdt(shared_lt[idx].p) != 0 {
            r = SHARED_IPCERR; /* if segment is resizable, then detach segment */
        }
        shared_lt[idx].p = ptr::null_mut();
    }
    shared_lt[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
    r2 = shared_demux(idx, SHARED_RDWRITE);
    return if r != 0 { r } else { r2 };
}

unsafe fn shared_check_locked_index(idx: usize) -> c_int /* verify that given idx is valid */ {
    let mut r: c_int = 0;

    if !shared_init_called
    /* delayed initialization */
    {
        (r = shared_init(0));
        if SHARED_OK != r {
            return r;
        }
    }
    if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
        return SHARED_BADARG;
    }
    if shared_lt[idx].p.is_null() {
        return SHARED_BADARG; /* NULL pointer, not attached ?? */
    }
    if 0 == shared_lt[idx].lkcnt {
        return SHARED_BADARG; /* not locked ?? */
    }
    if (SHARED_ID_0 != (shared_lt[idx].p).s.ID[0])
        || (SHARED_ID_1 != (shared_lt[idx].p).s.ID[1])
        || (BLOCK_SHARED != (shared_lt[idx].p).s.tflag)
    {
        /* invalid data in segment */
        return SHARED_BADARG;
    }
    return SHARED_OK;
}

unsafe fn shared_map(idx: usize) -> c_int /* map all tables for given idx, check for validity */ {
    let mut h: c_int = 0;
    let mut bp: *mut BLKHEAD = ptr::null_mut(); /* have to obtain excl. access before calling shared_map */

    if (idx < 0) || (idx >= shared_maxseg.try_into().unwrap()) {
        return SHARED_BADARG;
    }
    if SHARED_INVALID == shared_gt[idx].key {
        return SHARED_BADARG;
    }

    (h = shmget(shared_gt[idx].key, 1, shared_create_mode));
    if SHARED_INVALID == h {
        return SHARED_BADARG;
    }

    (bp = shmat(h, std::ptr::null(), 0) as *mut BLKHEAD);
    if (SHARED_INVALID as *const BLKHEAD) == bp {
        return SHARED_BADARG;
    }
    if (SHARED_ID_0 != bp.s.ID[0])
        || (SHARED_ID_1 != bp.s.ID[1])
        || (BLOCK_SHARED != bp.s.tflag)
        || (h != shared_gt[idx].handle)
    {
        shmdt(bp); /* invalid segment, detach everything */
        return SHARED_BADARG;
    }
    if shared_gt[idx].sem != semget(shared_gt[idx].semkey, 1, shared_create_mode)
    /* check if sema is still there */
    {
        shmdt(bp); /* cannot attach semaphore, detach everything */
        return SHARED_BADARG;
    }
    shared_lt[idx].p = bp; /* store pointer to shmem data */
    return SHARED_OK;
}

unsafe fn shared_validate(idx: usize, mode: c_int) -> c_int /* use intrnally inside crit.sect !!! */
{
    let mut r: c_int = 0;

    (r = shared_mux(idx, mode));
    if SHARED_OK != r {
        return r; /* idx checked by shared_mux */
    }
    if shared_lt[idx].p.is_null() {
        (r = shared_map(idx));
        if SHARED_OK != r {
            shared_demux(idx, mode);
            return r;
        }
    }
    if (SHARED_ID_0 != (shared_lt[idx].p).s.ID[0])
        || (SHARED_ID_1 != (shared_lt[idx].p).s.ID[1])
        || (BLOCK_SHARED != (shared_lt[idx].p).s.tflag)
    {
        shared_demux(idx, mode);
        return r;
    }
    return SHARED_OK;
}

unsafe fn shared_realloc(idx: usize, newsize: c_int) -> SHARED_P /* realloc shared memory segment */
{
    let mut h: c_int = 0;
    let mut key: c_int = 0;
    let mut i: c_int = 0;
    let mut r: c_int = 0;
    let mut bp: *mut BLKHEAD = ptr::null_mut();
    let mut transfersize: c_long = 0;

    r = SHARED_OK;
    if newsize < 0 {
        return ptr::null();
    }
    if shared_check_locked_index(idx) != 0 {
        return ptr::null();
    }
    if 0 == (shared_gt[idx].attr & SHARED_RESIZE) {
        return ptr::null();
    }
    if -1 != shared_lt[idx].lkcnt {
        return ptr::null(); /* check for RW lock */
    }
    if shared_adjust_size(shared_gt[idx].size) == shared_adjust_size(newsize) {
        shared_gt[idx].size = newsize;

        return (SHARED_P)((shared_lt[idx].p) + 1);
    }

    loop {
        if i >= shared_range {
            return ptr::null(); /* table full, signal error & exit */
        }
        key = shared_kbase + ((i + shared_get_hash(newsize, idx)) % shared_range);
        h = shmget(
            key,
            shared_adjust_size(newsize).try_into().unwrap(),
            IPC_CREAT | IPC_EXCL | shared_create_mode,
        );
        if SHARED_INVALID == h {
            i += 1;
            continue; /* segment already accupied */
        }
        bp = shmat(h, std::ptr::null(), 0) as *const BLKHEAD; /* try attach */
        if (SHARED_INVALID as *const BLKHEAD) == bp
        /* cannot attach, delete segment, try with another key */
        {
            shmctl(h, IPC_RMID, std::ptr::null_mut());
            i += 1;
            continue;
        }
        *bp = *(shared_lt[idx].p); /* copy header, then data */
        transfersize = if newsize < shared_gt[idx].size {
            newsize.into()
        } else {
            shared_gt[idx].size
        };
        if transfersize > 0 {
            memcpy(
                bp + 1,
                (shared_lt[idx].p) + 1,
                transfersize.try_into().unwrap(),
            );
        }
        if shmdt(shared_lt[idx].p) != 0 {
            r = SHARED_IPCERR; /* try to detach old segment */
        }
        if shmctl(shared_gt[idx].handle, IPC_RMID, std::ptr::null_mut()) != 0 {
            if SHARED_OK == r {
                r = SHARED_IPCERR; /* destroy old shared memory segment */
            }
        }
        shared_gt[idx].size = newsize; /* signal new size */
        shared_gt[idx].handle = h; /* signal new handle */
        shared_gt[idx].key = key; /* signal new key */
        shared_lt[idx].p = bp;
        break;
    }
    return (SHARED_P)(bp + 1);
}

unsafe fn shared_free(idx: usize) -> c_int /* detach segment, if last process & !PERSIST, destroy segment */
{
    let mut cnt: c_int = 0;
    let mut r: c_int = 0;
    let mut r2: c_int = 0;

    (r = shared_validate(idx, SHARED_RDWRITE | SHARED_WAIT));
    if SHARED_OK != r {
        return r;
    }

    (r = shared_detach_process(shared_gt[idx].sem));
    if SHARED_OK != r
    /* update number of processes using segment */
    {
        shared_demux(idx, SHARED_RDWRITE);
        return r;
    }
    shared_lt[idx].tcnt -= 1; /* update number of threads using segment */
    if shared_lt[idx].tcnt > 0 {
        return shared_demux(idx, SHARED_RDWRITE); /* if more threads are using segment we are done */
    }
    if shmdt(shared_lt[idx].p) != 0
    /* if, we are the last thread, try to detach segment */
    {
        shared_demux(idx, SHARED_RDWRITE);
        return SHARED_IPCERR;
    }
    shared_lt[idx].p = None; /* clear entry in local table */
    shared_lt[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
    (cnt = shared_process_count(shared_gt[idx].sem));
    if -1 == cnt
    /* get number of processes hanging on segment */
    {
        shared_demux(idx, SHARED_RDWRITE);
        return SHARED_IPCERR;
    }
    if (0 == cnt) && (0 == (shared_gt[idx].attr & SHARED_PERSIST)) {
        r = shared_destroy_entry(idx); /* no procs on seg, destroy it */
    }
    r2 = shared_demux(idx, SHARED_RDWRITE);
    return if r != 0 { r } else { r2 };
}

fn shared_lock(idx: usize, mode: c_int) -> SHARED_P /* lock given segment for exclusive access */ {
    let mut r: c_int = 0;

    if shared_mux(idx, mode) != 0 {
        return ptr::null(); /* idx checked by shared_mux */
    }
    if 0 != shared_lt[idx].lkcnt {
        /* are we already locked ?? */
        (r = shared_map(idx));
        if SHARED_OK != r {
            shared_demux(idx, mode);
            return ptr::null();
        }
    }
    if None == shared_lt[idx].p {
        /* stupid pointer ?? */
        (r = shared_map(idx));
        if SHARED_OK != r {
            shared_demux(idx, mode);
            return ptr::null();
        }
    }
    if (SHARED_ID_0 != (shared_lt[idx].p).s.ID[0])
        || (SHARED_ID_1 != (shared_lt[idx].p).s.ID[1])
        || (BLOCK_SHARED != (shared_lt[idx].p).s.tflag)
    {
        shared_demux(idx, mode);
        return ptr::null();
    }
    if mode & SHARED_RDWRITE {
        shared_lt[idx].lkcnt = -1;

        shared_gt[idx].nprocdebug += 1;
    } else {
        shared_lt[idx].lkcnt += 1;
    }
    shared_lt[idx].seekpos = 0; /* r/w pointer positioned at beg of block */
    return ((shared_lt[idx].p) + 1) as SHARED_P;
}

unsafe fn shared_unlock(idx: usize) -> c_int /* unlock given segment, assumes seg is locked !! */ {
    let mut r: c_int = 0;
    let mut r2: c_int = 0;
    let mut mode: c_int = 0;
    (r = shared_check_locked_index(idx));
    if SHARED_OK != r {
        return r;
    }
    if shared_lt[idx].lkcnt > 0 {
        shared_lt[idx].lkcnt -= 1; /* unlock read lock */
        mode = SHARED_RDONLY;
    } else {
        shared_lt[idx].lkcnt = 0; /* unlock write lock */
        shared_gt[idx].nprocdebug -= 1;
        mode = SHARED_RDWRITE;
    }
    if 0 == shared_lt[idx].lkcnt {
        if shared_gt[idx].attr & SHARED_RESIZE {
            if shmdt(shared_lt[idx].p) {
                r = SHARED_IPCERR; /* segment is resizable, then detach segment */
            }
            shared_lt[idx].p = None; /* signal detachment in local table */
        }
    }
    r2 = shared_demux(idx, mode); /* unlock segment, rest is only parameter checking */
    return if r != 0 { r } else { r2 };
}

/* API routines - support and info routines */

unsafe fn shared_attr(idx: usize) -> c_int /* get the attributes of the shared memory segment */ {
    let mut r: c_int = 0;

    if shared_check_locked_index(idx) {
        return SHARED_INVALID;
    }
    r = shared_gt[idx].attr;
    return r;
}

unsafe fn shared_set_attr(idx: usize, newattr: c_int) -> c_int /* get the attributes of the shared memory segment */
{
    let mut r: c_int = 0;

    if shared_check_locked_index(idx) {
        return SHARED_INVALID;
    }
    if -1 != shared_lt[idx].lkcnt {
        return SHARED_INVALID; /* ADDED - check for RW lock */
    }
    r = shared_gt[idx].attr;
    shared_gt[idx].attr = newattr;
    return r;
}

unsafe fn shared_set_debug(mode: c_int) -> c_int /* set/reset debug mode */ {
    let r: c_int = shared_debug;

    shared_debug = mode;
    return r;
}

unsafe fn shared_set_createmode(mode: c_int) -> c_int /* set/reset debug mode */ {
    let r: c_int = shared_create_mode;

    shared_create_mode = mode;
    return r;
}

unsafe fn shared_list(id: c_int) -> c_int {
    let mut i: c_int = 0;
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
    print!(" Idx    Key   Nproc   Size   Flags\n");
    print!("==============================================\n");
    for i in 0..shared_maxseg {
        if -1 != id {
            if i != id {
                continue;
            }
        }
        if SHARED_INVALID == shared_gt[i].key {
            continue; /* unused slot */
        }

        /* acquire exclusive access to segment, but do not wait */
        match shared_mux(i.try_into().unwrap(), SHARED_NOWAIT | SHARED_RDONLY) {
            SHARED_AGAIN => {
                print!(
                    "!{:3} {:08x} {:4}  {:8}",
                    i, shared_gt[i].key as c_ulong, shared_gt[i].nprocdebug, shared_gt[i].size
                );
                if SHARED_RESIZE & shared_gt[i].attr {
                    print!(" RESIZABLE");
                }
                if SHARED_PERSIST & shared_gt[i].attr {
                    print!(" PERSIST");
                }
                print!("\n");
                break;
            }
            SHARED_OK => {
                print!(
                    " {:3} {:08x} {:4}  {:8}",
                    i, shared_gt[i].key as c_ulong, shared_gt[i].nprocdebug, shared_gt[i].size
                );
                if SHARED_RESIZE & shared_gt[i].attr {
                    print!(" RESIZABLE");
                }
                if SHARED_PERSIST & shared_gt[i].attr {
                    print!(" PERSIST");
                }
                print!("\n");
                shared_demux(i.try_into().unwrap(), SHARED_RDONLY);
                break;
            }
            _ => continue,
        }
    }
    if shared_debug {
        print!(" done\n");
    }
    return r; /* table full */
}

unsafe fn shared_getaddr(id: c_int, address: &mut *mut c_char) -> c_int {
    let mut i: c_int = 0;
    let mut segname: [c_char; 10] = [0; 10];

    if shared_gt.is_empty() {
        return SHARED_NOTINIT; /* not initialized */
    }
    if shared_lt.is_empty() {
        return SHARED_NOTINIT; /* not initialized */
    }

    strcpy_safe(&mut segname, cs!("h"));
    snprintf(segname[1..].as_mut_ptr(), 9, c"{}".as_ptr(), id);

    if smem_open(segname.as_ptr(), 0, &mut i) != 0 {
        return SHARED_BADARG;
    }

    *address = (((shared_lt[i].p + 1) as *const DAL_SHM_SEGHEAD) + 1) as *const c_char;
    /*  smem_close(i); */
    return SHARED_OK;
}

unsafe fn shared_uncond_delete(id: c_int) -> c_int {
    let mut i: c_int = 0;
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
        if -1 != id {
            if i != id {
                continue;
            }
        }
        if shared_attach(i.try_into().unwrap()) != 0 {
            if -1 != id {
                print!("no such handle\n");
            }
            continue;
        }
        print!("handle {}:", i);
        if shared_lock(i.try_into().unwrap(), SHARED_RDWRITE | SHARED_NOWAIT).is_null() {
            print!(" cannot lock in RW mode, not deleted\n");
            continue;
        }
        if shared_set_attr(i.try_into().unwrap(), SHARED_RESIZE) >= SHARED_ERRBASE {
            print!(" cannot clear PERSIST attribute");
        }
        if shared_free(i.try_into().unwrap()) != 0 {
            print!(" delete failed\n");
        } else {
            print!(" deleted\n");
        }
    }
    if shared_debug {
        print!(" done\n");
    }
    return r; /* table full */
}

/************************* CFITSIO DRIVER FUNCTIONS ***************************/

pub(crate) fn smem_init() -> c_int {
    return 0;
}

pub(crate) fn smem_shutdown() -> c_int {
    if shared_init_called {
        shared_cleanup();
    }
    return 0;
}

pub(crate) fn smem_setoptions(mut option: c_int) -> c_int {
    option = 0;
    return 0;
}

pub(crate) fn smem_getoptions(options: &mut c_int) -> c_int {
    *options = 0;
    return 0;
}

pub(crate) fn smem_getversion(version: &mut c_int) -> c_int {
    *version = 10;
    return 0;
}

pub(crate) fn smem_open(filename: *const c_char, rwmode: c_int, driverhandle: &mut c_int) -> c_int {
    let mut h: c_int = 0;
    let mut nitems: c_int = 0;
    let mut r: c_int = 0;
    let mut sp: *mut DAL_SHM_SEGHEAD = ptr::null_mut();

    if filename.is_null() {
        return SHARED_NULPTR;
    }

    nitems = sscanf(filename, c"h%d".as_ptr(), &h);
    if 1 != nitems {
        return SHARED_BADARG;
    }

    (r = shared_attach(h.try_into().unwrap()));
    if SHARED_OK != r {
        return r;
    }

    (sp = shared_lock(
        h.try_into().unwrap(),
        if READWRITE == rwmode {
            SHARED_RDWRITE
        } else {
            SHARED_RDONLY
        },
    ) as *mut DAL_SHM_SEGHEAD);
    if sp.is_null() {
        shared_free(h.try_into().unwrap());
        return SHARED_BADARG;
    }

    if (h != sp.h) || (DAL_SHM_SEGHEAD_ID != sp.ID) {
        shared_unlock(h.try_into().unwrap());
        shared_free(h.try_into().unwrap());

        return SHARED_BADARG;
    }

    *driverhandle = h;
    return 0;
}

pub(crate) fn smem_create(filename: *const c_char, driverhandle: &mut c_int) -> c_int {
    let mut sp: *mut DAL_SHM_SEGHEAD = ptr::null_mut();
    let mut h: c_int = 0;
    let mut sz: c_int = 0;
    let mut nitems: c_int = 0;

    if filename.is_null() {
        return SHARED_NULPTR; /* currently ignored */
    }

    nitems = sscanf(filename, c"h%d".as_ptr(), &h);
    if 1 != nitems {
        return SHARED_BADARG;
    }

    sz = BL!() + size_of::<DAL_SHM_SEGHEAD>();
    h = shared_malloc(sz.into(), SHARED_RESIZE | SHARED_PERSIST, h);
    if SHARED_INVALID == h {
        return SHARED_NOMEM;
    }

    (sp = shared_lock(h.try_into().unwrap(), SHARED_RDWRITE) as *const DAL_SHM_SEGHEAD);
    if sp.is_null() {
        shared_free(h.try_into().unwrap());
        return SHARED_BADARG;
    }

    sp.ID = DAL_SHM_SEGHEAD_ID;
    sp.h = h;
    sp.size = sz;
    sp.nodeidx = -1;

    *driverhandle = h;

    return 0;
}

pub(crate) fn smem_close(driverhandle: usize) -> c_int {
    let mut r: c_int = 0;

    (r = shared_unlock(driverhandle));
    if SHARED_OK != r {
        return r;
    }
    return shared_free(driverhandle);
}

pub(crate) fn smem_remove(filename: *const c_char) -> c_int {
    let mut h: c_int = 0;
    let mut nitems: c_int = 0;
    let mut r: c_int = 0;

    if filename.is_null() {
        return SHARED_NULPTR;
    }
    nitems = sscanf(filename, c"h%d".as_ptr(), &h);
    if 1 != nitems {
        return SHARED_BADARG;
    }

    if 0 == shared_check_locked_index(h.try_into().unwrap())
    /* are we locked ? */
    {
        if -1 != shared_lt[h].lkcnt
        /* are we locked RO ? */
        {
            (r = shared_unlock(h.try_into().unwrap()));
            if SHARED_OK != r {
                return r; /* yes, so relock in RW */
            }
            if shared_lock(h.try_into().unwrap(), SHARED_RDWRITE).is_null() {
                return SHARED_BADARG;
            }
        }
    } else
    /* not locked */
    {
        (r = smem_open(filename, READWRITE, &h));
        if SHARED_OK != r {
            return r; /* so open in RW mode */
        }
    }

    shared_set_attr(h.try_into().unwrap(), SHARED_RESIZE); /* delete PERSIST attribute */
    return smem_close(h.try_into().unwrap()); /* detach segment (this will delete it) */
}

pub(crate) unsafe fn smem_size(driverhandle: usize, size: *mut LONGLONG) -> c_int {
    if size.is_null() {
        return SHARED_NULPTR;
    }
    if shared_check_locked_index(driverhandle) != 0 {
        return SHARED_INVALID;
    }
    *size = (shared_gt[driverhandle].size - size_of::<DAL_SHM_SEGHEAD>()) as LONGLONG;
    return 0;
}

pub(crate) fn smem_flush(driverhandle: usize) -> c_int {
    if shared_check_locked_index(driverhandle) != 0 {
        return SHARED_INVALID;
    }
    return 0;
}

pub(crate) fn smem_seek(driverhandle: usize, offset: LONGLONG) -> c_int {
    if offset < 0 {
        return SHARED_BADARG;
    }
    if shared_check_locked_index(driverhandle) != 0 {
        return SHARED_INVALID;
    }
    shared_lt[driverhandle].seekpos = offset;
    return 0;
}

pub(crate) fn smem_read(driverhandle: usize, buffer: *mut c_void, nbytes: usize) -> c_int {
    if buffer.is_null() {
        return SHARED_NULPTR;
    }
    if shared_check_locked_index(driverhandle) != 0 {
        return SHARED_INVALID;
    }
    if nbytes < 0 {
        return SHARED_BADARG;
    }
    if (shared_lt[driverhandle].seekpos + nbytes) > shared_gt[driverhandle].size {
        return SHARED_BADARG; /* read beyond EOF */
    }

    memcpy(
        buffer,
        ((((shared_lt[driverhandle].p + 1) as *const DAL_SHM_SEGHEAD) + 1) as *const c_char)
            + shared_lt[driverhandle].seekpos,
        nbytes,
    );

    shared_lt[driverhandle].seekpos += nbytes;
    return 0;
}

pub(crate) fn smem_write(driverhandle: usize, buffer: *mut c_void, nbytes: usize) -> c_int {
    if buffer.is_null() {
        return SHARED_NULPTR;
    }
    if shared_check_locked_index(driverhandle) != 0 {
        return SHARED_INVALID;
    }
    if -1 != shared_lt[driverhandle].lkcnt {
        return SHARED_INVALID; /* are we locked RW ? */
    }

    if nbytes < 0 {
        return SHARED_BADARG;
    }
    if (shared_lt[driverhandle].seekpos + nbytes) as c_ulong
        > (shared_gt[driverhandle].size - size_of::<DAL_SHM_SEGHEAD>()) as c_ulong
    {
        /* need to realloc shmem */
        if shared_realloc(
            driverhandle,
            shared_lt[driverhandle].seekpos + nbytes + size_of::<DAL_SHM_SEGHEAD>(),
        )
        .is_null()
        {
            return SHARED_NOMEM;
        }
    }

    memcpy(
        ((((shared_lt[driverhandle].p + 1) as *const DAL_SHM_SEGHEAD) + 1) as *const c_char)
            + shared_lt[driverhandle].seekpos,
        buffer,
        nbytes,
    );

    shared_lt[driverhandle].seekpos += nbytes;
    return 0;
}
