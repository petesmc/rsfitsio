/*  This file, cfileio.c, contains the low-level file access routines.     */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::ffi::CStr;
use std::sync::{Mutex, OnceLock};
use std::{cmp, mem, ptr};

use errno::{Errno, errno, set_errno};

use libc::{ERANGE, fclose, fgets, fopen, fprintf};

use crate::c_types::{FILE, c_char, c_int, c_long, c_void};
use bytemuck::{cast_mut, cast_slice, cast_slice_mut};
use cstr::cstr;

use crate::aliases::safer::fits_read_key_str;
use crate::aliases::{ffgisz_safe, ffpmsg_cstr};
use crate::drvrfile::{
    file_checkfile, file_close, file_compress_open, file_create, file_flush, file_getoptions,
    file_getversion, file_init, file_open, file_read, file_remove, file_seek, file_setoptions,
    file_shutdown, file_size, file_truncate, file_write, stream_close, stream_create, stream_flush,
    stream_open, stream_read, stream_seek, stream_size, stream_write,
};
use crate::drvrmem::{
    mem_close_comp_unsafe, mem_close_free_unsafe, mem_close_keep, mem_compress_open,
    mem_compress_openrw, mem_create, mem_create_comp_unsafe, mem_getoptions, mem_getversion,
    mem_init, mem_iraf_open, mem_openmem, mem_rawfile_open, mem_read_unsafe, mem_seek,
    mem_setoptions, mem_shutdown, mem_size, mem_truncate_unsafe, mem_write_unsafe, stdin_checkfile,
    stdin_open, stdout_close_unsafe,
};

#[cfg(feature = "shared_mem")]
use crate::drvrsmem::{
    smem_close, smem_create, smem_flush, smem_getoptions, smem_getversion, smem_init, smem_open,
    smem_read, smem_remove, smem_seek, smem_setoptions, smem_shutdown, smem_size, smem_write,
};
use crate::editcol::ffdcol_safe;
use crate::edithdu::ffcopy_safer;
use crate::eval_f::{ffcalc_safe, ffffrw_safer, fffrow_safe};
use crate::fitscore::{
    ffchdu, ffcmsg_safe, ffgcno_safe, ffgerr_safe, ffghdn_safe, ffghdt_safe, ffgidm_safe,
    ffgmsg_safe, ffgncl_safe, ffgnrw_safe, ffkeyn_safe, ffmahd_safe, ffmnhd_safe, ffmrhd_safe,
    ffpmsg_slice, ffpmsg_str, ffrhdu_safer, ffupch_safe, fits_strcasecmp, fits_strncasecmp,
};
use crate::getkey::{ffgcrd_safe, ffgkyl_safe, ffmaky_safe};
use crate::group::{fits_clean_url, fits_get_cwd, fits_path2url};
use crate::histo::{ffbinse, ffhist2e};
use crate::modkey::{ffdkey_safe, ffmkys_safe, ffmnam_safe};
use crate::putkey::ffphis_safe;
use crate::relibc::header::stdio::sscanf;
use crate::wrappers::*;
use crate::{FFLOCK, FFUNLOCK};
use crate::{bb, cs};
use crate::{buffers::*, raw_to_slice};
use crate::{fitsio::*, int_snprintf};
use crate::{fitsio2::*, slice_to_str};

pub const MAX_PREFIX_LEN: usize = 20; /* max length of file type prefix (e.g. 'http://') */
pub const MAX_DRIVERS: usize = 31; /* max number of file I/O drivers */

pub(crate) trait Driver {}

pub struct fitsdriver {
    /* structure containing pointers to I/O driver functions */
    pub prefix: [c_char; MAX_PREFIX_LEN],
    pub init: Option<fn() -> c_int>,
    pub shutdown: Option<fn() -> c_int>,
    pub setoptions: Option<fn(option: c_int) -> c_int>,
    pub getoptions: Option<fn(options: &mut c_int) -> c_int>,
    pub getversion: Option<fn(version: &mut c_int) -> c_int>,
    pub checkfile: Option<
        fn(
            urltype: &mut [c_char; MAX_PREFIX_LEN],
            infile: &mut [c_char; FLEN_FILENAME],
            outfile: &mut [c_char; FLEN_FILENAME],
        ) -> c_int,
    >,
    pub open: Option<fn(filename: &mut [c_char], rwmode: c_int, handle: &mut c_int) -> c_int>,
    pub create:
        Option<fn(filename: &mut [c_char; FLEN_FILENAME], drivehandle: &mut c_int) -> c_int>,
    pub truncate: Option<fn(drivehandle: c_int, size: usize) -> c_int>,
    pub close: fn(drivehandle: c_int) -> c_int,
    pub remove: Option<fn(filename: &[c_char]) -> c_int>,
    pub size: fn(drivehandle: c_int, size: &mut usize) -> c_int,
    pub flush: Option<fn(drivehandle: c_int) -> c_int>,
    pub seek: fn(drivehandle: c_int, offset: LONGLONG) -> c_int,
    pub read: fn(drivehandle: c_int, buffer: &mut [u8], nbytes: usize) -> c_int,
    pub write: fn(drivehandle: c_int, buffer: &[u8], nbyte: usize) -> c_int,
}

pub static DRIVER_TABLE: OnceLock<Vec<fitsdriver>> = OnceLock::new(); /* allocate driver tables */

/* this table of Fptr pointers is used by fits_already_open */
pub static mut FPTR_TABLE: [*mut FITSfile; NMAXFILES] =
    [ptr::null::<FITSfile>() as *mut FITSfile; NMAXFILES];

pub static NEED_TO_INITIALIZE: Mutex<bool> = Mutex::new(true); /* true if CFITSIO has not been initialized */

pub static STREAM_DRIVER: Mutex<c_int> = Mutex::new(0); /* number of currently defined I/O drivers */

/*--------------------------------------------------------------------------*/
pub(crate) fn fitsio_init_lock() -> c_int {
    0
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file in core memory.  This is a specialized version
/// of ffopen.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffomem(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - name of file to open                */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    buffptr: *const *const c_void,    /* I - address of memory pointer           */
    buffsize: *mut usize,             /* I - size of buffer, in bytes            */
    deltasize: usize,                 /* I - increment for future realloc's      */
    mem_realloc: unsafe extern "C" fn(p: *mut c_void, newsize: usize) -> *mut c_void, /* function       */
    status: *mut c_int, /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut driver: c_int = 0;
        let mut handle: c_int = 0;
        let mut hdutyp: c_int = 0;

        let mut movetotype: c_int = 0;
        let mut extvers: c_int = 0;
        let mut extnum: c_int = 0;
        let mut extname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut filesize: usize = 0;
        let mut urltype: [c_char; MAX_PREFIX_LEN] = [0; MAX_PREFIX_LEN];
        let mut infile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut outfile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut extspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut rowfilter: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut binspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut colspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut imagecolname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut rowexpress: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];

        let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let hdtype: [*const c_char; 3] = [
            cstr!(b"IMAGE").as_ptr(),
            cstr!("TABLE").as_ptr(),
            cstr!("BINTABLE").as_ptr(),
        ];

        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(name);

        if *status > 0 {
            return *status;
        }

        /* initialize null file pointer */
        let f_tmp = fptr.take();
        if let Some(f) = f_tmp {
            // WARNING: The c version doesn't null pointers after a close, so we have a dangling pointer.
            // We need to be careful with this, as it can cause double free errors.
            // Therefore, if this function is called with a Some(), then we will leak the pointer because
            // it's probably invalid.
            let _ = Box::into_raw(f);
        }

        if *NEED_TO_INITIALIZE.lock().unwrap() {
            /* this is called only once */
            *status = fits_init_cfitsio_safer();

            if *status > 0 {
                return *status;
            }
        }

        let mut url = 0;

        while name[url] == bb(b' ') {
            /* ignore leading spaces in the file spec */
            url += 1;
        }

        /* parse the input file specification */
        ffiurl_safer(
            &name[url..],
            urltype.as_mut_ptr(),
            infile.as_mut_ptr(),
            outfile.as_mut_ptr(),
            extspec.as_mut_ptr(),
            rowfilter.as_mut_ptr(),
            binspec.as_mut_ptr(),
            colspec.as_mut_ptr(),
            status,
        );

        strcpy_safe(&mut urltype, cs!("memkeep://")); /* URL type for pre-existing memory file */

        *status = urltype2driver(&urltype, &mut driver);

        if *status > 0 {
            ffpmsg_str("could not find driver for pre-existing memory file: (ffomem)");
            return *status;
        }

        /* call driver routine to open the memory file */
        let lock = FFLOCK(); /* lock this while searching for vacant handle */
        *status = mem_openmem(
            buffptr,
            buffsize.as_mut().unwrap(),
            deltasize,
            mem_realloc,
            &mut handle,
        );
        FFUNLOCK(lock);

        if *status > 0 {
            ffpmsg_str("failed to open pre-existing memory file: (ffomem)");
            return *status;
        }

        /* get initial file size */
        //let d = driverTable.lock().unwrap();
        let d = DRIVER_TABLE.get().unwrap();
        *status = (d[driver as usize].size)(handle, &mut filesize);

        if *status > 0 {
            (d[driver as usize].close)(handle); /* close the file */
            ffpmsg_str("failed get the size of the memory file: (ffomem)");
            return *status;
        }

        let Fptr = FITSfile::new(
            &d[driver as usize],
            handle,
            &name[url..],
            cs!("ffomem"),
            status,
        );
        if Fptr.is_err() {
            return *status;
        }

        let mut Fptr = Fptr.unwrap();

        /* initialize the ageindex array (relative age of the I/O buffers) */
        /* and initialize the bufrecnum array as being empty */
        for ii in 0..(NIOBUF as usize) {
            Fptr.ageindex[ii] = ii as c_int;
            Fptr.bufrecnum[ii] = -1;
        }

        /* store the parameters describing the file */
        Fptr.MAXHDU = 1000; /* initial size of headstart */
        Fptr.filehandle = handle; /* file handle */
        Fptr.driver = driver; /* driver number */
        strcpy(Fptr.filename, name[url..].as_ptr()); /* full input filename */
        Fptr.filesize = filesize as LONGLONG; /* physical file size */
        Fptr.logfilesize = filesize as LONGLONG; /* logical file size */
        Fptr.writemode = mode; /* read-write mode    */
        Fptr.datastart = DATA_UNDEFINED as LONGLONG; /* unknown start of data */
        Fptr.curbuf = -1; /* undefined current IO buffer */
        Fptr.open_count = 1; /* structure is currently used once */
        Fptr.validcode = VALIDSTRUC; /* flag denoting valid structure */
        Fptr.noextsyntax = 0; /* extended syntax can be used in filename */

        let f_fitsfile = Box::try_new(fitsfile {
            HDUposition: 0,
            Fptr,
        });

        if f_fitsfile.is_err() {
            let d = DRIVER_TABLE.get().unwrap();
            ((d[driver as usize]).close)(handle); /* close the file */
            ffpmsg_str("failed to allocate structure for following file: (ffomem)");
            ffpmsg_slice(&name[url..]);
            *status = MEMORY_ALLOCATION;
            return *status;
        }

        let mut f_fitsfile = f_fitsfile.unwrap();

        ffldrc(&mut f_fitsfile, 0, REPORT_EOF, status); /* load first record */

        fits_store_Fptr(&mut f_fitsfile.Fptr, status); /* store Fptr address */

        if ffrhdu_safer(&mut f_fitsfile, Some(&mut hdutyp), status) > 0 {
            /* determine HDU structure */
            ffpmsg_str("ffomem could not interpret primary array header of file: (ffomem)");
            ffpmsg_slice(&name[url..]);

            if *status == UNKNOWN_REC {
                ffpmsg_str("This does not look like a FITS file.");
            }

            ffclos_safer(f_fitsfile, status);
            *fptr = None; /* return null file pointer */
            return *status;
        }

        *fptr = Some(f_fitsfile);

        /* ---------------------------------------------------------- */
        /* move to desired extension, if specified as part of the URL */
        /* ---------------------------------------------------------- */

        imagecolname[0] = 0;
        rowexpress[0] = 0;

        if extspec[0] != 0 {
            /* parse the extension specifier into individual parameters */
            ffexts_safer(
                &extspec,
                &mut extnum,
                extname.as_mut_ptr(),
                &mut extvers,
                &mut movetotype,
                imagecolname.as_mut_ptr(),
                rowexpress.as_mut_ptr(),
                status,
            );

            if *status > 0 {
                return *status;
            }

            if extnum != 0 {
                ffmahd_safe(
                    (*fptr).as_mut().unwrap(),
                    extnum + 1,
                    Some(&mut hdutyp),
                    status,
                );
            } else if extname[0] != 0 {
                /* move to named extension, if specified */
                ffmnhd_safe(
                    (*fptr).as_mut().unwrap(),
                    movetotype,
                    &extname,
                    extvers,
                    status,
                );
            }

            if *status > 0 {
                ffpmsg_str("ffomem could not move to the specified extension:");
                if extnum > 0 {
                    int_snprintf!(
                        &mut errmsg,
                        FLEN_ERRMSG,
                        " extension number {} doesn't exist or couldn't be opened.",
                        extnum,
                    );
                    ffpmsg_slice(&errmsg);
                } else {
                    int_snprintf!(
                        &mut errmsg,
                        FLEN_ERRMSG,
                        " extension with EXTNAME = {},",
                        slice_to_str!(&extname),
                    );
                    ffpmsg_slice(&errmsg);

                    if extvers != 0 {
                        int_snprintf!(
                            &mut errmsg,
                            FLEN_ERRMSG,
                            "           and with EXTVERS = {},",
                            extvers,
                        );
                        ffpmsg_slice(&errmsg);
                    }

                    if movetotype != ANY_HDU {
                        int_snprintf!(
                            &mut errmsg,
                            FLEN_ERRMSG,
                            "           and with XTENSION = {},",
                            CStr::from_ptr(hdtype[movetotype as usize])
                                .to_str()
                                .unwrap(),
                        );
                        ffpmsg_slice(&errmsg);
                    }
                    ffpmsg_str(" doesn't exist or couldn't be opened.");
                }
                return *status;
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file on magnetic disk with either readonly or
/// read/write access.  The routine does not support CFITSIO's extended
/// filename syntax and simply uses the entire input 'name' string as
/// the name of the file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdkopn(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(name);

        if *status > 0 {
            return *status;
        }

        *status = OPEN_DISK_FILE;

        ffopen_safer(fptr, name, mode, status);

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access. and
/// move to the first HDU that contains 'interesting' data, if the primary
/// array contains a null image (i.e., NAXIS = 0).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdopn(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(name);

        if *status > 0 {
            return *status;
        }

        *status = SKIP_NULL_PRIMARY;

        ffopen_safer(fptr, name, mode, status);

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access. and
/// if the primary array contains a null image (i.e., NAXIS = 0) then attempt to
/// move to the first extension named in the extlist of extension names. If
/// none are found, then simply move to the 2nd extension.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffeopn(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    extlist: *mut c_char,             /* I - list of 'good' extensions to move to */
    hdutype: *mut c_int,              /* O - type of extension that is moved to  */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut hdunum: c_int = 0;
        let mut naxis: c_int = 0;
        let mut thdutype: c_int = 0;
        let mut gotext = false;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let hdutype = hdutype.as_mut();

        raw_to_slice!(name);

        if *status > 0 {
            return *status;
        }

        if ffopen_safer(fptr, name, mode, status) > 0 {
            return *status;
        }

        let f = (*fptr).as_mut().expect(NULL_MSG);

        ffghdn_safe(f, &mut hdunum);
        ffghdt_safe(f, &mut thdutype, status);

        if hdunum == 1 && thdutype == IMAGE_HDU {
            ffgidm_safe(f, &mut naxis, status);
        }

        /* We are in the "default" primary extension */
        /* look through the extension list */
        if (hdunum == 1) && (naxis == 0) {
            if !extlist.is_null() {
                raw_to_slice!(extlist); // We know its non-null here

                gotext = false;

                // HEAP ALLOCATION - Temporary
                let mut textlist = Vec::new();
                if textlist.try_reserve_exact(extlist.len()).is_err() {
                    *status = MEMORY_ALLOCATION;
                    return *status;
                } else {
                    textlist.resize(extlist.len(), 0);
                }

                strcpy_safe(&mut textlist, extlist);
                for ext_sub in textlist.split(|c| *c == bb(b' ')) {
                    ffmnhd_safe(f, ANY_HDU, ext_sub, 0, status);
                    if *status == 0 {
                        gotext = true;
                        break;
                    } else {
                        *status = 0;
                    }
                }
            }
            if !gotext {
                /* if all else fails, move to extension #2 and hope for the best */
                ffmahd_safe(f, 2, Some(&mut thdutype), status);
            }
        }

        if let Some(hdutype) = hdutype {
            ffghdt_safe(f, hdutype, status);
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access. and
/// move to the first HDU that contains 'interesting' table (not an image).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fftopn(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut hdutype: c_int = 0;

        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(name);

        if *status > 0 {
            return *status;
        }

        *status = SKIP_IMAGE;

        ffopen_safer(fptr, name, mode, status);

        let f = (*fptr).as_mut().expect(NULL_MSG);

        if ffghdt_safe(f, &mut hdutype, status) <= 0 && hdutype == IMAGE_HDU {
            *status = NOT_TABLE;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access. and
/// move to the first HDU that contains 'interesting' image (not an table).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffiopn(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut hdutype: c_int = 0;

        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(name);

        if *status > 0 {
            return *status;
        }

        *status = SKIP_TABLE;

        ffopen_safer(fptr, name, mode, status);

        let f = (*fptr).as_mut().expect(NULL_MSG);

        if ffghdt_safe(f, &mut hdutype, status) <= 0 && hdutype != IMAGE_HDU {
            *status = NOT_IMAGE;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access.
///
/// First test that the SONAME of fitsio.h used to build the CFITSIO library
/// is the same as was used in compiling the application program that
/// links to the library.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffopentest(
    soname: c_int, /* I - CFITSIO shared library version     */
    /*     application program (fitsio.h file) */
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(name);

        ffopentest_safe(soname, fptr, name, mode, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access.
///
/// First test that the SONAME of fitsio.h used to build the CFITSIO library
/// is the same as was used in compiling the application program that
/// links to the library.
pub fn ffopentest_safe(
    soname: c_int, /* I - CFITSIO shared library version  application program (fitsio.h file) */
    fptr: &mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: &[c_char], /* I - full name of file to open           */
    mode: c_int,   /* I - 0 = open readonly; 1 = read/write   */
    status: &mut c_int, /* IO - error status                       */
) -> c_int {
    unsafe {
        if soname != CFITSIO_SONAME as c_int {
            println!("\nERROR: Mismatch in the CFITSIO_SONAME value in the fitsio.h include file");
            println!(
                "that was used to build the CFITSIO library, and the value in the include file"
            );
            println!("that was used when compiling the application program:");
            println!("   Version used to build the CFITSIO library   = {CFITSIO_SONAME}");
            println!("   Version included by the application program = {soname}");
            print!("\nFix this by recompiling and then relinking this application program \n");
            println!("with the CFITSIO library.");

            *status = FILE_NOT_OPENED;
            return *status;
        }

        /* now call the normal file open routine */
        ffopen_safer(fptr, name, mode, status);
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffopen(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        raw_to_slice!(name);

        ffopen_safer(fptr, name, mode, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Open an existing FITS file with either readonly or read/write access.
pub(crate) unsafe fn ffopen_safer(
    fptr: &mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: &[c_char],                  /* I - full name of file to open           */
    mode: c_int,                      /* I - 0 = open readonly; 1 = read/write   */
    status: &mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut newptr: Option<Box<fitsfile>> = None;
        let mut driver = 0;
        let mut hdutyp = 0;
        let mut hdunum = 0;
        let mut slen = 0;
        let mut isopen = 0;
        let mut filesize = 0;
        let mut rownum: c_long = 0;
        let mut nrows: c_long = 0;
        let mut goodrows: c_long = 0;
        let mut extnum = 0;
        let mut extvers = 0;
        let mut handle = 0;
        let mut movetotype = 0;
        let mut tstatus = 0;
        let mut only_one = 0;
        let mut urltype: [c_char; MAX_PREFIX_LEN] = [0; MAX_PREFIX_LEN];
        let mut infile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut outfile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut origurltype: [c_char; MAX_PREFIX_LEN] = [0; MAX_PREFIX_LEN];
        let mut extspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut extname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut rowfilter: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut tblname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut imagecolname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut rowexpress: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut binspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut colspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut pixfilter: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut histfilename: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut filtfilename: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut compspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let wtcol: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
        let maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
        let binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
        let minin: [f64; 4] = [0.0; 4];
        let maxin: [f64; 4] = [0.0; 4];
        let binsizein: [f64; 4] = [0.0; 4];
        let mut weight: f64 = 0.0;
        let mut imagetype = 0;
        let mut naxis = 1;
        let mut haxis = 0;
        let mut recip = 0;
        let mut skip_null = false;
        let mut skip_image = false;
        let mut skip_table = false;
        let mut no_primary_data = false;
        let mut open_disk_file = false;
        let colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
        let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let hdtype: [*const c_char; 3] = [
            cstr!(b"IMAGE").as_ptr(),
            cstr!("TABLE").as_ptr(),
            cstr!("BINTABLE").as_ptr(),
        ];

        let rowselect = None;

        if *status > 0 {
            return *status;
        }

        if *status == SKIP_NULL_PRIMARY {
            /* this special status value is used as a flag by ffdopn to tell */
            /* ffopen to skip over a null primary array when opening the file. */

            skip_null = true;
            *status = 0;
        } else if *status == SKIP_IMAGE {
            /* this special status value is used as a flag by fftopn to tell */
            /* ffopen to move to 1st significant table when opening the file. */

            skip_image = true;
            *status = 0;
        } else if *status == SKIP_TABLE {
            /* this special status value is used as a flag by ffiopn to tell */
            /* ffopen to move to 1st significant image when opening the file. */

            skip_table = true;
            *status = 0;
        } else if *status == OPEN_DISK_FILE {
            /* this special status value is used as a flag by ffdkopn to tell */
            /* ffopen to not interpret the input filename using CFITSIO's    */
            /* extended filename syntax, and simply open the specified disk file */

            open_disk_file = true;
            *status = 0;
        }

        /* initialize null file pointer */
        let f_tmp = fptr.take();
        if let Some(f) = f_tmp {
            // WARNING: The c version doesn't null pointers after a close, so we have a dangling pointer.
            // We need to be careful with this, as it can cause double free errors.
            // Therefore, if this function is called with a Some(), then we will leak the pointer because
            // it's probably invalid.
            let _ = Box::into_raw(f);
        }

        let mut writecopy = false; /* have we made a write-able copy of the input file? */

        if *NEED_TO_INITIALIZE.lock().unwrap() {
            /* this is called only once */
            *status = fits_init_cfitsio_safer();
        }

        if *status > 0 {
            return *status;
        }

        let url = name;
        let mut ui = 0;

        while url[ui] == bb(b' ') {
            /* ignore leading spaces in the filename */
            ui += 1;
        }

        if url[ui] == 0 {
            ffpmsg_str("Name of file to open is blank. (ffopen)");
            *status = FILE_NOT_OPENED;
            return *status;
        }

        let url = &url[ui..];

        if open_disk_file {
            /* treat the input URL literally as the name of the file to open */
            /* and don't try to parse the URL using the extended filename syntax */

            if strlen_safe(url) > FLEN_FILENAME - 1 {
                ffpmsg_str("Name of file to open is too long. (ffopen)");
                *status = FILE_NOT_OPENED;
                return *status;
            }

            strcpy_safe(&mut infile, url);
            strcpy_safe(&mut urltype, cs!(b"file://"));
            outfile[0] = 0;
            extspec[0] = 0;
            binspec[0] = 0;
            colspec[0] = 0;
            rowfilter[0] = 0;
            pixfilter[0] = 0;
            compspec[0] = 0;
        } else {
            /* parse the input file specification */

            /* NOTE: This routine tests that all the strings do not */
            /* overflow the standard buffer sizes (FLEN_FILENAME, etc.) */
            /* therefore in general we do not have to worry about buffer */
            /* overflow of any of the returned strings. */

            /* call the newer version of this parsing routine that supports 'compspec' */

            ffifile2_safer(
                url,
                urltype.as_mut_ptr(),
                infile.as_mut_ptr(),
                outfile.as_mut_ptr(),
                extspec.as_mut_ptr(),
                rowfilter.as_mut_ptr(),
                binspec.as_mut_ptr(),
                colspec.as_mut_ptr(),
                pixfilter.as_mut_ptr(),
                compspec.as_mut_ptr(),
                status,
            );
        }

        if *status > 0 {
            ffpmsg_str("could not parse the input filename: (ffopen)");
            ffpmsg_slice(url);
            return *status;
        }

        imagecolname[0] = 0;
        rowexpress[0] = 0;

        if extspec[0] != 0 {
            slen = strlen_safe(&extspec);
            if extspec[slen - 1] == bb(b'#') {
                /* special symbol to mean only copy this extension */
                extspec[slen - 1] = 0;
                only_one = 1;
            }

            /* parse the extension specifier into individual parameters */
            ffexts_safer(
                &extspec,
                &mut extnum,
                extname.as_mut_ptr(),
                &mut extvers,
                &mut movetotype,
                imagecolname.as_mut_ptr(),
                rowexpress.as_mut_ptr(),
                status,
            );
            if *status > 0 {
                return *status;
            };
        }

        /*-------------------------------------------------------------------*/
        /* special cases:                                                    */
        /*-------------------------------------------------------------------*/

        histfilename[0] = 0;
        filtfilename[0] = 0;

        if outfile[0] != 0 && (binspec[0] != 0 || imagecolname[0] != 0 || pixfilter[0] != 0) {
            /* if binspec or imagecolumn are specified, then the  */
            /* output file name is intended for the final image,  */
            /* and not a copy of the input file.                  */

            strcpy_safe(&mut histfilename, &outfile);
            outfile[0] = 0;
        } else if outfile[0] != 0 && (rowfilter[0] != 0 || colspec[0] != 0) {
            /* if rowfilter or colspece are specified, then the    */
            /* output file name is intended for the filtered file  */
            /* and not a copy of the input file.                   */

            strcpy_safe(&mut filtfilename, &outfile);
            outfile[0] = 0;
        }

        /*-------------------------------------------------------------------*/
        /* check if this same file is already open, and if so, attach to it  */
        /*-------------------------------------------------------------------*/

        let lock = FFLOCK();

        if fits_already_open(
            fptr,
            url,
            &mut urltype,
            &mut infile,
            &mut extspec,
            &mut rowfilter,
            &mut binspec,
            &mut colspec,
            mode,
            open_disk_file,
            &mut isopen,
            status,
        ) > 0
        {
            FFUNLOCK(lock);
            return *status;
        }

        FFUNLOCK(lock);

        let mut move2hdu = false;
        if isopen != 0 {
            move2hdu = true;
        }

        if !move2hdu {
            /* get the driver number corresponding to this urltype */
            *status = urltype2driver(&urltype, &mut driver);

            if *status > 0 {
                ffpmsg_str("could not find driver for this file: (ffopen)");
                ffpmsg_slice(&urltype);
                ffpmsg_slice(url);
                return *status;
            }

            /*-------------------------------------------------------------------
              deal with all those messy special cases which may require that
              a different driver be used:
                  - is disk file compressed?
                  - are ftp:, gsiftp:, or http: files compressed?
                  - has user requested that a local copy be made of
                    the ftp or http file?
            -------------------------------------------------------------------*/

            let d = DRIVER_TABLE.get().unwrap();

            if let Some(checkfile) = (d[driver as usize]).checkfile {
                strcpy_safe(&mut origurltype, &urltype); /* Save the urltype */

                /* 'checkfile' may modify the urltype, infile and outfile strings */
                *status = checkfile(&mut urltype, &mut infile, &mut outfile);

                if *status != 0 {
                    ffpmsg_str("checkfile failed for this file: (ffopen)");
                    ffpmsg_slice(url);
                    return *status;
                }

                if strcmp_safe(&origurltype, &urltype) != 0 {
                    /* did driver changed on us? */
                    *status = urltype2driver(&urltype, &mut driver);
                    if *status > 0 {
                        ffpmsg_str("could not change driver for this file: (ffopen)");
                        ffpmsg_slice(url);
                        ffpmsg_slice(&urltype);
                        return *status;
                    };
                };
            }

            /* call appropriate driver to open the file */
            match (d[driver as usize]).open {
                Some(open) => {
                    let lock = FFLOCK(); /* lock this while searching for vacant handle */
                    *status = open(&mut infile, mode, &mut handle);
                    FFUNLOCK(lock);

                    if *status > 0 {
                        ffpmsg_str("failed to find or open the following file: (ffopen)");
                        ffpmsg_slice(url);
                        return *status;
                    };
                }
                None => {
                    ffpmsg_str("cannot open an existing file of this type: (ffopen)");
                    ffpmsg_slice(url);
                    *status = FILE_NOT_OPENED;
                    return *status;
                }
            };

            /* get initial file size */
            *status = ((d[driver as usize]).size)(handle, &mut filesize);

            if *status > 0 {
                ((d[driver as usize]).close)(handle); /* close the file */
                ffpmsg_str("failed get the size of the following file: (ffopen)");
                ffpmsg_slice(url);
                return *status;
            }

            let Fptr = FITSfile::new(&d[driver as usize], handle, url, cs!("ffopen"), status);
            if Fptr.is_err() {
                return *status;
            }
            let mut Fptr = Fptr.unwrap();

            /* initialize the ageindex array (relative age of the I/O buffers) */
            /* and initialize the bufrecnum array as being empty */
            for ii in 0..NIOBUF as usize {
                Fptr.ageindex[ii] = ii as c_int;
                Fptr.bufrecnum[ii] = -1;
            }

            /* store the parameters describing the file */
            Fptr.MAXHDU = 1000; /* initial size of headstart */
            Fptr.filehandle = handle; /* file handle */
            Fptr.driver = driver; /* driver number */
            strcpy(Fptr.filename, url.as_ptr()); /* full input filename */
            Fptr.filesize = filesize as LONGLONG; /* physical file size */
            Fptr.logfilesize = filesize as LONGLONG; /* logical file size */
            Fptr.writemode = mode; /* read-write mode    */
            Fptr.datastart = DATA_UNDEFINED as LONGLONG; /* unknown start of data */
            Fptr.curbuf = -1; /* undefined current IO buffer */
            Fptr.open_count = 1; /* structure is currently used once */
            Fptr.validcode = VALIDSTRUC; /* flag denoting valid structure */
            Fptr.only_one = only_one; /* flag denoting only copy single extension */
            Fptr.noextsyntax = if open_disk_file { 1 } else { 0 }; /* true if extended syntax is disabled */

            // HEAP ALLOCATION
            /* allocate fitsfile structure and initialize = 0 */
            let f_fitsfile = Box::try_new(fitsfile {
                HDUposition: 0,
                Fptr,
            });

            if f_fitsfile.is_err() {
                let d = DRIVER_TABLE.get().unwrap();
                ((d[driver as usize]).close)(handle); /* close the file */
                ffpmsg_str("failed to allocate structure for following file: (ffopen)");
                ffpmsg_slice(url);
                *status = MEMORY_ALLOCATION;
                return *status;
            }

            let mut f_fitsfile = f_fitsfile.unwrap();

            //drop(d); // To break mutex guard

            ffldrc(&mut f_fitsfile, 0, REPORT_EOF, status); /* load first record */

            fits_store_Fptr(&mut f_fitsfile.Fptr, status); /* store Fptr address */

            if ffrhdu_safer(&mut f_fitsfile, Some(&mut hdutyp), status) > 0 {
                /* determine HDU structure */
                ffpmsg_str("ffopen could not interpret primary array header of file: ");
                ffpmsg_slice(url);

                if *status == UNKNOWN_REC {
                    ffpmsg_str("This does not look like a FITS file.");
                }

                ffclos_safer(f_fitsfile, status);
                *fptr = None; /* return null file pointer */
                return *status;
            }

            *fptr = Some(f_fitsfile);

            /* ------------------------------------------------------------- */
            /* At this point, the input file has been opened. If outfile was */
            /* specified, then we have opened a copy of the file, not the    */
            /* original file so it is safe to modify it if necessary         */
            /* ------------------------------------------------------------- */

            if outfile[0] != 0 {
                writecopy = true;
            }
        }

        if extspec[0] != 0 {
            let f = (*fptr).as_mut().unwrap();

            if extnum != 0 {
                /* extension number was specified */

                ffmahd_safe(f, extnum + 1, Some(&mut hdutyp), status);
            } else if extname[0] != 0 {
                /* move to named extension, if specified */
                ffmnhd_safe(f, movetotype, &extname, extvers, status);
            }

            if *status > 0 {
                /* clean up after error */
                ffpmsg_str("ffopen could not move to the specified extension:");

                if extnum > 0 {
                    int_snprintf!(
                        &mut errmsg,
                        FLEN_ERRMSG,
                        " extension number {} doesn't exist or couldn't be opened.",
                        extnum,
                    );
                    ffpmsg_slice(&errmsg);
                } else {
                    int_snprintf!(
                        &mut errmsg,
                        FLEN_ERRMSG,
                        " extension with EXTNAME = {},",
                        slice_to_str!(&extname),
                    );
                    ffpmsg_slice(&errmsg);

                    if extvers != 0 {
                        int_snprintf!(
                            &mut errmsg,
                            FLEN_ERRMSG,
                            "           and with EXTVERS = {},",
                            extvers,
                        );
                        ffpmsg_slice(&errmsg);
                    }

                    if movetotype != ANY_HDU {
                        int_snprintf!(
                            &mut errmsg,
                            FLEN_ERRMSG,
                            "           and with XTENSION = {},",
                            CStr::from_ptr(hdtype[movetotype as usize])
                                .to_str()
                                .unwrap(),
                        );
                        ffpmsg_slice(&errmsg);
                    }
                    ffpmsg_str(" doesn't exist or couldn't be opened.");
                }

                let f_tmp = fptr.take().unwrap(); // Nulls pointer
                ffclos_safer(f_tmp, status);
                return *status;
            };
        } else if skip_null
            || skip_image
            || skip_table
            || (imagecolname[0] != 0 || colspec[0] != 0 || rowfilter[0] != 0 || binspec[0] != 0)
        {
            /* ------------------------------------------------------------------

            If no explicit extension specifier is given as part of the file
            name, and, if a) skip_null is true (set if ffopen is called by
            ffdopn) or b) skip_image or skip_table is true (set if ffopen is
            called by fftopn or ffdopn) or c) other file filters are
            specified, then CFITSIO will attempt to move to the first
            'interesting' HDU after opening an existing FITS file (or to
            first interesting table HDU if skip_image is true);

            An 'interesting' HDU is defined to be either an image with NAXIS
            > 0 (i.e., not a null array) or a table which has an EXTNAME
            value which does not contain any of the following strings:
               'GTI'  - Good Time Interval extension
               'OBSTABLE'  - used in Beppo SAX data files

            The main purpose for this is to allow CFITSIO to skip over a null
            primary and other non-interesting HDUs when opening an existing
            file, and move directly to the first extension that contains
            significant data.
            ------------------------------------------------------------------ */

            let f = (*fptr).as_mut().unwrap();

            ffghdn_safe(f, &mut hdunum);
            if hdunum == 1 {
                ffgidm_safe(f, &mut naxis, status);

                if naxis != 0 {
                    let mut naxes = vec![0; naxis as usize];
                    ffgisz_safe(f, naxis, &mut naxes, status);

                    for ii in 0..naxis as usize {
                        if naxes[ii] == 0 {
                            if ii == 0 {
                                /* NAXIS1=0 could be a random group indicator */
                                tstatus = 0;
                                ffmaky_safe(f, 2, status);

                                let mut group_val = 0;
                                if ffgkyl_safe(f, cs!("GROUPS"), &mut group_val, None, &mut tstatus)
                                    != 0
                                {
                                    no_primary_data = true; /* GROUPS keyword not found */
                                }
                            } else {
                                no_primary_data = true;
                            }
                        }
                    }
                } else {
                    no_primary_data = true;
                }

                if no_primary_data || skip_image {
                    /* skip primary array */

                    loop {
                        /* see if the next HDU is 'interesting' */
                        if ffmrhd_safe(f, 1, Some(&mut hdutyp), status) != 0 {
                            if *status == END_OF_FILE {
                                *status = 0; /* reset expected error */
                            }

                            /* didn't find an interesting HDU so move back to beginning */
                            ffmahd_safe(f, 1, Some(&mut hdutyp), status);
                            break;
                        }

                        if hdutyp == IMAGE_HDU && skip_image {
                            continue; /* skip images */
                        } else if hdutyp != IMAGE_HDU && skip_table {
                            continue; /* skip tables */
                        } else if hdutyp == IMAGE_HDU {
                            ffgidm_safe(f, &mut naxis, status);
                            if naxis > 0 {
                                break; /* found a non-null image */
                            };
                        } else {
                            tstatus = 0;
                            tblname[0] = 0;
                            fits_read_key_str(f, cs!(b"EXTNAME"), &mut tblname, None, &mut tstatus);

                            if (strstr_safe(&tblname, cs!(b"GTI")).is_none()
                                && strstr_safe(&tblname, cs!(b"gti")).is_none())
                                && fits_strncasecmp(&tblname, cs!(b"OBSTABLE"), 8) != 0
                            {
                                break; /* found an interesting table */
                            };
                        };
                    } /* end loop */
                }
            } /* end if (hdunum==1) */
        }

        if imagecolname[0] != 0 {
            /* ----------------------------------------------------------------- */
            /* we need to open an image contained in a single table cell         */
            /* First, determine which row of the table to use.                   */
            /* ----------------------------------------------------------------- */
            let f = (*fptr).as_mut().unwrap();

            if !isdigit_safe(rowexpress[0]) {
                /* is the row specification a number? */
                sscanf(rowexpress.as_ptr(), cstr!(b"%ld").as_ptr(), &mut rownum);
                if rownum < 1 {
                    ffpmsg_str("illegal rownum for image cell:");
                    ffpmsg_slice(&rowexpress);
                    ffpmsg_str("Could not open the following image in a table cell:");
                    ffpmsg_slice(&extspec);

                    let f_tmp = fptr.take().unwrap(); // Nulls fptr
                    ffclos_safer(f_tmp, status);
                    *status = BAD_ROW_NUM;
                    return *status;
                };
            } else if ffffrw_safer(f, &rowexpress, &mut rownum, status) > 0 {
                ffpmsg_str("Failed to find row matching this expression:");
                ffpmsg_slice(&rowexpress);
                ffpmsg_str("Could not open the following image in a table cell:");
                ffpmsg_slice(&extspec);
                let f_tmp = fptr.take().unwrap(); // Nulls fptr
                ffclos_safer(f_tmp, status);
                return *status;
            }
            if rownum == 0 {
                ffpmsg_str("row satisfying this expression doesn\'t exist::");
                ffpmsg_slice(&rowexpress);
                ffpmsg_str("Could not open the following image in a table cell:");
                ffpmsg_slice(&extspec);
                let f_tmp = fptr.take().unwrap(); // Nulls fptr
                ffclos_safer(f_tmp, status);
                *status = BAD_ROW_NUM;
                return *status;
            }

            /* determine the name of the new file to contain copy of the image */
            if histfilename[0] != 0 && (pixfilter[0]) == 0 {
                strcpy_safe(&mut outfile, &histfilename); /* the original outfile name */
            } else {
                strcpy_safe(&mut outfile, cs!(b"mem://_1")); /* create image file in memory */
            }

            /* Copy the image into new primary array and open it as the current */
            /* fptr.  This will close the table that contains the original image. */

            /* create new empty file to hold copy of the image */
            if ffinit_safer(&mut newptr, &outfile, status) > 0 {
                ffpmsg_str("failed to create file for copy of image in table cell:");
                ffpmsg_slice(&outfile);
                return *status;
            }

            if fits_copy_cell2image_safe(f, newptr.as_mut().unwrap(), &imagecolname, rownum, status)
                > 0
            {
                ffpmsg_str("Failed to copy table cell to new primary array:");
                ffpmsg_slice(&extspec);
                let f_tmp = fptr.take().unwrap(); // Nulls fptr
                ffclos_safer(f_tmp, status);
                return *status;
            }

            /* close the original file and set fptr to the new image */
            let f_tmp = fptr.take().unwrap(); // Nulls fptr
            ffclos_safer(f_tmp, status);

            *fptr = newptr; /* reset the pointer to the new table */

            writecopy = true; /* we are now dealing with a copy of the original file */

            /*  leave it up to calling routine to write any HISTORY keywords */
        }

        /* --------------------------------------------------------------------- */
        /* edit columns (and/or keywords) in the table, if specified in the URL  */
        /* --------------------------------------------------------------------- */

        if colspec[0] != 0 {
            /* the column specifier will modify the file, so make sure */
            /* we are already dealing with a copy, or else make a new copy */

            if !writecopy {
                /* Is the current file already a copy? */
                writecopy = fits_is_this_a_copy(urltype) != 0;
            }

            if !writecopy {
                if filtfilename[0] != 0 && outfile[0] == 0 {
                    strcpy_safe(&mut outfile, &filtfilename); /* the original outfile name */
                } else {
                    strcpy_safe(&mut outfile, cs!(b"mem://_1")); /* will create copy in memory */
                }

                writecopy = true;
            } else {
                let f = (*fptr).as_mut().unwrap();

                f.Fptr.writemode = READWRITE; /* we have write access */
                outfile[0] = 0;
            }

            if ffedit_columns(fptr, &outfile, &mut colspec, status) > 0 {
                ffpmsg_str("editing columns in input table failed (ffopen)");
                ffpmsg_str(" while trying to perform the following operation:");
                ffpmsg_slice(&colspec);
                let f_tmp = fptr.take().unwrap(); // Nulls fptr
                ffclos_safer(f_tmp, status);
                return *status;
            };
        }

        /* ------------------------------------------------------------------- */
        /* select rows from the table, if specified in the URL                 */
        /* or select a subimage (if this is an image HDU and not a table)      */
        /* ------------------------------------------------------------------- */

        if rowfilter[0] != 0 {
            let f = (*fptr).as_mut().unwrap();

            ffghdt_safe(f, &mut hdutyp, status); /* get type of HDU */
            if hdutyp == IMAGE_HDU {
                /* this is an image so 'rowfilter' is an image section specification */

                if filtfilename[0] != 0 && outfile[0] == 0 {
                    strcpy_safe(&mut outfile, &filtfilename); /* the original outfile name */
                } else if outfile[0] == 0 {
                    /* output file name not already defined? */
                    strcpy_safe(&mut outfile, cs!(b"mem://_2")); /* will create file in memory */
                }

                /* create new file containing the image section, plus a copy of */
                /* any other HDUs that exist in the input file.  This routine   */
                /* will close the original image file and return a pointer      */
                /* to the new file. */

                if fits_select_image_section_safer(fptr, &outfile, &rowfilter, status) > 0 {
                    ffpmsg_str("on-the-fly selection of image section failed (ffopen)");
                    ffpmsg_str(" while trying to use the following section filter:");
                    ffpmsg_slice(&rowfilter);

                    let f_tmp = fptr.take().unwrap(); // Nulls pointer
                    ffclos_safer(f_tmp, status);
                    return *status;
                };
            } else {
                /* this is a table HDU, so the rowfilter is really a row filter */

                if binspec[0] != 0 {
                    /*  since we are going to make a histogram of the selected rows,   */
                    /*  it would be a waste of time and memory to make a whole copy of */
                    /*  the selected rows.  Instead, just construct an array of TRUE   */
                    /*  or FALSE values that indicate which rows are to be included    */
                    /*  in the histogram and pass that to the histogram generating     */
                    /*  routine                                                        */

                    ffgnrw_safe(f, &mut nrows, status); /* get no. of rows */

                    // HEAP ALLOCATION - Temporary
                    let mut rowselect = Vec::new();
                    if rowselect.try_reserve_exact(nrows as usize).is_err() {
                        ffpmsg_str("failed to allocate memory for selected columns array (ffopen)");
                        ffpmsg_str(" while trying to select rows with the following filter:");
                        ffpmsg_slice(&rowfilter);

                        let f_tmp = fptr.take().unwrap(); // Nulls pointer
                        ffclos_safer(f_tmp, status);

                        *status = MEMORY_ALLOCATION;
                        return *status;
                    } else {
                        rowselect.resize(nrows as usize, 0);
                    }

                    if fffrow_safe(
                        f,
                        &rowfilter,
                        1,
                        nrows,
                        &mut goodrows,
                        &mut rowselect,
                        status,
                    ) > 0
                    {
                        ffpmsg_str("selection of rows in input table failed (ffopen)");
                        ffpmsg_str(" while trying to select rows with the following filter:");
                        ffpmsg_slice(&rowfilter);

                        let f_tmp = fptr.take().unwrap(); // Nulls pointer
                        ffclos_safer(f_tmp, status);

                        return *status;
                    };
                } else {
                    if !writecopy {
                        /* Is the current file already a copy? */
                        writecopy = fits_is_this_a_copy(urltype) != 0;
                    }

                    if !writecopy {
                        if filtfilename[0] != 0 && outfile[0] == 0 {
                            strcpy_safe(&mut outfile, &filtfilename); /* the original outfile name */
                        } else if outfile[0] == 0 {
                            /* output filename not already defined? */
                            strcpy_safe(&mut outfile, cs!(b"mem://_2")); /* will create copy in memory */
                        };
                    } else {
                        f.Fptr.writemode = READWRITE; /* we have write access */
                        outfile[0] = 0;
                    }

                    /* select rows in the table.  If a copy of the input file has */
                    /* not already been made, then this routine will make a copy */
                    /* and then close the input file, so that the modifications will */
                    /* only be made on the copy, not the original */

                    if ffselect_table(fptr, &outfile, &rowfilter, status) > 0 {
                        ffpmsg_str("on-the-fly selection of rows in input table failed (ffopen)");
                        ffpmsg_str(" while trying to select rows with the following filter:");
                        ffpmsg_slice(&rowfilter);
                        let f_tmp = fptr.take().unwrap(); // Nulls pointer
                        ffclos_safer(f_tmp, status);
                        return *status;
                    }

                    let f = (*fptr).as_mut().unwrap();

                    /* write history records */
                    ffphis_safe(
                    f,
                    cs!(b"CFITSIO used the following filtering expression to create this table:"),
                    status,
                );
                    ffphis_safe(f, name, status);
                } /* end of no binspec case */
            } /* end of table HDU case */
        } /* end of rowfilter exists case */

        /* ------------------------------------------------------------------- */
        /* make an image histogram by binning columns, if specified in the URL */
        /* ------------------------------------------------------------------- */

        if binspec[0] != 0 {
            let mut exprs = None; // This is a valid pointer to an empty value, i.e. is expecting a return value
            if histfilename[0] != 0 && (pixfilter[0]) == 0 {
                strcpy_safe(&mut outfile, &histfilename); /* the original outfile name */
            } else {
                strcpy_safe(&mut outfile, cs!(b"mem://_3")); /* create histogram in memory */
                /* if not already copied the file */
            }

            /* parse the binning specifier into individual parameters */
            ffbinse(
                &binspec,
                &mut imagetype,
                &mut haxis,
                &colname,
                &minin,
                &maxin,
                &binsizein,
                &minname,
                &maxname,
                &binname,
                &mut weight,
                &wtcol,
                &mut recip,
                Some(&mut exprs),
                status,
            );

            /* Create the histogram primary array and open it as the current fptr */
            /* This will close the table that was used to create the histogram. */

            match exprs {
                Some(e) => {
                    ffhist2e(
                        fptr,
                        &outfile,
                        imagetype,
                        haxis,
                        &colname,
                        Some(&[&e[0], &e[1], &e[2], &e[3]]),
                        &minin,
                        &maxin,
                        &binsizein,
                        Some(&minname),
                        Some(&maxname),
                        Some(&binname),
                        weight,
                        Some(&wtcol),
                        Some(&e[4]),
                        recip,
                        rowselect,
                        status,
                    );
                }
                None => {
                    ffhist2e(
                        fptr,
                        &outfile,
                        imagetype,
                        haxis,
                        &colname,
                        None,
                        &minin,
                        &maxin,
                        &binsizein,
                        Some(&minname),
                        Some(&maxname),
                        Some(&binname),
                        weight,
                        Some(&wtcol),
                        None,
                        recip,
                        rowselect,
                        status,
                    );
                }
            };

            let f = (*fptr).as_mut().unwrap();
            if *status > 0 {
                ffpmsg_str("on-the-fly histogramming of input table failed (ffopen)");
                ffpmsg_str(" while trying to execute the following histogram specification:");
                ffpmsg_slice(&binspec);

                let f_tmp = fptr.take().unwrap(); // Nulls pointer
                ffclos_safer(f_tmp, status);

                return *status;
            }

            /* write history records */
            ffphis_safe(
                f,
                cs!(b"CFITSIO used the following expression to create this histogram:"),
                status,
            );

            ffphis_safe(f, name, status);
        }

        if pixfilter[0] != 0 {
            let f = (*fptr).as_mut().unwrap();
            if histfilename[0] != 0 {
                strcpy_safe(&mut outfile, &histfilename); /* the original outfile name */
            } else {
                strcpy_safe(&mut outfile, cs!(b"mem://_4")); /* create in memory */
                /* if not already copied the file */
            }

            /* Ensure type of HDU is consistent with pixel filtering */
            ffghdt_safe(f, &mut hdutyp, status);

            if hdutyp == IMAGE_HDU {
                pixel_filter_helper(fptr, outfile, pixfilter, status);

                if *status > 0 {
                    ffpmsg_str("pixel filtering of input image failed (ffopen)");
                    ffpmsg_str(" while trying to execute the following:");
                    ffpmsg_slice(&pixfilter);

                    let f_tmp = fptr.take().unwrap(); // Nulls pointer
                    ffclos_safer(f_tmp, status);

                    return *status;
                }

                let f = (*fptr).as_mut().unwrap();

                /* write history records */
                ffphis_safe(
                    f,
                    cs!(b"CFITSIO used the following expression to create this image:"),
                    status,
                );

                ffphis_safe(f, name, status);
            } else {
                let f = (*fptr).as_mut().unwrap();
                ffpmsg_str("cannot use pixel filter on non-IMAGE HDU");
                ffpmsg_slice(&pixfilter);

                let f_tmp = fptr.take().unwrap(); // Nulls pointer
                ffclos_safer(f_tmp, status);
                *status = NOT_IMAGE;
                return *status;
            };
        }

        /* parse and save image compression specification, if given */
        if compspec[0] != 0 {
            let f = (*fptr).as_mut().unwrap();
            ffparsecompspec(f, &compspec, status);
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Reopen an existing FITS file with either readonly or read/write access.
///
/// The reopened file shares the same FITSfile structure but may point to a
/// different HDU within the file.
/// SATEFY: This is in no ways safe, multiple fptrs sharing the same underlying data
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffreopen(
    openfptr: *mut fitsfile,     /* I - FITS file pointer to open file  */
    newfptr: *mut *mut fitsfile, /* O - pointer to new re opened file   */
    status: *mut c_int,          /* IO - error status                   */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        /* check that the open file pointer is valid */
        if openfptr.is_null() {
            *status = NULL_INPUT_PTR;
            return *status;
        }

        let openfptr = openfptr.as_mut().expect(NULL_MSG);

        if openfptr.Fptr.validcode != VALIDSTRUC {
            /* check magic value */
            *status = BAD_FILEPTR;
            return *status;
        }

        /* allocate fitsfile structure and initialize = 0 */
        // HEAP ALLOCATION
        let mut n = Box::new(fitsfile {
            HDUposition: 0, /* set initial position to primary array */
            Fptr: Box::from_raw(&mut *openfptr.Fptr as *mut FITSfile), /* both point to the same structure */ // TODO this is very unsafe!
        });

        n.Fptr.open_count += 1; /* increment the file usage counter */

        *newfptr = Box::into_raw(n);

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// store the new Fptr address for future use by fits_already_open
pub(crate) unsafe fn fits_store_Fptr(
    Fptr: &mut FITSfile, /* O - FITS file pointer               */
    status: &mut c_int,  /* IO - error status                   */
) -> c_int {
    unsafe {
        if *status > 0 {
            return *status;
        }

        let lock = FFLOCK();
        for ii in 0..NMAXFILES {
            if FPTR_TABLE[ii].is_null() {
                FPTR_TABLE[ii] = Fptr;
                break;
            };
        }
        FFUNLOCK(lock);

        *status
    }
}

/*--------------------------------------------------------------------------*/
///  clear the Fptr address from the Fptr Table  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_clear_Fptr(
    Fptr: *mut FITSfile, /* O - FITS file pointer               */
    status: *mut c_int,  /* IO - error status                   */
) -> c_int {
    unsafe {
        let lock = FFLOCK();
        for ii in 0..NMAXFILES {
            if std::ptr::eq(FPTR_TABLE[ii], Fptr) {
                FPTR_TABLE[ii] = ptr::null_mut();
                break;
            };
        }
        FFUNLOCK(lock);
        *status
    }
}

/*--------------------------------------------------------------------------*/
///  clear the Fptr address from the Fptr Table  
fn fits_clear_Fptr_safer(
    Fptr: &mut FITSfile, /* O - FITS file pointer               */
    status: &mut c_int,  /* IO - error status                   */
) -> c_int {
    let lock = FFLOCK();
    for ii in 0..NMAXFILES {
        unsafe {
            if FPTR_TABLE[ii] == Fptr {
                FPTR_TABLE[ii] = ptr::null_mut();
                break;
            };
        }
    }
    FFUNLOCK(lock);
    *status
}

/*--------------------------------------------------------------------------*/
/// Check if the file to be opened is already open.  If so, then attach to it.
///
/// the input strings must not exceed the standard lengths
/// of FLEN_FILENAME, MAX_PREFIX_LEN, etc.
///
///      
/// this function was changed so that for files of access method FILE://
/// the file paths are compared using standard URL syntax and absolute
/// paths (as opposed to relative paths). This eliminates some instances
/// where a file is already opened but it is not realized because it
/// was opened with another file path. For instance, if the CWD is
/// /a/b/c and I open /a/b/c/foo.fits then open ./foo.fits the previous
/// version of this function would not have reconized that the two files
/// were the same. This version does recognize that the two files are
/// the same.
pub(crate) unsafe fn fits_already_open(
    fptr: &mut Option<Box<fitsfile>>, /* I/O - FITS file pointer       */
    url: &[c_char],
    urltype: &mut [c_char],
    infile: &mut [c_char],
    extspec: &mut [c_char],
    rowfilter: &mut [c_char],
    binspec: &mut [c_char],
    colspec: &mut [c_char],
    mode: c_int,        /* I - 0 = open readonly; 1 = read/write   */
    noextsyn: bool,     /* I - 0 = ext syntax may be used; 1 = ext syntax disabled */
    isopen: &mut c_int, /* O - 1 = file is already open            */
    status: &mut c_int, /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut oldFptr: &mut FITSfile;
        let mut iMatch = None;
        let mut oldurltype: [c_char; MAX_PREFIX_LEN] = [0; MAX_PREFIX_LEN];
        let mut oldinfile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut oldextspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut oldoutfile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut oldrowfilter: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut oldbinspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut oldcolspec: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let cwd: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let tmpStr: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut tmpinfile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];

        *isopen = 0;

        /*  When opening a file with readonly access then we simply let
            the operating system open the file again, instead of using the CFITSIO
            trick of attaching to the previously opened file.  This is required
            if CFITSIO is running in a multi-threaded environment, because 2 different
            threads cannot share the same FITSfile pointer.

            If the file is opened/reopened with write access, then the file MUST
            only be physically opened once..
        */

        if mode == 0 {
            return *status;
        }

        strcpy_safe(&mut tmpinfile, infile);
        if fits_strcasecmp(urltype, cs!(b"FILE://")) == 0
            && standardize_path(&mut tmpinfile, status) != 0
        {
            return *status;
        };
        for ii in 0..NMAXFILES {
            /* check every buffer */

            if !FPTR_TABLE[ii].is_null() {
                oldFptr = FPTR_TABLE[ii].as_mut().expect(NULL_MSG);

                if oldFptr.noextsyntax != 0 {
                    /* old urltype must be "file://" */
                    if fits_strcasecmp(urltype, cs!(b"FILE://")) == 0 {
                        /* compare tmpinfile to adjusted oldFptr->filename */

                        let ofn = cast_slice(CStr::from_ptr(oldFptr.filename).to_bytes_with_nul());

                        /* This shouldn't be possible, but check anyway */
                        if strlen_safe(ofn) > FLEN_FILENAME - 1 {
                            ffpmsg_str("Name of old file is too long. (fits_already_open)");
                            *status = FILE_NOT_OPENED;
                            return *status;
                        }

                        let ofn = cast_slice(CStr::from_ptr(oldFptr.filename).to_bytes_with_nul());
                        strcpy_safe(&mut oldinfile, ofn);
                        if standardize_path(&mut oldinfile, status) != 0 {
                            return *status;
                        }

                        if strcmp_safe(&tmpinfile, &oldinfile) == 0 {
                            /* if infile is not noextsyn, must check that it is not
                            using filters of any kind */
                            if noextsyn || (rowfilter[0] == 0 && binspec[0] == 0 && colspec[0] == 0)
                            {
                                if mode == READWRITE as c_int
                                    && oldFptr.writemode == READONLY as c_int
                                {
                                    /*
                                      cannot assume that a file previously opened with READONLY
                                      can now be written to (e.g., files on CDROM, or over the
                                      the network, or STDIN), so return with an error.
                                    */

                                    ffpmsg_str(
                                        "cannot reopen file READWRITE when previously opened READONLY",
                                    );
                                    ffpmsg_slice(url);
                                    *status = FILE_NOT_OPENED;
                                    return *status;
                                }
                                iMatch = Some(ii);
                            };
                        };
                    };

                    /* end if old file has disabled extended syntax */
                } else {
                    let filename = cast_slice(CStr::from_ptr(oldFptr.filename).to_bytes_with_nul());

                    ffiurl_safer(
                        filename,
                        oldurltype.as_mut_ptr(),
                        oldinfile.as_mut_ptr(),
                        oldoutfile.as_mut_ptr(),
                        oldextspec.as_mut_ptr(),
                        oldrowfilter.as_mut_ptr(),
                        oldbinspec.as_mut_ptr(),
                        oldcolspec.as_mut_ptr(),
                        status,
                    );

                    if *status > 0 {
                        ffpmsg_str("could not parse the previously opened filename: (ffopen)");
                        ffpmsg_cstr(CStr::from_ptr(oldFptr.filename));
                        return *status;
                    }

                    if fits_strcasecmp(&oldurltype, cs!(b"FILE://")) == 0
                        && standardize_path(&mut oldinfile, status) != 0
                    {
                        return *status;
                    };

                    if strcmp_safe(urltype, &oldurltype) == 0
                        && strcmp_safe(&tmpinfile, &oldinfile) == 0
                    {
                        /* identical type of file and root file name */

                        if (rowfilter[0] == 0
                        && oldrowfilter[0] == 0
                        && binspec[0] == 0
                        && oldbinspec[0] == 0
                        && colspec[0] == 0
                        && oldcolspec[0] == 0)

                     /* no filtering or binning specs for either file, so */
                     /* this is a case where the same file is being reopened. */
                     /* It doesn't matter if the extensions are different */


                        || (strcmp_safe(rowfilter, &oldrowfilter) == 0
                            && strcmp_safe(binspec, &oldbinspec) == 0
                            && strcmp_safe(colspec, &oldcolspec) == 0
                            && strcmp_safe(extspec, &oldextspec) == 0)
                        /* filtering specs are given and are identical, and */
                        /* the same extension is specified */
                        {
                            if mode == READWRITE as c_int && oldFptr.writemode == READONLY as c_int
                            {
                                /*
                                  cannot assume that a file previously opened with READONLY
                                  can now be written to (e.g., files on CDROM, or over the
                                  the network, or STDIN), so return with an error.
                                */

                                ffpmsg_str(
                                    "cannot reopen file READWRITE when previously opened READONLY",
                                );
                                ffpmsg_slice(url);
                                *status = FILE_NOT_OPENED;
                                return *status;
                            }
                            iMatch = Some(ii);
                        }
                    }
                } /* end if old file recognizes extended syntax */
            } /* end if old fptr exists */
        } /* end loop over NMAXFILES */

        if iMatch.is_some() {
            let iMatch = iMatch.unwrap();
            let oldFptr = FPTR_TABLE[iMatch];

            // HEAP ALLOCATION
            let f = Box::try_new(fitsfile {
                HDUposition: 0,               /* set initial position */
                Fptr: Box::from_raw(oldFptr), /* point to the structure */
            });

            if f.is_err() {
                ffpmsg_str("failed to allocate structure for following file: (ffopen)");
                ffpmsg_slice(url);
                *status = MEMORY_ALLOCATION;
                return *status;
            }

            let mut f = f.unwrap();

            (f.Fptr.open_count) += 1; /* increment usage counter */

            *fptr = Some(f);

            if binspec[0] != 0 {
                /* if binning specified, don't move */
                extspec[0] = 0;
            }

            /* all the filtering has already been applied, so ignore */
            rowfilter[0] = 0;
            binspec[0] = 0;
            colspec[0] = 0;
            *isopen = 1;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Utility function for common operation in fits_already_open
/// fullpath:  I/O string to be standardized. Assume len = FLEN_FILENAME
fn standardize_path(fullpath: &mut [c_char], status: &mut c_int) -> c_int {
    let mut tmpPath: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
    let mut cwd: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];

    if fits_path2url(fullpath, FLEN_FILENAME, &mut tmpPath, status) > 0 {
        return *status;
    }

    if tmpPath[0] != bb(b'/') {
        fits_get_cwd(&mut cwd, status);
        if strlen_safe(&cwd) + strlen_safe(&tmpPath) + 1 > FLEN_FILENAME - 1 {
            ffpmsg_str("Tile name is too long. (standardize_path)");
            *status = FILE_NOT_OPENED;
            return *status;
        }
        strcat_safe(&mut cwd, cs!("/"));
        strcat_safe(&mut cwd, &tmpPath);
        fits_clean_url(&cwd, &mut tmpPath, status);
    }

    strcpy_safe(fullpath, &tmpPath);

    *status
}

/*--------------------------------------------------------------------------*/
/// specialized routine that returns 1 if the file is known to be a temporary
/// copy of the originally opened file.  Otherwise it returns 0.
pub(crate) fn fits_is_this_a_copy(urltype: [c_char; 20] /* I - type of file */) -> c_int {
    let mut iscopy: c_int = 0;

    if strncmp_safe(&urltype, cs!("mem"), 3) == 0 {
        iscopy = 1; /* file copy is in memory */
    } else if strncmp_safe(&urltype, cs!("compress"), 8) == 0 {
        iscopy = 1; /* compressed diskfile that is uncompressed in memory */
    } else if strncmp_safe(&urltype, cs!("http"), 4) == 0 {
        iscopy = 1; /* copied file using http protocol */
    } else if strncmp_safe(&urltype, cs!("ftp"), 3) == 0 {
        iscopy = 1; /* copied file using ftp protocol */
    } else if strncmp_safe(&urltype, cs!("gsiftp"), 6) == 0 {
        iscopy = 1; /* copied file using gsiftp protocol */
    } else if strncmp_safe(&urltype, cs!("stdin"), 5) == 0 {
        //NOTE looks like a bug in original code
        iscopy = 1; /* piped stdin has been copied to memory */
    } else {
        iscopy = 0; /* file is not known to be a copy */
    }

    iscopy
}

/*--------------------------------------------------------------------------*/
/// Look for the closing single quote character in the input string
fn find_quote(string: &[c_char]) -> Option<usize> {
    let mut i = 0;
    let len = string.len();

    while i < len && string[i] != 0 {
        if string[i] == bb(b'\'') {
            /* found the closing quote */
            return Some(i + 1); /* set pointer to next char */
        }
        i += 1;
    }

    None /* opps, didn't find the closing character */
}

/*--------------------------------------------------------------------------*/
/*
  Find matching delimiter, respecting quoting and (potentially nested) parentheses

  char *string - null-terminated string to be searched for delimiter
  char delim - single delimiter to search for (one of '")]} )

  returns: pointer to character after delimiter, or 0 if not found
*/
fn fits_find_match_delim(string: &mut [c_char], delim: c_char) -> Option<usize> {
    if string.is_empty() {
        return None;
    }

    match delim as u8 {
        b'\'' => find_quote(string),
        b'"' => find_doublequote(string),
        b'}' => find_curlybracket(string),
        b']' => find_bracket(string),
        b')' => find_paren(string),
        _ => None, // Invalid delimiter, return failure
    }
}

/*--------------------------------------------------------------------------*/
/// Look for the closing double quote character in the input string
fn find_doublequote(string: &[c_char]) -> Option<usize> {
    let mut i = 0;
    let len = string.len();

    while i < len && string[i] != 0 {
        if string[i] == bb(b'"') {
            /* found the closing quote */
            return Some(i + 1); /* set pointer to next char */
        }
        i += 1;
    }

    None /* opps, didn't find the closing character */
}

/*--------------------------------------------------------------------------*/
/// look for the closing parenthesis character in the input string
fn find_paren(string: &[c_char]) -> Option<usize> {
    let mut i = 0;
    let len = string.len();

    while i < len && string[i] != 0 {
        if string[i] == bb(b')') {
            /* found the closing parens */
            return Some(i + 1); /* set pointer to next char */
        } else if string[i] == bb(b'(') {
            /* found another level of parens */
            i += 1;
            let p = find_paren(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'[') {
            /* found another level of parens */
            i += 1;
            let p = find_bracket(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'{') {
            /* found another level of parens */
            i += 1;
            let p = find_curlybracket(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'"') {
            /* found another level of parens */
            i += 1;
            let p = find_doublequote(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'\'') {
            /* found another level of parens */
            i += 1;
            let p = find_quote(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else {
            i += 1;
        }
    }

    None /* opps, didn't find the closing character */
}

/*--------------------------------------------------------------------------*/
/// look for the closing bracket character in the input string
fn find_bracket(string: &[c_char]) -> Option<usize> {
    let mut i = 0;
    let len = string.len();

    while i < len && string[i] != 0 {
        if string[i] == bb(b']') {
            /* found the closing bracket */
            return Some(i + 1); /* set pointer to next char */
        } else if string[i] == bb(b'(') {
            /* found another level of parens */
            i += 1;
            let p = find_paren(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'[') {
            /* found another level of parens */
            i += 1;
            let p = find_bracket(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'{') {
            /* found another level of parens */
            i += 1;
            let p = find_curlybracket(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'"') {
            /* found another level of parens */
            i += 1;
            let p = find_doublequote(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'\'') {
            /* found another level of parens */
            i += 1;
            let p = find_quote(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else {
            i += 1;
        }
    }

    None /* opps, didn't find the closing character */
}

/*--------------------------------------------------------------------------*/
/// look for the closing curly bracket character in the input string
fn find_curlybracket(string: &[c_char]) -> Option<usize> {
    let mut i = 0;
    let len = string.len();

    while i < len && string[i] != 0 {
        if string[i] == bb(b'}') {
            /* found the closing curly bracket */
            return Some(i + 1); /* set pointer to next char */
        } else if string[i] == bb(b'(') {
            /* found another level of parens */
            i += 1;
            let p = find_paren(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'[') {
            /* found another level of parens */
            i += 1;
            let p = find_bracket(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'{') {
            /* found another level of parens */
            i += 1;
            let p = find_curlybracket(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'"') {
            /* found another level of parens */
            i += 1;
            let p = find_doublequote(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'\'') {
            /* found another level of parens */
            i += 1;
            let p = find_quote(&string[i..]);
            match p {
                None => return None,
                Some(x) => i += x,
            }
        } else {
            i += 1;
        }
    }

    None /* opps, didn't find the closing character */
}

/*--------------------------------------------------------------------------*/
/// replace commas with semicolons, unless the comma is within a quoted or bracketed expression
fn comma2semicolon(string: &mut [c_char]) -> c_int {
    let mut i = 0;
    let len = string.len();

    while i < len && string[i] != 0 {
        if string[i] == bb(b',') {
            /* found a comma */
            string[i] = bb(b';');
            i += 1;
        } else if string[i] == bb(b'(') {
            /* found another level of parens */
            i += 1;
            let p = find_paren(&string[i..]);

            match p {
                None => return 1,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'[') {
            i += 1;
            let p = find_bracket(&string[i..]);
            match p {
                None => return 1,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'{') {
            i += 1;
            let p = find_curlybracket(&string[i..]);
            match p {
                None => return 1,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'"') {
            i += 1;
            let p = find_doublequote(&string[i..]);
            match p {
                None => return 1,
                Some(x) => i += x,
            }
        } else if string[i] == bb(b'\'') {
            i += 1;
            let p = find_quote(&string[i..]);
            match p {
                None => return 1,
                Some(x) => i += x,
            }
        } else {
            i += 1;
        }
    }

    0 /* reached end of string */
}

/*--------------------------------------------------------------------------*/
/// modify columns in a table and/or header keywords in the HDU
pub(crate) unsafe fn ffedit_columns(
    infptr: &mut Option<Box<fitsfile>>, /* IO - pointer to input table; on output it  */
    /*      points to the new selected rows table */
    outfile: &[c_char],  /* I - name for output file */
    expr: &mut [c_char], /* I - column edit expression    */
    status: &mut c_int,
) -> c_int {
    unsafe {
        let mut newptr: Option<Box<fitsfile>> = None;
        let ii: c_int = 0;
        let mut hdunum: c_int = 0;
        let mut slen: c_int = 0;
        let mut colnum: c_int = -1;
        let mut testnum: c_int = 0;
        let mut deletecol: c_int = 0 as c_int;
        let mut savecol: c_int = 0 as c_int;
        let mut numcols: c_int = 0 as c_int;

        let mut tstatus: c_int = 0 as c_int;

        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut colname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut oldname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut colformat: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

        let mut testname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

        let mut colindex: Vec<c_int> = Vec::new();

        let fptr: &mut Box<fitsfile> = (infptr).as_mut().unwrap();

        if outfile[0] != 0 {
            /* create new empty file in to hold the selected rows */
            if ffinit_safer(&mut newptr, outfile, status) > 0 {
                ffpmsg_str("failed to create file for copy (ffedit_columns)");
                return *status;
            }

            let mut newptr = newptr.expect(NULL_MSG);

            ffghdn_safe(fptr, &mut hdunum); /* current HDU number in input file */

            /* copy all HDUs to the output copy, if the 'only_one' flag is not set */
            if fptr.Fptr.only_one == 0 {
                let mut ii = 1;
                loop {
                    if ffmahd_safe(fptr, ii, None, status) > 0 {
                        break;
                    }

                    ffcopy_safer(fptr, &mut newptr, 0, status);
                    ii += 1;
                }

                if *status == END_OF_FILE {
                    *status = 0; /* got the expected EOF error; reset = 0  */
                } else if *status > 0 {
                    ffclos_safer(newptr, status);
                    ffpmsg_str("failed to copy all HDUs from input file (ffedit_columns)");
                    return *status;
                }
            } else {
                /* only copy the primary array and the designated table extension */
                ffmahd_safe(fptr, 1, None, status);
                ffcopy_safer(fptr, &mut newptr, 0, status);
                ffmahd_safe(fptr, hdunum, None, status);
                ffcopy_safer(fptr, &mut newptr, 0, status);
                if *status > 0 {
                    ffclos_safer(newptr, status);
                    ffpmsg_str("failed to copy all HDUs from input file (ffedit_columns)");
                    return *status;
                }
                hdunum = 2;
            }

            /* close the original file and return ptr to the new image */
            let f_tmp = infptr.take().unwrap();
            ffclos_safer(f_tmp, status);

            *infptr = Some(newptr); /* reset the pointer to the new table */

            let fptr: &mut Box<fitsfile> = (infptr).as_mut().unwrap();

            /* move back to the selected table HDU */
            if ffmahd_safe(fptr, hdunum, None, status) > 0 {
                ffpmsg_str("failed to copy the input file (ffedit_columns)");
                return *status;
            }
        }

        let fptr: &mut Box<fitsfile> = (infptr).as_mut().unwrap();

        /* remove the "col " from the beginning of the column edit expression */
        let mut cptr = &mut expr[4..]; // &expr

        while cptr[0] == bb(b' ') {
            cptr = &mut cptr[1..]; /* skip leading white space */
        }

        /* Check if need to import expression from a file */
        let mut file_expr = ptr::null_mut();
        if cptr[0] == bb(b'@') {
            if ffimport_file_safer(&cptr[1..], &mut file_expr, status) != 0 {
                return *status;
            }

            let _len = CStr::from_ptr(file_expr).to_bytes_with_nul().len();

            cptr = slice::from_raw_parts_mut(file_expr, _len);
            while cptr[0] == bb(b' ') {
                cptr = &mut cptr[1..]; /* skip leading white space... again */
            }
        }

        tstatus = 0;
        ffgncl_safe(fptr, &mut numcols, &mut tstatus); /* get initial # of cols */

        /* as of July 2012, the CFITSIO column filter syntax was modified */
        /* so that commas may be used to separate clauses, as well as semi-colons. */
        /* This was done because users cannot enter the semi-colon in the HEASARC's */
        /* Hera on-line data processing system for computer security reasons.  */
        /* Therefore, we must convert those commas back to semi-colons here, but we */
        /* must not convert any columns that occur within parenthesies.  */

        if comma2semicolon(cptr) != 0 {
            ffpmsg_str("parsing error in column filter expression");
            ffpmsg_slice(cptr);

            *status = PARSE_SYNTAX_ERR;
            return *status;
        }

        /* parse expression and get first clause, if more than 1 */
        let mut tmp_clause = None;
        let mut cptr_idx = 0;

        loop {
            slen = fits_get_token2_safe(
                cptr,
                &mut cptr_idx,
                cast_slice(cstr!(";").to_bytes_with_nul()),
                &mut tmp_clause,
                None,
                status,
            );

            if slen <= 0 {
                break;
            }

            let _clause_len = tmp_clause.as_deref().unwrap().len();
            let clause = tmp_clause.as_deref_mut().unwrap();
            if cptr[0] == bb(b';') {
                cptr = &mut cptr[1..];
            }
            clause[slen as usize] = 0;

            if clause[0] == bb(b'!') || clause[0] == bb(b'-') {
                let mut clause1 = 1;
                let mut clen = if clause[clause1] != 0 {
                    strlen_safe(&clause[clause1..])
                } else {
                    0
                };
                /* ===================================== */
                /* Case I. delete this column or keyword */
                /* ===================================== */

                /* Case Ia. delete column names with 0-or-more wildcard
                        -COLNAME+ - delete repeated columns with exact name
                    -COLNAM*+ - delete columns matching patterns
                */
                if *status == 0
                    && clen > 1
                    && clause[clause1] != bb(b'#')
                    && clause[clause1 + clen - 1] == bb(b'+')
                {
                    clause[clause1 + clen - 1] = 0;
                    clen -= 1;

                    /* Note that this is a delete 0 or more specification,
                    which means that no matching columns is not an error. */
                    loop {
                        let mut status_del = 0;

                        /* Have to set status=0 so we can reset the search at
                        start column.  Because we are deleting columns on
                        the fly here, we have to reset the search every
                        time. The only penalty here is execution time
                        because leaving *status == COL_NOT_UNIQUE is merely
                        an optimization for tables assuming the tables do
                        not change from one call to the next. (an
                        assumption broken in this loop) */
                        *status = 0;
                        ffgcno_safe(
                            fptr,
                            CASEINSEN as c_int,
                            &clause[clause1..],
                            &mut colnum,
                            status,
                        );
                        /* ffgcno returns COL_NOT_UNIQUE if there are multiple columns,
                        and COL_NOT_FOUND after the last column is found, and
                        COL_NOT_FOUND if no matches were found */
                        if *status != 0 && *status != COL_NOT_UNIQUE {
                            break;
                        }

                        if ffdcol_safe(fptr, colnum, &mut status_del) > 0 {
                            ffpmsg_str("failed to delete column in input file:");
                            ffpmsg_slice(clause);

                            *status = status_del;
                            return *status;
                        }
                        deletecol = 1; /* set flag that at least one col was deleted */
                        numcols -= 1;
                        if *status != COL_NOT_UNIQUE {
                            break;
                        }
                    }

                    *status = 0; /* No matches are still successful */
                    colnum = -1; /* Ignore the column we found */

                /* Case Ib. delete column names with wildcard or not
                        -COLNAME  - deleted exact column
                    -COLNAM*  - delete first column that matches pattern
                   Note no leading '#'
                */
                } else if clause[clause1] != 0
                    && clause[clause1] != bb(b'#')
                    && ((ffgcno_safe(
                        fptr,
                        CASEINSEN as c_int,
                        &clause[clause1..],
                        &mut colnum,
                        status,
                    ) <= 0)
                        || *status == COL_NOT_UNIQUE)
                {
                    /* a column with this name exists, so try to delete it */
                    *status = 0; /* Clear potential status=COL_NOT_UNIQUE */
                    if ffdcol_safe(fptr, colnum, status) > 0 {
                        ffpmsg_str("failed to delete column in input file:");
                        ffpmsg_slice(clause);

                        return *status;
                    }
                    deletecol = 1; /* set flag that at least one col was deleted */
                    numcols -= 1;
                    colnum = -1;
                }
                /* Case Ic. delete keyword(s)
                        -KEYNAME,#KEYNAME  - delete exact keyword (first match)
                    -KEYNAM*,#KEYNAM*  - delete first matching keyword
                    -KEYNAME+,-#KEYNAME+ - delete 0-or-more exact matches of exact keyword
                    -KEYNAM*+,-#KEYNAM*+ - delete 0-or-more wildcard matches
                   Note the preceding # is optional if no conflicting column name exists
                   and that wildcard patterns are described in "colfilter" section of
                   documentation.
                */
                else {
                    let mut delall = false;
                    let mut haswild = false;
                    ffcmsg_safe(); /* clear previous error message from ffgcno */
                    /* try deleting a keyword with this name */
                    *status = 0;
                    /* skip past leading '#' if any */
                    if clause[clause1] == bb(b'#') {
                        clause1 += 1;
                    }
                    clen = strlen_safe(&clause[clause1..]);

                    /* Repeat deletion of keyword if requested with trailing '+' */
                    if clen > 1 && clause[clause1 + clen - 1] == bb(b'+') {
                        delall = true;
                        clause[clause1 + clen - 1] = 0;
                    }
                    /* Determine if this pattern has wildcards */
                    if strchr_safe(&clause[clause1..], bb(b'?')).is_some()
                        || strchr_safe(&clause[clause1..], bb(b'*')).is_some()
                        || strchr_safe(&clause[clause1..], bb(b'#')).is_some()
                    {
                        haswild = true;
                    }

                    if haswild {
                        /* ffdkey() behaves differently if the pattern has a wildcard:
                        it only checks from the "current" header position to the end, and doesn't
                        check before the "current" header position.  Therefore, for the
                        case of wildcards we will have to reset to the beginning. */
                        ffmaky_safe(fptr, 1, status); /* reset pointer to beginning of header */
                    }

                    /* Single or repeated deletions until done */
                    loop {
                        if ffdkey_safe(fptr, &clause[clause1..], status) > 0 {
                            if delall && *status == KEY_NO_EXIST {
                                /* Found last wildcard item. Stop deleting */
                                ffcmsg_safe();
                                *status = 0;
                                delall = false; /* Force end of this loop */
                            } else {
                                /* This was not a wildcard deletion, or it resulted in
                                another kind of error */
                                ffpmsg_str("column or keyword to be deleted does not exist:");
                                ffpmsg_slice(&clause[clause1..]);
                                return *status;
                            }
                        }
                        if !delall {
                            break;
                        }
                    } /* end do{} */
                }
            } else {
                /* ===================================================== */
                /* Case II:
                this is either a column name, (case 1)

                    or a new column name followed by double = ("==") followed
                    by the old name which is to be renamed. (case 2A)

                    or a column or keyword name followed by a single "=" and a
                calculation expression (case 2B) */
                /* ===================================================== */
                let mut cptr2 = 0; //clause
                let mut tstbuff = None;
                let mut ptr_index = 0;
                slen = fits_get_token2_safe(
                    &clause[cptr2..],
                    &mut ptr_index,
                    cs!("( ="),
                    &mut tstbuff,
                    None,
                    status,
                );

                cptr2 += ptr_index;

                if slen == 0 || *status != 0 {
                    ffpmsg_str("error: column or keyword name is blank (ffedit_columns):");
                    ffpmsg_slice(clause);
                    if *status == 0 {
                        *status = URL_PARSE_ERROR;
                    }
                    return *status;
                }

                let tstbuff = tstbuff.unwrap(); // Since slen != 0, can assume this is valid unwrap

                if strlen_safe(&tstbuff) > FLEN_VALUE - 1 {
                    ffpmsg_str("error: column or keyword name is too long (ffedit_columns):");
                    ffpmsg_slice(clause);
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                strcpy_safe(&mut colname, &tstbuff);
                drop(tstbuff);

                /* If this is a keyword of the form
                     #KEYWORD#
                   then transform to the form
                     #KEYWORDn
                   where n is the previously used column number
                */
                if colname[0] == bb(b'#')
                    && strstr_safe(&colname[1..], cs!("#")) == Some(strlen_safe(&colname) - 1)
                {
                    if colnum <= 0 {
                        ffpmsg_str("The keyword name:");
                        ffpmsg_slice(&colname);
                        ffpmsg_str("is invalid unless a column has been previously");
                        ffpmsg_str("created or editted by a calculator command");
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }
                    colname[strlen_safe(&colname) - 1] = 0;
                    /* Make keyword name and put it in oldname */
                    ffkeyn_safe(&colname[1..], colnum, &mut oldname, status);
                    if *status != 0 {
                        return *status;
                    }
                    /* Re-copy back into colname */
                    strcpy_safe(&mut colname[1..], &oldname);
                } else if strstr_safe(&colname, cs!("#")) == Some(strlen_safe(&colname) - 1) {
                    /*  colname is of the form "NAME#";  if
                          a) colnum is defined, and
                          b) a column with literal name "NAME#" does not exist, and
                          c) a keyword with name "NAMEn" (where n=colnum) exists, then
                        transfrom the colname string to "NAMEn", otherwise
                        do nothing.
                    */
                    if colnum > 0 {
                        /* colnum must be defined */
                        tstatus = 0;
                        ffgcno_safe(
                            fptr,
                            CASEINSEN as c_int,
                            &colname,
                            &mut testnum,
                            &mut tstatus,
                        );
                        if tstatus != 0 && tstatus != COL_NOT_UNIQUE {
                            /* OK, column doesn't exist, now see if keyword exists */
                            ffcmsg_safe(); /* clear previous error message from ffgcno */
                            strcpy_safe(&mut testname, &colname);
                            testname[strlen_safe(&testname) - 1] = 0;
                            /* Make keyword name and put it in oldname */
                            ffkeyn_safe(&testname, colnum, &mut oldname, status);
                            if *status != 0 {
                                return *status;
                            }

                            tstatus = 0;
                            if ffgcrd_safe(fptr, &oldname, &mut card, &mut tstatus) == 0 {
                                /* Keyword does exist; copy real name back into colname */
                                strcpy_safe(&mut colname, &oldname);
                            }
                        }
                    }
                }

                /* if we encountered an opening parenthesis, then we need to */
                /* find the closing parenthesis, and concatinate the 2 strings */
                /* This supports expressions like:
                    [col #EXTNAME(Extension name)="GTI"]
                */
                if clause[cptr2] == bb(b'(') {
                    let mut tstbuff = None;
                    let mut ptr_index = 0;
                    if fits_get_token2_safe(
                        &clause[cptr2..],
                        &mut ptr_index,
                        cs!(")"),
                        &mut tstbuff,
                        None,
                        status,
                    ) == 0
                    {
                        strcat_safe(&mut colname, cs!(")"));
                    } else {
                        cptr2 += ptr_index;
                        let tstbuff = tstbuff.unwrap();
                        if (strlen_safe(&tstbuff) + strlen_safe(&colname) + 1) > FLEN_VALUE - 1 {
                            ffpmsg_str("error: column name is too long (ffedit_columns):");
                            *status = URL_PARSE_ERROR;
                            return *status;
                        }
                        strcat_safe(&mut colname, &tstbuff);
                        strcat_safe(&mut colname, cs!(")"));
                        drop(tstbuff);
                    }
                    cptr2 += 1;
                }

                while clause[cptr2] == bb(b' ') {
                    cptr2 += 1; /* skip white space */
                }

                if clause[cptr2] != bb(b'=') {
                    /* ------------------------------------ */
                    /* case 1 - simply the name of a column */
                    /* ------------------------------------ */

                    /* look for matching column */
                    ffgcno_safe(fptr, CASEINSEN as c_int, &colname, &mut testnum, status);

                    while *status == COL_NOT_UNIQUE {
                        /* the column name contained wild cards, and it */
                        /* matches more than one column in the table. */

                        colnum = testnum;

                        /* keep this column in the output file */
                        savecol = 1;

                        if colindex.is_empty() {
                            colindex = vec![0; 999];
                        }

                        colindex[(colnum - 1) as usize] = 1; /* flag this column number */

                        /* look for other matching column names */
                        ffgcno_safe(fptr, CASEINSEN as c_int, &colname, &mut testnum, status);

                        if *status == COL_NOT_FOUND {
                            *status = 999; /* temporary status flag value */
                        }
                    }

                    if *status <= 0 {
                        colnum = testnum;

                        /* keep this column in the output file */
                        savecol = 1;

                        if colindex.is_empty() {
                            colindex = vec![0; 999];
                        }

                        colindex[(colnum - 1) as usize] = 1; /* flag this column number */
                    } else if *status == 999 {
                        /* this special flag value does not represent an error */
                        *status = 0;
                    } else {
                        ffpmsg_str("Syntax error in columns specifier in input URL:");
                        ffpmsg_slice(&clause[cptr2..]);
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }
                } else {
                    /* ----------------------------------------------- */
                    /* case 2 where the token ends with an equals sign */
                    /* ----------------------------------------------- */

                    cptr2 += 1; /* skip over the first '=' */

                    if clause[cptr2] == bb(b'=') {
                        /*................................................. */
                        /*  Case A:  rename a column or keyword;  syntax is
                        "new_name == old_name"  */
                        /*................................................. */

                        cptr2 += 1; /* skip the 2nd '=' */
                        while clause[cptr2] == bb(b' ') {
                            cptr2 += 1; /* skip white space */
                        }

                        let mut tstbuff = None;
                        let mut ptr_index = 0;
                        if fits_get_token2_safe(
                            &clause[cptr2..],
                            &mut ptr_index,
                            cs!(" "),
                            &mut tstbuff,
                            None,
                            status,
                        ) == 0
                        {
                            oldname[0] = 0;
                        } else {
                            cptr2 += ptr_index;
                            let tstbuff = tstbuff.unwrap();

                            if strlen_safe(&tstbuff) > FLEN_VALUE - 1 {
                                ffpmsg_str(
                                    "error: column name syntax is too long (ffedit_columns):",
                                );
                                *status = URL_PARSE_ERROR;
                                return *status;
                            }
                            strcpy_safe(&mut oldname, &tstbuff);
                        }

                        /* get column number of the existing column */
                        if ffgcno_safe(fptr, CASEINSEN as c_int, &oldname, &mut colnum, status) <= 0
                        {
                            /* modify the TTYPEn keyword value with the new name */
                            ffkeyn_safe(cs!("TTYPE"), colnum, &mut keyname, status);

                            if ffmkys_safe(fptr, &keyname, &colname, None, status) > 0 {
                                ffpmsg_str("failed to rename column in input file");
                                ffpmsg_str(" oldname =");
                                ffpmsg_slice(&oldname);
                                ffpmsg_str(" newname =");
                                ffpmsg_slice(&colname);
                                return *status;
                            }
                            /* keep this column in the output file */
                            savecol = 1;
                            if colindex.is_empty() {
                                colindex = vec![0; 999];
                            }

                            colindex[(colnum - 1) as usize] = 1; /* flag this column number */
                        } else {
                            /* try renaming a keyword */
                            ffcmsg_safe(); /* clear error message stack */
                            *status = 0;
                            if ffmnam_safe(fptr, &oldname, &colname, status) > 0 {
                                ffpmsg_str("column or keyword to be renamed does not exist:");
                                ffpmsg_slice(clause);
                                return *status;
                            }
                        }
                    } else {
                        /*...................................................... */
                        /* Case B: */
                        /* this must be a general column/keyword calc expression */
                        /* "name = expression" or "colname(TFORM) = expression" */
                        /*...................................................... */

                        /* parse the name and TFORM values, if present */
                        colformat[0] = 0;
                        let mut cptr3 = &colname[..];

                        let mut tstbuff = None;
                        let mut ptr_index = 0;
                        if fits_get_token2_safe(
                            cptr3,
                            &mut ptr_index,
                            cs!("("),
                            &mut tstbuff,
                            None,
                            status,
                        ) == 0
                        {
                            oldname[0] = 0;
                        } else {
                            cptr3 = &cptr3[ptr_index..];
                            let tstbuff = tstbuff.unwrap();

                            if strlen_safe(&tstbuff) > FLEN_VALUE - 1 {
                                ffpmsg_str("column expression is too long (ffedit_columns)");
                                *status = URL_PARSE_ERROR;
                                return *status;
                            }
                            strcpy_safe(&mut oldname, &tstbuff);
                            drop(tstbuff);
                        }

                        if cptr3[0] == bb(b'(') {
                            cptr3 = &cptr3[1..]; /* skip the '(' */

                            let mut tstbuff = None;
                            let mut ptr_index = 0;
                            if fits_get_token2_safe(
                                cptr3,
                                &mut ptr_index,
                                cs!(")"),
                                &mut tstbuff,
                                None,
                                status,
                            ) == 0
                            {
                                colformat[0] = 0;
                            } else {
                                cptr3 = &cptr3[ptr_index..];
                                let tstbuff = tstbuff.unwrap();
                                if strlen_safe(&tstbuff) > FLEN_VALUE - 1 {
                                    ffpmsg_str("column expression is too long (ffedit_columns)");
                                    *status = URL_PARSE_ERROR;
                                    return *status;
                                }
                                strcpy_safe(&mut colformat, &tstbuff);
                                drop(tstbuff);
                            }
                        }

                        /* calculate values for the column or keyword */
                        /*   cptr2 = the expression to be calculated */
                        /*   oldname = name of the column or keyword */
                        /*   colformat = column format, or keyword comment string */

                        // WARNING / SAFETY / TODO
                        let same_ftpr = &mut *(fptr.as_mut() as *mut _);

                        if ffcalc_safe(
                            fptr,
                            &clause[cptr2..],
                            same_ftpr,
                            &oldname,
                            &colformat,
                            status,
                        ) > 0
                        {
                            ffpmsg_str("Unable to calculate expression");
                            return *status;
                        }

                        /* test if this is a column and not a keyword */
                        tstatus = 0;
                        ffgcno_safe(
                            fptr,
                            CASEINSEN as c_int,
                            &oldname,
                            &mut testnum,
                            &mut tstatus,
                        );
                        if tstatus == 0 {
                            /* keep this column in the output file */
                            colnum = testnum;
                            savecol = 1;

                            if colindex.is_empty() {
                                colindex = vec![0; 999];
                            }

                            colindex[(colnum - 1) as usize] = 1;
                            if colnum > numcols {
                                numcols += 1;
                            }
                        } else {
                            ffcmsg_safe(); /* clear the error message stack */
                        }
                    }
                }
            }
            //clause = NULL;
        }

        let clause = tmp_clause.as_deref_mut().unwrap();

        if savecol != 0 && deletecol == 0 {
            /* need to delete all but the specified columns */
            let mut ii = numcols;
            while ii > 0 {
                if colindex[(ii - 1) as usize] == 0 {
                    /* delete this column */

                    if ffdcol_safe(fptr, ii, status) > 0 {
                        ffpmsg_str("failed to delete column in input file:");
                        ffpmsg_slice(clause); // SAFETY: This is a hack and probably isn't safe. TODO
                        return *status;
                    }
                }
                ii -= 1;
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Copy a table cell of a given row and column into an image extension.
/// The output file must already have been created.  A new image
/// extension will be created in that file.
///
/// This routine was written by Craig Markwardt, GSFC
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_copy_cell2image(
    fptr: *mut fitsfile,                  /* I - point to input table */
    newptr: *mut fitsfile, /* O - existing output file; new image HDU will be appended to it */
    colname: *const [c_char; FLEN_VALUE], /* I - column name / number containing the image*/
    rownum: c_long,        /* I - number of the row containing the image */
    status: *mut c_int,    /* IO - error status */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let newptr = newptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let colname = &*colname;

        fits_copy_cell2image_safe(fptr, newptr, colname, rownum, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Copy a table cell of a given row and column into an image extension.
/// The output file must already have been created.  A new image
/// extension will be created in that file.
///
/// This routine was written by Craig Markwardt, GSFC
pub(crate) fn fits_copy_cell2image_safe(
    fptr: &mut fitsfile,   /* I - point to input table */
    newptr: &mut fitsfile, /* O - existing output file; new image HDU will be appended to it */
    colname: &[c_char],    /* I - column name / number containing the image*/
    rownum: c_long,        /* I - number of the row containing the image */
    status: &mut c_int,    /* IO - error status */
) -> c_int {
    todo!()
}

/*--------------------------------------------------------------------------*/
/// Copy an image extension into a table cell at a given row and
/// column.  The table must have already been created.  If the "colname"
/// column exists, it will be used, otherwise a new column will be created
/// in the table.
///
/// The "copykeyflag" parameter controls which keywords to copy from the
/// input image to the output table header (with any appropriate translation).
///
/// copykeyflag = 0  -- no keywords will be copied
/// copykeyflag = 1  -- essentially all keywords will be copied
/// copykeyflag = 2  -- copy only the WCS related keywords
///
/// This routine was written by Craig Markwardt, GSFC
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_copy_image2cell(
    fptr: *mut fitsfile,   /* I - pointer to input image extension */
    newptr: *mut fitsfile, /* I - pointer to output table */
    colname: *mut c_char,  /* I - name of column containing the image    */
    rownum: c_long,        /* I - number of the row containing the image */
    copykeyflag: c_int,    /* I - controls which keywords to copy */
    status: *mut c_int,    /* IO - error status */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// copies an image section from the input file to a new output file.
/// Any HDUs preceding or following the image are also copied to the
/// output file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_select_image_section(
    fptr: *mut Option<Box<fitsfile>>, /* IO - pointer to input image; on output it  */
    /*      points to the new subimage */
    outfile: *const c_char, /* I - name for output file        */
    expr: *const c_char,    /* I - Image section expression    */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(outfile);
        raw_to_slice!(expr);

        fits_select_image_section_safer(fptr, outfile, expr, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copies an image section from the input file to a new output file.
/// Any HDUs preceding or following the image are also copied to the
/// output file.
pub(crate) fn fits_select_image_section_safer(
    fptr: &mut Option<Box<fitsfile>>, /* IO - pointer to input image; on output it  */
    /*      points to the new subimage */
    outfile: &[c_char], /* I - name for output file        */
    expr: &[c_char],    /* I - Image section expression    */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// copies an image section from the input file to a new output HDU
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_copy_image_section(
    fptr: *mut fitsfile,   /* I - pointer to input image */
    newptr: *mut fitsfile, /* I - pointer to output image */
    expr: *const c_char,   /* I - Image section expression    */
    status: *mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Parse the input image section specification string, returning
/// the  min, max and increment values.
/// Typical string =   "1:512:2"  or "1:512"
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_get_section_range(
    ptr: *mut *mut c_char,
    secmin: *mut c_long,
    secmax: *mut c_long,
    incre: *mut c_long,
    status: *mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
pub(crate) fn ffselect_table(
    fptr: &mut Option<Box<fitsfile>>, /* IO - pointer to input table; on output it  */
    /*      points to the new selected rows table */
    outfile: &[c_char; FLEN_FILENAME], /* I - name for output file */
    rowfilter: &[c_char; FLEN_FILENAME], /* I - Boolean expression    */
    status: &mut c_int,
) -> c_int {
    todo!()
}

/*--------------------------------------------------------------------------*/
/// Parse the image compression specification that was give in square brackets
/// following the output FITS file name, as in these examples:
///
///   myfile.fits[compress]  - default Rice compression, row by row
///   myfile.fits[compress TYPE] -  the first letter of TYPE defines the
///                                 compression algorithm:
///                                  R = Rice
///                                  G = GZIP
///                                  H = HCOMPRESS
///                                  HS = HCOMPRESS (with smoothing)
///                  B - BZIP2
///                                  P = PLIO
///
///   myfile.fits[compress TYPE 100,100] - the numbers give the dimensions
///                                        of the compression tiles.  Default
///                                        is NAXIS1, 1, 1, ...
///
///   other optional parameters may be specified following a semi-colon
///
///   myfile.fits[compress; q 8.0]          q specifies the floating point
///   mufile.fits[compress TYPE; q -.0002]        quantization level;
///   myfile.fits[compress TYPE 100,100; q 10, s 25]  s specifies the HCOMPRESS
///                                                    integer scaling parameter
///
/// The compression parameters are saved in the fptr->Fptr structure for use
/// when writing FITS images.
pub(crate) fn ffparsecompspec(
    fptr: &mut fitsfile, /* I - FITS file pointer               */
    compspec: &[c_char], /* I - image compression specification */
    status: *mut c_int,  /* IO - error status                       */
) -> c_int {
    todo!()
}

/*--------------------------------------------------------------------------*/
/// Create and initialize a new FITS file on disk.  This routine differs
/// from ffinit in that the input 'name' is literally taken as the name
/// of the disk file to be created, and it does not support CFITSIO's
/// extended filename syntax.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdkinit(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - name of file to create              */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(name);

        /* initialize null file pointer */
        let f_tmp = fptr.take();
        if let Some(f) = f_tmp {
            // WARNING: The c version doesn't null pointers after a close, so we have a dangling pointer.
            // We need to be careful with this, as it can cause double free errors.
            // Therefore, if this function is called with a Some(), then we will leak the pointer because
            // it's probably invalid.
            let _ = Box::into_raw(f);
        }

        /* regardless of the value of *status */
        if *status > 0 {
            return *status;
        }

        *status = CREATE_DISK_FILE;

        ffinit_safer(fptr, name, status);

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Create and initialize a new FITS file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffinit(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: *const c_char,              /* I - name of file to create              */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(name);

        ffinit_safer(fptr, name, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Create and initialize a new FITS file.
pub unsafe fn ffinit_safer(
    fptr: &mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    name: &[c_char],                  /* I - name of file to create              */
    status: &mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let ii: c_int = 0;
        let mut driver: c_int = 0;
        let slen: c_int = 0;
        let mut clobber: c_int = 0 as c_int;

        let mut urltype: [c_char; MAX_PREFIX_LEN] = [0; MAX_PREFIX_LEN];
        let mut outfile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut tmplfile: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        let mut compspec: [c_char; 80] = [0; 80];
        let mut handle: c_int = 0;
        let mut create_disk_file: c_int = 0 as c_int;

        /* initialize null file pointer */
        let f_tmp = fptr.take();
        if let Some(f) = f_tmp {
            // WARNING: The c version doesn't null pointers after a close, so we have a dangling pointer.
            // We need to be careful with this, as it can cause double free errors.
            // Therefore, if this function is called with a Some(), then we will leak the pointer because
            // it's probably invalid.
            let _ = Box::into_raw(f);
        }

        /* regardless of the value of *status */
        if *status > 0 {
            return *status;
        }

        if *status == CREATE_DISK_FILE {
            create_disk_file = 1;
            *status = 0;
        }

        if *NEED_TO_INITIALIZE.lock().unwrap() {
            /* this is called only once */
            *status = fits_init_cfitsio_safer();
        }

        if *status > 0 {
            return *status;
        }

        let url = name;

        let mut j = 0;
        while url[j] == bb(b' ') {
            /* ignore leading spaces in the filename */
            j += 1;
        }
        if url[j] == 0 {
            ffpmsg_str("Name of file to create is blank. (ffinit)");
            *status = FILE_NOT_CREATED;
            return *status;
        }

        if create_disk_file > 0 {
            if strlen_safe(&url[j..]) > FLEN_FILENAME - 1 {
                ffpmsg_str("Filename is too long. (ffinit)");
                *status = FILE_NOT_CREATED;
                return *status;
            }

            strcpy_safe(&mut outfile, &url[j..]);
            strcpy_safe(&mut urltype, cs!("file://"));
            tmplfile[0] = 0;
            compspec[0] = 0;
        } else {
            /* check for clobber symbol, i.e,  overwrite existing file */
            if url[j] == bb(b'!') {
                clobber = TRUE as c_int;
                j += 1;
            } else {
                clobber = FALSE as c_int;
            }
            /* parse the output file specification */
            /* this routine checks that the strings will not overflow */

            ffourl(
                &url[j..],
                urltype.as_mut_ptr(),
                outfile.as_mut_ptr(),
                tmplfile.as_mut_ptr(),
                compspec.as_mut_ptr(),
                status,
            );

            if *status > 0 {
                ffpmsg_str("could not parse the output filename: (ffinit)");
                ffpmsg_slice(&url[j..]);
                return *status;
            }
        }

        let url = &url[j..];

        /* find which driver corresponds to the urltype */
        *status = urltype2driver(&urltype, &mut driver);

        if *status > 0 {
            ffpmsg_str("could not find driver for this file: (ffinit)");
            ffpmsg_slice(url);
            return *status;
        }

        /* delete pre-existing file, if asked to do so */
        if clobber > 0 {
            {
                //let d = driverTable.lock().unwrap();
                let d = DRIVER_TABLE.get().unwrap();
                if d[driver as usize].remove.is_some() {
                    (d[driver as usize].remove.unwrap())(&outfile);
                }
            }
        }

        {
            //let d = driverTable.lock().unwrap();
            let d = DRIVER_TABLE.get().unwrap();

            /* call appropriate driver to create the file */
            if d[driver as usize].create.is_some() {
                let lock = FFLOCK(); /* lock this while searching for vacant handle */
                *status = (d[driver as usize].create.unwrap())(&mut outfile, &mut handle);
                FFUNLOCK(lock);

                if *status > 0 {
                    ffpmsg_str("failed to create new file (already exists?):");
                    ffpmsg_slice(url);
                    return *status;
                }
            } else {
                ffpmsg_str("cannot create a new file of this type: (ffinit)");
                ffpmsg_slice(url);
                *status = FILE_NOT_CREATED;
                return *status;
            }
        }

        let d = DRIVER_TABLE.get().unwrap();

        /* allocate fitsfile structure and initialize = 0 */
        let Fptr = FITSfile::new(&d[driver as usize], handle, url, cs!("ffinit"), status);
        if Fptr.is_err() {
            return *status;
        }
        let mut Fptr = Fptr.unwrap();

        /* initialize the ageindex array (relative age of the I/O buffers) */
        /* and initialize the bufrecnum array as being empty */
        let mut ii = 0;
        while ii < NIOBUF as usize {
            Fptr.ageindex[ii] = ii as c_int;
            Fptr.bufrecnum[ii] = -1;
            ii += 1;
        }

        /* store the parameters describing the file */
        Fptr.MAXHDU = 1000; /* initial size of headstart */
        Fptr.filehandle = handle; /* store the file pointer */
        Fptr.driver = driver; /*  driver number         */
        strcpy(Fptr.filename, url.as_ptr()); /* full input filename    */
        Fptr.filesize = 0; /* physical file size     */
        Fptr.logfilesize = 0; /* logical file size      */
        Fptr.writemode = 1; /* read-write mode        */
        Fptr.datastart = DATA_UNDEFINED as LONGLONG; /* unknown start of data  */
        Fptr.curbuf = -1; /* undefined current IO buffer   */
        Fptr.open_count = 1; /* structure is currently used once */
        Fptr.validcode = VALIDSTRUC; /* flag denoting valid structure */
        Fptr.noextsyntax = create_disk_file; /* true if extended syntax is disabled */

        // HEAP ALLOCATION
        /* allocate fitsfile structure and initialize = 0 */
        let f_fitsfile = Box::try_new(fitsfile {
            HDUposition: 0,
            Fptr,
        });

        if f_fitsfile.is_err() {
            let d = DRIVER_TABLE.get().unwrap();
            ((d[driver as usize]).close)(handle); /* close the file */
            ffpmsg_str("failed to allocate structure for following file: (ffopen)");
            ffpmsg_slice(url);
            *status = MEMORY_ALLOCATION;
            return *status;
        }

        let mut f_fitsfile = f_fitsfile.unwrap();

        ffldrc(&mut f_fitsfile, 0, IGNORE_EOF, status); /* initialize first record */

        fits_store_Fptr(&mut f_fitsfile.Fptr, status); /* store Fptr address */

        /* if template file was given, use it to define structure of new file */

        if tmplfile[0] > 0 {
            ffoptplt(&mut f_fitsfile, &tmplfile, status);
        }

        /* parse and save image compression specification, if given */
        if compspec[0] > 0 {
            ffparsecompspec(&mut f_fitsfile, &compspec, status);
        }

        *fptr = Some(f_fitsfile);

        *status /* successful return */
    }
}

/*--------------------------------------------------------------------------*/
/// ffimem == fits_create_memfile
/// Create and initialize a new FITS file in memory
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffimem(
    fptr: *mut *mut fitsfile,  /* O - FITS file pointer                   */
    buffptr: *mut *mut c_void, /* I - address of memory pointer           */
    buffsize: *mut usize,      /* I - size of buffer, in bytes            */
    deltasize: usize,          /* I - increment for future realloc's      */
    mem_realloc: unsafe extern "C" fn(p: *mut c_void, newsize: usize) -> *mut c_void, /* function       */
    status: *mut c_int, /* IO - error status                       */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Initialize anything that is required before using the CFITSIO routines
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_init_cfitsio() -> c_int {
    fits_init_cfitsio_safer()
}

/*--------------------------------------------------------------------------*/
/// Initialize anything that is required before using the CFITSIO routines
pub(crate) fn fits_init_cfitsio_safer() -> c_int {
    pub union u_tag {
        pub ival: i16,
        pub cval: [c_char; 2],
    }

    fitsio_init_lock();

    let lock = FFLOCK(); /* lockout other threads while executing this critical */
    /* section of code  */

    if !*NEED_TO_INITIALIZE.lock().unwrap() {
        /* already initialized? */
        FFUNLOCK(lock);
        return 0;
    }

    let mut u_ival: i16 = 1;
    let u_cval: &mut [c_char; 2] = cast_mut(&mut u_ival);

    // Because of the FFLOCK mutex, we can do this in two stages.
    let mut drivers = Vec::with_capacity(MAX_DRIVERS);
    let d = &mut drivers;

    /*   test for correct byteswapping.   */
    if (BYTESWAPPED && u_cval[0] != 1) || (!BYTESWAPPED && u_cval[1] != 1) {
        println!("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        println!(" Byteswapping is not being done correctly on this system.");
        println!(" Check the CFITSIO_MACHINE and BYTESWAPPED definitions in fitsio2.h");
        println!(" Please report this problem to the CFITSIO developers.");
        println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        FFUNLOCK(lock);
        return 1;
    }

    /*  test that LONGLONG is an 8 byte integer */
    if std::mem::size_of::<LONGLONG>() != 8 {
        println!("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        println!(" CFITSIO did not find an 8-byte long integer data type.");
        println!("   sizeof(LONGLONG) = {}", std::mem::size_of::<LONGLONG>());
        println!(" Please report this problem to the CFITSIO developers.");
        println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        FFUNLOCK(lock);
        return 1;
    }

    /* register the standard I/O drivers that are always available */

    /* 1--------------------disk file driver-----------------------*/
    let status = fits_register_driver(
        d,
        c"file://",
        Some(file_init),
        Some(file_shutdown),
        Some(file_setoptions),
        Some(file_getoptions),
        Some(file_getversion),
        Some(file_checkfile),
        Some(file_open),
        Some(file_create),
        Some(file_truncate),
        file_close,
        Some(file_remove),
        file_size,
        Some(file_flush),
        file_seek,
        file_read,
        file_write,
    );

    if status != 0 {
        ffpmsg_str("failed to register the file:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 2------------ output temporary memory file driver ----------------*/
    let status = fits_register_driver(
        d,
        c"mem://",
        Some(mem_init),
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None, /* checkfile not needed */
        None, /* open function not allowed */
        Some(mem_create),
        Some(mem_truncate_unsafe),
        mem_close_free_unsafe,
        None, /* remove function not required */
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the mem:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 3--------------input pre-existing memory file driver----------------*/
    let status = fits_register_driver(
        d,
        c"memkeep://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None, /* checkfile not needed */
        None, /* file open driver function is not used */
        None, /* create function not allowed */
        Some(mem_truncate_unsafe),
        mem_close_keep,
        None, /* remove function not required */
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the memkeep:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* ==================== SHARED MEMORY DRIVER SECTION ======================= */

    #[cfg(feature = "shared_mem")]
    {
        /* 22--------------------shared memory driver-----------------------*/
        let status = fits_register_driver(
            d,
            c"shmem://",
            Some(smem_init),
            Some(smem_shutdown),
            Some(smem_setoptions),
            Some(smem_getoptions),
            Some(smem_getversion),
            None, /* checkfile not needed */
            Some(smem_open),
            Some(smem_create),
            None, /* truncate file not supported yet */
            smem_close,
            Some(smem_remove),
            smem_size,
            Some(smem_flush),
            smem_seek,
            Some(smem_read),
            Some(smem_write),
        );

        if status != 0 {
            ffpmsg_str("failed to register the shmem:// driver (init_cfitsio)");
            FFUNLOCK(lock);
            return status;
        }
    }

    /* 4-------------------stdin stream driver----------------------*/
    /*  the stdin stream is copied to memory then opened in memory */
    let status = fits_register_driver(
        d,
        c"stdin://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        Some(stdin_checkfile),
        Some(stdin_open),
        None, /* create function not allowed */
        Some(mem_truncate_unsafe),
        mem_close_free_unsafe,
        None, /* remove function not required */
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the stdin:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 5-------------------stdin file stream driver----------------------*/
    /*  the stdin stream is copied to a disk file then the disk file is opened */
    let status = fits_register_driver(
        d,
        c"stdinfile://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None, /* checkfile not needed */
        Some(stdin_open),
        None, /* create function not allowed */
        Some(file_truncate),
        file_close,
        Some(file_remove),
        file_size,
        Some(file_flush),
        file_seek,
        file_read,
        file_write,
    );

    if status != 0 {
        ffpmsg_str("failed to register the stdinfile:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 6-----------------------stdout stream driver------------------*/
    let status = fits_register_driver(
        d,
        c"stdout://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None,
        None,
        Some(mem_create),
        Some(mem_truncate_unsafe),
        stdout_close_unsafe,
        None,
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the stdout:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 7------------------iraf disk file to memory driver -----------*/
    let status = fits_register_driver(
        d,
        c"irafmem://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None,
        Some(mem_iraf_open),
        None,
        Some(mem_truncate_unsafe),
        mem_close_free_unsafe,
        None,
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the irafmem:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 8------------------raw binary file to memory driver -----------*/
    let status = fits_register_driver(
        d,
        c"rawfile://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None,
        Some(mem_rawfile_open),
        None,
        Some(mem_truncate_unsafe),
        mem_close_free_unsafe,
        None,
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the rawfile:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 9------------------compressed disk file to memory driver -----------*/
    let status = fits_register_driver(
        d,
        c"compress://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None,
        Some(mem_compress_open),
        None,
        Some(mem_truncate_unsafe),
        mem_close_free_unsafe,
        None,
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the compress:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 10------------------compressed disk file to memory driver -----------*/
    let status = fits_register_driver(
        d,
        c"compressmem://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None,
        Some(mem_compress_openrw),
        None,
        Some(mem_truncate_unsafe),
        mem_close_free_unsafe,
        None, /* remove function not required */
        mem_size,
        None, /* flush function not required */
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the compressmem:// driver (init_cfitsio)"); /* checkfile not needed */
        FFUNLOCK(lock);
        return status; /* no file truncate function */
    }

    /* 11------------------compressed disk file to disk file driver -------*/
    let status = fits_register_driver(
        d,
        c"compressfile://",
        None,
        Some(file_shutdown),
        Some(file_setoptions),
        Some(file_getoptions),
        Some(file_getversion),
        None,
        Some(file_compress_open),
        Some(file_create),
        None,
        file_close,
        Some(file_remove),
        file_size,
        Some(file_flush),
        file_seek,
        file_read,
        file_write,
    );

    if status != 0 {
        ffpmsg_str("failed to register the compressfile:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    /* 12---create file in memory, then compress it to disk file on close--*/
    let status = fits_register_driver(
        d,
        c"compressoutfile://",
        None,
        Some(mem_shutdown),
        Some(mem_setoptions),
        Some(mem_getoptions),
        Some(mem_getversion),
        None,
        None,
        Some(mem_create_comp_unsafe),
        Some(mem_truncate_unsafe),
        mem_close_comp_unsafe,
        Some(file_remove),
        mem_size,
        None,
        mem_seek,
        mem_read_unsafe,
        mem_write_unsafe,
    );

    if status != 0 {
        ffpmsg_str("failed to register the compressoutfile:// driver (init_cfitsio)"); /* 13--------------------root driver-----------------------*/
        FFUNLOCK(lock);
        return status; /* checkfile not needed */
    }

    /* Register Optional drivers */

    /* 24---------------stdin and stdout stream driver-------------------*/
    let status = fits_register_driver(
        d,
        c"stream://",
        None,
        None,
        None,
        None,
        None,
        None,
        Some(stream_open),
        Some(stream_create),
        None, /* no stream truncate function */
        stream_close,
        None, /* no stream remove */
        stream_size,
        Some(stream_flush),
        stream_seek,
        stream_read,
        stream_write,
    );

    if status != 0 {
        ffpmsg_str("failed to register the stream:// driver (init_cfitsio)");
        FFUNLOCK(lock);
        return status;
    }

    let dd = DRIVER_TABLE.set(drivers);

    let l = NEED_TO_INITIALIZE.lock();
    *l.unwrap() = false;

    FFUNLOCK(lock);
    status
}

/*--------------------------------------------------------------------------*/
/// register all the functions needed to support an I/O driver
pub(crate) fn fits_register_driver(
    drivers: &mut Vec<fitsdriver>,
    prefix: &CStr,
    init: Option<fn() -> c_int>,
    shutdown: Option<fn() -> c_int>,
    setoptions: Option<fn(option: c_int) -> c_int>,
    getoptions: Option<fn(options: &mut c_int) -> c_int>,
    getversion: Option<fn(version: &mut c_int) -> c_int>,
    checkfile: Option<
        fn(
            urltype: &mut [c_char; MAX_PREFIX_LEN],
            infile: &mut [c_char; FLEN_FILENAME],
            outfile: &mut [c_char; FLEN_FILENAME],
        ) -> c_int,
    >,
    open: Option<fn(filename: &mut [c_char], rwmode: c_int, handle: &mut c_int) -> c_int>,
    create: Option<fn(filename: &mut [c_char; FLEN_FILENAME], drivehandle: &mut c_int) -> c_int>,
    truncate: Option<fn(drivehandle: c_int, size: usize) -> c_int>,
    close: fn(drivehandle: c_int) -> c_int,
    remove: Option<fn(filename: &[c_char]) -> c_int>,
    size: fn(drivehandle: c_int, sizex: &mut usize) -> c_int,
    flush: Option<fn(drivehandle: c_int) -> c_int>,
    seek: fn(drivehandle: c_int, offset: LONGLONG) -> c_int,
    read: fn(drivehandle: c_int, buffer: &mut [u8], nbytes: usize) -> c_int,
    write: fn(drivehandle: c_int, buffer: &[u8], nbyte: usize) -> c_int,
) -> c_int {
    let status; /* increment the number of drivers */

    if drivers.len() + 1 > MAX_DRIVERS {
        return TOO_MANY_DRIVERS;
    }

    if prefix.is_empty() {
        return BAD_URL_PREFIX;
    }

    if init.is_some() {
        status = (init.unwrap())(); /* initialize the driver */
        if status != 0 {
            return status;
        };
    }

    /*  fill in data in table */
    let mut pf = [0; MAX_PREFIX_LEN];
    strncpy_safe(
        &mut pf,
        cast_slice(prefix.to_bytes_with_nul()),
        MAX_PREFIX_LEN,
    );
    pf[MAX_PREFIX_LEN - 1] = 0;

    let dd = fitsdriver {
        prefix: pf,
        init,
        shutdown,
        setoptions,
        getoptions,
        getversion,
        checkfile,
        open,
        create,
        truncate,
        close,
        remove,
        size,
        flush,
        seek,
        read,
        write,
    };

    drivers.push(dd); /* increment the number of drivers */

    0
}

/*--------------------------------------------------------------------------*/
/// fits_parse_input_url
/// parse the input URL into its basic components.
/// This routine does not support the pixfilter or compspec components.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffiurl(
    url: *const c_char,      /* input filename */
    urltype: *mut c_char,    /* e.g., 'file://', 'http://', 'mem://' */
    infilex: *mut c_char,    /* root filename (may be complete path) */
    outfile: *mut c_char,    /* optional output file name            */
    extspec: *mut c_char,    /* extension spec: +n or [extname, extver]  */
    rowfilterx: *mut c_char, /* boolean row filter expression */
    binspec: *mut c_char,    /* histogram binning specifier   */
    colspec: *mut c_char,    /* column or keyword modifier expression */
    status: *mut c_int,
) -> c_int {
    unsafe {
        raw_to_slice!(url);

        let status = status.as_mut().expect(NULL_MSG);

        ffiurl_safer(
            url, urltype, infilex, outfile, extspec, rowfilterx, binspec, colspec, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// fits_parse_input_url
/// parse the input URL into its basic components.
/// This routine does not support the pixfilter or compspec components.
pub(crate) fn ffiurl_safer(
    url: &[c_char],          /* input filename */
    urltype: *mut c_char,    /* e.g., 'file://', 'http://', 'mem://' */
    infilex: *mut c_char,    /* root filename (may be complete path) */
    outfile: *mut c_char,    /* optional output file name            */
    extspec: *mut c_char,    /* extension spec: +n or [extname, extver]  */
    rowfilterx: *mut c_char, /* boolean row filter expression */
    binspec: *mut c_char,    /* histogram binning specifier   */
    colspec: *mut c_char,    /* column or keyword modifier expression */
    status: &mut c_int,
) -> c_int {
    ffifile2_safer(
        url,
        urltype,
        infilex,
        outfile,
        extspec,
        rowfilterx,
        binspec,
        colspec,
        ptr::null_mut(),
        ptr::null_mut(),
        status,
    )
}

/*--------------------------------------------------------------------------*/
/// fits_parse_input_filename
/// parse the input URL into its basic components.
/// This routine does not support the compspec component.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffifile(
    url: *const c_char,      /* input filename */
    urltype: *mut c_char,    /* e.g., 'file://', 'http://', 'mem://' */
    infilex: *mut c_char,    /* root filename (may be complete path) */
    outfile: *mut c_char,    /* optional output file name            */
    extspec: *mut c_char,    /* extension spec: +n or [extname, extver]  */
    rowfilterx: *mut c_char, /* boolean row filter expression */
    binspec: *mut c_char,    /* histogram binning specifier   */
    colspec: *mut c_char,    /* column or keyword modifier expression */
    pixfilter: *mut c_char,  /* pixel filter expression */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        raw_to_slice!(url);

        ffifile2_safer(
            url,
            urltype,
            infilex,
            outfile,
            extspec,
            rowfilterx,
            binspec,
            colspec,
            pixfilter,
            ptr::null_mut(),
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// fits_parse_input_filename
/// parse the input URL into its basic components.
/// This routine is big and ugly and should be redesigned someday!
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffifile2(
    url: *const c_char,      /* input filename */
    urltype: *mut c_char,    /* e.g., 'file://', 'http://', 'mem://' */
    infilex: *mut c_char,    /* root filename (may be complete path) */
    outfile: *mut c_char,    /* optional output file name            */
    extspec: *mut c_char,    /* extension spec: +n or [extname, extver]  */
    rowfilterx: *mut c_char, /* boolean row filter expression */
    binspec: *mut c_char,    /* histogram binning specifier   */
    colspec: *mut c_char,    /* column or keyword modifier expression */
    pixfilter: *mut c_char,  /* pixel filter expression */
    compspec: *mut c_char,   /* image compression specification */
    status: *mut c_int,
) -> c_int {
    unsafe {
        raw_to_slice!(url);

        let status = status.as_mut().expect(NULL_MSG);

        ffifile2_safer(
            url, urltype, infilex, outfile, extspec, rowfilterx, binspec, colspec, pixfilter,
            compspec, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// fits_parse_input_filename
/// parse the input URL into its basic components.
/// This routine is big and ugly and should be redesigned someday!
pub(crate) fn ffifile2_safer(
    url: &[c_char],          /* input filename */
    urltype: *mut c_char,    /* e.g., 'file://', 'http://', 'mem://' */
    infilex: *mut c_char,    /* root filename (may be complete path) */
    outfile: *mut c_char,    /* optional output file name            */
    extspec: *mut c_char,    /* extension spec: +n or [extname, extver]  */
    rowfilterx: *mut c_char, /* boolean row filter expression */
    binspec: *mut c_char,    /* histogram binning specifier   */
    colspec: *mut c_char,    /* column or keyword modifier expression */
    pixfilter: *mut c_char,  /* pixel filter expression */
    compspec: *mut c_char,   /* image compression specification */
    status: &mut c_int,
) -> c_int {
    unsafe {
        let mut infilelen = 0;
        let mut plus_ext = 0;
        let collen = 0;
        let mut ptr2: *mut c_char = ptr::null_mut();
        let mut ptr3: *mut c_char = ptr::null_mut();
        let mut ptr4: *mut c_char = ptr::null_mut();
        let mut hasAt = 0;
        let mut hasDot = 0;
        let mut hasOper = 0;
        let mut followingOper = 0;
        let mut spaceTerm = 0;
        let mut rowFilter = 0;
        let mut colStart = 0;
        let mut binStart = 0;
        let mut pixStart = 0;
        let mut compStart = 0;

        if *status > 0 {
            return *status;
        }

        /* Initialize null strings */
        if !infilex.is_null() {
            *infilex = 0;
        }
        if !urltype.is_null() {
            *urltype = 0;
        }
        if !outfile.is_null() {
            *outfile = 0;
        }
        if !extspec.is_null() {
            *extspec = 0;
        }
        if !binspec.is_null() {
            *binspec = 0;
        }
        if !colspec.is_null() {
            *colspec = 0;
        }
        if !rowfilterx.is_null() {
            *rowfilterx = 0;
        }
        if !pixfilter.is_null() {
            *pixfilter = 0;
        }
        if !compspec.is_null() {
            *compspec = 0;
        }

        let mut slen = strlen_safe(url);

        if slen == 0 {
            /* blank filename ?? */
            return *status;
        }

        // TEMP HEAP ALLOCATION
        /* allocate memory for 3 strings, each as long as the input url */
        let mut holding_tank = Vec::new();
        if holding_tank.try_reserve_exact(3 * (slen + 1)).is_err() {
            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            holding_tank.resize(3 * (slen + 1), 0);
        }

        // Split the 1 vec into 3 mutables slices
        let (infile, rowfilter) = holding_tank.split_at_mut(slen + 1);
        let (rowfilter, tmpstr) = rowfilter.split_at_mut(slen + 1);

        let mut ptr1 = url;
        let mut ptr1_index = 0;

        /* -------------------------------------------------------- */
        /*  get urltype (e.g., file://, ftp://, http://, etc.)  */
        /* --------------------------------------------------------- */

        if ptr1[0] == bb(b'-')
            && (ptr1[1] == 0 || ptr1[1] == bb(b' ') || ptr1[1] == bb(b'[') || ptr1[1] == bb(b'('))
        {
            /* "-" means read file from stdin. Also support "- ",        */
            /* "-[extname]" and '-(outfile.fits)" but exclude disk file  */
            /* names that begin with a minus sign, e.g., "-55d33m.fits"  */

            if !urltype.is_null() {
                strcat(urltype, cstr!(b"stdin://").as_ptr());
            }

            ptr1_index += 1;
        } else if fits_strncasecmp(ptr1, cs!(b"stdin"), 5) == 0 {
            if !urltype.is_null() {
                strcat(urltype, cstr!(b"stdin://").as_ptr());
            }
            ptr1_index += 5;
        } else {
            let mut ptr2 = strstr_safe(ptr1, cs!(b"://"));
            let ptr3 = strstr_safe(ptr1, cs!(b"("));
            if ptr3.is_some() && (ptr3.unwrap() < ptr2.unwrap()) {
                /* the urltype follows a '(' character, so it must apply */
                /* to the output file, and is not the urltype of the input file */
                ptr2 = None; /* so reset pointer to zero */
            }

            if let Some(ptr2) = ptr2 {
                /* copy the explicit urltype string */
                if (ptr2 + 3) >= MAX_PREFIX_LEN {
                    ffpmsg_str("Name of urltype is too long.");
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                if !urltype.is_null() {
                    strncat(urltype, ptr1.as_ptr(), ptr2 + 3);
                }
                ptr1_index += 3;
            } else if strncmp_safe(ptr1, cs!(b"ftp:"), 4) == 0 {
                /* the 2 //'s are optional */
                if !urltype.is_null() {
                    strcat(urltype, cstr!(b"ftp://").as_ptr());
                }
                ptr1_index += 4;
            } else if strncmp_safe(ptr1, cs!(b"gsiftp:"), 7) == 0 {
                /* the 2 //'s are optional */
                if !urltype.is_null() {
                    strcat(urltype, cstr!(b"gsiftp://").as_ptr());
                }
                ptr1_index += 7;
            } else if strncmp_safe(ptr1, cs!(b"http:"), 5) == 0 {
                /* the 2 //'s are optional */
                if !urltype.is_null() {
                    strcat(urltype, cstr!(b"http://").as_ptr());
                }
                ptr1_index += 5;
            } else if strncmp_safe(ptr1, cs!(b"mem:"), 4) == 0 {
                /* the 2 //'s are optional */
                if !urltype.is_null() {
                    strcat(urltype, cstr!(b"mem://").as_ptr());
                }
                ptr1_index += 4;
            } else if strncmp_safe(ptr1, cs!(b"shmem:"), 6) == 0 {
                /* the 2 //'s are optional */
                if !urltype.is_null() {
                    strcat(urltype, cstr!(b"shmem://").as_ptr());
                }
                ptr1_index += 6;
            } else if strncmp_safe(ptr1, cs!(b"file:"), 5) == 0 {
                /* the 2 //'s are optional */
                if !urltype.is_null() {
                    strcat(urltype, cstr!(b"file://").as_ptr());
                }
                ptr1_index += 5;
            } else {
                /* assume file driver    */
                if !urltype.is_null() {
                    strcat(urltype, cstr!(b"file://").as_ptr());
                };
            };
        }

        /* ----------------------------------------------------------
        If this is a http:// type file, then the cgi file name could
        include the '[' character, which should not be interpreted
        as part of CFITSIO's Extended File Name Syntax.  Test for this
        case by seeing if the last character is a ']' or ')'.  If it
        is not, then just treat the whole input string as the file name
        and do not attempt to interprete the name using the extended
        filename syntax.
        ----------------------------------------------------------- */

        // Advance slice
        ptr1 = &ptr1[ptr1_index..];
        ptr1_index = 0;

        if !urltype.is_null() && strncmp(urltype, cstr!(b"http://").as_ptr(), 7) == 0 {
            /* test for opening parenthesis or bracket in the file name */
            if strchr_safe(ptr1, bb(b'(')).is_some() || strchr_safe(ptr1, bb(b'[')).is_some() {
                slen = strlen_safe(ptr1);

                let mut ptr3_index = slen - 1;

                while ptr1[ptr3_index] == bb(b' ') {
                    /* ignore trailing blanks */
                    ptr3_index -= 1;
                }

                if ptr1[ptr3_index] != bb(b']') && ptr1[ptr3_index] != bb(b')') {
                    /* name doesn't end with a ']' or ')' so don't try */
                    /* to parse this unusual string (may be cgi string)  */
                    if !infilex.is_null() {
                        if strlen_safe(ptr1) > FLEN_FILENAME - 1 {
                            ffpmsg_str("Name of file is too long.");
                            *status = URL_PARSE_ERROR;
                            return *status;
                        }
                        strcpy(infilex, ptr1.as_ptr());
                    }

                    return *status;
                };
            };
        }
        /* ----------------------------------------------------------
        Look for VMS style filenames like:
             disk:[directory.subdirectory]filename.ext, or
                  [directory.subdirectory]filename.ext

        Check if the first character is a '[' and urltype != stdin
        or if there is a ':[' string in the remaining url string. If
        so, then need to move past this bracket character before
        search for the opening bracket of a filter specification.
        ----------------------------------------------------------- */
        let mut tmptr = ptr1; /* this bracket encloses a VMS directory name */
        if ptr1[0] == bb(b'[') {
            if url[0] != bb(b'-') {
                tmptr = &ptr1[1..]; /* these 2 chars are part of the VMS disk and directory */
            };
        } else {
            let tmp_index = strstr_safe(ptr1, cs!(b":["));
            if let Some(tmp_index) = tmp_index {
                tmptr = &ptr1[(tmp_index + 2)..];
            } else {
                tmptr = ptr1;
            };
        }

        /* ------------------------ */
        /*  get the input file name */
        /* ------------------------ */
        ptr2 = strchr(tmptr.as_ptr(), bb(b'(') as c_int); /* search for opening parenthesis ( */
        ptr3 = strchr(tmptr.as_ptr(), bb(b'[') as c_int); /* search for opening bracket [ */
        if !ptr2.is_null() {
            ptr4 = strchr(ptr2, bb(b')') as c_int); /* search for closing parenthesis ) */
            while !ptr4.is_null() && !ptr2.is_null() {
                loop {
                    ptr4 = ptr4.offset(1); /* find next non-blank char after ')' */
                    if *ptr4 != bb(b' ') {
                        break;
                    }; /* simple case: no [ or ( in the file name */
                } /* no bracket, so () enclose output file name */
                if *ptr4 == 0 || *ptr4 == bb(b'[') {
                    break; /* () enclose output name before bracket */
                } /* search for closing ) */
                ptr2 = strchr(ptr2.offset(1), bb(b'(') as c_int); /* error, no closing ) */
                ptr4 = strchr(ptr4, bb(b')') as c_int);
            }
        }

        if std::ptr::eq(ptr2, ptr3) {
            /* simple case: no [ or ( in the file name */
            strcat_safe(infile, ptr1);
        } else if ptr3.is_null() || (!ptr2.is_null() && (ptr2 < ptr3)) {
            let t = ptr2.offset_from(ptr1.as_ptr());
            assert!(t >= 0, "Pointer math is bad mkay");

            strncat_safe(infile, ptr1, t as usize);

            ptr2 = ptr2.add(1);

            let p1 = strchr(ptr2, bb(b')') as c_int);
            if p1.is_null() {
                *status = URL_PARSE_ERROR;
                return *status;
            }

            if !outfile.is_null() {
                let t = p1.offset_from(ptr2);
                assert!(t >= 0, "Pointer math is bad mkay");
                let t = t as usize;

                if t > FLEN_FILENAME - 1 {
                    *status = URL_PARSE_ERROR;
                    return *status;
                }
                strncat(outfile, ptr2, t);
            }
            /* the opening [ could have been part of output name,    */
            /*      e.g., file(out[compress])[3][#row > 5]           */
            /* so search again for opening bracket following the closing ) */
            ptr3 = strchr(p1, bb(b'[') as c_int); /*   bracket comes first, so there is no output name */
        } else {
            let t = ptr3.offset_from(ptr1.as_ptr());
            assert!(t >= 0, "Pointer math is bad mkay");
            let t = t as usize;
            strncat_safe(infile, ptr1, t);
        }

        /* strip off any trailing blanks in the names */
        let mut slen = strlen_safe(infile) as isize;

        slen -= 1;
        while slen > 0 && infile[slen as usize] == bb(b' ') {
            infile[slen as usize] = 0;
            slen -= 1;
        }

        if !outfile.is_null() {
            slen = strlen(outfile) as isize;
            slen -= 1;
            while slen > 0 && *outfile.offset(slen) == bb(b' ') {
                *outfile.offset(slen) = 0;
                slen -= 1;
            }
        }
        /* --------------------------------------------- */
        /* check if this is an IRAF file (.imh extension */
        /* --------------------------------------------- */
        let ptr4 = strstr_safe(infile, cs!(b".imh"));
        /* did the infile name end with cstr!(b".imh").as_ptr() ? */
        if let Some(p4) = ptr4 {
            if infile[p4] == 0 && !urltype.is_null() {
                strcpy(urltype, cstr!(b"irafmem://").as_ptr());
            };
        }

        /* --------------------------------------------- */
        /* check if the 'filename+n' convention has been */
        /* used to specifiy which HDU number to open     */
        /* --------------------------------------------- */
        let jj = strlen_safe(infile) as isize; /* search backwards for '+' sign */
        let mut ii: isize = jj as isize - 1;
        while ii >= 0 {
            if infile[ii as usize] == bb(b'+') {
                break; /* limit extension numbers to 5 digits */
            }; /* pointer to start of sequence */
            ii -= 1;
        }
        /* are all the chars digits? */

        if ii > 0 && (jj - ii) < 7 {
            infilelen = ii; /* yes, the '+n' convention was used.  Copy */
            ii += 1; /* the digits to the output extspec string. */
            let ptr1 = &infile[(ii as usize)..]; /* delete the extension number */
            for ii in ii..jj {
                if !isdigit_safe(infile[ii as usize]) {
                    break;
                };
            }
            if ii == jj {
                plus_ext = 1;
                if !extspec.is_null() {
                    if jj - infilelen > FLEN_FILENAME as isize - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }
                    strncpy(extspec, ptr1.as_ptr(), (jj - infilelen) as usize);
                }
                infile[infilelen as usize] = 0;
            };
        }
        /* -------------------------------------------------------------------- */
        /* if '*' was given for the output name expand it to the root file name */
        /* -------------------------------------------------------------------- */
        if !outfile.is_null() && *outfile.offset(0) == bb(b'*') {
            /* scan input name backwards to the first '/' character */
            let mut ii = jj - 1;
            while ii >= 0 {
                if infile[ii as usize] == bb(b'/') || ii == 0 {
                    if strlen_safe(&infile[((ii + 1) as usize)..]) > FLEN_FILENAME - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }
                    strcpy(outfile, infile[((ii + 1) as usize)..].as_ptr());
                    break;
                };
                ii -= 1
            }
        }
        /* ------------------------------------------ */
        /* copy strings from local copy to the output */
        /* ------------------------------------------ */
        if !infilex.is_null() {
            if strlen_safe(infile) > FLEN_FILENAME - 1 {
                *status = URL_PARSE_ERROR;
                return *status;
            }
            strcpy(infilex, infile.as_ptr());
        }
        /* ---------------------------------------------------------- */
        /* if no '[' character in the input string, then we are done. */
        /* ---------------------------------------------------------- */
        if ptr3.is_null() {
            return *status;
        }
        /* ------------------------------------------- */
        /* see if [ extension specification ] is given */
        /* ------------------------------------------- */
        if plus_ext == 0 {
            /* extension no. not already specified?  Then      */
            /* first brackets must enclose extension name or # */
            /* or it encloses a image subsection specification */
            /* or a raw binary image specifier */
            /* or a image compression specifier */
            /* Or, the extension specification may have been */
            /* omitted and we have to guess what the user intended */

            let mut ptr1 = ptr3.add(1); /* pointer to first char after the [ */
            let ptr2 = strchr(ptr1, bb(b']') as c_int); /* search for closing ] */
            if ptr2.is_null() {
                ffpmsg_str("input file URL is missing closing bracket ']'"); /* error, no closing ] */
                *status = URL_PARSE_ERROR;
                return *status;
            }
            /* ---------------------------------------------- */
            /* First, test if this is a rawfile specifier     */
            /* which looks something like: '[ib512,512:2880]' */
            /* Test if first character is b,i,j,d,r,f, or u,  */
            /* and optional second character is b or l,       */
            /* followed by one or more digits,                */
            /* finally followed by a ',', ':', or ']'         */
            /* ---------------------------------------------- */
            if *ptr1 == bb(b'b')
                || *ptr1 == bb(b'B')
                || *ptr1 == bb(b'i')
                || *ptr1 == bb(b'I')
                || *ptr1 == bb(b'j')
                || *ptr1 == bb(b'J')
                || *ptr1 == bb(b'd')
                || *ptr1 == bb(b'D')
                || *ptr1 == bb(b'r')
                || *ptr1 == bb(b'R')
                || *ptr1 == bb(b'f')
                || *ptr1 == bb(b'F')
                || *ptr1 == bb(b'u')
                || *ptr1 == bb(b'U')
            {
                /* next optional character may be a b or l (for Big or Little) */
                ptr1 = ptr1.offset(1); /* must have at least 1 digit */
                if *ptr1 == bb(b'b') || *ptr1 == bb(b'B') || *ptr1 == bb(b'l') || *ptr1 == bb(b'L')
                {
                    ptr1 = ptr1.offset(1); /* skip over digits */
                } /* OK, this looks like a rawfile specifier */

                if isdigit_safe(*ptr1) {
                    while isdigit_safe(*ptr1) {
                        ptr1 = ptr1.offset(1); /* append the raw array specifier to infilex */
                    } /* find the closing ] char */
                    if *ptr1 == bb(b',') || *ptr1 == bb(b':') || *ptr1 == bb(b']') {
                        if !urltype.is_null() {
                            if !strstr(urltype, cstr!(b"stdin").as_ptr()).is_null() {
                                strcpy(urltype, cstr!(b"rawstdin://").as_ptr()); /* terminate string after the ] */
                            } else {
                                strcpy(urltype, cstr!(b"rawfile://").as_ptr()); /* the 0 ext number is implicit */
                            }; /* search for another [ char */
                        } /* copy any remaining characters into rowfilterx  */
                        if !infilex.is_null() {
                            if strlen(infilex) + strlen(ptr3) > FLEN_FILENAME - 1 {
                                *status = URL_PARSE_ERROR; /* overwrite the ] with null terminator */
                                return *status; /* finished parsing, so return */
                            } /* end of rawfile specifier test */
                            strcat(infilex, ptr3);
                            let ptr1 = strchr(infilex, bb(b']') as c_int);
                            if !ptr1.is_null() {
                                *ptr1.offset(1) = 0;
                            };
                        }
                        if !extspec.is_null() {
                            strcpy(extspec, cast_slice(b"0\0").as_ptr());
                        }
                        let tmptr = strchr(ptr2.offset(1), bb(b'[') as c_int);
                        if !tmptr.is_null() && !rowfilterx.is_null() {
                            if strlen(rowfilterx) + strlen(tmptr.offset(1)) > FLEN_FILENAME - 1 {
                                *status = URL_PARSE_ERROR;
                                return *status;
                            }
                            strcat(rowfilterx, tmptr.offset(1));
                            let tmptr = strchr(rowfilterx, bb(b']') as c_int);
                            if !tmptr.is_null() {
                                *tmptr = 0;
                            };
                        }
                        return *status;
                    };
                };
            }
            /* -------------------------------------------------------- */
            /* Not a rawfile, so next, test if this is an image section */
            /* i.e., an integer followed by a ':' or a '*' or '-*'      */
            /* -------------------------------------------------------- */
            ptr1 = ptr3.offset(1); /* reset pointer to first char after the [ */
            let mut tmptr = ptr1; /* skip leading blanks */
            while *tmptr == bb(b' ') {
                tmptr = tmptr.offset(1); /* skip over leading digits */
            } /* this is an image section specifier */
            while isdigit_safe(*tmptr) {
                tmptr = tmptr.offset(1);
            }
            if *tmptr == bb(b':') || *tmptr == bb(b'*') || *tmptr == bb(b'-') {
                strcat(rowfilter.as_mut_ptr(), ptr3);
                /*
                don't want to assume 0 extension any more; may imply an image extension.
                    if (extspec)
                       strcpy(extspec, "0");
                */
            } else {
                /* -----------------------------------------------------------------
                Not an image section or rawfile spec so may be an extension spec.

                Examples of valid extension specifiers:
                   [3]                - 3rd extension; 0 = primary array
                   [events]           - events extension
                   [events, 2]        - events extension, with EXTVER = 2
                   [events,2]         - spaces are optional
                   [events, 3, b]     - same as above, plus XTENSION = 'BINTABLE'
                   [PICS; colName(12)] - an image in row 12 of the colName column
                                             in the PICS table extension
                   [PICS; colName(exposure > 1000)] - as above, but find image in
                                 first row with with exposure column value > 1000.
                   [Rate Table] - extension name can contain spaces!
                   [Rate Table;colName(exposure>1000)]

                Examples of other types of specifiers (Not extension specifiers)

                   [bin]  !!! this is ambiguous, and can't be distinguished from
                              a valid extension specifier
                   [bini X=1:512:16]  (also binb, binj, binr, and bind are allowed)
                   [binr (X,Y) = 5]
                   [bin @binfilter.txt]

                   [col Time;rate]
                   [col PI=PHA * 1.1]
                   [col -Time; status]

                   [X > 5]
                   [X>5]
                   [@filter.txt]
                   [StatusCol]  !!! this is ambiguous, and can't be distinguished
                              from a valid extension specifier
                   [StatusCol==0]
                   [StatusCol || x>6]
                   [gtifilter()]
                   [regfilter(cstr!(b"region.reg").as_ptr())]

                   [compress Rice]

                There will always be some ambiguity between an extension name and
                a boolean row filtering expression, (as in a couple of the above
                examples).  If there is any doubt, the expression should be treated
                as an extension specification;  The user can always add an explicit
                expression specifier to override this interpretation.

                The following decision logic will be used:

                1) locate the first token, terminated with a space, comma,
                   semi-colon, or closing bracket.

                2) the token is not part of an extension specifier if any of
                   the following is true:

                   - if the token begins with '@' and contains a '.'
                   - if the token contains an operator: = > < || &&
                   - if the token begins with cstr!(b"gtifilter(" or "regfilter(").as_ptr()
                   - if the token is terminated by a space and is followed by
                      additional characters (not a ']')  AND any of the following:
                        - the token is 'col'
                        - the token is 3 or 4 chars long and begins with 'bin'
                        - the second token begins with an operator:
                            ! = < > | & + - * / %


                3) otherwise, the string is assumed to be an extension specifier

                ----------------------------------------------------------------- */
                tmptr = ptr1; /* test for leading @ symbol */
                while *tmptr == bb(b' ') {
                    tmptr = tmptr.offset(1); /* parse the first token of the expression */
                } /* a space char? */
                hasAt = 0; /* skip spaces */
                hasDot = 0; /* is this the end? */
                hasOper = 0; /* 1st token is terminated by space */
                followingOper = 0; /* test if this is a column or binning specifier */
                spaceTerm = 0; /* check if next character is an operator */
                rowFilter = 0; /* test if this is NOT an extension specifier */
                colStart = 0; /* this is (probably) not an extension specifier */
                binStart = 0; /* so copy all chars to filter spec string */
                pixStart = 0; /* this appears to be a legit extension specifier */
                compStart = 0; /* copy the extension specification */
                if *tmptr == bb(b'@') {
                    hasAt = 1; /* copy any remaining chars to filter spec string */
                } /* end of  if (!plus_ext)     */

                let _tp = cast_slice(CStr::from_ptr(tmptr).to_bytes_with_nul());

                if fits_strncasecmp(_tp, cs!(b"col "), 4) == 0 {
                    colStart = 1;
                }
                if fits_strncasecmp(_tp, cs!(b"bin"), 3) == 0 {
                    binStart = 1;
                }
                if fits_strncasecmp(_tp, cs!(b"pix"), 3) == 0 {
                    pixStart = 1;
                }
                if fits_strncasecmp(_tp, cs!(b"compress "), 9) == 0
                    || fits_strncasecmp(_tp, cs!(b"compress]"), 9) == 0
                {
                    compStart = 1;
                }
                if fits_strncasecmp(_tp, cs!(b"gtifilter("), 10) == 0
                    || fits_strncasecmp(_tp, cs!(b"regfilter("), 10) == 0
                {
                    rowFilter = 1;
                } else {
                    let mut ii = 0;
                    let t = ptr2.offset_from(ptr1);
                    while ii < t + 1 {
                        if *tmptr == bb(b'.') {
                            hasDot = 1;
                        } else if *tmptr == bb(b'=')
                            || *tmptr == bb(b'>')
                            || *tmptr == bb(b'<')
                            || (*tmptr == bb(b'|') && *tmptr.offset(1) == bb(b'|'))
                            || (*tmptr == bb(b'&') && *tmptr.offset(1) == bb(b'&'))
                        {
                            hasOper = 1;
                        } else if *tmptr == bb(b',') || *tmptr == bb(b';') || *tmptr == bb(b']') {
                            break;
                        } else if *tmptr == bb(b' ') {
                            while *tmptr == bb(b' ') {
                                tmptr = tmptr.offset(1);
                            }
                            if *tmptr == bb(b']') {
                                break;
                            }
                            spaceTerm = 1;
                            if colStart != 0 || (ii <= 4 && (binStart != 0 || pixStart != 0)) {
                                rowFilter = 1;
                            } else if *tmptr == bb(b'=')
                                || *tmptr == bb(b'>')
                                || *tmptr == bb(b'<')
                                || *tmptr == bb(b'|')
                                || *tmptr == bb(b'&')
                                || *tmptr == bb(b'!')
                                || *tmptr == bb(b'+')
                                || *tmptr == bb(b'-')
                                || *tmptr == bb(b'*')
                                || *tmptr == bb(b'/')
                                || *tmptr == bb(b'%')
                            {
                                followingOper = 1;
                            }
                            break;
                        };
                        {
                            ii += 1;
                            tmptr = tmptr.offset(1)
                        }
                    }
                }
                if rowFilter != 0
                    || (pixStart != 0 && spaceTerm != 0)
                    || (hasAt != 0 && hasDot != 0)
                    || hasOper != 0
                    || compStart != 0
                    || (spaceTerm != 0 && followingOper != 0)
                {
                    strcat(rowfilter.as_mut_ptr(), ptr3);
                } else {
                    if !extspec.is_null() {
                        let t = ptr2.offset_from(ptr1);
                        assert!(t >= 0, "Pointer math is bad mkay");
                        if t as usize > FLEN_FILENAME - 1 {
                            *status = URL_PARSE_ERROR;
                            return *status;
                        }

                        strncat(extspec, ptr1, t as usize);
                    }
                    strcat(rowfilter.as_mut_ptr(), ptr2.offset(1));
                };
            };
        } else {
            /* ------------------------------------------------------------------ */
            /* already have extension, so this must be a filter spec of some sort */
            /* ------------------------------------------------------------------ */
            strcat(rowfilter.as_mut_ptr(), ptr3);
        }

        /* strip off any trailing blanks from filter */
        slen = strlen_safe(rowfilter) as isize;
        slen -= 1;
        while slen >= 0 && rowfilter[slen as usize] == bb(b' ') {
            rowfilter[slen as usize] = 0;
        }

        if rowfilter[0] == 0 {
            return *status; /* nothing left to parse */
        }

        /* ------------------------------------------------ */
        /* does the filter contain a binning specification? */
        /* ------------------------------------------------ */

        let mut ptr1 = strstr_safe(rowfilter, cs!(b"[bin")); /* search for "[bin" */

        if ptr1.is_none() {
            ptr1 = strstr_safe(rowfilter, cs!(b"[BIN")); /* search for "[BIN" */
        }

        if ptr1.is_none() {
            ptr1 = strstr_safe(rowfilter, cs!(b"[Bin")); /* search for "[Bin" */
        }

        if ptr1.is_some() {
            let p1 = ptr1.unwrap();
            let mut p2 = p1 + 4; /* end of the '[bin' string */
            if rowfilter[p2] == bb(b'b')
                || rowfilter[p2] == bb(b'i')
                || rowfilter[p2] == bb(b'j')
                || rowfilter[p2] == bb(b'r')
                || rowfilter[p2] == bb(b'd')
            {
                p2 += 1; /* skip the datatype code letter */
            }
            if rowfilter[p2] != bb(b' ') && rowfilter[p2] != bb(b']') {
                ptr1 = None; /* bin string must be followed by space or ] */
            };
        }

        if ptr1.is_some() {
            let p1: usize = ptr1.unwrap();

            /* found the binning string */
            if !binspec.is_null() {
                if strlen_safe(&rowfilter[(p1 + 1)..]) > FLEN_FILENAME - 1 {
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                strcpy(binspec, rowfilter[(p1 + 1)..].as_ptr());
                let _bs = CStr::from_ptr(binspec);
                let _bs_len = _bs.to_bytes_with_nul().len();
                let binspec = slice::from_raw_parts_mut(binspec, _bs_len);

                let ptr2 = fits_find_match_delim(binspec, bb(b']'));
                if let Some(mut p2) = ptr2 {
                    p2 -= 1;
                    binspec[p2] = 0;

                    p2 -= 1;
                    if binspec[p2] == bb(b' ') {
                        binspec[p2] = 0;
                    };
                } else {
                    ffpmsg_str("input file URL is missing closing bracket ']'");
                    ffpmsg_slice(rowfilter);
                    *status = URL_PARSE_ERROR;
                    return *status;
                };
            }

            /* delete the binning spec from the row filter string */
            let ptr2 = fits_find_match_delim(&mut rowfilter[(p1 + 1)..], bb(b']'));
            if let Some(mut p2) = ptr2 {
                p2 = p2 + p1 + 1;

                strcpy_safe(tmpstr, &rowfilter[p2..]); /* copy any chars after the binspec */
                strcpy_safe(&mut rowfilter[p1..], tmpstr); /* overwrite binspec */
            } else {
                ffpmsg_str("input file URL is missing closing bracket ']'");
                ffpmsg_slice(rowfilter);
                *status = URL_PARSE_ERROR; /* error, no closing ] */
                return *status;
            };
        }

        /* --------------------------------------------------------- */
        /* does the filter contain a column selection specification? */
        /* --------------------------------------------------------- */

        let mut ptr1 = strstr_safe(rowfilter, cs!(b"[col "));
        if ptr1.is_none() {
            ptr1 = strstr_safe(rowfilter, cs!(b"[COL "));
            if ptr1.is_none() {
                ptr1 = strstr_safe(rowfilter, cs!(b"[Col "));
            };
        }

        hasAt = 0;

        while ptr1.is_some() {
            let p1 = ptr1.unwrap();

            /* find the end of the column specifier */
            let mut p2 = p1 + 5;
            /* Scan past any whitespace and check for @filename */
            while rowfilter[p2] == bb(b' ') {
                p2 += 1; /* error, no closing ] */
            } /* start of a literal string */

            if rowfilter[p2] == bb(b'@') {
                hasAt = 1; /* find closing quote */
            } /* error, no closing ] */

            while rowfilter[p2] != bb(b']') {
                if rowfilter[p2] == 0 {
                    ffpmsg_str("input file URL is missing closing bracket ']'"); /* set of nested square brackets */
                    /* find closing bracket */
                    *status = URL_PARSE_ERROR; /* error, no closing ] */
                    return *status; /* continue search for the closing bracket character */
                } /* copy the column specifier to output string */

                if rowfilter[p2] == bb(b'\\') {
                    let _ptr2 = strchr_safe(&rowfilter[(p2 + 1)..], bb(b'\\')); /* Pre-existing colspec, append with ";" */
                    if _ptr2.is_none() {
                        ffpmsg_str(
                            "literal string in input file URL is missing closing single quote",
                        );
                        *status = URL_PARSE_ERROR;
                        return *status;
                    } else {
                        p2 += _ptr2.unwrap();
                    }
                }
                if rowfilter[p2] == bb(b'[') {
                    let _ptr2 = strchr_safe(&rowfilter[(p2 + 1)..], bb(b']'));
                    if _ptr2.is_none() {
                        ffpmsg_str("nested brackets in input file URL is missing closing bracket");
                        *status = URL_PARSE_ERROR;
                        return *status;
                    } else {
                        p2 += _ptr2.unwrap();
                    }
                }

                p2 += 1;
            }

            let collen = p2 - p1 - 1;

            if !colspec.is_null() {
                /* copy the column specifier to output string */

                if collen + strlen(colspec) > FLEN_FILENAME - 1 {
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                if *colspec == 0 {
                    strncpy(colspec, rowfilter[(p1 + 1)..].as_ptr(), collen);
                    *colspec.add(collen) = 0;
                } else {
                    strcat(colspec, cstr!(b";").as_ptr());
                    strncat(colspec, rowfilter[(p1 + 5)..].as_ptr(), collen - 4);
                    /* Note that strncat always null-terminates the destination string */

                    /* Special error checking here.  We can't allow there to be a
                    col @filename.txt includes if there are multiple col expressions */
                    if hasAt != 0 {
                        ffpmsg_str("input URL multiple column filter cannot use @filename.txt");
                        *status = URL_PARSE_ERROR;
                        return *status;
                    };
                }

                let mut collen = strlen(colspec) as isize;
                collen -= 1;
                while *colspec.offset(collen) == bb(b' ') {
                    *colspec.offset(collen) = 0; /* strip trailing blanks */
                }
            }

            /* delete the column selection spec from the row filter string */
            strcpy_safe(tmpstr, &rowfilter[(p2 + 1)..]); /* copy any chars after the colspec */
            strcpy_safe(&mut rowfilter[p1..], tmpstr); /* overwrite binspec */

            /* Check for additional column specifiers */
            ptr1 = strstr_safe(rowfilter, cs!(b"[col "));
            if ptr1.is_none() {
                ptr1 = strstr_safe(rowfilter, cs!(b"[COL "));
            }
            if ptr1.is_none() {
                ptr1 = strstr_safe(rowfilter, cs!(b"[Col "));
            };
        }

        /* --------------------------------------------------------- */
        /* does the filter contain a pixel filter specification?     */
        /* --------------------------------------------------------- */

        let mut ptr1 = strstr_safe(rowfilter, cs!(b"[pix"));
        if ptr1.is_none() {
            ptr1 = strstr_safe(rowfilter, cs!(b"[PIX"));
            if ptr1.is_none() {
                ptr1 = strstr_safe(rowfilter, cs!(b"[Pix"));
            };
        }

        let mut p2 = 0;
        if ptr1.is_some() {
            let p1 = ptr1.unwrap();
            p2 = p1 + 4; /* end of the '[pix' string */

            if rowfilter[p2] == bb(b'b')
                || rowfilter[p2] == bb(b'i')
                || rowfilter[p2] == bb(b'j')
                || rowfilter[p2] == bb(b'B')
                || rowfilter[p2] == bb(b'I')
                || rowfilter[p2] == bb(b'J')
                || rowfilter[p2] == bb(b'r')
                || rowfilter[p2] == bb(b'd')
                || rowfilter[p2] == bb(b'R')
                || rowfilter[p2] == bb(b'D')
            {
                p2 += 1 /* skip the datatype code letter */
            }
            if rowfilter[p2] == bb(b'1') {
                p2 += 1
            }
            if rowfilter[p2] != bb(b' ') {
                ptr1 = None;
            };
        }

        if ptr1.is_some() {
            let p1 = ptr1.unwrap();

            while rowfilter[p2] != bb(b']') {
                if rowfilter[p2] == 0 {
                    ffpmsg_str("input file URL is missing closing bracket ']'"); /* copy the column specifier to output string */

                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                if rowfilter[p2] == bb(b'\\') {
                    let _ptr2 = strchr_safe(&rowfilter[(p2 + 1)..], bb(b'\\'));
                    if let Some(_ptr2) = _ptr2 {
                        p2 += _ptr2;
                    } else {
                        ffpmsg_str(
                            "literal string in input file URL is missing closing single quote",
                        );
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }
                }

                if rowfilter[p2] == bb(b'[') {
                    let _ptr2 = strchr_safe(&rowfilter[(p2 + 1)..], bb(b']'));
                    if _ptr2.is_none() {
                        ffpmsg_str("nested brackets in input file URL is missing closing bracket");
                        *status = URL_PARSE_ERROR;
                        return *status;
                    } else {
                        p2 += _ptr2.unwrap();
                    }
                }

                p2 += 1;
            }

            let mut collen = p2 - p1 - 1;

            if !pixfilter.is_null() {
                if collen as usize > FLEN_FILENAME - 1 {
                    *status = URL_PARSE_ERROR;
                    return *status;
                }
                strncpy(pixfilter, rowfilter[(p1 + 1)..].as_ptr(), collen);
                *pixfilter.add(collen) = 0;

                collen -= 1;
                while *pixfilter.add(collen) == bb(b' ') {
                    *pixfilter.add(collen) = 0;
                }
            }
            /* delete the pixel filter from the row filter string */
            strcpy_safe(tmpstr, &rowfilter[(p2 + 1)..]); /* copy any chars after the pixel filter */
            strcpy_safe(&mut rowfilter[p1..], tmpstr); /* overwrite binspec */
        }

        /* ------------------------------------------------------------ */
        /* does the filter contain an image compression specification?  */
        /* ------------------------------------------------------------ */

        let mut ptr1 = strstr_safe(rowfilter, cs!(b"[compress")); /* end of the '[compress' string */

        if ptr1.is_some() {
            let p2 = ptr1.unwrap() + 9; /* compress string must be followed by space or ] */
            if rowfilter[p2] != bb(b' ') && rowfilter[p2] != bb(b']') {
                ptr1 = None;
            };
        }

        if ptr1.is_some() {
            let p1 = ptr1.unwrap();

            /* found the compress string */
            if !compspec.is_null() {
                if strlen_safe(&rowfilter[(p1 + 1)..]) > FLEN_FILENAME - 1 {
                    *status = URL_PARSE_ERROR; /* delete trailing spaces */
                    return *status; /* error, no closing ] */
                }

                strcpy(compspec, rowfilter[(p1 + 1)..].as_ptr());

                let _ptr2 = strchr(compspec, bb(b']') as c_int);
                if !ptr2.is_null() {
                    *ptr2 = 0;
                    if *{
                        ptr2 = ptr2.offset(-1);
                        ptr2
                    } == bb(b' ')
                    {
                        *ptr2 = 0;
                    };
                } else {
                    ffpmsg_str("input file URL is missing closing bracket ']'");
                    ffpmsg_slice(rowfilter);
                    *status = URL_PARSE_ERROR;
                    return *status;
                };
            }

            /* delete the compression spec from the row filter string */
            let ptr2 = strchr_safe(&rowfilter[p1..], bb(b']'));
            strcpy_safe(tmpstr, &rowfilter[(p1 + ptr2.unwrap() + 1)..]); /* copy any chars after the binspec */
            strcpy_safe(&mut rowfilter[p1..], tmpstr); /* overwrite binspec */
        }

        /* copy the remaining string to the rowfilter output... should only */
        /* contain a rowfilter expression of the form cstr!(b"[expr]").as_ptr()
         */
        if !rowfilterx.is_null() && rowfilter[0] != 0 {
            hasAt = 0;

            /* Check for multiple expressions, which would appear as cstr!(b"[expr][expr]...").as_ptr() */
            let mut p1 = 0; // rowfilter;
            let mut p2 = strstr_safe(rowfilter, cs!(b"][")).unwrap() - p1;

            while (rowfilter[p1] == bb(b'[')) && p2 > 2 {
                /* Advance past any white space */
                let mut p3 = p1 + 1;
                while rowfilter[p3] == bb(b' ') {
                    p3 += 1;
                }

                /* Check for @filename.txt */
                if rowfilter[p3] == bb(b'@') {
                    hasAt = 1;
                }

                /* Add expression of the form cstr!(b"((expr))&&").as_ptr(), note the addition of 6 characters */
                if (strlen(rowfilterx) + (p2 - p1) + 6) > FLEN_FILENAME - 1 {
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                /* Special error checking here.  We can't allow there to be a
                @filename.txt includes if there are multiple row expressions */
                if *rowfilterx != 0 && hasAt != 0 {
                    ffpmsg_str("input URL multiple row filter cannot use @filename.txt");
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                /* Append the expression */
                strcat(rowfilterx, cstr!(b"((").as_ptr());
                strncat(rowfilterx, rowfilter[(p1 + 1)..].as_ptr(), p2 - p1 - 1);
                /* Note that strncat always null-terminates the destination string */
                strcat(rowfilterx, cstr!(b"))&&").as_ptr());

                /* Advance to next expression */
                p1 = p2 + 1;

                p2 = strstr_safe(rowfilter, cs!(b"][")).unwrap() - p1;
            }

            /* At final iteration, ptr1 points to beginning [ and ptr2 to ending ] */
            let p2 = strlen_safe(rowfilter) - 1;
            if rowfilter[p1] == bb(b'[') && rowfilter[p2] == bb(b']') {
                /* Check for @include in final position */
                let mut p3 = p1 + 1;
                while rowfilter[p3] == bb(b' ') {
                    p3 += 1;
                }

                if rowfilter[p3] == bb(b'@') {
                    hasAt = 1;
                }

                /* Check for overflow; add extra 4 characters if we have pre-existing expression */
                if strlen(rowfilterx) + (p2 - p1 + (if *rowfilterx != 0 { 4 } else { 0 }))
                    > FLEN_FILENAME - 1
                {
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                /* Special error checking here.  We can't allow there to be a
                @filename.txt includes if there are multiple row expressions */
                if *rowfilterx != 0 && hasAt != 0 {
                    ffpmsg_str("input URL multiple row filter cannot use @filename.txt");

                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                if *rowfilterx != 0 {
                    /* A pre-existing row filter: we bracket by ((expr)) to be sure */
                    strcat(rowfilterx, cstr!(b"((").as_ptr());

                    strncat(rowfilterx, rowfilter[(p1 + 1)..].as_ptr(), p2 - p1 - 1);
                    strcat(rowfilterx, cstr!(b"))").as_ptr());
                } else {
                    /* We have only one filter, so just copy the expression alone.
                    This will be the most typical case */

                    strncat(rowfilterx, rowfilter[(p1 + 1)..].as_ptr(), p2 - p1 - 1);
                };
            } else {
                ffpmsg_str("input file URL lacks valid row filter expression");
                *status = URL_PARSE_ERROR;
            };
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// test if the input file specifier is an existing file on disk
/// If the specified file can't be found, it then searches for a
/// compressed version of the file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffexist(
    infile: *const c_char, /* I - input filename or URL */
    exists: *mut c_int,    /* O -  2 = a compressed version of file exists */
    /*      1 = yes, disk file exists               */
    /*      0 = no, disk file could not be found    */
    /*     -1 = infile is not a disk file (could    */
    /*   be a http, ftp, gsiftp, smem, or stdin file) */
    status: *mut c_int, /* I/O  status  */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// parse the input URL, returning the root name (filetype://basename).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffrtnm(
    url: *const c_char,
    rootname: *mut c_char,
    status: *mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// parse the output URL into its basic components.
pub(crate) unsafe fn ffourl(
    url: &[c_char],        /* I - full input URL   */
    urltype: *mut c_char,  /* O - url type         */
    outfile: *mut c_char,  /* O - base file name   */
    tpltfile: *mut c_char, /* O - template file name, if any */
    compspec: *mut c_char, /* O - compression specification, if any */
    status: *mut c_int,
) -> c_int {
    unsafe {
        if *status > 0 {
            return *status;
        }
        if !urltype.is_null() {
            *urltype = 0;
        }
        if !outfile.is_null() {
            *outfile = 0;
        }
        if !tpltfile.is_null() {
            *tpltfile = 0;
        }
        if !compspec.is_null() {
            *compspec = 0;
        }

        let mut ptr1 = url; // url

        while ptr1[0] == bb(b' ') {
            /* ignore leading blanks */
            ptr1 = &ptr1[1..];
        }

        if ((ptr1[0] == bb(b'-')) && (ptr1[1] == 0 || ptr1[1] == bb(b' ')))
            || strcmp_safe(ptr1, cs!("stdout")) == 0
            || strcmp_safe(ptr1, cs!("STDOUT")) == 0
        /* "-" means write to stdout;  also support "- "            */
        /* but exclude disk file names that begin with a minus sign */
        /* e.g., "-55d33m.fits"   */
        {
            if !urltype.is_null() {
                strcpy(urltype, cstr!("stdout://").as_ptr());
            }
        } else {
            /* not writing to stdout */
            /*  get urltype (e.g., file://, ftp://, http://, etc.)  */

            let ptr2 = strstr_safe(ptr1, cs!("://"));

            if let Some(ptr2) = ptr2 {
                /* copy the explicit urltype string */

                if !urltype.is_null() {
                    if ptr2 + 3 > MAX_PREFIX_LEN - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }

                    strncat(urltype, ptr1.as_ptr(), ptr2 + 3);
                }

                ptr1 = &ptr1[(ptr2 + 3)..];
            } else {
                /* assume file driver    */

                if !urltype.is_null() {
                    strcat(urltype, cstr!("file://").as_ptr());
                }
            }

            /* look for template file name, enclosed in parenthesis */
            let ptr2 = strchr_safe(ptr1, bb(b'('));

            /* look for image compression parameters, enclosed in sq. brackets */
            let ptr3 = strchr_safe(ptr1, bb(b'['));

            if !outfile.is_null() {
                if let Some(ptr2) = ptr2 {
                    /* template file was specified  */
                    if ptr2 > FLEN_FILENAME - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }

                    strncat(outfile, ptr1.as_ptr(), ptr2);
                } else if let Some(ptr3) = ptr3 {
                    /* compression was specified  */
                    if ptr3 > FLEN_FILENAME - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }

                    strncat(outfile, ptr1.as_ptr(), ptr3);
                } else {
                    /* no template file or compression */
                    if strlen_safe(ptr1) > FLEN_FILENAME - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }

                    strcpy(outfile, ptr1.as_ptr());
                }
            }

            if let Some(mut ptr2) = ptr2 {
                /* template file was specified  */

                ptr2 += 1;

                let tmp_ptr1 = strchr_safe(&ptr1[ptr2..], bb(b')')); /* search for closing ) */

                if tmp_ptr1.is_none() {
                    *status = URL_PARSE_ERROR; /* error, no closing ) */
                    return *status;
                }

                let tmp_ptr1 = tmp_ptr1.unwrap();
                ptr1 = &ptr1[(ptr2 + tmp_ptr1)..];

                if !tpltfile.is_null() {
                    if tmp_ptr1 > FLEN_FILENAME - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }
                    strncat(tpltfile, ptr1[ptr2..].as_ptr(), tmp_ptr1);
                }
            }

            if let Some(mut ptr3) = ptr3 {
                /* compression was specified  */

                ptr3 += 1;

                let tmp_ptr1 = strchr_safe(&ptr1[ptr3..], bb(b']')); /* search for closing ] */

                if tmp_ptr1.is_none() {
                    *status = URL_PARSE_ERROR; /* error, no closing ] */
                    return *status;
                }

                let tmp_ptr1 = tmp_ptr1.unwrap();
                ptr1 = &ptr1[(ptr3 + tmp_ptr1)..];

                if !compspec.is_null() {
                    if tmp_ptr1 > FLEN_FILENAME - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }

                    strncat(compspec, ptr1[ptr3..].as_ptr(), tmp_ptr1);
                }
            }

            /* check if a .gz compressed output file is to be created */
            /* by seeing if the filename ends in '.gz'   */
            if !urltype.is_null() && !outfile.is_null() {
                raw_to_slice!(urltype);
                raw_to_slice!(outfile);

                if strcmp_safe(urltype, cs!("file://")) == 0 {
                    let ptr1 = strstr_safe(outfile, cs!(".gz"));
                    if let Some(mut ptr1) = ptr1 {
                        /* make sure the ".gz" is at the end of the file name */
                        ptr1 += 3;
                        if outfile[ptr1] == 0 || outfile[ptr1] == bb(b' ') {
                            strcpy(
                                urltype.as_ptr() as *mut _,
                                cstr!("compressoutfile://").as_ptr(),
                            );
                        }
                    }
                }
            }
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Parse the input extension specification string, returning either the
/// extension number or the values of the EXTNAME, EXTVERS, and XTENSION
/// keywords in desired extension. Also return the name of the column containing
/// an image, and an expression to be used to determine which row to use,
/// if present.
///
/// DANGER: Don't know size of extname, imagecolname, rowexpress
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffexts(
    extspec: *const c_char,
    extnum: *mut c_int,
    extname: *mut c_char,
    extvers: *mut c_int,
    hdutype: *mut c_int,
    imagecolname: *mut c_char,
    rowexpress: *mut c_char,
    status: *mut c_int,
) -> c_int {
    unsafe {
        let extnum = extnum.as_mut().expect(NULL_MSG);
        let extvers = extvers.as_mut().expect(NULL_MSG);
        let hdutype = hdutype.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(extspec);

        // Don't assign, using this to test for null pointer.
        extname.as_mut().expect(NULL_MSG);
        imagecolname.as_mut().expect(NULL_MSG);
        rowexpress.as_mut().expect(NULL_MSG);

        ffexts_safer(
            extspec,
            extnum,
            extname,
            extvers,
            hdutype,
            imagecolname,
            rowexpress,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Parse the input extension specification string, returning either the
/// extension number or the values of the EXTNAME, EXTVERS, and XTENSION
/// keywords in desired extension. Also return the name of the column containing
/// an image, and an expression to be used to determine which row to use,
/// if present.
///
/// DANGER: Don't know size of extname, imagecolname, rowexpress
pub(crate) unsafe fn ffexts_safer(
    extspec: &[c_char],
    extnum: &mut c_int,
    extname: *mut c_char,
    extvers: &mut c_int,
    hdutype: &mut c_int,
    imagecolname: *mut c_char,
    rowexpress: *mut c_char,
    status: &mut c_int,
) -> c_int {
    unsafe {
        let mut slen: usize = 0;
        let mut nvals: c_int = 0;
        let mut notint: c_int = 1;
        let mut tmpname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

        *extnum = 0;
        *extvers = 0;
        *hdutype = ANY_HDU;

        *extname = 0;
        *imagecolname = 0;
        *rowexpress = 0;

        // Use these as intermediate variables then copy into the FFI variables
        let mut _extname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut _imagecolname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut _rowexpress: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];

        if *status > 0 {
            return *status;
        }

        let mut ptr1 = 0; // ptr to extspec

        while extspec[ptr1] == bb(b' ') {
            /* skip over any leading blanks */
            ptr1 += 1;
        }

        /* is the extension specification a number? */
        if isdigit_safe(extspec[ptr1]) {
            notint = 0; /* looks like extname may actually be the ext. number */

            set_errno(Errno(0)); /* reset this prior to calling strtol */

            let mut loc = 0;
            *extnum = strtol_safe(&extspec[ptr1..], &mut loc, 10) as c_int; /* read the string as an integer */

            let mut loc = ptr1 + loc;
            while extspec[loc] == bb(b' ') {
                /* skip over trailing blanks */
                loc += 1;
            }

            /* check for read error, or junk following the integer */
            if (extspec[loc] != 0 && extspec[loc] != bb(b';')) || (errno().0 == ERANGE) {
                *extnum = 0;
                notint = 1; /* no, extname was not a simple integer after all */
                set_errno(Errno(0)); /* reset error condition flag if it was set */
            }

            if *extnum < 0 || *extnum > 99999 {
                *extnum = 0; /* this is not a reasonable extension number */
                ffpmsg_str("specified extension number is out of range:");
                ffpmsg_slice(extspec);
                *status = URL_PARSE_ERROR;
                return *status;
            }
        }

        /*  This logic was too simple, and failed on extnames like '1000TEMP'
            where it would try to move to the 1000th extension

            if (isdigit((int) extspec[ptr1]))
            {
                sscanf(ptr1, "%d", extnum);
                if (*extnum < 0 || *extnum > 9999)
                {
                    *extnum = 0;
                    ffpmsg("specified extension number is out of range:");
                    ffpmsg_slice(extspec);
                    *status = URL_PARSE_ERROR;
                    return *status;
                }
            }
        */

        if notint != 0 {
            /* not a number, so EXTNAME must be specified, followed by */
            /* optional EXTVERS and XTENSION  values */

            /* don't use space char as end indicator, because there */
            /* may be imbedded spaces in the EXTNAME value */
            slen = strcspn_safe(&extspec[ptr1..], cs!(",:;")); /* length of EXTNAME */

            if slen > FLEN_VALUE - 1 {
                *status = URL_PARSE_ERROR;
                return *status;
            }

            strncat_safe(&mut _extname, &extspec[ptr1..], slen); /* EXTNAME value */

            /* now remove any trailing blanks */
            while slen > 0 && (_extname[slen - 1]) == bb(b' ') {
                (_extname[slen - 1]) = 0;
                slen -= 1;
            }

            ptr1 += slen;
            slen = strspn_safe(&extspec[ptr1..], cs!(" ,:")); /* skip delimiter characters */
            ptr1 += slen;

            slen = strcspn_safe(&extspec[ptr1..], cs!(" ,:;")); /* length of EXTVERS */
            if slen != 0 {
                nvals = sscanf(extspec[ptr1..].as_ptr(), cstr!("%d").as_ptr(), extvers); /* EXTVERS value */
                if nvals != 1 {
                    ffpmsg_str("illegal EXTVER value in input URL:");
                    ffpmsg_slice(extspec);
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                ptr1 += slen;
                slen = strspn_safe(&extspec[ptr1..], cs!(" ,:")); /* skip delimiter characters */
                ptr1 += slen;

                slen = strcspn_safe(&extspec[ptr1..], cs!(";")); /* length of HDUTYPE */
                if slen != 0 {
                    if extspec[ptr1] == bb(b'b') || extspec[ptr1] == bb(b'B') {
                        *hdutype = BINARY_TBL;
                    } else if extspec[ptr1] == bb(b't')
                        || extspec[ptr1] == bb(b'T')
                        || extspec[ptr1] == bb(b'a')
                        || extspec[ptr1] == bb(b'A')
                    {
                        *hdutype = ASCII_TBL;
                    } else if extspec[ptr1] == bb(b'i') || extspec[ptr1] == bb(b'I') {
                        *hdutype = IMAGE_HDU;
                    } else {
                        ffpmsg_str("unknown type of HDU in input URL:");
                        ffpmsg_slice(extspec);
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }
                }
            } else {
                strcpy_safe(&mut tmpname, &_extname);
                ffupch_safe(&mut tmpname);
                if strcmp_safe(&tmpname, cs!("PRIMARY")) == 0
                    || strcmp_safe(&tmpname, cs!("P")) == 0
                {
                    _extname[0] = 0; /* return extnum = 0 */
                }
            }
        }

        let ptr1 = strchr_safe(&extspec[ptr1..], bb(b';'));

        if let Some(mut ptr1) = ptr1 {
            /* an image is to be opened; the image is contained in a single */
            /* cell of a binary table.  A column name and an expression to  */
            /* determine which row to use has been entered.                 */

            ptr1 += 1; /* skip over the ';' delimiter */
            while extspec[ptr1] == bb(b' ') {
                /* skip over any leading blanks */
                ptr1 += 1;
            }

            let ptr2 = strchr_safe(&extspec[ptr1..], bb(b'('));
            match ptr2 {
                None => {
                    ffpmsg_str("illegal specification of image in table cell in input URL:");
                    ffpmsg_str(" did not find a row expression enclosed in ( )");
                    ffpmsg_slice(extspec);
                    *status = URL_PARSE_ERROR;
                    return *status;
                }
                Some(mut ptr2) => {
                    ptr2 += ptr1;
                    if ptr2 - ptr1 > FLEN_FILENAME - 1 {
                        *status = URL_PARSE_ERROR;
                        return *status;
                    }

                    strncat_safe(&mut _imagecolname, &extspec[ptr1..], ptr2 - ptr1); /* copy column name */

                    ptr2 += 1; /* skip over the '(' delimiter */
                    while extspec[ptr2] == bb(b' ') {
                        /* skip over any leading blanks */
                        ptr2 += 1;
                    }

                    let ptr1 = strchr_safe(&extspec[ptr2..], bb(b')'));
                    match ptr1 {
                        None => {
                            ffpmsg_str(
                                "illegal specification of image in table cell in input URL:",
                            );
                            ffpmsg_str(" missing closing ')' character in row expression");
                            ffpmsg_slice(extspec);
                            *status = URL_PARSE_ERROR;
                            return *status;
                        }
                        Some(ptr1) => {
                            if ptr1 - ptr2 > FLEN_FILENAME - 1 {
                                *status = URL_PARSE_ERROR;
                                return *status;
                            }

                            strncat_safe(&mut _rowexpress, &extspec[ptr2..], ptr1 - ptr2);
                            /* row expression */
                        }
                    }
                }
            }
        }

        strcpy(extname, _extname.as_ptr());
        strcpy(imagecolname, _imagecolname.as_ptr());
        strcpy(rowexpress, _rowexpress.as_ptr());

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Parse the input url string and return the number of the extension that
/// CFITSIO would automatically move to if CFITSIO were to open this input URL.
/// The extension numbers are one's based, so 1 = the primary array, 2 = the
/// first extension, etc.
///
/// The extension number that gets returned is determined by the following
/// algorithm:
///
/// 1. If the input URL includes a binning specification (e.g.
///    'myfile.fits[3][bin X,Y]') then the returned extension number
///    will always = 1, since CFITSIO would create a temporary primary
///    image on the fly in this case.  The same is true if an image
///    within a single cell of a binary table is opened.
///
/// 2.  Else if the input URL specifies an extension number (e.g.,
///     'myfile.fits[3]' or 'myfile.fits+3') then the specified extension
///     number (+ 1) is returned.  
///
/// 3.  Else if the extension name is specified in brackets
///     (e.g., this 'myfile.fits[EVENTS]') then the file will be opened and searched
///     for the extension number.  If the input URL is '-'  (reading from the stdin
///     file stream) this is not possible and an error will be returned.
///
/// 4.  Else if the URL does not specify an extension (e.g. 'myfile.fits') then
///     a special extension number = -99 will be returned to signal that no
///     extension was specified.  This feature is mainly for compatibility with
///     existing FTOOLS software.  CFITSIO would open the primary array by default
///     (extension_num = 1) in this case.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffextn(
    url: *const c_char,        /* I - input filename/URL  */
    extension_num: *mut c_int, /* O - returned extension number */
    status: *mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// return the prefix string associated with the driver in use by the
/// fitsfile pointer fptr
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffurlt(
    fptr: *mut fitsfile,
    urlType: *mut c_char,
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_ref().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let urlType = slice::from_raw_parts_mut(urlType, MAX_PREFIX_LEN);

        ffurlt_safe(fptr, urlType, status)
    }
}

/*--------------------------------------------------------------------------*/
/// return the prefix string associated with the driver in use by the
/// fitsfile pointer fptr
pub(crate) fn ffurlt_safe(fptr: &fitsfile, urlType: &mut [c_char], status: &mut c_int) -> c_int {
    let d = DRIVER_TABLE.get().unwrap();

    strcpy_safe(urlType, &(d[fptr.Fptr.driver as usize]).prefix);

    *status
}

/*--------------------------------------------------------------------------*/
/// Read and concatenate all the lines from the given text file.  User
/// must free the pointer returned in contents.  Pointer is guaranteed
/// to hold 2 characters more than the length of the text... allows the
/// calling routine to append (or prepend) a newline (or quotes?) without
/// reallocating memory.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffimport_file(
    filename: *const c_char,    /* Text file to read                   */
    contents: *mut *mut c_char, /* Pointer to pointer to hold file     */
    status: *mut c_int,         /* CFITSIO error code                  */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        raw_to_slice!(filename);

        ffimport_file_safer(filename, contents, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read and concatenate all the lines from the given text file.  User
/// must free the pointer returned in contents.  Pointer is guaranteed
/// to hold 2 characters more than the length of the text... allows the
/// calling routine to append (or prepend) a newline (or quotes?) without
/// reallocating memory.
pub(crate) unsafe fn ffimport_file_safer(
    filename: &[c_char],        /* Text file to read                   */
    contents: *mut *mut c_char, /* Pointer to pointer to hold file     */
    status: &mut c_int,         /* CFITSIO error code                  */
) -> c_int {
    unsafe {
        let mut eoline = true;
        let mut line: [c_char; 256] = [0; 256];

        if *status > 0 {
            return *status;
        }

        let mut totalLen = 0;
        let mut allocLen = 1024;

        // HEAP ALLOCATION
        let mut lines = Vec::new();
        if lines.try_reserve_exact(allocLen).is_err() {
            ffpmsg_str("Couldn't allocate memory to hold ASCII file contents.");
            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            lines.resize(allocLen, 0);
        }
        lines[0] = 0;

        let aFile = fopen(filename.as_ptr(), cstr!("r").as_ptr());
        if aFile.is_null() {
            int_snprintf!(
                &mut line,
                256,
                "Could not open ASCII file {}.",
                slice_to_str!(&filename),
            );
            ffpmsg_slice(&line);

            *status = FILE_NOT_OPENED;
            return *status;
        }

        while !fgets(line.as_mut_ptr(), 256, aFile).is_null() {
            let mut llen = strlen_safe(&line);
            if eoline && (llen > 1) && (line[0] == bb(b'/') && line[1] == bb(b'/')) {
                continue; /* skip comment lines begging with // */
            }

            eoline = false;

            /* replace CR and newline chars at end of line with nulls */
            if (llen > 0) && (line[llen - 1] == bb(b'\n') || line[llen - 1] == bb(b'\r')) {
                llen -= 1;
                line[llen] = 0;
                eoline = true; /* found an end of line character */

                if (llen > 0) && (line[llen - 1] == bb(b'\n') || line[llen - 1] == bb(b'\r')) {
                    llen -= 1;
                    line[llen] = 0;
                }
            }

            if totalLen + llen + 3 >= allocLen {
                allocLen += 256;

                if lines.try_reserve_exact(256).is_err() {
                    ffpmsg_str("Couldn't allocate memory to hold ASCII file contents.");
                    *status = MEMORY_ALLOCATION;
                    break;
                } else {
                    lines.resize(allocLen, 0);
                }
            }
            strcpy_safe(&mut lines[totalLen..], &line);
            totalLen += llen;

            if eoline {
                strcpy_safe(&mut lines[totalLen..], cs!(" ")); /* add a space between lines */
                totalLen += 1;
            }
        }
        fclose(aFile);

        // HEAP ALLOCATION
        let (v, _, _) = lines.into_raw_parts();
        *contents = v;
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// parse off the next token, delimited by a character in 'delimiter',
/// from the input ptr string;  increment *ptr to the end of the token.
/// Returns the length of the token, not including the delimiter char;
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_get_token(
    ptr: *mut *mut c_char,
    delimiter: *const c_char,
    token: *mut c_char,
    isanumber: *mut c_int, /* O - is this token a number? */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// parse off the next token, delimited by a character in 'delimiter',
/// from the input ptr string;  increment *ptr to the end of the token.
/// Returns the length of the token, not including the delimiter char;
///
/// This routine allocates the *token string;  the calling routine must free it
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_get_token2(
    ptr: *mut *mut c_char,
    delimiter: *const c_char,
    token: *mut *mut c_char,
    isanumber: *mut c_int, /* O - is this token a number? */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let isanumber = isanumber.as_mut();
        let status = status.as_mut().expect(NULL_MSG);
        raw_to_slice!(delimiter);

        let mut ptr_index = 0;
        let mut token_vec = None;

        let p: &[c_char] = cast_slice(CStr::from_ptr(*ptr).to_bytes_with_nul());

        let len = fits_get_token2_safe(
            p,
            &mut ptr_index,
            delimiter,
            &mut token_vec,
            isanumber,
            status,
        );

        //HEAP ALLOCATION
        if let Some(mut token_vec) = token_vec {
            *token = token_vec.as_mut_ptr();
            mem::forget(token_vec);
        } else {
            *token = ptr::null_mut();
        }

        *ptr = *ptr.add(ptr_index);

        len
    }
}

/*---------------------------------------------------------------------------*/
/// A sequence of calls to fits_split_names will split the input string
/// into name tokens.  The string typically contains a list of file or
/// column names.  The names must be delimited by a comma and/or spaces.
/// This routine ignores spaces and commas that occur within parentheses,
/// brackets, or curly brackets.  It also strips any leading and trailing
/// blanks from the returned name.
///
/// This routine is similar to the ANSI C 'strtok' function:
///
/// The first call to fits_split_names has a non-null input string.
/// It finds the first name in the string and terminates it by
/// overwriting the next character of the string with a '\0' and returns
/// a pointer to the name.  Each subsequent call, indicated by a NULL
/// value of the input string, returns the next name, searching from
/// just past the end of the previous name.  It returns NULL when no
/// further names are found.
///
/// The following line illustrates how a string would be split into 3 names:
///  myfile[1][bin (x,y)=4], file2.fits  file3.fits
///  ^^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^^  ^^^^^^^^^^
///    1st name               2nd name    3rd name
///
///
/// NOTE:  This routine is not thread-safe.  
/// This routine is simply provided as a utility routine for other external
/// software. It is not used by any CFITSIO routine.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_split_names(
    list: *const c_char, /* I   - input list of names */
) -> *mut c_char {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// compare input URL with list of known drivers, returning the
/// matching driver numberL.
pub(crate) fn urltype2driver(urltype: &[c_char], driver: &mut c_int) -> c_int {
    let d = DRIVER_TABLE.get().unwrap();
    let mut ii = (d.len() - 1) as isize;

    /* find matching driver; search most recent drivers first */
    while ii >= 0 {
        if 0 == strcmp_safe(&d[ii as usize].prefix, urltype) {
            *driver = ii as c_int;
            return 0;
        }
        ii -= 1;
    }

    NO_MATCHING_DRIVER
}

/*--------------------------------------------------------------------------*/
/// close the FITS file by completing the current HDU, flushing it to disk,
/// then calling the system dependent routine to physically close the FITS file
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffclos(
    fptr: Option<Box<fitsfile>>, /* I - FITS file pointer */
    status: *mut c_int,          /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        match fptr {
            None => {
                *status = NULL_INPUT_PTR;
                return *status;
            }
            Some(fptr) => {
                return ffclos_safer(fptr, status);
            }
        };
    }
}

/*--------------------------------------------------------------------------*/
/// close the FITS file by completing the current HDU, flushing it to disk,
/// then calling the system dependent routine to physically close the FITS file
pub unsafe fn ffclos_safer(mut fptr: Box<fitsfile>, status: &mut c_int) -> c_int {
    unsafe {
        let mut tstatus = NO_CLOSE_ERROR;
        let mut zerostatus = 0;

        if fptr.Fptr.validcode != VALIDSTRUC {
            /* BAD_FILEPTR */
            *status = BAD_FILEPTR;
            return *status;
        }

        /* close and flush the current HDU */
        if *status > 0 {
            ffchdu(&mut fptr, &mut tstatus); /* turn off the error message from ffchdu */
        } else {
            ffchdu(&mut fptr, status);
        }

        (fptr.Fptr.open_count) -= 1; /* decrement usage counter */

        if fptr.Fptr.open_count == 0 {
            /* if no other files use structure */

            ffflsh_safe(&mut fptr, true, status); /* flush and disassociate IO buffers */
            //let d = driverTable.lock().unwrap();
            let d = DRIVER_TABLE.get().unwrap();

            /* call driver function to actually close the file */
            if ((d[fptr.Fptr.driver as usize]).close)(fptr.Fptr.filehandle) != 0 && *status <= 0 {
                *status = FILE_NOT_CLOSED; /* report if no previous error */
                ffpmsg_str("failed to close the following file: (ffclos)");
                ffpmsg_cstr(CStr::from_ptr(fptr.Fptr.filename));
            };

            fits_clear_Fptr_safer(&mut fptr.Fptr, status); /* clear Fptr address */

            drop(fptr);
        } else {
            /*
               to minimize the fallout from any previous error (e.g., trying to
               open a non-existent extension in a already opened file),
               always call ffflsh with status = 0.
            */

            /* just flush the buffers, don't disassociate them */
            if *status > 0 {
                ffflsh_safe(&mut fptr, false, &mut zerostatus);
            } else {
                ffflsh_safe(&mut fptr, false, status);
            }

            // WARNING: In the C version of this code, the inner FITSfile can be shared many times
            // and as such, given the open_count > 0, the inner FITSfile is not dropped.
            // To replicate this, we need to deconstruct the fitsfile struct and
            // drop the inner FITSfile struct, but not the outer one.

            let fitsfile { HDUposition, Fptr } = *fptr;
            let _ = Box::into_raw(Fptr); // WARNING: Dangling pointer
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// close and DELETE the FITS file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdelt(
    mut fptr: *mut fitsfile, /* I - FITS file pointer */
    status: *mut c_int,      /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if fptr.is_null() {
            *status = NULL_INPUT_PTR;
            return *status;
        }

        let mut boxed_fptr = Some(Box::from_raw(fptr));

        let result = ffdelt_safer(&mut (boxed_fptr), status);

        if result != 0 {
            fptr = ptr::null_mut(); /* set to null to  avoid dangling pointer */
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// close and DELETE the FITS file.
pub(crate) unsafe fn ffdelt_safer(
    fptr: &mut Option<Box<fitsfile>>, /* I - FITS file pointer */
    status: &mut c_int,               /* IO - error status     */
) -> c_int {
    unsafe {
        let mut basename: Vec<c_char> = Vec::new();
        let mut slen: c_int = 0;
        let mut tstatus = NO_CLOSE_ERROR;
        let mut zerostatus: c_int = 0;

        let local_fptr = fptr.as_deref_mut().unwrap();

        if local_fptr.Fptr.validcode != VALIDSTRUC {
            /* check for magic value */
            *status = BAD_FILEPTR;
            return *status;
        }

        if *status > 0 {
            ffchdu(local_fptr, &mut tstatus); /* turn off the error message from ffchdu */
        } else {
            ffchdu(local_fptr, status);
        }

        ffflsh_safe(local_fptr, true, status); /* flush and disassociate IO buffers */

        /* call driver function to actually close the file */
        //let d = driverTable.lock().unwrap();
        let d = DRIVER_TABLE.get().unwrap();
        if (d[local_fptr.Fptr.driver as usize].close)(local_fptr.Fptr.filehandle) != 0
            && *status <= 0
        {
            *status = FILE_NOT_CLOSED; /* report error if no previous error */

            ffpmsg_str("failed to close the following file: (ffdelt)");
            ffpmsg_cstr(CStr::from_ptr(local_fptr.Fptr.filename));
        }

        /* call driver function to actually delete the file */
        if let Some(remove) = d[local_fptr.Fptr.driver as usize].remove {
            /* parse the input URL to get the base filename */
            slen = strlen(local_fptr.Fptr.filename) as c_int;
            if basename.try_reserve_exact((slen + 1) as usize).is_err() {
                *status = MEMORY_ALLOCATION;
                return *status;
            } else {
                basename.resize((slen + 1) as usize, 0)
            }

            let url = CStr::from_ptr(local_fptr.Fptr.filename).to_bytes_with_nul();

            ffiurl_safer(
                cast_slice(url),
                ptr::null_mut(),
                basename.as_mut_ptr(),
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
                &mut zerostatus,
            );

            if remove(&basename) != 0 {
                ffpmsg_str("failed to delete the following file: (ffdelt)");
                ffpmsg_cstr(CStr::from_ptr(local_fptr.Fptr.filename));
                if (*status) == 0 {
                    *status = FILE_NOT_CLOSED;
                }
            }
        }

        fits_clear_Fptr_safer(&mut local_fptr.Fptr, status); /* clear Fptr address */
        local_fptr.Fptr.filename = ptr::null_mut();
        local_fptr.Fptr.validcode = 0; /* magic value to indicate invalid fptr */

        *fptr = None;

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// low level routine to truncate a file to a new smaller size.
pub(crate) fn fftrun(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    filesize: LONGLONG,  /* I - size to truncate the file   */
    status: &mut c_int,  /* O - error status                */
) -> c_int {
    //let d = driverTable.lock().unwrap();
    let d = DRIVER_TABLE.get().unwrap();
    let truncate = d[fptr.Fptr.driver as usize].truncate;

    if let Some(truncate) = truncate {
        ffflsh_safe(fptr, false, status); /* flush all the buffers first */
        fptr.Fptr.filesize = filesize;
        fptr.Fptr.io_pos = filesize;
        fptr.Fptr.logfilesize = filesize;
        fptr.Fptr.bytepos = filesize;
        ffbfeof(fptr, status); /* eliminate any buffers beyond current EOF */
        *status = truncate(fptr.Fptr.filehandle, filesize as usize);
        *status
    } else {
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// low level routine to flush internal file buffers to the file.
pub(crate) fn ffflushx(fptr: &mut FITSfile, /* I - FITS file pointer                  */) -> c_int {
    let d = DRIVER_TABLE.get().unwrap();
    if let Some(flush) = (d[fptr.driver as usize]).flush {
        flush(fptr.filehandle)
    } else {
        0 /* no flush function defined for this driver */
    }
}

/*--------------------------------------------------------------------------*/
/// low level routine to seek to a position in a file.
pub(crate) fn ffseek(
    fptr: &mut FITSfile, /* I - FITS file pointer              */
    position: LONGLONG,  /* I - byte position to seek to       */
) -> c_int {
    let d = DRIVER_TABLE.get().unwrap();
    // SAFETY: This is just a hack so I can make other functions 'safe'
    // TODO: Fix
    (d[fptr.driver as usize].seek)(fptr.filehandle, position)
}

/*--------------------------------------------------------------------------*/
/// low level routine to write bytes to a file.
pub(crate) fn ffwrite(
    fptr: &mut FITSfile, /* I - FITS file pointer              */
    nbytes: c_long,      /* I - number of bytes to write       */
    buffer: &[u8],       /* I - buffer to write                */
    status: &mut c_int,  /* O - error status                   */
) -> c_int {
    {
        let d = DRIVER_TABLE.get().unwrap();
        // SAFETY: This is just a hack so I can make other functions 'safe'
        // TODO: Fix
        unsafe {
            if (d[fptr.driver as usize].write)(fptr.filehandle, buffer, nbytes as usize) > 0 {
                ffpmsg_str("Error writing data buffer to file:");
                ffpmsg_cstr(CStr::from_ptr(fptr.filename));

                *status = WRITE_ERROR;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// low level routine to write bytes to a file.
pub(crate) fn ffwrite_int(
    fptr: &mut FITSfile, /* I - FITS file pointer              */
    nbytes: usize,       /* I - number of bytes to write       */
    nbuff: usize,        /* I - buffer offset to write to      */
    status: &mut c_int,  /* O - error status                   */
) -> c_int {
    {
        let d = DRIVER_TABLE.get().unwrap();
        // SAFETY: This is just a hack so I can make other functions 'safe'
        // TODO: Fix
        let buffer = cast_slice(&fptr.iobuffer[(nbuff * IOBUFLEN as usize)..]);
        unsafe {
            if (d[fptr.driver as usize].write)(fptr.filehandle, buffer, nbytes) > 0 {
                ffpmsg_str("Error writing data buffer to file:");
                ffpmsg_cstr(CStr::from_ptr(fptr.filename));

                *status = WRITE_ERROR;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// low level routine to read bytes from a file.
pub(crate) fn ffread(
    fptr: &FITSfile,    /* I - FITS file pointer              */
    nbytes: c_long,     /* I - number of bytes to read        */
    buffer: &mut [u8],  /* O - buffer to read into            */
    status: &mut c_int, /* O - error status                   */
) -> c_int {
    let d = DRIVER_TABLE.get().unwrap();
    let readstatus = (d[fptr.driver as usize].read)(fptr.filehandle, buffer, nbytes as usize);

    if readstatus == END_OF_FILE {
        *status = END_OF_FILE;
    } else if readstatus > 0 {
        ffpmsg_str("Error reading data buffer from file:");

        // SAFETY: This is just a hack so I can make other functions 'safe'
        // TODO: Fix
        unsafe {
            ffpmsg_cstr(CStr::from_ptr(fptr.filename));
        }

        *status = READ_ERROR;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// low level routine to read bytes from a file.
pub(crate) fn ffread_int(
    fptr: &mut FITSfile, /* I - FITS file pointer              */
    nbytes: usize,       /* I - number of bytes to read        */
    nbuff: usize,        /* I - buffer offset to read into            */
    status: &mut c_int,  /* O - error status                   */
) -> c_int {
    let d = DRIVER_TABLE.get().unwrap();
    let buffer = cast_slice_mut(&mut fptr.iobuffer[(nbuff * IOBUFLEN as usize)..]);
    let readstatus = (d[fptr.driver as usize].read)(fptr.filehandle, buffer, nbytes);

    if readstatus == END_OF_FILE {
        *status = END_OF_FILE;
    } else if readstatus > 0 {
        ffpmsg_str("Error reading data buffer from file:");

        // SAFETY: This is just a hack so I can make other functions 'safe'
        // TODO: Fix
        unsafe {
            ffpmsg_cstr(CStr::from_ptr(fptr.filename));
        }

        *status = READ_ERROR;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Create and initialize a new FITS file  based on a template file.
/// Uses C fopen and fgets functions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fftplt(
    fptr: *mut Option<Box<fitsfile>>, /* O - FITS file pointer                   */
    filename: *const c_char,          /* I - name of file to create              */
    tempname: *const c_char,          /* I - name of template file               */
    status: *mut c_int,               /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(filename);
        raw_to_slice!(tempname);

        /* initialize null file pointer */
        let f_tmp = fptr.take();
        if let Some(f) = f_tmp {
            // WARNING: The c version doesn't null pointers after a close, so we have a dangling pointer.
            // We need to be careful with this, as it can cause double free errors.
            // Therefore, if this function is called with a Some(), then we will leak the pointer because
            // it's probably invalid.
            let _ = Box::into_raw(f);
        }

        /* regardless of the value of *status */
        if *status > 0 {
            return *status;
        }

        if ffinit_safer(fptr, filename, status) != 0 {
            /* create empty file */
            return *status;
        }

        let f = (*fptr).as_mut().expect(NULL_MSG);

        ffoptplt(f, tempname, status); /* open and use template */

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// open template file and use it to create new file
pub(crate) unsafe fn ffoptplt(
    fptr: &mut fitsfile, /* O - FITS file pointer                   */
    tempname: &[c_char], /* I - name of template file               */
    status: &mut c_int,  /* IO - error status                       */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Print out report of cfitsio error status and messages on the error stack.
/// Uses C FILE stream.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffrprt(stream: *mut FILE, status: c_int) {
    unsafe {
        let mut status_str: [c_char; FLEN_STATUS] = [0; FLEN_STATUS];
        let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

        if status > 0 {
            ffgerr_safe(status, &mut status_str); /* get the error description */
            fprintf(
                stream,
                cstr!("\nFITSIO status = %d: %s\n").as_ptr(),
                status,
                status_str.as_ptr(),
            );

            while ffgmsg_safe(&mut errmsg) > 0 {
                /* get error stack messages */
                fprintf(stream, cstr!("%s\n").as_ptr(), errmsg.as_ptr());
            }
        }
    }
}

/*--------------------------------------------------------------------------*/
pub fn pixel_filter_helper(
    fptr: &mut Option<Box<fitsfile>>, /* IO - pointer to input image; on output it  */
    /*      points to the new image */
    outfile: [c_char; FLEN_FILENAME], /* I - name for output file        */
    pixfilter: [c_char; FLEN_FILENAME], /* I - Image filter expression    */
    status: *mut c_int,
) -> c_int {
    todo!()
}

/*-------------------------------------------------------------------*/
/// Wrapper function for global initialization of curl library.
/// This is NOT THREAD-SAFE
pub(crate) fn ffihtps() {
    todo!();
}

/*-------------------------------------------------------------------*/
/// Wrapper function for global cleanup of curl library.
/// This is NOT THREAD-SAFE
pub(crate) fn ffchtps() {
    todo!();
}

/*-------------------------------------------------------------------*/
/// Turn libcurl's verbose output on (1) or off (0).
/// This is NOT THREAD-SAFE
pub(crate) fn ffvhtps(flag: c_int) {
    todo!();
}

/*-------------------------------------------------------------------*/
/// Display download status bar (to stderr), where applicable.
/// This is NOT THREAD-SAFE
pub(crate) fn ffshdwn(flag: c_int) {
    todo!();
}

/*-------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtmo() -> c_int {
    todo!();
}

/*-------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffstmo(sec: c_int, status: *mut c_int) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// parse off the next token, delimited by a character in 'delimiter',
/// from the input ptr string;  increment *ptr to the end of the token.
/// Returns the length of the token, not including the delimiter char;
///
/// This routine allocates the *token string;  the calling routine must free it
pub(crate) fn fits_get_token2_safe(
    ptr: &[c_char],
    ptr_index: &mut usize,
    delimiter: &[c_char],
    token: &mut Option<Vec<c_char>>,
    isanumber: Option<&mut c_int>, /* O - is this token a number? */
    status: &mut c_int,
) -> c_int {
    let mut tval: [c_char; 73] = [0; 73];
    let mut slen = 0;
    let mut dval: f64 = 0.0;
    let loc = 0;

    *ptr_index = 0; // Ensure starts at beginning of ptr

    if *status != 0 {
        return 0;
    }

    while ptr[*ptr_index] == bb(b' ') {
        /* skip over leading blanks */
        *ptr_index += 1;
    }

    let p_tmp = &ptr[*ptr_index..];
    slen = strcspn_safe(p_tmp, delimiter); /* length of next token */
    if slen != 0 {
        // HEAP ALLOCATION
        let mut t: Vec<c_char> = Vec::new();
        if t.try_reserve_exact(slen + 1).is_err() {
            ffpmsg_str("Couldn't allocate memory to hold token string (fits_get_token2).");
            *status = MEMORY_ALLOCATION;
            return 0;
        } else {
            t.resize(slen + 1, 0);
        }

        let p_tmp = &ptr[*ptr_index..];
        strncat_safe(&mut t, p_tmp, slen); /* copy token */
        *ptr_index += 1; /* skip over the token */

        if let Some(isanumber) = isanumber {
            /* check if token is a number */

            *isanumber = 1;

            let mut loc = 0;
            let r;
            if (strchr_safe(&t, bb(b'D'))).is_some() {
                strncpy_safe(&mut tval, &t, 72);
                tval[72] = 0;

                /*  The C language does not support a 'D'; replace with 'E' */
                let tmp_loc: Option<usize> = strchr_safe(&tval, bb(b'D'));
                if let Some(tmp_loc) = tmp_loc {
                    tval[tmp_loc] = bb(b'E');
                }

                r = &tval[..];
                dval = strtod_safe(r, &mut loc);
            } else {
                r = t.as_slice();
                dval = strtod_safe(r, &mut loc);
            }

            /* check for read error, or junk following the value */
            if r[loc] != 0 && r[loc] != bb(b' ') {
                *isanumber = 0;
            }

            if errno().0 == ERANGE {
                *isanumber = 0;
            }
        }

        *token = Some(t);
    }

    slen as c_int
}
