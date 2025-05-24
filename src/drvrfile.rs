/*  This file, drvrfile.c contains driver routines for disk files.         */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use std::env;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Mutex;

use crate::c_types::{c_char, c_int};

use bytemuck::cast_slice;
use cstr::cstr;

#[cfg(target_family = "unix")]
use pwd::Passwd;

use crate::cfileio::MAX_PREFIX_LEN;
use crate::fitscore::{ffpmsg_slice, ffpmsg_str};

use crate::fitsio::*;
use crate::fitsio2::*;
use crate::group::{fits_get_cwd, fits_relurl2url};
use crate::wrappers::*;
use crate::{bb, cs};

const IO_SEEK: c_int = 0; /* last file I/O operation was a seek */
const IO_READ: c_int = 1; /* last file I/O operation was a read */
const IO_WRITE: c_int = 2; /* last file I/O operation was a write */

static FILE_OUTFILE: Mutex<[c_char; FLEN_FILENAME]> = Mutex::new([0; FLEN_FILENAME]);

#[derive(Default)]
pub(crate) struct diskdriver {
    /* structure containing disk file structure */
    fileptr: Option<File>,
    currentpos: usize,
    last_io_op: c_int,
}

//pub static mut handleTable: Lazy<Mutex<[diskdriver; NMAXFILES]>> = Lazy::new(|| Mutex::new([diskdriver::default(); NMAXFILES]));
/* allocate diskfile handle tables */
pub(crate) static HANDLE_TABLE: Mutex<Vec<diskdriver>> = Mutex::new(Vec::new());

/*--------------------------------------------------------------------------*/
pub(crate) fn file_init() -> c_int {
    let mut h = HANDLE_TABLE.lock().unwrap();

    h.clear();
    h.resize_with(NMAXFILES, ||
         /* initialize all empty slots in table */
            diskdriver {
                fileptr: None,
                currentpos: 0,
                last_io_op: 0,
            });

    0
}

/*--------------------------------------------------------------------------*/
pub(crate) fn file_setoptions(options: c_int) -> c_int {
    /* do something with the options argument, to stop compiler warning */
    0
}

/*--------------------------------------------------------------------------*/
pub(crate) fn file_getoptions(options: &mut c_int) -> c_int {
    *options = 0;
    0
}

/*--------------------------------------------------------------------------*/
pub(crate) fn file_getversion(version: &mut c_int) -> c_int {
    *version = 10;
    0
}

/*--------------------------------------------------------------------------*/
pub(crate) fn file_shutdown() -> c_int {
    0
}

/*--------------------------------------------------------------------------*/
pub(crate) fn file_open(filename: &mut [c_char], rwmode: c_int, handle: &mut c_int) -> c_int {
    let mut copyhandle = 0;
    let mut ii = 0;
    let mut status = 0;
    let mut recbuf: [u8; IOBUFLEN as usize] = [0; IOBUFLEN as usize];
    let nread: usize = 0;
    let mut diskfile: Option<File> = None;
    /*
       if an output filename has been specified as part of the input
       file, as in "inputfile.fits(outputfile.fit)" then we have to
       create the output file, copy the input to it, then reopen the
       the new copy.
    */

    let mut ofile = FILE_OUTFILE.lock().unwrap();

    if ofile[0] != 0 {
        /* open the original file, with readonly access */
        status = file_openfile(filename, READONLY, &mut diskfile);
        if status != 0 {
            ofile[0] = 0;
            return status;
        }

        /* create the output file */
        status = file_create(&mut ofile, handle);

        if status != 0 {
            ffpmsg_str("Unable to create output file for copy of input file:");
            ffpmsg_slice(&*ofile);
            ofile[0] = 0;
            return status;
        }

        /* copy the file from input to output */
        let mut nread = diskfile.as_mut().unwrap().read(&mut recbuf).unwrap();
        while 0 != nread {
            status = file_write(*handle, &recbuf, nread);
            if status != 0 {
                ofile[0] = 0;
                return status;
            };
            nread = diskfile.as_mut().unwrap().read(&mut recbuf).unwrap();
        }

        /* close both files */
        copyhandle = *handle;
        file_close(*handle);
        *handle = copyhandle; /* reuse the old file handle */

        /* reopen the new copy, with correct rwmode */
        status = file_openfile(&*ofile, rwmode, &mut diskfile);
        ofile[0] = 0;
    } else {
        //let h = handleTable.lock().unwrap();
        let h = HANDLE_TABLE.lock().unwrap();
        *handle = -1;
        ii = 0;
        while ii < NMAXFILES {
            /* find empty slot in table */
            if (h[ii]).fileptr.is_none() {
                *handle = ii as c_int;
                break;
            };
            ii += 1
        }

        if *handle == -1 {
            return TOO_MANY_FILES; /* too many files opened */
        }

        /*open the file */
        status = file_openfile(filename, rwmode, &mut diskfile);
    }

    //let mut h = handleTable.lock().unwrap();
    let mut h = HANDLE_TABLE.lock().unwrap();
    h[*handle as usize].fileptr = diskfile;
    h[*handle as usize].currentpos = 0;
    h[*handle as usize].last_io_op = IO_SEEK;
    status
}

/*--------------------------------------------------------------------------*/
/// lowest level routine to physically open a disk file
pub(crate) fn file_openfile(
    filename: &[c_char],
    rwmode: c_int,
    diskfile: &mut Option<File>,
) -> c_int {
    let mode: [c_char; 4] = [0; 4];

    let file: Result<File, std::io::Error>;

    if CFITSIO_MACHINE == ALPHAVMS || CFITSIO_MACHINE == VAXVMS {
        todo!();
    }

    #[cfg(target_family = "windows")]
    {
        file = File::options().read(true).write(rwmode == READWRITE).open(
            CStr::from_bytes_until_nul(cast_slice(filename))
                .unwrap()
                .to_str()
                .unwrap(),
        );
    }

    #[cfg(target_family = "unix")]
    {
        let mut tempname: [c_char; 1024] = [0; 1024];
        let mut user: [c_char; 80] = [0; 80];
        let mut ii = 0;

        /* support the ~user/file.fits or ~/file.fits filenames in UNIX */

        if filename[0] == bb(b'~') {
            if filename[1] == bb(b'/') {
                let home_dir = std::env::var("HOME");

                match home_dir {
                    Ok(home_dir) => {
                        let h = CString::new(home_dir).unwrap();
                        let home_dir = cast_slice(h.as_bytes_with_nul());
                        if strlen_safe(home_dir) + strlen_safe(&filename[1..]) > 1023 {
                            return FILE_NOT_OPENED;
                        }
                        strcpy_safe(&mut tempname, home_dir);
                        strcat_safe(&mut tempname, &filename[1..]);
                    }
                    Err(_) => {
                        if strlen_safe(filename) > 1023 {
                            return FILE_NOT_OPENED;
                        }
                        strcpy_safe(&mut tempname, filename);
                    }
                }
            } else {
                /* copy user name */
                let mut ci = 1;
                while filename[ci] != 0 && (filename[ci] != bb(b'/')) {
                    user[ii] = filename[ci];
                    ci += 1;
                    ii += 1;
                }

                user[ii] = 0;

                /* get structure that includes name of user's home directory */
                let pwd = Passwd::from_name(
                    CStr::from_bytes_with_nul(cast_slice(&user))
                        .unwrap()
                        .to_str()
                        .unwrap(),
                )
                .unwrap()
                .unwrap();
                let pw_dir_cstr = CString::new(pwd.dir).unwrap();
                let pw_dir = cast_slice(pw_dir_cstr.to_bytes_with_nul());

                /* copy user's home directory */
                if strlen_safe(pw_dir) + strlen_safe(&filename[ci..]) > 1023 {
                    return FILE_NOT_OPENED;
                }

                strcpy_safe(&mut tempname, pw_dir);
                strcat_safe(&mut tempname, &filename[ci..]);
            }
            file = File::options().read(true).write(rwmode == READWRITE).open(
                CStr::from_bytes_until_nul(cast_slice(&tempname))
                    .unwrap()
                    .to_str()
                    .unwrap(),
            );
        } else {
            /* don't need to expand the input file name */
            file = File::options().read(true).write(rwmode == READWRITE).open(
                CStr::from_bytes_until_nul(cast_slice(filename))
                    .unwrap()
                    .to_str()
                    .unwrap(),
            );
        }
    }

    if (file).is_err() {
        /* couldn't open file */
        return FILE_NOT_OPENED;
    } else {
        *diskfile = Some(file.unwrap());
    }

    0
}

/*--------------------------------------------------------------------------*/
pub(crate) fn file_create(filename: &mut [c_char; FLEN_FILENAME], handle: &mut c_int) -> c_int {
    let mut mode: [c_char; 4] = [0; 4]; /* make sure the CWD ends with slash */
    let mut status = 0;
    let mut rootlen = 0;
    let mut rootlen2 = 0;
    let mut slen = 0;
    let mut cwd: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
    let absURL: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
    let mut rootstring: [c_char; 256] = [0; 256];
    let mut rootstring2: [c_char; 256] = [0; 256];
    let mut username: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
    let mut userroot: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
    let mut userroot2: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];

    let cptr = env::var("HERA_DATA_DIRECTORY");

    if let Ok(cptr) = cptr {
        let cptr_str = CString::new(cptr).unwrap();
        let hera_dir = cast_slice(cptr_str.as_bytes_with_nul());

        /* This environment variable is defined in the Hera data analysis environment. */
        /* It specifies the root directory path to the users data directories.  */
        /* CFITSIO will verify that the path to the file that is to be created */
        /* is within this root directory + the user's home directory name. */

        if strlen_safe(hera_dir) > 200 {
            /* guard against possible string overflows */
            return FILE_NOT_CREATED;
        }

        /* environment variable has the form "path/one/;/path/two/" where the */
        /* second path is optional */

        strcpy_safe(&mut rootstring, hera_dir);
        let cpos = strchr_safe(&rootstring, bb(b';'));
        if let Some(mut cpos) = cpos {
            rootstring[cpos] = 0;
            cpos += 1;
            strcpy_safe(&mut rootstring2, &rootstring[cpos..]);
        } else {
            rootstring2[0] = 0;
        }

        /* Get the current working directory */
        fits_get_cwd(&mut cwd, &mut status);
        slen = strlen_safe(&cwd);
        if (slen < FLEN_FILENAME) && cwd[slen - 1] != bb(b'/') {
            /* make sure the CWD ends with slash */
            strcat_safe(&mut cwd, cs!(b"/"));
        }

        /* check that CWD string matches the rootstring */
        rootlen = strlen_safe(&rootstring);
        if strncmp_safe(&rootstring, &cwd, rootlen) != 0 {
            ffpmsg_str("invalid CWD: does not match root data directory");
            return FILE_NOT_CREATED;
        } else {
            /* get the user name from CWD (it follows the root string) */
            strncpy_safe(&mut username, &cwd[rootlen..], 50); /* limit length of user name */
            username[50] = 0;

            let cpos = strchr_safe(&username, bb(b'/'));
            match cpos {
                None => {
                    ffpmsg_str("invalid CWD: not equal to root data directory + username");
                    return FILE_NOT_CREATED;
                }
                Some(cpos) => {
                    username[cpos + 1] = 0; /* truncate user name string */

                    /* construct full user root name */
                    strcpy_safe(&mut userroot, &rootstring);
                    strcat_safe(&mut userroot, &username);
                    rootlen = strlen_safe(&userroot);

                    /* construct alternate full user root name */
                    strcpy_safe(&mut userroot2, &rootstring2);
                    strcat_safe(&mut userroot2, &username);
                    rootlen2 = strlen_safe(&userroot2);

                    /* convert the input filename to absolute path relative to the CWD */
                    fits_relurl2url(cwd, filename.as_mut_ptr(), absURL, &mut status);

                    /*
                    printf("username = %s\n", username);
                    printf("userroot = %s\n", userroot);
                    printf("userroot2 = %s\n", userroot2);
                    printf("filename = %s\n", filename);
                    printf("ABS = %s\n", absURL);
                    */
                    /* check that CWD string matches the rootstring or alternate root string */
                    if strncmp_safe(&userroot, &absURL, rootlen) != 0
                        && strncmp_safe(&userroot2, &absURL, rootlen2) != 0
                    {
                        ffpmsg_str("invalid filename: path not within user directory");
                        return FILE_NOT_CREATED;
                    }
                }
            }
            /* if we got here, then the input filename appears to be valid */
        }
    }

    *handle = -1;
    let mut ii = 0;
    {
        //let h = handleTable.lock().unwrap();
        let h = HANDLE_TABLE.lock().unwrap();
        while ii < NMAXFILES {
            /* find empty slot in table */
            if (h[ii]).fileptr.is_none() {
                *handle = ii as c_int;
                break;
            }
            ii += 1;
        }
    }

    if *handle == -1 {
        return TOO_MANY_FILES; /* too many files opened */
    }

    strcpy_safe(&mut mode, cs!(b"w+b")); /* create new file with read-write */

    /* does file already exist? */
    let diskfile = File::options()
        .read(true)
        .write(false)
        .truncate(false)
        .open(
            CStr::from_bytes_until_nul(cast_slice(filename))
                .unwrap()
                .to_str()
                .unwrap(),
        );

    if diskfile.is_ok() {
        /* close file and exit with error */
        return FILE_NOT_CREATED;
    }

    if CFITSIO_MACHINE == ALPHAVMS || CFITSIO_MACHINE == VAXVMS {
        todo!();
    }

    let diskfile = File::options()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(
            CStr::from_bytes_until_nul(cast_slice(filename))
                .unwrap()
                .to_str()
                .unwrap(),
        );

    match diskfile {
        Err(x) => return FILE_NOT_CREATED, /* couldn't create file */
        Ok(x) => {
            //let mut h = handleTable.lock().unwrap();

            let mut h = HANDLE_TABLE.lock().unwrap();
            (h[ii]).fileptr = Some(x);
            (h[ii]).currentpos = 0;
            (h[ii]).last_io_op = IO_SEEK;
        }
    }

    0
}

/*--------------------------------------------------------------------------*/
/// truncate the diskfile to a new smaller size
pub(crate) fn file_truncate(handle: c_int, filesize: usize) -> c_int {
    if HAVE_FTRUNCATE {
        let handle = handle as usize;
        let mut h = HANDLE_TABLE.lock().unwrap();

        if let Some(f) = &mut h[handle].fileptr {
            f.set_len(filesize as u64).unwrap();
            f.seek(SeekFrom::Start(filesize as u64)).unwrap();

            h[handle].currentpos = filesize;
            h[handle].last_io_op = IO_SEEK;
        }
    }

    0
}

/*--------------------------------------------------------------------------*/
/// return the size of the file in bytes
pub(crate) fn file_size(handle: c_int, filesize: &mut usize) -> c_int {
    //let h = handleTable.lock().unwrap();
    let mut h = HANDLE_TABLE.lock().unwrap();

    if let Some(diskfile) = &mut (h[handle as usize].fileptr) {
        let position1 = diskfile.stream_position().unwrap(); /* save current postion */
        let position2 = diskfile.seek(SeekFrom::End(0)).unwrap(); /* get file size */
        diskfile.seek(SeekFrom::Start(position1)).unwrap(); /* seek back to original pos */

        *filesize = position2 as usize;
    }

    0
}

/*--------------------------------------------------------------------------*/
/// close the file
pub(crate) fn file_close(handle: c_int) -> c_int {
    //let mut h = handleTable.lock().unwrap();
    let mut h = HANDLE_TABLE.lock().unwrap();

    if let Some(f) = &mut h[handle as usize].fileptr {
        if f.sync_all().is_err() {
            return FILE_NOT_CLOSED;
        }
    }

    (h[handle as usize]).fileptr = None; // Implicitly drop the file

    0
}

/*--------------------------------------------------------------------------*/
/// delete the file from disk
pub(crate) fn file_remove(filename: &[c_char]) -> c_int {
    let filename = CStr::from_bytes_until_nul(cast_slice(filename)).unwrap();
    let filename = filename.to_str().unwrap();

    let _ = std::fs::remove_file(filename);

    0
}

/*--------------------------------------------------------------------------*/
/// flush the file
pub(crate) fn file_flush(handle: c_int) -> c_int {
    let mut h = HANDLE_TABLE.lock().unwrap();

    if let Some(f) = &mut h[handle as usize].fileptr {
        if f.flush().is_err() {
            return WRITE_ERROR;
        }
    }

    /* The flush operation is not supposed to move the internal */
    /* file pointer, but it does on some Windows-95 compilers and */
    /* perhaps others, so seek to original position to be sure. */
    /* This seek will do no harm on other systems.   */

    if CFITSIO_MACHINE == IBMPC && file_seek(handle, h[handle as usize].currentpos as _) != 0 {
        return SEEK_ERROR;
    }

    0
}

/*--------------------------------------------------------------------------*/
/// seek to position relative to start of the file
pub(crate) fn file_seek(handle: c_int, offset: LONGLONG) -> c_int {
    //let mut h = handleTable.lock().unwrap();
    let mut h = HANDLE_TABLE.lock().unwrap();

    file_seek_internal(&mut h[handle as usize], offset as u64)
}

/*--------------------------------------------------------------------------*/
/// seek to position relative to start of the file
fn file_seek_internal(handle: &mut diskdriver, offset: u64) -> c_int {
    //let mut h = handleTable.lock().unwrap();

    if let Some(f) = &mut handle.fileptr {
        if f.seek(SeekFrom::Start(offset)).is_err() {
            return SEEK_ERROR;
        } else {
            handle.currentpos = offset as usize;
        }
    }

    0
}

/*--------------------------------------------------------------------------*/
/// read bytes from the current position in the file
pub(crate) fn file_read(hdl: c_int, buffer: &mut [u8], nbytes: usize) -> c_int {
    //let h = handleTable.lock().unwrap();
    let mut h = HANDLE_TABLE.lock().unwrap();
    let handle = &mut h[hdl as usize];

    if handle.last_io_op == IO_WRITE {
        let currentpos = handle.currentpos as u64;
        if file_seek_internal(handle, currentpos) != 0 {
            return SEEK_ERROR;
        };
    }

    let nread = handle
        .fileptr
        .as_mut()
        .unwrap()
        .read(&mut buffer[..nbytes])
        .unwrap();

    if nread == 1 {
        /* some editors will add a single end-of-file character to a file */
        /* Ignore it if the character is a zero, 10, or 32 */

        if buffer[0] == 0 || buffer[0] == 10 || buffer[0] == 32 {
            return END_OF_FILE;
        } else {
            return READ_ERROR;
        };
    } else if nread != nbytes {
        return READ_ERROR;
    }

    (handle).currentpos += nbytes;
    (handle).last_io_op = IO_READ;

    0
}

/*--------------------------------------------------------------------------*/
/// write bytes at the current position in the file
pub(crate) fn file_write(hdl: c_int, buffer: &[u8], nbytes: usize) -> c_int {
    //let h = handleTable.lock().unwrap();
    let mut h = HANDLE_TABLE.lock().unwrap();
    let handle = &mut h[hdl as usize];

    if handle.last_io_op == IO_READ {
        let offset = handle.currentpos as u64;
        if file_seek_internal(handle, offset) != 0 {
            return SEEK_ERROR;
        };
    }

    let stream = handle.fileptr.as_mut().unwrap();
    if stream.write_all(&buffer[..nbytes]).is_err() {
        return WRITE_ERROR;
    }

    (handle).currentpos += nbytes;
    (handle).last_io_op = IO_WRITE;

    0
}

/*--------------------------------------------------------------------------*/
/// This routine opens the compressed diskfile by creating a new uncompressed
/// file then opening it.  The input file name (the name of the compressed
/// file) gets replaced with the name of the uncompressed file, which is
/// initially stored in the global file_outfile string.   file_outfile
/// then gets set to a null string.
pub(crate) fn file_compress_open(filename: &mut [c_char], rwmode: c_int, hdl: &mut c_int) -> c_int {
    let mut indiskfile: Option<File> = None;
    let outdiskfile: Option<File> = None;

    /* open the compressed disk file */
    let status = file_openfile(filename, READONLY, &mut indiskfile);
    if status != 0 {
        ffpmsg_str("failed to open compressed disk file (file_compress_open)");
        ffpmsg_slice(filename);
        return status;
    }

    /* name of the output uncompressed file is stored in the */
    /* global variable called 'file_outfile'.                */
    let mut ofile = FILE_OUTFILE.lock().unwrap();
    let mut cptr = ofile.as_slice();
    if cptr[0] == b'!' as c_char {
        /* clobber any existing file with the same name */
        cptr = &ofile[1..];
        let cptr = CStr::from_bytes_until_nul(cast_slice(cptr)).unwrap();
        let cptr = cptr.to_str().unwrap();
        let _ = std::fs::remove_file(cptr);
    } else {
        let tmp_outdiskfile = File::options().read(true).write(false).open(
            CStr::from_bytes_until_nul(cast_slice(&*ofile))
                .unwrap()
                .to_str()
                .unwrap(),
        ); /* does file already exist? */

        match tmp_outdiskfile {
            Ok(f) => {
                ffpmsg_str("uncompressed file already exists: (file_compress_open)");
                ffpmsg_slice(&*ofile);
                ofile[0] = 0;
                return FILE_NOT_CREATED;
            }
            Err(e) => {}
        }
    }

    let tmp_outdiskfile = File::options()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(
            CStr::from_bytes_until_nul(cast_slice(cptr))
                .unwrap()
                .to_str()
                .unwrap(),
        );
    /* create new file */

    if tmp_outdiskfile.is_err() {
        ffpmsg_str("could not create uncompressed file: (file_compress_open)");
        ffpmsg_slice(&*ofile);
        ofile[0] = 0;
        return FILE_NOT_CREATED;
    }

    /* uncompress file into another file */
    todo!();
    //uncompress2file(filename, indiskfile, outdiskfile, &mut status);

    if status != 0 {
        ffpmsg_str("error in file_compress_open: failed to uncompressed file:");
        ffpmsg_slice(filename);
        ffpmsg_str(" into new output file:");
        ffpmsg_slice(&*ofile);
        ofile[0] = 0;
        return status;
    }
    strcpy_safe(filename, cptr); /* switch the names */
    ofile[0] = 0;

    status = file_open(filename, rwmode, hdl);

    status
}

/*--------------------------------------------------------------------------*/
/// Test if the disk file is compressed.  Returns 1 if compressed, 0 if not.
/// This may modify the filename string by appending a compression suffix.
/// WARNING: Need to hope that enough memory was allocated for extending the filename
/// Logic here is very odd.
pub(crate) fn file_is_compressed(filename: &mut [c_char; FLEN_FILENAME]) -> c_int {
    let mut diskfile = None;
    let mut buffer: [u8; 2] = [0; 2];
    let mut tmpfilename: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];

    /* Open file.  Try various suffix combinations */

    if file_openfile(filename, 0, &mut diskfile) != 0 {
        if strlen_safe(filename) > FLEN_FILENAME - 5 {
            return 0;
        }
        strcpy_safe(&mut tmpfilename, filename);
        strcat_safe(filename, cs!(b".gz"));

        if file_openfile(filename, 0, &mut diskfile) != 0 {
            strcpy_safe(filename, &tmpfilename);
            strcat_safe(filename, cs!(b".Z"));

            if file_openfile(filename, 0, &mut diskfile) != 0 {
                strcpy_safe(filename, &tmpfilename);
                strcat_safe(filename, cs!(b".z")); /* it's often lower case on CDROMs */

                if file_openfile(filename, 0, &mut diskfile) != 0 {
                    strcpy_safe(filename, &tmpfilename);
                    strcat_safe(filename, cs!(b".zip"));

                    if file_openfile(filename, 0, &mut diskfile) != 0 {
                        strcpy_safe(filename, &tmpfilename);
                        strcat_safe(filename, cs!(b"-z")); /* VMS suffix */

                        if file_openfile(filename, 0, &mut diskfile) != 0 {
                            strcpy_safe(filename, &tmpfilename);
                            strcat_safe(filename, cs!(b"-gz")); /* VMS suffix */

                            if file_openfile(filename, 0, &mut diskfile) != 0 {
                                strcpy_safe(filename, &tmpfilename); /* restore original name */
                                return 0; /* file not found */
                            };
                        };
                    };
                };
            };
        };
    }

    /* read the first 2 bytes of the file */
    if diskfile.unwrap().read(&mut buffer).unwrap() != 2 {
        /* read 2 bytes */
        return 0;
    }

    /*
    see if the 2 bytes have the magic values for a compressed file
    1 = this is a compressed file
    0 = not a compressed file
    */
    match buffer[..2] {
        [37, 13] => 1,     /* GZIP  */
        [b'P', b'K'] => 1, /* PKZIP */
        [37, 36] => 1,     /* PACK  */
        [37, 35] => 1,     /* LZW   */
        [b'B', b'Z'] => 1, /* BZip2 */
        [37, 40] => 1,     /* LZH   */
        _ => 0,
    }
}

/*--------------------------------------------------------------------------*/
pub(crate) fn file_checkfile(
    urltype: &mut [c_char; MAX_PREFIX_LEN],
    infile: &mut [c_char; FLEN_FILENAME],
    outfile: &mut [c_char; FLEN_FILENAME],
) -> c_int {
    let mut ofile = FILE_OUTFILE.lock().unwrap();

    /* special case: if file:// driver, check if the file is compressed */
    if file_is_compressed(infile) != 0 {
        /* if output file has been specified, save the name for future use: */
        /* This is the name of the uncompressed file to be created on disk. */
        if strlen_safe(outfile) != 0 {
            if strncmp_safe(outfile, cs!(b"mem:"), 4) == 0 {
                /* uncompress the file in memory, with READ and WRITE access */
                strcpy_safe(urltype, cs!(b"compressmem://")); /* use special driver */
                ofile[0] = 0;
            } else {
                strcpy_safe(urltype, cs!(b"compressfile://")); /* use special driver */

                /* don't copy the "file://" prefix, if present.  */
                if strncmp_safe(outfile, cs!(b"file://"), 7) == 0 {
                    strcpy_safe(&mut *ofile, &outfile[7..]);
                } else {
                    strcpy_safe(&mut *ofile, outfile);
                };
            };
        } else {
            /* uncompress the file in memory */
            strcpy_safe(urltype, cs!(b"compress://")); /* use special driver */
            ofile[0] = 0; /* no output file was specified */
        }
    } else {
        /* an ordinary, uncompressed FITS file on disk */

        /* save the output file name for later use when opening the file. */
        /* In this case, the file to be opened will be opened READONLY,   */
        /* and copied to this newly created output file.  The original file */
        /* will be closed, and the copy will be opened by CFITSIO for     */
        /* subsequent processing (possibly with READWRITE access).        */
        if strlen_safe(outfile) != 0 {
            ofile[0] = 0;
            strncat_safe(&mut *ofile, outfile, FLEN_FILENAME - 1);
        }
    }

    0
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/****  driver routines for stream//: device (stdin or stdout)  ********/

/*--------------------------------------------------------------------------*/
/// read from stdin
pub(crate) fn stream_open(filename: &mut [c_char], rwmode: c_int, handle: &mut c_int) -> c_int {
    *handle = 1; /*  1 = stdin */

    0
}

/*--------------------------------------------------------------------------*/
///  write to stdout
pub(crate) fn stream_create(filename: &mut [c_char; FLEN_FILENAME], handle: &mut c_int) -> c_int {
    *handle = 2; /*  2 = stdout */

    0
}

/*--------------------------------------------------------------------------*/
/// return the size of the file in bytes
pub(crate) fn stream_size(handle: c_int, filesize: &mut usize) -> c_int {
    /* this operation is not supported in a stream; return large value */
    *filesize = LONG_MAX as usize;
    0
}

/*--------------------------------------------------------------------------*/
/// don't have to close stdin or stdout
pub(crate) fn stream_close(handle: c_int) -> c_int {
    0
}

/*--------------------------------------------------------------------------*/
/// flush the file
pub(crate) fn stream_flush(handle: c_int) -> c_int {
    if handle == 2 {
        let _ = std::io::stdout().flush();
    }

    0
}
/*--------------------------------------------------------------------------*/
/// seeking is not allowed in a stream
pub(crate) fn stream_seek(handle: c_int, offset: LONGLONG) -> c_int {
    1
}

/*--------------------------------------------------------------------------*/
/// reading from stdin stream
pub(crate) fn stream_read(hdl: c_int, buffer: &mut [u8], nbytes: usize) -> c_int {
    if hdl != 1 {
        return 1; /* can only read from stdin */
    }

    let mut stdin_hdl = std::io::stdin();

    let nread = stdin_hdl.read(buffer);

    match nread {
        Ok(nread) => {
            if nread != nbytes {
                return END_OF_FILE;
            }
        }
        Err(_) => {
            return READ_ERROR;
        }
    }

    0
}

/*--------------------------------------------------------------------------*/
///  write bytes at the current position in the file
pub(crate) fn stream_write(hdl: c_int, buffer: &[u8], nbytes: usize) -> c_int {
    if hdl != 2 {
        return 1; /* can only write to stdout */
    }

    let mut stdout_hdl = std::io::stdout();
    let nwrite = stdout_hdl.write(buffer);

    match nwrite {
        Ok(nwrite) => {
            if nwrite != nbytes {
                return WRITE_ERROR;
            }
        }
        Err(_) => {
            return WRITE_ERROR;
        }
    }

    0
}
