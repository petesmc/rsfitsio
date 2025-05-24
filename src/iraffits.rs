/*------------------------------------------------------------------------*/
/*                                                                        */
/*  These routines have been modified by William Pence for use by CFITSIO */
/*        The original files were provided by Doug Mink                   */
/*------------------------------------------------------------------------*/

/* File imhfile.c
* August 6, 1998
* By Doug Mink, based on Mike VanHilst's readiraf.c

* Module:      imhfile.c (IRAF .imh image file reading and writing)
* Purpose:     Read and write IRAF image files (and translate headers)
* Subroutine:  irafrhead (filename, lfhead, fitsheader, lihead)
*              Read IRAF image header
* Subroutine:  irafrimage (fitsheader)
*              Read IRAF image pixels (call after irafrhead)
* Subroutine:	same_path (pixname, hdrname)
*		Put filename and header path together
* Subroutine:	iraf2fits (hdrname, irafheader, nbiraf, nbfits)
*		Convert IRAF image header to FITS image header
* Subroutine:  irafgeti4 (irafheader, offset)
*		Get 4-byte integer from arbitrary part of IRAF header
* Subroutine:  irafgetc2 (irafheader, offset)
*		Get character string from arbitrary part of IRAF v.1 header
* Subroutine:  irafgetc (irafheader, offset)
*		Get character string from arbitrary part of IRAF header
* Subroutine:  iraf2str (irafstring, nchar)
* 		Convert 2-byte/char IRAF string to 1-byte/char string
* Subroutine:	irafswap (bitpix,string,nbytes)
*		Swap bytes in string in place, with FITS bits/pixel code
* Subroutine:	irafswap2 (string,nbytes)
*		Swap bytes in string in place
* Subroutine	irafswap4 (string,nbytes)
*		Reverse bytes of Integer*4 or Real*4 vector in place
* Subroutine	irafswap8 (string,nbytes)
*		Reverse bytes of Real*8 vector in place


* Copyright:   2000 Smithsonian Astrophysical Observatory
*              You may do anything you like with this file except remove
*              this copyright.  The Smithsonian Astrophysical Observatory
*              makes no representations about the suitability of this
*              software for any purpose.  It is provided "as is" without
*              express or implied warranty.
*/

use crate::c_types::{c_char, c_int, c_long};
use bytemuck::{cast_slice, cast_slice_mut};
use cstr::cstr;
use std::ffi::CStr;
use std::fs::File;
use std::io::{Read, Seek};
use std::sync::Mutex;
use std::{
    cmp,
    io::{Error, ErrorKind},
    ptr,
};

use crate::fitscore::{ffpmsg_slice, ffpmsg_str};
use crate::fitsio::NULL_MSG;
use crate::int_snprintf;
use crate::{
    BL, bb, cs,
    fitsio::FLEN_ERRMSG,
    raw_to_slice,
    wrappers::{
        atof_safe, atoi_safe, strchr_safe, strcmp_safe, strcpy_safe, strlen_safe, strncat_safe,
        strncmp_safe, strncpy_safe,
    },
};

const FILE_NOT_OPENED: c_int = 104;

/* Parameters from iraf/lib/imhdr.h for IRAF version 1 images */
const SZ_IMPIXFILE: usize = 79; /* name of pixel storage file */
const SZ_IMHDRFILE: usize = 79; /* length of header storage file */
const SZ_IMTITLE: usize = 79; /* image title string */
const LEN_IMHDR: usize = 2052; /* length of std header */

/* Parameters from iraf/lib/imhdr.h for IRAF version 2 images */
const SZ_IM2PIXFILE: usize = 255; /* name of pixel storage file */
const SZ_IM2HDRFILE: usize = 255; /* name of header storage file */
const SZ_IM2TITLE: usize = 383; /* image title string */
const LEN_IM2HDR: usize = 2046; /* length of std header */

/* Offsets into header in bytes for parameters in IRAF version 1 images */
const IM_HDRLEN: usize = 12; /* Length of header in 4-byte ints */
const IM_PIXTYPE: usize = 16; /* Datatype of the pixels */
const IM_NDIM: usize = 20; /* Number of dimensions */
const IM_LEN: usize = 24; /* Length (as stored) */
const IM_PHYSLEN: usize = 52; /* Physical length (as stored) */
const IM_PIXOFF: usize = 88; /* Offset of the pixels */
const IM_CTIME: usize = 108; /* Time of image creation */
const IM_MTIME: usize = 112; /* Time of last modification */
const IM_LIMTIME: usize = 116; /* Time of min,max computation */
const IM_MAX: usize = 120; /* Maximum pixel value */
const IM_MIN: usize = 124; /* Maximum pixel value */
const IM_PIXFILE: usize = 412; /* Name of pixel storage file */
const IM_HDRFILE: usize = 572; /* Name of header storage file */
const IM_TITLE: usize = 732; /* Image name string */

/* Offsets into header in bytes for parameters in IRAF version 2 images */
const IM2_HDRLEN: usize = 6; /* Length of header in 4-byte ints */
const IM2_PIXTYPE: usize = 10; /* Datatype of the pixels */
const IM2_SWAPPED: usize = 14; /* Pixels are byte swapped */
const IM2_NDIM: usize = 18; /* Number of dimensions */
const IM2_LEN: usize = 22; /* Length (as stored) */
const IM2_PHYSLEN: usize = 50; /* Physical length (as stored) */
const IM2_PIXOFF: usize = 86; /* Offset of the pixels */
const IM2_CTIME: usize = 106; /* Time of image creation */
const IM2_MTIME: usize = 110; /* Time of last modification */
const IM2_LIMTIME: usize = 114; /* Time of min,max computation */
const IM2_MAX: usize = 118; /* Maximum pixel value */
const IM2_MIN: usize = 122; /* Maximum pixel value */
const IM2_PIXFILE: usize = 126; /* Name of pixel storage file */
const IM2_HDRFILE: usize = 382; /* Name of header storage file */
const IM2_TITLE: usize = 638; /* Image name string */

/* Codes from iraf/unix/hlib/iraf.h */
const TY_CHAR: c_int = 2;
const TY_SHORT: c_int = 3;
const TY_INT: c_int = 4;
const TY_LONG: c_int = 5;
const TY_REAL: c_int = 6;
const TY_DOUBLE: c_int = 7;
const TY_COMPLEX: c_int = 8;
const TY_POINTER: c_int = 9;
const TY_STRUCT: c_int = 10;
const TY_USHORT: c_int = 11;
const TY_UBYTE: c_int = 12;

const LEN_PIXHDR: usize = 1024;
const MAXINT: c_int = 2147483647; /* Biggest number that can fit in long */

static SWAPHEAD: Mutex<bool> = Mutex::new(false); /* =1 to swap data bytes of IRAF header values */
static SWAPDATA: Mutex<bool> = Mutex::new(false); /* =1 to swap bytes in IRAF data pixels */

/*--------------------------------------------------------------------------*/
/// Delete the iraf .imh header file and the associated .pix data file
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_delete_iraf_file(
    filename: *const c_char, /* name of input file      */
    status: *mut c_int,      /* IO - error status       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        raw_to_slice!(filename);

        fits_delete_iraf_file_safe(filename, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Delete the iraf .imh header file and the associated .pix data file
pub(crate) fn fits_delete_iraf_file_safe(
    filename: &[c_char], /* name of input file      */
    status: &mut c_int,  /* IO - error status       */
) -> c_int {
    let mut lenirafhead: usize = 0;
    let mut pixfilename: [c_char; SZ_IM2PIXFILE + 1] = [0; SZ_IM2PIXFILE + 1];

    /* read IRAF header into dynamically created char array (free it later!) */
    let irafheader = irafrdhead(filename, &mut lenirafhead);

    if irafheader.is_err() {
        *status = FILE_NOT_OPENED;
        return *status;
    }

    let mut irafheader = irafheader.unwrap();

    getirafpixname(filename, &mut irafheader, &mut pixfilename, status);

    /* don't need the IRAF header any more */
    // free(irafheader);

    if *status > 0 {
        return *status;
    }

    let _ = std::fs::remove_file(
        CStr::from_bytes_with_nul(cast_slice(filename))
            .unwrap()
            .to_str()
            .unwrap(),
    );

    let _ = std::fs::remove_file(
        CStr::from_bytes_with_nul(cast_slice(&pixfilename))
            .unwrap()
            .to_str()
            .unwrap(),
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// Driver routine that reads an IRAF image into memory, also converting
/// it into FITS format.
pub(crate) unsafe fn iraf2mem(
    filename: &[c_char],       /* name of input file                 */
    buffptr: *mut *mut c_char, /* O - memory pointer (initially ptr::null_mut())    */
    buffsize: &mut usize,      /* O - size of mem buffer, in bytes        */
    filesize: &mut usize,      /* O - size of FITS file, in bytes         */
    status: &mut c_int,        /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut lenirafhead: usize = 0;

        *buffptr = ptr::null_mut();
        *buffsize = 0;
        *filesize = 0;

        let mut b_ptr = Vec::new();

        /* read IRAF header into dynamically created char array (free it later!) */
        let irafheader = irafrdhead(filename, &mut lenirafhead);

        if irafheader.is_err() {
            *status = FILE_NOT_OPENED;
            return *status;
        }

        let mut irafheader = irafheader.unwrap();

        /* convert IRAF header to FITS header in memory */
        let tmp_buffer = iraftofits(
            filename,
            &mut irafheader,
            lenirafhead,
            &mut b_ptr,
            buffsize,
            filesize,
            status,
        );

        /* don't need the IRAF header any more */
        // free(irafheader);

        if *status > 0 {
            return *status;
        }

        *filesize = (((*filesize - 1) / BL!()) + 1) * BL!(); /* multiple of 2880 */

        /* append the image data onto the FITS header */
        irafrdimage(&mut b_ptr, buffsize, filesize, status);

        let (raw_ptr, _, _) = b_ptr.into_raw_parts();
        *buffptr = raw_ptr;

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// irafrdhead  (was irafrhead in D. Mink's original code)
/// Open and read the iraf .imh file.
/// The imhdr format is defined in iraf/lib/imhdr.h, some of which defines or mimicked, above.
fn irafrdhead(
    filename: &[c_char], /* Name of IRAF header file */
    lihead: &mut usize,  /* Length of IRAF image header in bytes (returned) */
) -> Result<Vec<c_char>, Error> {
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    *lihead = 0;

    /* open the image header file */
    let fd = File::options().read(true).write(false).open(
        CStr::from_bytes_with_nul(cast_slice(filename))
            .unwrap()
            .to_str()
            .unwrap(),
    );

    if fd.is_err() {
        ffpmsg_str("unable to open IRAF header file:");
        ffpmsg_slice(filename);
        return Err(Error::from(ErrorKind::NotFound));
    }

    let mut fd = fd.unwrap(); //Already checked non-null above

    /* Find size of image header file */
    let pos = fd.seek(std::io::SeekFrom::End(0));
    if pos.is_err() {
        /* move to end of the file */
        ffpmsg_str("IRAFRHEAD: cannot seek in file:");
        ffpmsg_slice(filename);
        return Err(Error::from(ErrorKind::Unsupported));
    }

    let nbhead = pos.unwrap(); /* position = size of file */

    if fd.seek(std::io::SeekFrom::Start(0)).is_err() {
        /* move back to beginning */
        ffpmsg_str("IRAFRHEAD: cannot seek to beginning of file:");
        ffpmsg_slice(filename);
        return Err(Error::from(ErrorKind::Unsupported));
    }

    /* allocate initial sized buffer */
    // HEAP ALLOCATION
    let nihead = nbhead;
    let mut irafheader: Vec<c_char> = Vec::new();
    if irafheader.try_reserve_exact(nihead as usize).is_err() {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "IRAFRHEAD Cannot allocate {}-byte header",
            nihead,
        );
        ffpmsg_slice(&errmsg);
        ffpmsg_slice(filename);
        return Err(Error::from(ErrorKind::OutOfMemory));
    } else {
        irafheader.resize(nihead as usize, 0);
    }

    *lihead = nihead as usize;

    /* Read IRAF header */
    let nbr = fd
        .read(cast_slice_mut(&mut irafheader[..(nbhead as usize)]))
        .unwrap();

    drop(fd); // Close the file

    /* Reject if header less than minimum length */
    if nbr < LEN_PIXHDR {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "IRAFRHEAD header file: {} / {} bytes read.",
            nbr,
            LEN_PIXHDR,
        );
        ffpmsg_slice(&errmsg);
        ffpmsg_slice(filename);

        return Err(Error::from(ErrorKind::InvalidData));
    }

    Ok(irafheader)
}

/*--------------------------------------------------------------------------*/
fn irafrdimage(
    buffptr: &mut Vec<c_char>, /* FITS image header (filled) */
    buffsize: &mut usize,      /* allocated size of the buffer */
    filesize: &mut usize,      /* actual size of the FITS file */
    status: &mut c_int,
) -> c_int {
    let mut nax: c_int = 1;
    let mut naxis1: c_int = 1;
    let mut naxis2: c_int = 1;
    let mut naxis3: c_int = 1;
    let mut naxis4: c_int = 1;
    let mut npaxis1: c_int = 1;
    let mut npaxis2: c_int = 0;
    let mut bitpix: c_int = 0;
    let mut bytepix: c_int = 0;
    let i: c_int = 0;
    let mut nbr: usize = 0;
    let mut nbimage: usize = 0;
    let mut nbaxis: usize = 0;
    let mut nbl: usize = 0;
    let mut nbdiff: c_long = 0;

    let mut imhver: c_int = 0;
    let mut lpixhead: c_int = 0;
    let mut pixname: [c_char; SZ_IM2PIXFILE + 1] = [0; SZ_IM2PIXFILE + 1];
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut newfilesize: usize = 0;

    let fitsheader = buffptr; /* pointer to start of header */

    /* Convert pixel file name to character string */
    hgets(fitsheader, cs!("PIXFILE"), SZ_IM2PIXFILE, &mut pixname);
    hgeti4(fitsheader, cs!("PIXOFF"), &mut lpixhead);

    let lpixhead = lpixhead as usize;

    /* Open pixel file, ignoring machine name if present */
    let bang = strchr_safe(&pixname, bb(b'!'));
    let fd = if let Some(bang) = bang {
        File::options().read(true).write(false).open(
            CStr::from_bytes_with_nul(cast_slice(&pixname[bang + 1..]))
                .unwrap()
                .to_str()
                .unwrap(),
        )
    } else {
        File::options().read(true).write(false).open(
            CStr::from_bytes_with_nul(cast_slice(&pixname))
                .unwrap()
                .to_str()
                .unwrap(),
        )
    };

    /* Print error message and exit if pixel file is not found */
    if fd.is_err() {
        ffpmsg_str("IRAFRIMAGE: Cannot open IRAF pixel file:");
        ffpmsg_slice(&pixname);
        *status = FILE_NOT_OPENED;
        return *status;
    }

    let mut fd = fd.unwrap(); //Already checked non-null above

    /* Read pixel header */
    // HEAP ALLOCATION
    let mut pixheader: Vec<c_char> = Vec::new();
    if pixheader.try_reserve_exact(lpixhead as usize).is_err() {
        ffpmsg_str("IRAFRIMAGE: Cannot alloc memory for pixel header");
        ffpmsg_slice(&pixname);
        drop(fd);
        *status = FILE_NOT_OPENED;
        return *status;
    } else {
        pixheader.resize(lpixhead as usize, 0);
    }

    let r = fd.read(cast_slice_mut(&mut pixheader[..lpixhead])).unwrap();

    /* Check size of pixel header */
    if nbr < lpixhead {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "IRAF pixel file: {} / {} bytes read.",
            nbr,
            LEN_PIXHDR,
        );
        ffpmsg_slice(&errmsg);

        drop(fd);
        *status = FILE_NOT_OPENED;
        return *status;
    }

    /* check pixel header magic word */
    imhver = pix_version(&mut pixheader);
    if imhver < 1 {
        ffpmsg_str("File not valid IRAF pixel file:");
        ffpmsg_slice(&pixname);

        drop(fd);
        *status = FILE_NOT_OPENED;
        return *status;
    }

    /* Find number of bytes to read */
    hgeti4(fitsheader, cs!("NAXIS"), &mut nax);
    hgeti4(fitsheader, cs!("NAXIS1"), &mut naxis1);
    hgeti4(fitsheader, cs!("NPAXIS1"), &mut npaxis1);
    if nax > 1 {
        hgeti4(fitsheader, cs!("NAXIS2"), &mut naxis2);
        hgeti4(fitsheader, cs!("NPAXIS2"), &mut npaxis2);
    }
    if nax > 2 {
        hgeti4(fitsheader, cs!("NAXIS3"), &mut naxis3);
    }
    if nax > 3 {
        hgeti4(fitsheader, cs!("NAXIS4"), &mut naxis4);
    }

    hgeti4(fitsheader, cs!("BITPIX"), &mut bitpix);
    if bitpix < 0 {
        bytepix = -bitpix / 8;
    } else {
        bytepix = bitpix / 8;
    }

    nbimage = (naxis1 * naxis2 * naxis3 * naxis4 * bytepix) as usize; // Will be positive because of above

    newfilesize = *filesize + nbimage; /* header + data */
    newfilesize = (((newfilesize - 1) / BL!()) + 1) * BL!();

    if newfilesize > *buffsize {
        /* need to allocate more memory? */
        if fitsheader
            .try_reserve_exact(newfilesize - fitsheader.len())
            .is_err()
        {
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "IRAFRIMAGE Cannot allocate {}-byte image buffer",
                (*filesize) as c_int,
            );
            ffpmsg_slice(&errmsg);
            ffpmsg_slice(&pixname);
            drop(fd);
            *status = FILE_NOT_OPENED;
            return *status;
        } else {
            fitsheader.resize(newfilesize, 0);
        }
    }

    // *buffptr = fitsheader; No need since using Vec<>'s now
    *buffsize = newfilesize;

    *filesize = newfilesize;

    /* Read IRAF image all at once if physical and image dimensions are the same */
    if npaxis1 == naxis1 {
        let image = &mut fitsheader[*filesize..];

        nbr = fd.read(cast_slice_mut(&mut image[..nbimage])).unwrap();

    /* Read IRAF image one line at a time if physical and image dimensions differ */
    } else {
        nbdiff = ((npaxis1 - naxis1) * bytepix) as c_long;
        nbaxis = (naxis1 * bytepix) as usize;
        let mut linebuff = &mut fitsheader[*filesize..];
        nbr = 0;
        if naxis2 == 1 && naxis3 > 1 {
            naxis2 = naxis3;
        }

        for i in 0..naxis2 {
            nbl = fd.read(cast_slice_mut(&mut linebuff[..nbaxis])).unwrap();
            nbr += nbl;

            let _ = fd.seek(std::io::SeekFrom::Start(0)).unwrap();

            linebuff = &mut linebuff[nbaxis..];
        }
    }

    drop(fd);

    /* Check size of image */
    if nbr < nbimage {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "IRAF pixel file: {} / {} bytes read.",
            nbr,
            nbimage,
        );
        ffpmsg_slice(&errmsg);
        ffpmsg_slice(&pixname);
        *status = FILE_NOT_OPENED;
        return *status;
    }

    /* Byte-reverse image, if necessary */
    if *SWAPDATA.lock().unwrap() {
        let image = &mut fitsheader[*filesize..];
        irafswap(bitpix, image, nbimage);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Return IRAF image format version number from magic word in IRAF header
fn head_version(irafheader: &[c_char] /* IRAF image header from file */) -> c_int {
    /* Check header file magic word */
    if irafncmp(irafheader, cs!("imhdr"), 5) != 0 {
        if strncmp_safe(irafheader, cs!("imhv2"), 5) != 0 {
            0
        } else {
            2
        }
    } else {
        1
    }
}

/*--------------------------------------------------------------------------*/
/// Return IRAF image format version number from magic word in IRAF pixel file
fn pix_version(irafheader: &[c_char] /* IRAF image header from file */) -> c_int {
    /* Check pixel file header magic word */
    if irafncmp(irafheader, cs!("impix"), 5) != 0 {
        if strncmp_safe(irafheader, cs!("impv2"), 5) != 0 {
            0
        } else {
            2
        }
    } else {
        1
    }
}

/*--------------------------------------------------------------------------*/
/// Verify that file is valid IRAF imhdr or impix by checking first 5 chars
/// Returns 0 on success, 1 on failure
fn irafncmp(
    irafheader: &[c_char], /* IRAF image header from file */
    teststring: &[c_char], /* C character string to compare */
    nc: usize,             /* Number of characters to compate */
) -> c_int {
    let line = iraf2str(irafheader, nc);

    if line.is_none() {
        return 1;
    }

    let line = line.unwrap();

    if strncmp_safe(&line, teststring, nc) == 0 {
        0
    } else {
        1
    }
}
/*--------------------------------------------------------------------------*/
/// Convert IRAF image header to FITS image header, returning FITS header
fn iraftofits(
    hdrname: &[c_char],        /* IRAF header file name (may be path) */
    irafheader: &[c_char],     /* IRAF image header */
    nbiraf: usize,             /* Number of bytes in IRAF header */
    buffptr: &mut Vec<c_char>, /* pointer to the FITS header  */
    nbfits: &mut usize,        /* allocated size of the FITS header buffer */
    fitssize: &mut usize,      /* Number of bytes in FITS header (returned) */
    /*  = number of bytes to the end of the END keyword */
    status: &mut c_int,
) -> c_int {
    let mut lstr: usize = 0;
    let i: c_int = 0;
    let mut j: usize = 0;
    let k: c_int = 0;
    let ib: c_int = 0;
    let mut nax: c_int = 0;
    let mut nbits: c_int = 0;
    let mut newpixname: Option<Vec<c_char>>;
    let mut nblock: usize = 0;
    let mut nlines: usize = 0;
    let mut fp: *mut c_char;
    let mut endline: [c_char; 81] = [0; 81];
    let mut irafchar: c_char = 0;
    let mut fitsline: [c_char; 81] = [0; 81];
    let mut pixtype: c_int = 0;
    let mut imhver: c_int = 0;
    let mut n: c_int = 0;
    let mut imu: usize = 0;
    let pixoff: usize = 0;
    let mut impixoff: usize = 0;
    let mut imndim: usize = 0;
    let mut imlen: usize = 0;
    let mut imphyslen: usize = 0;
    let mut impixtype: usize = 0;
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    /* Set up last line of FITS header */
    strncpy_safe(&mut endline, cs!("END"), 3);
    for i in 3..80 {
        endline[i] = bb(b' ');
    }
    endline[80] = 0;

    /* Check header magic word */
    imhver = head_version(irafheader);
    if imhver < 1 {
        ffpmsg_str("File not valid IRAF image header");
        ffpmsg_slice(hdrname);
        *status = FILE_NOT_OPENED;
        return *status;
    }
    if imhver == 2 {
        nlines = 24 + ((nbiraf - LEN_IM2HDR) / 81);
        imndim = IM2_NDIM;
        imlen = IM2_LEN;
        imphyslen = IM2_PHYSLEN;
        impixtype = IM2_PIXTYPE;
        impixoff = IM2_PIXOFF;
    /*	imtime = IM2_MTIME; */
    /*	immax = IM2_MAX;  */
    /*	immin = IM2_MIN; */
    } else {
        nlines = 24 + ((nbiraf - LEN_IMHDR) / 162);
        imndim = IM_NDIM;
        imlen = IM_LEN;
        imphyslen = IM_PHYSLEN;
        impixtype = IM_PIXTYPE;
        impixoff = IM_PIXOFF;
        /*	imtime = IM_MTIME; */
        /*	immax = IM_MAX; */
        /*	immin = IM_MIN; */
    }

    /*  Initialize FITS header */
    nblock = (nlines * 80) / BL!();
    *nbfits = (nblock + 5) * BL!() + 4;

    // HEAP ALLOCATION
    let mut fitsheader = Vec::new();
    if fitsheader.try_reserve_exact(*nbfits).is_err() {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "IRAF2FITS Cannot allocate {}-byte FITS header",
            (*nbfits) as c_int,
        );
        ffpmsg_slice(hdrname);
        *status = FILE_NOT_OPENED;
        return *status;
    } else {
        fitsheader.resize(*nbfits, 0);
    }

    let mut fhead = 0; // &fitsheader

    *buffptr = fitsheader;

    let fitsheader = buffptr;
    strncpy_safe(fitsheader, &endline, 80);
    hputl(fitsheader, cs!("SIMPLE"), 1);
    fhead += 80;

    /*  check if the IRAF file is in big endian (sun) format (= 0) or not. */
    /*  This is done by checking the 4 byte integer in the header that     */
    /*  represents the iraf pixel type.  This 4-byte word is guaranteed to */
    /*  have the least sig byte != 0 and the most sig byte = 0,  so if the */
    /*  first byte of the word != 0, then the file in little endian format */
    /*  like on an Alpha machine.                                          */

    let mut swaphead = SWAPHEAD.lock().unwrap();
    let mut swapdata = SWAPDATA.lock().unwrap();

    *swaphead = isirafswapped(irafheader, impixtype);
    if imhver == 1 {
        *swapdata = *swaphead; /* vers 1 data has same swapness as header */
    } else {
        *swapdata = irafgeti4(irafheader, IM2_SWAPPED) != 0;
    }

    /*  Set pixel size in FITS header */
    pixtype = irafgeti4(irafheader, impixtype);
    match pixtype {
        TY_CHAR => {
            nbits = 8;
        }
        TY_UBYTE => {
            nbits = 8;
        }
        TY_SHORT => {
            nbits = 16;
        }
        TY_USHORT => {
            nbits = -16;
        }
        TY_INT => {}
        TY_LONG => {
            nbits = 32;
        }
        TY_REAL => {
            nbits = -32;
        }
        TY_DOUBLE => {
            nbits = -64;
        }
        _ => {
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "Unsupported IRAF data type: {}",
                pixtype,
            );
            ffpmsg_slice(&errmsg);
            ffpmsg_slice(hdrname);
            *status = FILE_NOT_OPENED;
            return *status;
        }
    }

    hputi4(fitsheader, cs!("BITPIX"), nbits);
    hputcom(fitsheader, cs!("BITPIX"), cs!("IRAF .imh pixel type"));
    fhead += 80;

    /*  Set image dimensions in FITS header */
    nax = irafgeti4(irafheader, imndim);
    hputi4(fitsheader, cs!("NAXIS"), nax);
    hputcom(fitsheader, cs!("NAXIS"), cs!("IRAF .imh naxis"));
    fhead += 80;

    n = irafgeti4(irafheader, imlen);
    hputi4(fitsheader, cs!("NAXIS1"), n);
    hputcom(fitsheader, cs!("NAXIS1"), cs!("IRAF .imh image naxis[1]"));
    fhead += 80;

    if nax > 1 {
        n = irafgeti4(irafheader, imlen + 4);
        hputi4(fitsheader, cs!("NAXIS2"), n);
        hputcom(fitsheader, cs!("NAXIS2"), cs!("IRAF .imh image naxis[2]"));
        fhead += 80;
    }
    if nax > 2 {
        n = irafgeti4(irafheader, imlen + 8);
        hputi4(fitsheader, cs!("NAXIS3"), n);
        hputcom(fitsheader, cs!("NAXIS3"), cs!("IRAF .imh image naxis[3]"));
        fhead += 80;
    }
    if nax > 3 {
        n = irafgeti4(irafheader, imlen + 12);
        hputi4(fitsheader, cs!("NAXIS4"), n);
        hputcom(fitsheader, cs!("NAXIS4"), cs!("IRAF .imh image naxis[4]"));
        fhead += 80;
    }

    /* Set object name in FITS header */
    let objname: Option<Vec<c_char>> = if imhver == 2 {
        irafgetc(irafheader, IM2_TITLE, SZ_IM2TITLE)
    } else {
        irafgetc2(irafheader, IM_TITLE, SZ_IMTITLE)
    };

    let mut objname = objname.unwrap();

    lstr = strlen_safe(&objname);
    if lstr < 8 {
        objname.resize(8, 0);
        for i in lstr..8 {
            objname[i] = bb(b' ');
        }
        objname[8] = 0;
    }
    hputs(fitsheader, cs!("OBJECT"), &objname);
    hputcom(fitsheader, cs!("OBJECT"), cs!("IRAF .imh title"));

    fhead += 80;

    /* Save physical axis lengths so image file can be read */
    n = irafgeti4(irafheader, imphyslen);
    hputi4(fitsheader, cs!("NPAXIS1"), n);
    hputcom(
        fitsheader,
        cs!("NPAXIS1"),
        cs!("IRAF .imh physical naxis[1]"),
    );
    fhead += 80;
    if nax > 1 {
        n = irafgeti4(irafheader, imphyslen + 4);
        hputi4(fitsheader, cs!("NPAXIS2"), n);
        hputcom(
            fitsheader,
            cs!("NPAXIS2"),
            cs!("IRAF .imh physical naxis[2]"),
        );
        fhead += 80;
    }
    if nax > 2 {
        n = irafgeti4(irafheader, imphyslen + 8);
        hputi4(fitsheader, cs!("NPAXIS3"), n);
        hputcom(
            fitsheader,
            cs!("NPAXIS3"),
            cs!("IRAF .imh physical naxis[3]"),
        );
        fhead += 80;
    }
    if nax > 3 {
        n = irafgeti4(irafheader, imphyslen + 12);
        hputi4(fitsheader, cs!("NPAXIS4"), n);
        hputcom(
            fitsheader,
            cs!("NPAXIS4"),
            cs!("IRAF .imh physical naxis[4]"),
        );
        fhead += 80;
    }

    /* Save image header filename in header */
    hputs(fitsheader, cs!("IMHFILE"), hdrname);
    hputcom(fitsheader, cs!("IMHFILE"), cs!("IRAF header file name"));
    fhead += 80;

    /* Save image pixel file pathname in header */
    let pixname: Option<Vec<c_char>> = if imhver == 2 {
        irafgetc(irafheader, IM2_PIXFILE, SZ_IM2PIXFILE)
    } else {
        irafgetc2(irafheader, IM_PIXFILE, SZ_IMPIXFILE)
    };

    let mut pixname = pixname.unwrap();

    if strncmp_safe(&pixname, cs!("HDR"), 3) == 0 {
        newpixname = same_path(&pixname, hdrname);
        if let Some(newpixname) = newpixname {
            pixname = newpixname;
        }
    }
    if strchr_safe(&pixname, bb(b'/')).is_none() && strchr_safe(&pixname, bb(b'$')).is_none() {
        newpixname = same_path(&pixname, hdrname);
        if let Some(newpixname) = newpixname {
            pixname = newpixname;
        }
    }

    let bang = strchr_safe(&pixname, bb(b'!'));

    if let Some(bang) = bang {
        hputs(fitsheader, cs!("PIXFILE"), &pixname[(bang + 1)..]);
    } else {
        hputs(fitsheader, cs!("PIXFILE"), &pixname);
    }

    hputcom(fitsheader, cs!("PIXFILE"), cs!("IRAF .pix pixel file"));
    fhead += 80;

    /* Save image offset from star of pixel file */
    let mut pixoff = irafgeti4(irafheader, impixoff);
    pixoff = (pixoff - 1) * 2;
    hputi4(fitsheader, cs!("PIXOFF"), pixoff);
    hputcom(
        fitsheader,
        cs!("PIXOFF"),
        cs!("IRAF .pix pixel offset (Do not change!)"),
    );
    fhead += 80;

    /* Save IRAF file format version in header */
    hputi4(fitsheader, cs!("IMHVER"), imhver);
    hputcom(
        fitsheader,
        cs!("IMHVER"),
        cs!("IRAF .imh format version (1 or 2)"),
    );
    fhead += 80;

    /* Save flag as to whether to swap IRAF data for this file and machine */
    if *swapdata {
        hputl(fitsheader, cs!("PIXSWAP"), 1);
    } else {
        hputl(fitsheader, cs!("PIXSWAP"), 0);
    }
    hputcom(
        fitsheader,
        cs!("PIXSWAP"),
        cs!("IRAF pixels, FITS byte orders differ if T"),
    );
    fhead += 80;

    /* Add user portion of IRAF header to FITS header */
    fitsline[80] = 0;
    if imhver == 2 {
        imu = LEN_IM2HDR;
        let chead = irafheader;
        j = 0;
        for k in 0..80 {
            fitsline[k] = bb(b' ');
        }
        for i in imu..nbiraf {
            irafchar = chead[i];
            if irafchar == 0 {
                break;
            } else if irafchar == 10 {
                strncpy_safe(&mut fitsheader[fhead..], &fitsline, 80);
                /* fprintf (stderr,"%80s\n",fitsline); */
                if strncmp_safe(&fitsline, cs!("OBJECT "), 7) != 0 {
                    fhead += 80;
                }
                for k in 0..80 {
                    fitsline[k] = bb(b' ');
                }
                j = 0;
            } else {
                if j > 80 {
                    if strncmp_safe(&fitsline, cs!("OBJECT "), 7) != 0 {
                        strncpy_safe(&mut fitsheader[fhead..], &fitsline, 80);
                        /* fprintf (stderr,"%80s\n",fitsline); */
                        j = 9;
                        fhead += 80;
                    }
                    for k in 0..80 {
                        fitsline[k] = bb(b' ');
                    }
                }
                if irafchar > 32 && irafchar < 127 {
                    fitsline[j] = irafchar;
                }
                j += 1;
            }
        }
    } else {
        imu = LEN_IMHDR;
        let chead = irafheader;
        let ib: usize = if *swaphead { 0 } else { 1 };

        for k in 0..80 {
            fitsline[k] = bb(b' ');
        }
        let mut j = 0;
        for i in (imu..nbiraf).step_by(2) {
            irafchar = chead[i + ib];
            if irafchar == 0 {
                break;
            } else if irafchar == 10 {
                if strncmp_safe(&fitsline, cs!("OBJECT "), 7) != 0 {
                    strncpy_safe(&mut fitsheader[fhead..], &fitsline, 80);
                    fhead += 80;
                }
                /* fprintf (stderr,"%80s\n",fitsline); */
                j = 0;
                for k in 0..80 {
                    fitsline[k] = bb(b' ');
                }
            } else {
                if j > 80 {
                    if strncmp_safe(&fitsline, cs!("OBJECT "), 7) != 0 {
                        strncpy_safe(&mut fitsheader[fhead..], &fitsline, 80);
                        j = 9;
                        fhead += 80;
                    }
                    /* fprintf (stderr,"%80s\n",fitsline); */
                    for k in 0..80 {
                        fitsline[k] = bb(b' ');
                    }
                }
                if irafchar > 32 && irafchar < 127 {
                    fitsline[j] = irafchar;
                }
                j += 1;
            }
        }
    }

    /* Add END to last line */
    strncpy_safe(&mut fitsheader[fhead..], &endline, 80);

    /* Find end of last 2880-byte block of header */
    fhead = ksearch(fitsheader, cs!("END")).unwrap() + 80;
    nblock = *nbfits / BL!();
    let fhead1 = nblock * BL!(); //&fitsheader
    *fitssize = fhead; /* no. of bytes to end of END keyword */

    /* Pad rest of header with spaces */
    strncpy_safe(&mut endline, cs!("   "), 3);
    for fp in (fhead..fhead1).step_by(80) {
        strncpy_safe(&mut fitsheader[fp..], &endline, 80);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// get the IRAF pixel file name
fn getirafpixname(
    hdrname: &[c_char],         /* IRAF header file name (may be path) */
    irafheader: &[c_char],      /* IRAF image header */
    pixfilename: &mut [c_char], /* IRAF pixel file name */
    status: &mut c_int,
) -> c_int {
    let mut newpixname = None;

    /* Check header magic word */
    let imhver = head_version(irafheader);
    if imhver < 1 {
        ffpmsg_str("File not valid IRAF image header");
        ffpmsg_slice(hdrname);
        *status = FILE_NOT_OPENED;
        return *status;
    }

    /* get image pixel file pathname in header */
    let mut pixname = if imhver == 2 {
        irafgetc(irafheader, IM2_PIXFILE, SZ_IM2PIXFILE)
    } else {
        irafgetc2(irafheader, IM_PIXFILE, SZ_IMPIXFILE)
    }
    .unwrap();

    if strncmp_safe(&pixname, cs!("HDR"), 3) == 0 {
        newpixname = same_path(&pixname, hdrname);
        if let Some(newpixname) = newpixname {
            pixname = newpixname;
        }
    }

    if strchr_safe(&pixname, bb(b'/')).is_none() && strchr_safe(&pixname, bb(b'$')).is_none() {
        newpixname = same_path(&pixname, hdrname);
        if let Some(newpixname) = newpixname {
            pixname = newpixname;
        }
    }

    let bang = strchr_safe(&pixname, bb(b'!'));
    if let Some(bang) = bang {
        strcpy_safe(pixfilename, &pixname[(bang + 1)..]);
    } else {
        strcpy_safe(pixfilename, &pixname);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/* Put filename and header path together */
fn same_path(
    pixname: &[c_char], /* IRAF pixel file pathname */
    hdrname: &[c_char], /* IRAF image header file pathname */
) -> Option<Vec<c_char>> {
    /*  WDP - 10/16/2007 - increased allocation to avoid possible overflow */
    /*    newpixname =  calloc (SZ_IM2PIXFILE, sizeof (char)); */

    let mut newpixname = Vec::new();
    if newpixname.try_reserve_exact(2 * SZ_IM2PIXFILE + 1).is_err() {
        ffpmsg_str("iraffits same_path: Cannot alloc memory for newpixname");
        return None;
    } else {
        newpixname.resize(2 * SZ_IM2PIXFILE + 1, 0)
    }

    /* Pixel file is in same directory as header */
    if strncmp_safe(pixname, cs!("HDR$"), 4) == 0 {
        strncpy_safe(&mut newpixname, hdrname, SZ_IM2PIXFILE);

        /* find the end of the pathname */
        let mut len = strlen_safe(&newpixname);
        if cfg!(vms) {
            while (len > 0) && (newpixname[len - 1] != bb(b'/')) {
                len -= 1;
            }
        } else {
            while (len > 0)
                && (newpixname[len - 1] != bb(b']'))
                && (newpixname[len - 1] != bb(b':'))
            {
                len -= 1;
            }
        }

        /* add name */
        newpixname[len] = 0;
        strncat_safe(&mut newpixname, &pixname[4..], SZ_IM2PIXFILE);
    }
    /* Bare pixel file with no path is assumed to be same as HDR$filename */
    else if strchr_safe(pixname, bb(b'/')).is_none() && strchr_safe(pixname, bb(b'$')).is_none() {
        strncpy_safe(&mut newpixname, hdrname, SZ_IM2PIXFILE);

        /* find the end of the pathname */
        let mut len = strlen_safe(&newpixname);
        if cfg!(vms) {
            while (len > 0) && (newpixname[len - 1] != bb(b'/')) {
                len -= 1;
            }
        } else {
            while (len > 0)
                && (newpixname[len - 1] != bb(b']'))
                && (newpixname[len - 1] != bb(b':'))
            {
                len -= 1;
            }
        }

        /* add name */
        newpixname[len] = 0;
        strncat_safe(&mut newpixname, pixname, SZ_IM2PIXFILE);
    }
    /* Pixel file has same name as header file, but with .pix extension */
    else if strncmp_safe(pixname, cs!("HDR"), 3) == 0 {
        /* load entire header name string into name buffer */
        strncpy_safe(&mut newpixname, hdrname, SZ_IM2PIXFILE);
        let len = strlen_safe(&newpixname);
        newpixname[len - 3] = bb(b'p');
        newpixname[len - 2] = bb(b'i');
        newpixname[len - 1] = bb(b'x');
    }

    Some(newpixname)
}

/*--------------------------------------------------------------------------*/
///  check if the IRAF file is in big endian (sun) format (= 0) or not
///  This is done by checking the 4 byte integer in the header that
///  represents the iraf pixel type.  This 4-byte word is guaranteed to
///  have the least sig byte != 0 and the most sig byte = 0,  so if the
///  first byte of the word != 0, then the file in little endian format
///  like on an Alpha machine.                                          
fn isirafswapped(
    irafheader: &[c_char], /* IRAF image header */
    offset: usize,         /* Number of bytes to skip before number */
) -> bool {
    let mut swapped = false;

    swapped = irafheader[offset] != 0;

    swapped
}

/*--------------------------------------------------------------------------*/
fn irafgeti4(
    irafheader: &[c_char], /* IRAF image header */
    offset: usize,         /* Number of bytes to skip before number */
) -> c_int {
    let cheader = irafheader;

    let mut itemp: [c_int; 1] = [0; 1];
    let ctemp: &mut [c_char] = cast_slice_mut(&mut itemp);

    let swaphead = SWAPHEAD.lock().unwrap();

    if machswap() != *swaphead {
        ctemp[3] = cheader[offset];
        ctemp[2] = cheader[offset + 1];
        ctemp[1] = cheader[offset + 2];
        ctemp[0] = cheader[offset + 3];
    } else {
        ctemp[0] = cheader[offset];
        ctemp[1] = cheader[offset + 1];
        ctemp[2] = cheader[offset + 2];
        ctemp[3] = cheader[offset + 3];
    }

    itemp[0]
}

/*--------------------------------------------------------------------------*/
/// IRAFGETC2 -- Get character string from arbitrary part of v.1 IRAF header
fn irafgetc2(
    irafheader: &[c_char], /* IRAF image header */
    offset: usize,         /* Number of bytes to skip before string */
    nc: usize,             /* Maximum number of characters in string */
) -> Option<Vec<c_char>> {
    let irafstring = irafgetc(irafheader, offset, 2 * (nc + 1));

    iraf2str(&mut irafstring.unwrap(), nc)
}

/*--------------------------------------------------------------------------*/
/// IRAFGETC -- Get character string from arbitrary part of IRAF header
fn irafgetc(
    irafheader: &[c_char], /* IRAF image header */
    offset: usize,         /* Number of bytes to skip before string */
    nc: usize,             /* Maximum number of characters in string */
) -> Option<Vec<c_char>> {
    let cheader = irafheader;

    let mut ctemp = Vec::new();

    if ctemp.try_reserve_exact(nc + 1).is_err() {
        ffpmsg_str("IRAFGETC Cannot allocate memory for string variable");
        return None;
    }

    for i in 0..nc {
        ctemp.push(cheader[offset + i]);
        if ctemp[i] > 0 && ctemp[i] < 32 {
            ctemp[i] = bb(b' ');
        }
    }

    Some(ctemp)
}

/*--------------------------------------------------------------------------*/
/// Convert IRAF 2-byte/char string to 1-byte/char string
fn iraf2str(
    irafstring: &[c_char], /* IRAF 2-byte/character string */
    nchar: usize,          /* Number of characters in string */
) -> Option<Vec<c_char>> {
    let mut string = Vec::new();

    if string.try_reserve_exact(nchar + 1).is_err() {
        ffpmsg_str("IRAF2STR Cannot allocate memory for string variable");
        return None;
    }

    /* the chars are in bytes 1, 3, 5, ... if bigendian format (SUN) */
    /* else in bytes 0, 2, 4, ... if little endian format (Alpha)    */

    let mut j = if irafstring[0] != 0 { 0 } else { 1 };

    /* Convert appropriate byte of input to output character */
    for i in 0..(nchar) {
        string.push(irafstring[j]);
        j += 2;
    }

    Some(string)
}

/*--------------------------------------------------------------------------*/
/// IRAFSWAP -- Reverse bytes of any type of vector in place
fn irafswap(
    bitpix: c_int, /* Number of bits per pixel */
    /*  16 = short, -16 = unsigned short, 32 = int */
    /* -32 = float, -64 = double */
    string: &mut [c_char], /* Address of starting point of bytes to swap */
    nbytes: usize,         /* Number of bytes to swap */
) {
    match bitpix {
        16 => {
            if nbytes < 2 {
                return;
            }
            irafswap2(string, nbytes);
        }
        32 => {
            if nbytes < 4 {
                return;
            }
            irafswap4(string, nbytes);
        }
        -16 => {
            if nbytes < 2 {
                return;
            }
            irafswap2(string, nbytes);
        }
        -32 => {
            if nbytes < 4 {
                return;
            }
            irafswap4(string, nbytes);
        }
        -64 => {
            if nbytes < 8 {
                return;
            }
            irafswap8(string, nbytes);
        }
        _ => (),
    }
}

/*--------------------------------------------------------------------------*/
///IRAFSWAP2 -- Swap bytes in string in place
fn irafswap2(
    string: &mut [c_char], /* Address of starting point of bytes to swap */
    nbytes: usize,         /* Number of bytes to swap */
) {
    let slast = nbytes;
    let mut sbyte = 0;
    while sbyte < slast {
        string.swap(sbyte, sbyte + 1);
        sbyte += 2;
    }
}

/*--------------------------------------------------------------------------*/
/// IRAFSWAP4 -- Reverse bytes of Integer*4 or Real*4 vector in place
fn irafswap4(
    string: &mut [c_char], /* Address of Integer*4 or Real*4 vector */
    nbytes: usize,         /* Number of bytes to reverse */
) {
    let slast = nbytes;
    let mut sbyte = 0;
    while sbyte < slast {
        let temp3 = string[sbyte];
        let temp2 = string[sbyte + 1];
        let temp1 = string[sbyte + 2];
        let temp0 = string[sbyte + 3];
        string[sbyte] = temp0;
        string[sbyte + 1] = temp1;
        string[sbyte + 2] = temp2;
        string[sbyte + 3] = temp3;
        sbyte += 4;
    }
}

/*--------------------------------------------------------------------------*/
/// IRAFSWAP8 -- Reverse bytes of Real*8 vector in place
fn irafswap8(
    string: &mut [c_char], /* Address of Real*8 vector */
    nbytes: usize,         /* Number of bytes to reverse */
) {
    let slast = nbytes;
    let mut sbyte = 0;
    while sbyte < slast {
        let temp7 = string[sbyte];
        let temp6 = string[sbyte + 1];
        let temp5 = string[sbyte + 2];
        let temp4 = string[sbyte + 3];
        let temp3 = string[sbyte + 4];
        let temp2 = string[sbyte + 5];
        let temp1 = string[sbyte + 6];
        let temp0 = string[sbyte + 7];
        string[sbyte] = temp0;
        string[sbyte + 1] = temp1;
        string[sbyte + 2] = temp2;
        string[sbyte + 3] = temp3;
        string[sbyte + 4] = temp4;
        string[sbyte + 5] = temp5;
        string[sbyte + 6] = temp6;
        string[sbyte + 7] = temp7;
        sbyte += 8;
    }
}

/*--------------------------------------------------------------------------*/
fn machswap() -> bool {
    let mut itest: [c_int; 1] = [0; 1];
    let ctest: &[c_char] = cast_slice_mut(&mut itest);

    ctest[0] != 0
}

/*--------------------------------------------------------------------------*/
/*             the following routines were originally in hget.c             */
/*--------------------------------------------------------------------------*/

static LHEAD0: usize = 0;

/*--------------------------------------------------------------------------*/
/// Extract long value for variable from FITS header string
fn hgeti4(
    hstring: &mut [c_char], /* character string containing FITS header information
                            in the format <keyword>= <value> {/ <comment>} */
    keyword: &[c_char], /* character string containing the name of the keyword
                        the value of which is returned.  hget searches for a
                        line beginning with this string.  if "[n]" is present,
                        the n'th token in the value is returned.
                        (the first 8 characters must be unique) */
    ival: &mut c_int,
) -> c_int {
    let mut dval: f64 = 0.0;
    let mut minint: c_int = 0;
    let mut val: [c_char; 30] = [0; 30];

    /* Get value and comment from header string */
    let value = hgetc(hstring, keyword);

    /* Translate value from ASCII to binary */
    if let Some(value) = value {
        minint = c_int::MIN;
        if strlen_safe(&value) > 29 {
            return 0;
        }
        strcpy_safe(&mut val, &value);
        dval = atof_safe(&val);
        if dval + 0.001 > MAXINT as f64 {
            *ival = MAXINT;
        } else if dval >= 0.0 {
            *ival = (dval + 0.001) as c_int;
        } else if dval - 0.001 < minint as f64 {
            *ival = minint;
        } else {
            *ival = (dval - 0.001) as c_int;
        }
        1
    } else {
        0
    }
}

/*-------------------------------------------------------------------*/
/// Extract string value for variable from FITS header string
fn hgets(
    hstring: &mut [c_char], /* character string containing FITS header information
                            in the format <keyword>= <value> {/ <comment>} */
    keyword: &[c_char], /* character string containing the name of the keyword
                        the value of which is returned.  hget searches for a
                        line beginning with this string.  if "[n]" is present,
                        the n'th token in the value is returned.
                        (the first 8 characters must be unique) */
    lstr: usize,        /* Size of str in characters */
    str: &mut [c_char], /* String (returned) */
) -> c_int {
    let mut lval: usize = 0;

    /* Get value and comment from header string */
    let value = hgetc(hstring, keyword);

    if let Some(value) = value {
        lval = strlen_safe(&value);
        if lval < lstr {
            strcpy_safe(str, &value);
        } else if lstr > 1 {
            strncpy_safe(str, &value, lstr - 1);
            str[lstr - 1] = 0;
        } else {
            str[0] = value[0];
        }
        1
    } else {
        0
    }
}

/*-------------------------------------------------------------------*/
/// Extract character value for variable from FITS header string
fn hgetc(
    hstring: &mut [c_char], /* character string containing FITS header information
                            in the format <keyword>= <value> {/ <comment>} */
    keyword0: &[c_char], /* character string containing the name of the keyword
                         the value of which is returned.  hget searches for a
                         line beginning with this string.  if "[n]" is present,
                         the n'th token in the value is returned.
                         (the first 8 characters must be unique) */
) -> Option<Vec<c_char>> {
    let mut cval: [c_char; 80] = [0; 80];
    let value: Option<usize>;
    let mut squot: [c_char; 2] = [0; 2];
    let mut dquot: [c_char; 2] = [0; 2];
    let mut lbracket: [c_char; 2] = [0; 2];
    let mut rbracket: [c_char; 2] = [0; 2];
    let mut slash: [c_char; 2] = [0; 2];
    let mut comma: [c_char; 2] = [0; 2];
    let mut keyword: [c_char; 81] = [0; 81]; /* large for ESO hierarchical keywords */
    let mut line: [c_char; 100] = [0; 100];

    squot[0] = 39;
    squot[1] = 0;
    dquot[0] = 34;
    dquot[1] = 0;
    lbracket[0] = 91;
    lbracket[1] = 0;
    comma[0] = 44;
    comma[1] = 0;
    rbracket[0] = 93;
    rbracket[1] = 0;
    slash[0] = 47;
    slash[1] = 0;

    /* Find length of variable name */
    let key_len = keyword.len();
    strncpy_safe(&mut keyword, keyword0, key_len - 1);
    keyword[80] = 0;
    let mut brack1 = strsrch(&keyword, &lbracket);
    if brack1.is_none() {
        brack1 = strsrch(&keyword, &comma);
    }

    if let Some(mut brack1) = brack1 {
        keyword[brack1] = 0;
        brack1 += 1;
    }

    /* Search header string for variable name */
    let vpos = ksearch(hstring, &keyword);

    /* Exit if not found */
    vpos?;

    let vpos = vpos.unwrap();

    /* Initialize line to nulls */
    for i in 0..100 {
        line[i] = 0;
    }

    /* In standard FITS, data lasts until 80th character */

    /* Extract entry for this variable from the header */
    strncpy_safe(&mut line, &hstring[vpos..], 80);

    /* check for quoted value */
    let mut q1 = strsrch(&line, &squot);
    let mut c1 = strsrch(&line, &slash);
    let mut q2 = None;

    if let Some(q1_tmp) = q1 {
        if let Some(c1_tmp) = c1 {
            if q1_tmp < c1_tmp {
                q2 = strsrch(&line[(q1_tmp + 1)..], &squot);
            } else {
                q1 = None;
            }
        } else {
            q2 = strsrch(&line[(q1_tmp + 1)..], &squot);
        }
    } else {
        q1 = strsrch(&line, &dquot);
        if let Some(q1_tmp) = q1 {
            if let Some(c1_tmp) = c1 {
                if q1_tmp < c1_tmp {
                    q2 = strsrch(&line[(q1_tmp + 1)..], &dquot);
                } else {
                    q1 = None;
                }
            } else {
                q2 = strsrch(&line[(q1_tmp + 1)..], &dquot);
            }
        } else {
            q1 = None;
            q2 = Some(10); // &line
        }
    }

    /* Extract value and remove excess spaces */
    let mut v1;
    let mut v2;
    if let Some(q1) = q1 {
        q2?;
        let q2 = q2.unwrap();

        v1 = q1 + 1;
        v2 = q2;
        c1 = strsrch(&line[q2..], cs!("/"));
    } else {
        v1 = strsrch(&line, cs!("=")).expect("Shouldn't be null, not handled in original code") + 1;
        c1 = strsrch(&line, cs!("/"));
        if let Some(c1) = c1 {
            v2 = c1;
        } else {
            v2 = 79; //&line
        }
    }

    /* Ignore leading spaces */
    while line[v1] == bb(b' ') && v1 < v2 {
        v1 += 1;
    }

    /* Drop trailing spaces */
    line[v2] = 0;
    v2 -= 1;
    while line[v2] == bb(b' ') && v2 > v1 {
        line[v2] = 0;
        v2 -= 1;
    }

    if strcmp_safe(&line[v1..], cs!("-0")) == 0 {
        v1 += 1;
    }
    strcpy_safe(&mut cval, &line[v1..]);

    /* If keyword has brackets, extract appropriate token from value */
    if let Some(brack1) = brack1 {
        let brack2 = strsrch(&keyword[brack1..], &rbracket);
        if let Some(brack2) = brack2 {
            keyword[brack2] = 0;
        }

        let ipar = atoi_safe(&keyword[brack1..]);
        if ipar > 0 {
            let mut cpar = None;
            let v_tok = &line[v1..];

            let mut tokens = v_tok.split(|f| (*f == b' ' as c_char) || (*f == 0));

            for i in 1..=ipar {
                cpar = tokens.next();
            }

            if let Some(cpar) = cpar {
                let len = cpar.len();
                if len > 0 {
                    cval[0..len].copy_from_slice(cpar);
                    cval[len] = 0; // Null termiante
                }
            } else {
                value = None;
            }
        }
    }

    Some(Vec::from(cval))
}

/*-------------------------------------------------------------------*/
/// Find beginning of fillable blank line before FITS header keyword line
fn blsearch(
    /* Find entry for keyword keyword in FITS header string hstring.
    (the keyword may have a maximum of eight letters)
    ptr::null_mut() is returned if the keyword is not found */
    hstring: &mut [c_char], /* character string containing fits-style header
                            information in the format <keyword>= <value> {/ <comment>}
                            the default is that each entry is 80 characters long;
                            however, lines may be of arbitrary length terminated by
                            nulls, carriage returns or linefeeds, if packed is true.  */
    keyword: &[c_char], /* character string containing the name of the variable
                        to be returned.  ksearch searches for a line beginning
                        with this string.  The string may be a character
                        literal or a character variable terminated by a null
                        or '$'.  it is truncated to 8 characters. */
) -> Option<usize> {
    let mut lhstr: usize = 0;

    /* Search header string for variable name */
    if LHEAD0 != 0 {
        lhstr = LHEAD0;
    } else {
        lhstr = 0;
        while lhstr < 57600 && hstring[lhstr] != 0 {
            lhstr += 1;
        }
    }

    let headlast = lhstr;
    let mut headnext = 0;
    let mut pval = None;
    while headnext < headlast {
        let nleft = headlast - headnext;
        let loc = strnsrch(&hstring[headnext..], keyword, nleft);

        /* Exit if keyword is not found */
        if loc.is_none() {
            break;
        }

        let loc = loc.unwrap();

        let icol = (loc) % 80;
        let lkey = strlen_safe(keyword);
        let nextchar = hstring[loc + lkey];

        /* If this is not in the first 8 characters of a line, keep searching */
        if icol > 7 {
            headnext = loc + 1;

        /* If parameter name in header is longer, keep searching */
        } else if nextchar != 61 && nextchar > 32 && nextchar < 127 {
            headnext = loc + 1;

        /* If preceeding characters in line are not blanks, keep searching */
        } else {
            let line = loc - icol;
            for lc in line..loc {
                if hstring[lc] != bb(b' ') {
                    headnext = loc + 1;
                }
            }

            /* Return pointer to start of line if match */
            if loc >= headnext {
                pval = Some(line);
                break;
            }
        }
    }

    /* Return ptr::null_mut() if keyword is found at start of FITS header string */
    pval?;

    let pval = pval.unwrap();

    /* Return ptr::null_mut() if found the first keyword in the header */
    if pval == 0 {
        return None;
    }

    /* Find last nonblank line before requested keyword */
    let mut bval = pval - 80;
    while strncmp_safe(&hstring[bval..], cs!("        "), 8) == 0 {
        bval -= 80;
    }
    bval += 80;

    /* Return pointer to calling program if blank lines found */
    if bval < pval { Some(bval) } else { None }
}

/*-------------------------------------------------------------------*/
/// Find FITS header line containing specified keyword
fn ksearch(
    /* Find entry for keyword keyword in FITS header string hstring.
    (the keyword may have a maximum of eight letters)
    ptr::null_mut() is returned if the keyword is not found */
    hstring: &mut [c_char], /* character string containing fits-style header
                            information in the format <keyword>= <value> {/ <comment>}
                            the default is that each entry is 80 characters long;
                            however, lines may be of arbitrary length terminated by
                            nulls, carriage returns or linefeeds, if packed is true.  */
    keyword: &[c_char], /* character string containing the name of the variable
                        to be returned.  ksearch searches for a line beginning
                        with this string.  The string may be a character
                        literal or a character variable terminated by a null
                        or '$'.  it is truncated to 8 characters. */
) -> Option<usize> {
    let mut icol: usize = 0;
    let mut lkey: usize = 0;
    let mut nleft: usize = 0;
    let mut lhstr: usize = 0;

    /* Search header string for variable name */
    if LHEAD0 != 0 {
        lhstr = LHEAD0;
    } else {
        lhstr = 0;
        while lhstr < 57600 && hstring[lhstr] != 0 {
            lhstr += 1;
        }
    }

    let headlast = lhstr;
    let mut headnext = 0;
    let mut pval = None;
    while headnext < headlast {
        nleft = headlast - headnext;
        let loc = strnsrch(&hstring[headnext..], keyword, nleft);

        /* Exit if keyword is not found */
        if loc.is_none() {
            break;
        }

        let loc = loc.unwrap();

        icol = (loc) % 80;
        lkey = strlen_safe(keyword);
        let nextchar = hstring[loc + lkey];

        /* If this is not in the first 8 characters of a line, keep searching */
        if icol > 7 {
            headnext = loc + 1;

        /* If parameter name in header is longer, keep searching */
        } else if nextchar != 61 && nextchar > 32 && nextchar < 127 {
            headnext = loc + 1;

        /* If preceeding characters in line are not blanks, keep searching */
        } else {
            let line = loc - icol;
            for lc in line..loc {
                if hstring[lc] != bb(b' ') {
                    headnext = loc + 1;
                }
            }

            /* Return pointer to start of line if match */
            if loc >= headnext {
                pval = Some(line);
                break;
            }
        }
    }

    /* Return pointer to calling program */
    pval
}

/*-------------------------------------------------------------------*/
/// Find string s2 within null-terminated string s1
fn strsrch(
    s1: &[c_char], /* String to search */
    s2: &[c_char], /* String to look for */
) -> Option<usize> {
    let ls1 = strlen_safe(s1);
    strnsrch(s1, s2, ls1)
}

/*-------------------------------------------------------------------*/
/// Find string s2 within string s1
fn strnsrch(
    s1: &[c_char], /* String to search */
    s2: &[c_char], /* String to look for */
    ls1: usize,    /* Length of string being searched */
) -> Option<usize> {
    let mut cfirst: c_char = 0;
    let mut clast: c_char = 0;
    let mut i: usize = 0;

    /* Return null string if either pointer is ptr::null_mut() */
    if s1.is_empty() {
        return None;
    }

    /* A zero-length pattern is found in any string */
    let ls2 = strlen_safe(s2);
    if ls2 == 0 {
        return Some(0);
    }

    /* Only a zero-length string can be found in a zero-length string */
    if ls1 == 0 {
        return None;
    }

    cfirst = s2[0];
    clast = s2[ls2 - 1];
    let s1e = ls1 - ls2 + 1; // Index in s1
    let mut s = 0; // Index in S1

    while s < s1e {
        /* Search for first character in pattern string */
        if s1[s] == cfirst {
            /* If single character search, return */
            if ls2 == 1 {
                return Some(s);
            }

            /* Search for last character in pattern string if first found */
            if s1[s + ls2 - 1] == clast {
                /* If two-character search, return */
                if ls2 == 2 {
                    return Some(s);
                }

                /* If 3 or more characters, check for rest of search string */
                i = 1;
                while i < ls2 && s1[s + i] == s2[i] {
                    i += 1;
                }

                /* If entire string matches, return */
                if i >= ls2 {
                    return Some(s);
                }
            }
        }
        s += 1;
    }
    None
}

/*-------------------------------------------------------------------*/
/*             the following routines were originally in hget.c      */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/// HPUTI4 - Set int keyword = ival in FITS header string
fn hputi4(
    hstring: &mut [c_char], /* character string containing FITS-style header
                            information in the format
                            <keyword>= <value> {/ <comment>}
                            each entry is padded with spaces to 80 characters */

    keyword: &[c_char], /* character string containing the name of the variable
                        to be returned.  hput searches for a line beginning
                        with this string, and if there isn't one, creates one.
                           The first 8 characters of keyword must be unique. */
    ival: c_int, /* int number */
) {
    let mut value: [c_char; 30] = [0; 30];

    /* Translate value from binary to ASCII */
    int_snprintf!(&mut value, 30, "{}", ival);

    /* Put value into header string */
    hputc(hstring, keyword, &value);

    /* Return to calling program */
}

/*-------------------------------------------------------------------*/
/// HPUTL - Set keyword = F if lval=0, else T, in FITS header string
fn hputl(
    hstring: &mut [c_char], /* FITS header */
    keyword: &[c_char],     /* Keyword name */
    lval: c_int,            /* logical variable (0=false, else true) */
) {
    let mut value: [c_char; 8] = [0; 8];

    /* Translate value from binary to ASCII */
    if lval != 0 {
        strcpy_safe(&mut value, cs!("T"));
    } else {
        strcpy_safe(&mut value, cs!("F"));
    }

    /* Put value into header string */
    hputc(hstring, keyword, &value);

    /* Return to calling program */
}

/*-------------------------------------------------------------------*/
///  HPUTS - Set character string keyword = 'cval' in FITS header string
fn hputs(
    hstring: &mut [c_char], /* FITS header */
    keyword: &[c_char],     /* Keyword name */
    cval: &[c_char],        /* character string containing the value for variable
                            keyword.  trailing and leading blanks are removed.  */
) {
    let squot = 39;
    let mut value: [c_char; 70] = [0; 70];

    /*  find length of variable string */

    let mut lcval = strlen_safe(cval);
    if lcval > 67 {
        lcval = 67;
    }

    /* Put quotes around string */
    value[0] = squot;
    strncpy_safe(&mut value[1..], cval, lcval);
    value[lcval + 1] = squot;
    value[lcval + 2] = 0;

    /* Put value into header string */
    hputc(hstring, keyword, &value);

    /* Return to calling program */
}

/*---------------------------------------------------------------------*/
/// HPUTC - Set character string keyword = value in FITS header string
fn hputc(
    hstring: &mut [c_char],
    keyword: &[c_char],
    value: &[c_char], /* character string containing the value for variable
                      keyword.  trailing and leading blanks are removed.  */
) {
    let squot: c_char = 39 as c_int as c_char;
    let mut line: [c_char; 100] = [0; 100];
    let mut newcom: [c_char; 50] = [0; 50];
    let mut blank: [c_char; 80] = [0; 80];

    let mut vp = 0;
    let mut v1 = 0;
    let mut v2 = 0;
    let mut q2 = 0;
    let c1;

    let mut lcom: usize = 0;
    let mut lc: usize = 0;

    for i in 0..80 {
        blank[i] = bb(b' ');
    }

    /*  find length of keyword and value */
    let lkeyword = strlen_safe(keyword);
    let lval = strlen_safe(value);

    /*  If COMMENT or HISTORY, always add it just before the END */
    if lkeyword == 7
        && (strncmp_safe(keyword, cs!("COMMENT"), 7) == 0
            || strncmp_safe(keyword, cs!("HISTORY"), 7) == 0)
    {
        /* Find end of header */
        let v1_tmp = ksearch(hstring, cs!("END"));
        if v1_tmp.is_none() {
            return; /* UB in original code */
        } else {
            v1 = v1_tmp.unwrap();
        }
        v2 = v1 + 80;

        /* Move END down one line */
        let (h1, h2) = hstring.split_at_mut(v2);
        //strncpy_safe(&mut hstring[v2..], &hstring[v1..], 80);
        strncpy_safe(h2, &h1[v1..], 80);

        /* Insert keyword */
        strncpy_safe(&mut hstring[v1..], keyword, 7);

        /* Pad with spaces */
        for vp in (v1 + lkeyword)..v2 {
            hstring[vp] = bb(b' ');
        }

        /* Insert comment */
        strncpy_safe(&mut hstring[(v1 + 9)..], value, lval);
        return;
    }

    /* Otherwise search for keyword */
    let v1_tmp = ksearch(hstring, keyword);

    /*  If parameter is not found, find a place to put it */
    if v1_tmp.is_none() {
        /* First look for blank lines before END */
        let v1_tmp = blsearch(hstring, cs!("END"));

        /*  Otherwise, create a space for it at the end of the header */
        if v1_tmp.is_none() {
            let ve_tmp = ksearch(hstring, cs!("END"));
            if ve_tmp.is_none() {
                return; // UB in original code
            }
            v1 = ve_tmp.unwrap();
            v2 = v1 + 80;

            let (h1, h2) = hstring.split_at_mut(v2);
            //strncpy_safe(&mut hstring[v2..], &hstring[v1..], 80);
            strncpy_safe(h2, &h1[v1..], 80);
        } else {
            v1 = v1_tmp.unwrap();
            v2 = v1 + 80;
        }
        lcom = 0;
        newcom[0] = 0;
    } else {
        /*  Otherwise, extract the entry for this keyword from the header */

        strncpy_safe(&mut line, &hstring[v1..], 80);
        line[80] = 0;
        v2 = v1 + 80;

        /*  check for quoted value */
        let q1 = strchr_safe(&line, squot);
        if let Some(q1) = q1 {
            let q2_tmp = strchr_safe(&line[(q1 + 1)..], squot);
            if let Some(x) = q2_tmp {
                q2 = x;
            } else {
                return; // UB in original code
            }
        } else {
            q2 = 0;
        }

        /*  extract comment and remove trailing spaces */

        c1 = strchr_safe(&line[q2..], bb(b'/'));
        if let Some(c1) = c1 {
            lcom = 80 - c1;
            strncpy_safe(&mut newcom, &line[(c1 + 1)..], lcom);
            vp = lcom - 1;

            while vp > 0 {
                vp -= 1;
                if newcom[vp] == bb(b' ') {
                    newcom[vp] = 0;
                }
            }
            lcom = strlen_safe(&newcom);
        } else {
            newcom[0] = 0;
            lcom = 0;
        }
    }

    /* Fill new entry with spaces */
    for vp in v1..v2 {
        hstring[vp] = bb(b' ');
    }

    /*  Copy keyword to new entry */
    strncpy_safe(&mut hstring[v1..], keyword, lkeyword);

    /*  Add parameter value in the appropriate place */
    vp = v1 + 8;
    hstring[vp] = bb(b'=');
    vp = v1 + 9;
    hstring[vp] = bb(b' ');
    vp += 1;
    if value[0] == squot {
        strncpy_safe(&mut hstring[vp..], value, lval);
        if lval + 12 > 31 {
            lc = lval + 12;
        } else {
            lc = 30;
        }
    } else {
        vp = v1 + 30 - lval;
        strncpy_safe(&mut hstring[vp..], value, lval);
        lc = 30;
    }

    /* Add comment in the appropriate place */
    if lcom > 0 {
        if lc + 2 + lcom > 80 {
            lcom = 78 - lc;
        }
        vp = v1 + lc + 2; /* Jul 16 1997: was vp = v1 + lc * 2 */
        hstring[vp] = bb(b'/');
        vp += 1;
        strncpy_safe(&mut hstring[vp..], &newcom, lcom);
        for v in (vp + lcom)..v2 {
            hstring[v] = bb(b' ');
        }
    }
}

/*-------------------------------------------------------------------*/
/// HPUTCOM - Set comment for keyword or on line in FITS header string
fn hputcom(hstring: &mut [c_char], keyword: &[c_char], comment: &[c_char]) {
    let mut squot: c_char = 0;
    let mut line: [c_char; 100] = [0; 100];

    let v1;
    let v2;
    let mut c0 = 0;
    let c1;
    let q2;

    squot = 39;

    /*  Find length of variable name */
    let lkeyword = strlen_safe(keyword);

    /*  If COMMENT or HISTORY, always add it just before the END */
    if lkeyword == 7
        && (strncmp_safe(keyword, cs!("COMMENT"), 7) == 0
            || strncmp_safe(keyword, cs!("HISTORY"), 7) == 0)
    {
        /* Find end of header */
        let v1_tmp = ksearch(hstring, cs!("END"));
        if v1_tmp.is_none() {
            return; /* UB in original code */
        } else {
            v1 = v1_tmp.unwrap();
        }

        v2 = v1 + 80;

        let (h1, h2) = hstring.split_at_mut(v2);
        //strncpy_safe(&mut hstring[v2..], &hstring[v1..], 80);
        strncpy_safe(h2, &h1[v1..], 80);

        /*  blank out new line and insert keyword */
        for vp in v1..v2 {
            hstring[vp] = bb(b' ');
        }

        strncpy_safe(&mut hstring[v1..], keyword, lkeyword);
    } else {
        /* search header string for variable name */

        let v1_tmp = ksearch(hstring, keyword);
        if v1_tmp.is_none() {
            return; /* if parameter is not found, return without doing anything */
        } else {
            v1 = v1_tmp.unwrap();
        }

        v2 = v1 + 80;

        /* otherwise, extract entry for this variable from the header */
        strncpy_safe(&mut line, &hstring[v1..], 80);

        /* check for quoted value */
        let q1 = strchr_safe(&line, squot);
        q2 = match q1 {
            Some(x) => strchr_safe(&line[(x + 1)..], squot),
            None => None,
        };

        if let Some(q2) = q2 {
            if q2 < 31 {
                c0 = v1 + 31;
            } else {
                c0 = v1 + q2 + 2; /* allan: 1997-09-30, was c0=q2+2 */
            }
        } else {
            c0 = v1 + 31;
        }

        strncpy_safe(&mut line[c0..], cs!("/ "), 2);
    }

    /* create new entry */
    let mut lcom = strlen_safe(comment);

    if lcom > 0 {
        c1 = c0 + 2;
        if c1 + lcom > v2 {
            lcom = v2 - c1;
        }
        strncpy_safe(&mut hstring[c1..], comment, lcom);
    }
}

#[cfg(test)]
mod tests {

    use std::ffi::CString;

    use super::*;

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_read_iraf_header() {
        let iraf_imh_file = "test_resources/pix.imh";

        let mut header_len = 0;
        let z = irafrdhead(
            CString::new(iraf_imh_file).unwrap().to_bytes_with_nul(),
            &mut header_len,
        );

        assert!(z.is_ok());
        assert_eq!(header_len, 5288);

        let z = z.unwrap();

        let header_version = head_version(&z);
        let pix_version = pix_version(&z);

        assert_eq!(header_version, 2);
        assert_eq!(pix_version, 2);
    }
}
