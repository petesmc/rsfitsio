use std::{io::Read, ptr};

use crate::c_types::{FILE, c_char, c_int, c_uint, c_ulong, c_void};
use libc::{feof, ferror, fread, free, fwrite, malloc};

use libz_rs_sys::{
    Z_BEST_SPEED, Z_BUF_ERROR, Z_DEFAULT_STRATEGY, Z_DEFLATED, Z_FINISH, Z_NO_FLUSH, Z_OK,
    Z_STREAM_END, Z_STREAM_ERROR, deflate, deflateEnd, deflateInit2_, inflate, inflateEnd,
    inflateInit2_, uInt, uLong, voidpf, z_stream, z_streamp, zlibVersion,
};

use crate::fitsio::{DATA_COMPRESSION_ERR, DATA_DECOMPRESSION_ERR, MEMORY_ALLOCATION};

const GZBUFSIZE: usize = 115200; /* 40 FITS blocks */
const BUFFINCR: usize = 28800; /* 10 FITS blocks */

pub(crate) unsafe fn inflateInit2(strm: z_streamp, windowBits: c_int) -> c_int {
    unsafe {
        inflateInit2_(
            strm,
            windowBits,
            zlibVersion(),
            std::mem::size_of::<z_stream>() as c_int,
        )
    }
}

pub(crate) unsafe fn deflateInit2(
    strm: z_streamp,
    level: c_int,
    method: c_int,
    windowBits: c_int,
    memLevel: c_int,
    strategy: c_int,
) -> c_int {
    unsafe {
        deflateInit2_(
            strm,
            level,
            method,
            windowBits,
            memLevel,
            strategy,
            zlibVersion(),
            std::mem::size_of::<z_stream>() as c_int,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Uncompress the disk file into memory.  Fill whatever amount of memory has
/// already been allocated, then realloc more memory, using the supplied
/// input function, if necessary.
pub(crate) unsafe fn uncompress2mem<T: Read>(
    filename: &[c_char],   /* name of input file                 */
    diskfile: &mut T,      /* I - file pointer                        */
    buffptr: *mut *mut u8, /* IO - memory pointer                     */
    buffsize: &mut usize,  /* IO - size of buffer, in bytes           */
    mem_realloc: Option<unsafe extern "C" fn(p: *mut c_void, newsize: usize) -> *mut c_void>, /* function     */
    filesize: &mut usize, /* O - size of file, in bytes              */
    status: &mut c_int,   /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut err: c_int = 0;
        let mut d_stream: z_stream; /* decompression stream */
        /* Input args buffptr and buffsize may refer to a block of memory
        larger than the 2^32 4 byte limit.  If so, must be broken
        up into "pages" when assigned to d_stream.
        (d_stream.avail_out is a uInt type, which might be smaller
        than buffsize's size_t type.)
        */
        let nPages: libz_rs_sys::uLong = (*buffsize as uLong) / (c_uint::MAX as uLong);
        let mut iPage: uLong = 0;
        let outbuffsize: uInt = if nPages > 0 {
            c_uint::MAX
        } else {
            *buffsize as uInt
        };

        if *status > 0 {
            return *status;
        }

        /* Allocate memory to hold compressed bytes read from the file. */
        let mut filebuff = Vec::new();
        if filebuff.try_reserve_exact(GZBUFSIZE).is_err() {
            *status = MEMORY_ALLOCATION;
            return *status; /* memory error */
        } else {
            filebuff.resize(GZBUFSIZE, 0);
        }

        d_stream = z_stream {
            next_in: ptr::null_mut(),
            avail_in: Default::default(),
            total_in: Default::default(),
            next_out: *buffptr,
            avail_out: outbuffsize,
            total_out: Default::default(),
            msg: ptr::null_mut(),
            state: ptr::null_mut(),
            zalloc: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, u32, u32) -> *mut c_void>,
            >(std::ptr::null_mut()),
            zfree: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, *mut c_void)>,
            >(std::ptr::null_mut()),
            opaque: ptr::null_mut() as voidpf,
            data_type: Default::default(),
            adler: Default::default(),
            reserved: Default::default(),
        };

        /* Initialize the decompression.  The argument (15+16) tells the
        decompressor that we are to use the gzip algorithm */

        err = inflateInit2(&mut d_stream, 15 + 16);
        if err != Z_OK {
            *status = DATA_DECOMPRESSION_ERR;
            return *status;
        }

        /* loop through the file, reading a buffer and uncompressing it */
        loop {
            let len = diskfile.read(&mut filebuff[..GZBUFSIZE]);
            if len.is_err() {
                inflateEnd(&mut d_stream);
                *status = DATA_DECOMPRESSION_ERR;
                return *status;
            }

            let len = len.unwrap();

            if len == 0 {
                break; /* no more data */
            }

            d_stream.next_in = filebuff.as_mut_ptr();
            d_stream.avail_in = len as uInt;

            loop {
                /* uncompress as much of the input as will fit in the output */
                err = inflate(&mut d_stream, Z_NO_FLUSH);

                if err == Z_STREAM_END {
                    /* We reached the end of the input */
                    break;
                } else if err == Z_OK || err == Z_BUF_ERROR {
                    /* Z_BUF_ERROR means need more input data to make progress */
                    if d_stream.avail_in == 0 {
                        break; /* need more input */
                    }

                    /* need more space in output buffer */
                    /* First check if more memory is available above the 4Gb limit in the originally input buffptr array */
                    if iPage < nPages {
                        iPage += 1;
                        d_stream.next_out =
                            (*buffptr).add((iPage * (c_uint::MAX as uLong)) as usize);
                        if iPage < nPages {
                            d_stream.avail_out = c_uint::MAX;
                        } else {
                            d_stream.avail_out =
                                ((*buffsize as uLong) % (c_uint::MAX as uLong)) as uInt;
                        }
                    } else if let Some(mem_realloc) = mem_realloc {
                        *buffptr =
                            mem_realloc(*buffptr as *mut c_void, *buffsize + BUFFINCR) as *mut _;
                        if (*buffptr).is_null() {
                            inflateEnd(&mut d_stream);
                            *status = DATA_DECOMPRESSION_ERR;
                            return *status; /* memory allocation failed */
                        }

                        d_stream.avail_out = BUFFINCR as uInt;
                        d_stream.next_out = (*buffptr).add(*buffsize);
                        *buffsize += BUFFINCR;
                    } else {
                        /* error: no realloc function available */
                        inflateEnd(&mut d_stream);
                        *status = DATA_DECOMPRESSION_ERR;
                        return *status;
                    }
                } else {
                    /* some other error */
                    inflateEnd(&mut d_stream);
                    *status = DATA_DECOMPRESSION_ERR;
                    return *status;
                }
            }

            /*
            // Should not happen since break occurs higher
            if feof(diskfile) {
                break;
            }
            */

            /*
            These settings for next_out and avail_out appear to be redundant,
            as the inflate() function should already be re-setting these.
            For case where *buffsize < 4Gb this did not matter, but for
            > 4Gb it would produce the wrong value in the avail_out assignment.
            (C. Gordon Jul 2016)
            d_stream.next_out = (unsigned char*) (*buffptr + d_stream.total_out);
            d_stream.avail_out = *buffsize - d_stream.total_out;
            */
        }

        /* Set the output file size to be the total output data */
        *filesize = d_stream.total_out as usize;

        err = inflateEnd(&mut d_stream); /* End the decompression */
        if err != Z_OK {
            *status = DATA_DECOMPRESSION_ERR;
            return *status;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
///Uncompress the file in memory into memory.  Fill whatever amount of memory has
///already been allocated, then realloc more memory, using the supplied
///input function, if necessary.
pub(crate) unsafe fn uncompress2mem_from_mem(
    inmemptr: &[c_char],   /* I - memory pointer to compressed bytes */
    inmemsize: usize,      /* I - size of input compressed file      */
    buffptr: *mut *mut u8, /* IO - memory pointer                      */
    buffsize: &mut usize,  /* IO - size of buffer, in bytes           */
    mem_realloc: Option<unsafe extern "C" fn(p: *mut c_void, newsize: usize) -> *mut c_void>, /* function     */
    filesize: Option<&mut usize>, /* O - size of file, in bytes              */
    status: &mut c_int,           /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut err: c_int = 0;
        let mut d_stream: z_stream; /* decompression stream */

        if *status > 0 {
            return *status;
        }

        d_stream = z_stream {
            next_in: ptr::null_mut(),
            avail_in: Default::default(),
            total_in: Default::default(),
            next_out: ptr::null_mut(),
            avail_out: Default::default(),
            total_out: Default::default(),
            msg: ptr::null_mut(),
            state: ptr::null_mut(),
            zalloc: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, u32, u32) -> *mut c_void>,
            >(std::ptr::null_mut()),
            zfree: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, *mut c_void)>,
            >(std::ptr::null_mut()),
            opaque: ptr::null_mut() as voidpf,
            data_type: Default::default(),
            adler: Default::default(),
            reserved: Default::default(),
        };

        /* Initialize the decompression.  The argument (15+16) tells the
        decompressor that we are to use the gzip algorithm */
        err = inflateInit2(&mut d_stream, 15 + 16);
        if err != Z_OK {
            *status = DATA_DECOMPRESSION_ERR;
            return *status;
        }

        d_stream.next_in = inmemptr.as_ptr() as *mut u8; // Yes convert from const to mut
        d_stream.avail_in = inmemsize as uInt;

        d_stream.next_out = *buffptr;
        d_stream.avail_out = *buffsize as uInt;

        loop {
            /* uncompress as much of the input as will fit in the output */
            err = inflate(&mut d_stream, Z_NO_FLUSH);

            if err == Z_STREAM_END {
                /* We reached the end of the input */
                break;
            } else if err == Z_OK || err == Z_BUF_ERROR {
                /* need more space in output buffer */
                /* Z_BUF_ERROR means need more input data to make progress */

                if let Some(mem_realloc) = mem_realloc {
                    *buffptr = mem_realloc(*buffptr as *mut c_void, *buffsize + BUFFINCR) as *mut _;
                    if (*buffptr).is_null() {
                        inflateEnd(&mut d_stream);
                        *status = DATA_DECOMPRESSION_ERR;
                        return *status; /* memory allocation failed */
                    }

                    d_stream.avail_out = BUFFINCR as uInt;
                    d_stream.next_out = (*buffptr).add(*buffsize);
                    *buffsize += BUFFINCR;
                } else {
                    /* error: no realloc function available */
                    inflateEnd(&mut d_stream);
                    *status = DATA_DECOMPRESSION_ERR;
                    return *status;
                }
            } else {
                /* some other error */
                inflateEnd(&mut d_stream);
                *status = DATA_DECOMPRESSION_ERR;
                return *status;
            }
        }

        /* Set the output file size to be the total output data */
        if let Some(filesize) = filesize {
            *filesize = d_stream.total_out as usize;
        }

        /* End the decompression */
        err = inflateEnd(&mut d_stream);

        if err != Z_OK {
            *status = DATA_DECOMPRESSION_ERR;
            return *status;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Uncompress the file into another file.
pub(crate) unsafe fn uncompress2file(
    filename: &[c_char],    /* name of input file                  */
    indiskfile: &mut FILE,  /* I - input file pointer                */
    outdiskfile: &mut FILE, /* I - output file pointer               */
    status: &mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut err: c_int = 0;

        let mut bytes_out: c_ulong = 0;
        //char *infilebuff, *outfilebuff;
        let mut d_stream: z_stream; /* decompression stream */

        if *status > 0 {
            return *status;
        }

        /* Allocate buffers to hold compressed and uncompressed */
        let infilebuff = malloc(GZBUFSIZE) as *mut u8;
        if infilebuff.is_null() {
            *status = MEMORY_ALLOCATION;
            return *status; /* memory error */
        }

        let outfilebuff = malloc(GZBUFSIZE) as *mut u8;
        if outfilebuff.is_null() {
            *status = MEMORY_ALLOCATION;
            return *status; /* memory error */
        }

        d_stream = z_stream {
            next_in: ptr::null_mut(),
            avail_in: Default::default(),
            total_in: Default::default(),
            next_out: outfilebuff,
            avail_out: GZBUFSIZE as uInt,
            total_out: Default::default(),
            msg: ptr::null_mut(),
            state: ptr::null_mut(),
            zalloc: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, u32, u32) -> *mut c_void>,
            >(std::ptr::null_mut()),
            zfree: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, *mut c_void)>,
            >(std::ptr::null_mut()),
            opaque: ptr::null_mut() as voidpf,
            data_type: Default::default(),
            adler: Default::default(),
            reserved: Default::default(),
        };

        /* Initialize the decompression.  The argument (15+16) tells the
        decompressor that we are to use the gzip algorithm */

        err = inflateInit2(&mut d_stream, 15 + 16);
        if err != Z_OK {
            *status = DATA_DECOMPRESSION_ERR;
            return *status;
        }

        /* loop through the file, reading a buffer and uncompressing it */
        loop {
            let len = fread(infilebuff as *mut c_void, 1, GZBUFSIZE, indiskfile);
            if ferror(indiskfile) != 0 {
                inflateEnd(&mut d_stream);
                free(infilebuff as *mut c_void);
                free(outfilebuff as *mut c_void);
                *status = DATA_DECOMPRESSION_ERR;
                return *status;
            }

            if len == 0 {
                break; /* no more data */
            }

            d_stream.next_in = infilebuff;
            d_stream.avail_in = len as uInt;

            loop {
                /* uncompress as much of the input as will fit in the output */
                err = inflate(&mut d_stream, Z_NO_FLUSH);

                if err == Z_STREAM_END {
                    /* We reached the end of the input */
                    break;
                } else if err == Z_OK || err == Z_BUF_ERROR {
                    /* Z_BUF_ERROR means need more input data to make progress */

                    if d_stream.avail_in == 0 {
                        break; /* need more input */
                    }

                    /* flush out the full output buffer */
                    if fwrite(outfilebuff as *const c_void, 1, GZBUFSIZE, outdiskfile) != GZBUFSIZE
                    {
                        inflateEnd(&mut d_stream);
                        free(infilebuff as *mut c_void);
                        free(outfilebuff as *mut c_void);
                        *status = DATA_DECOMPRESSION_ERR;
                        return *status;
                    }
                    bytes_out += GZBUFSIZE as c_ulong;
                    d_stream.next_out = outfilebuff;
                    d_stream.avail_out = GZBUFSIZE as _;
                } else {
                    /* some other error */
                    inflateEnd(&mut d_stream);
                    free(infilebuff as *mut c_void);
                    free(outfilebuff as *mut c_void);
                    *status = DATA_DECOMPRESSION_ERR;
                    return *status;
                }
            }

            if feof(indiskfile) != 0 {
                break;
            }
        }

        /* write out any remaining bytes in the buffer */
        if d_stream.total_out > bytes_out
            && fwrite(
                outfilebuff as *const c_void,
                1,
                (d_stream.total_out - bytes_out) as usize,
                outdiskfile,
            ) != (d_stream.total_out - bytes_out) as usize
        {
            inflateEnd(&mut d_stream);
            free(infilebuff as *mut c_void);
            free(outfilebuff as *mut c_void);
            *status = DATA_DECOMPRESSION_ERR;
            return *status;
        }

        free(infilebuff as *mut c_void); /* free temporary output data buffer */
        free(outfilebuff as *mut c_void); /* free temporary output data buffer */

        err = inflateEnd(&mut d_stream); /* End the decompression */
        if err != Z_OK {
            *status = DATA_DECOMPRESSION_ERR;
            return *status;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Compress the file into memory.  Fill whatever amount of memory has
/// already been allocated, then realloc more memory, using the supplied
/// input function, if necessary.
pub(crate) unsafe fn compress2mem_from_mem(
    inmemptr: &[c_char],   /* I - memory pointer to uncompressed bytes */
    inmemsize: usize,      /* I - size of input uncompressed file      */
    buffptr: *mut *mut u8, /* IO - memory pointer for compressed file    */
    buffsize: &mut usize,  /* IO - size of buffer, in bytes           */
    mem_realloc: Option<unsafe extern "C" fn(p: *mut c_void, newsize: usize) -> *mut c_void>, /* function     */
    filesize: Option<&mut usize>, /* O - size of file, in bytes              */
    status: &mut c_int,           /* IO - error status                       */
) -> c_int {
    unsafe {
        let mut err: c_int = 0;
        let mut c_stream: z_stream; /* compression stream */

        if *status > 0 {
            return *status;
        }

        c_stream = z_stream {
            next_in: ptr::null_mut(),
            avail_in: Default::default(),
            total_in: Default::default(),
            next_out: ptr::null_mut(),
            avail_out: Default::default(),
            total_out: Default::default(),
            msg: ptr::null_mut(),
            state: ptr::null_mut(),
            zalloc: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, u32, u32) -> *mut c_void>,
            >(std::ptr::null_mut()),
            zfree: std::mem::transmute::<
                *mut u32,
                Option<unsafe extern "C" fn(*mut c_void, *mut c_void)>,
            >(std::ptr::null_mut()),
            opaque: ptr::null_mut() as voidpf,
            data_type: Default::default(),
            adler: Default::default(),
            reserved: Default::default(),
        };

        /* Initialize the compression.  The argument (15+16) tells the
        compressor that we are to use the gzip algorythm.
        Also use Z_BEST_SPEED for maximum speed with very minor loss
        in compression factor. */
        err = deflateInit2(
            &mut c_stream,
            Z_BEST_SPEED,
            Z_DEFLATED,
            15 + 16,
            8,
            Z_DEFAULT_STRATEGY,
        );

        if err != Z_OK {
            *status = DATA_COMPRESSION_ERR;
            return *status;
        }

        c_stream.next_in = inmemptr.as_ptr() as *mut u8; // Yes convert const to mut
        c_stream.avail_in = inmemsize as uInt;

        c_stream.next_out = *buffptr;
        c_stream.avail_out = *buffsize as uInt;

        loop {
            /* compress as much of the input as will fit in the output */
            err = deflate(&mut c_stream, Z_FINISH);

            if err == Z_STREAM_END {
                /* We reached the end of the input */
                break;
            } else if err == Z_OK {
                /* need more space in output buffer */

                if let Some(mem_realloc) = mem_realloc {
                    *buffptr = mem_realloc(*buffptr as *mut c_void, *buffsize + BUFFINCR) as *mut _;
                    if (*buffptr).is_null() {
                        deflateEnd(&mut c_stream);
                        *status = DATA_COMPRESSION_ERR;
                        return *status; /* memory allocation failed */
                    }

                    c_stream.avail_out = BUFFINCR as uInt;
                    c_stream.next_out = (*buffptr).add(*buffsize);
                    *buffsize += BUFFINCR;
                } else {
                    /* error: no realloc function available */
                    deflateEnd(&mut c_stream);
                    *status = DATA_COMPRESSION_ERR;
                    return *status;
                }
            } else {
                /* some other error */
                deflateEnd(&mut c_stream);
                *status = DATA_COMPRESSION_ERR;
                return *status;
            }
        }

        /* Set the output file size to be the total output data */
        if let Some(filesize) = filesize {
            *filesize = c_stream.total_out as usize;
        }

        /* End the compression */
        err = deflateEnd(&mut c_stream);

        if err != Z_OK {
            *status = DATA_COMPRESSION_ERR;
            return *status;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Compress the memory file into disk file.
pub(crate) unsafe fn compress2file_from_mem(
    inmemptr: &[c_char], /* I - memory pointer to uncompressed bytes */
    inmemsize: usize,    /* I - size of input uncompressed file      */
    outdiskfile: &mut FILE,
    filesize: Option<&mut usize>, /* O - size of file, in bytes              */
    status: &mut c_int,
) -> c_int {
    unsafe {
        let mut err: c_int = 0;
        let mut flushflag: c_int = 0;

        let mut nPages: uLong = 1;
        let mut nBytesToFile: uInt = 0;

        let mut c_stream: z_stream; /* compression stream */

        if *status > 0 {
            return *status;
        }

        /* Allocate buffer to hold compressed bytes */
        let outfilebuff = malloc(GZBUFSIZE) as *mut u8;
        if outfilebuff.is_null() {
            *status = MEMORY_ALLOCATION;
            return *status; /* memory error */
        }

        c_stream = z_stream {
            next_in: ptr::null_mut(),
            avail_in: Default::default(),
            total_in: Default::default(),
            next_out: ptr::null_mut(),
            avail_out: Default::default(),
            total_out: Default::default(),
            msg: ptr::null_mut(),
            state: ptr::null_mut(),
            zalloc: std::mem::transmute::<
                *mut u32,
                std::option::Option<unsafe extern "C" fn(*mut c_void, u32, u32) -> *mut c_void>,
            >(std::ptr::null_mut()),
            zfree: std::mem::transmute::<
                *mut u32,
                std::option::Option<unsafe extern "C" fn(*mut c_void, *mut c_void)>,
            >(std::ptr::null_mut()),
            opaque: ptr::null_mut() as voidpf,
            data_type: Default::default(),
            adler: Default::default(),
            reserved: Default::default(),
        };

        /* Initialize the compression.  The argument (15+16) tells the
        compressor that we are to use the gzip algorythm.
        Also use Z_BEST_SPEED for maximum speed with very minor loss
        in compression factor. */
        err = deflateInit2(
            &mut c_stream,
            Z_BEST_SPEED,
            Z_DEFLATED,
            15 + 16,
            8,
            Z_DEFAULT_STRATEGY,
        );

        if err != Z_OK {
            *status = DATA_COMPRESSION_ERR;
            return *status;
        }

        // ####
        if inmemsize > 0 {
            nPages = 1 + (inmemsize as uLong - 1) / (uInt::MAX as uLong);
        }

        /*
        c_stream.next_in = inmemptr.as_ptr() as *mut u8; // Yes convert const to mut
        c_stream.avail_in = inmemsize as uInt;

        c_stream.next_out = outfilebuff;
        c_stream.avail_out = GZBUFSIZE as uInt;

        loop {
            /* compress as much of the input as will fit in the output */
            err = deflate(&mut c_stream, Z_FINISH);

            if err == Z_STREAM_END {
                /* We reached the end of the input */
                break;
            } else if err == Z_OK {
                /* need more space in output buffer */

                /* flush out the full output buffer */
                if fwrite(outfilebuff as *const c_void, 1, GZBUFSIZE, outdiskfile) != GZBUFSIZE {
                    deflateEnd(&mut c_stream);
                    free(outfilebuff as *mut c_void);
                    *status = DATA_COMPRESSION_ERR;
                    return *status;
                }
                bytes_out += GZBUFSIZE as c_ulong;
                c_stream.next_out = outfilebuff;
                c_stream.avail_out = GZBUFSIZE as uInt;
            } else {
                /* some other error */
                deflateEnd(&mut c_stream);
                free(outfilebuff as *mut c_void);
                *status = DATA_COMPRESSION_ERR;
                return *status;
            }
        }
        */

        let iPage: uLong = 0;
        for iPage in 0..nPages {
            // SAFETY: Converted a const pointer to a mutable pointer, don't know why it needs to be.
            c_stream.next_in =
                (inmemptr.as_ptr() as *mut u8).add((iPage * uInt::MAX as uLong) as usize);
            c_stream.avail_in = if iPage == nPages - 1 {
                (inmemsize as uLong - iPage * uInt::MAX as uLong) as uInt
            } else {
                uInt::MAX
            };

            flushflag = if iPage < nPages - 1 {
                Z_NO_FLUSH
            } else {
                Z_FINISH
            };
            loop {
                c_stream.next_out = outfilebuff;
                c_stream.avail_out = GZBUFSIZE as uInt;

                /* compress as much of the input as will fit in the output */
                err = deflate(&mut c_stream, flushflag);

                if err == Z_STREAM_ERROR {
                    deflateEnd(&mut c_stream);
                    free(outfilebuff as *mut c_void);
                    *status = DATA_COMPRESSION_ERR;
                    return *status;
                } else {
                    /* c_stream.avail_out will be 0 unless we've reached the end of the avail_in
                    stream.  When that happens avail_out MAY also be 0, if by chance the output
                    buffer fills up just as the input stream ends.  That's OK though, as it will
                    execute just one more do/while where the deflate call won't actually do
                    anything.  */
                    nBytesToFile = GZBUFSIZE as uInt - c_stream.avail_out;
                    if nBytesToFile != 0
                        && fwrite(
                            outfilebuff as *const c_void,
                            1,
                            nBytesToFile as usize,
                            outdiskfile,
                        ) != nBytesToFile as usize
                    {
                        deflateEnd(&mut c_stream);
                        free(outfilebuff as *mut c_void);
                        *status = DATA_COMPRESSION_ERR;
                        return *status;
                    }
                }
                if c_stream.avail_out != 0 {
                    break;
                }
            }
        }

        free(outfilebuff as *mut c_void); /* free temporary output data buffer */

        /* Set the output file size to be the total output data */
        if let Some(filesize) = filesize {
            *filesize = c_stream.total_out as usize;
        }

        /* End the compression */
        err = deflateEnd(&mut c_stream);

        if err != Z_OK {
            *status = DATA_COMPRESSION_ERR;
            return *status;
        }

        *status
    }
}
