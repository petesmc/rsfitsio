//! This is an old implementation of the LZW decompression algorithm from the
//! Unix compress utility.

use bytemuck::{cast_slice, cast_slice_mut};

use crate::c_types::{FILE, c_char, c_int, c_long, c_uchar, c_uint, c_ulong, c_ushort, c_void};
use libc::{EOF, fread, fwrite, memcpy};

use crate::{
    aliases::{ffpmsg_slice, ffpmsg_str},
    fitsio::DATA_DECOMPRESSION_ERR,
    wrappers::strncat_safe,
};

/* Return codes from gzip */
const OK: c_int = 0;
const ERROR: c_int = 1;
const COMPRESSED: c_int = 1;
const DEFLATED: c_int = 8;
const INBUFSIZ: usize = 0x8000; /* input buffer size */
const INBUF_EXTRA: usize = 64; /* required by unlzw() */
const OUTBUFSIZ: usize = 16384; /* output buffer size */
const OUTBUF_EXTRA: usize = 2048; /* required by unlzw() */
const DIST_BUFSIZE: usize = 0x8000; /* buffer for distances, see trees.c */
const WSIZE: usize = 0x8000; /* window size--must be a power of two, and */

const LZW_MAGIC: [u8; 2] = [0x1f, 0x9d]; /* magic header for LZW files */

const BITS: c_int = 16; /* code bits */
const INIT_BITS: c_int = 9; /* initial number of bits/code */
const BIT_MASK: c_int = 0x1f; /* mask for 'n_bits' bits */
const BLOCK_MODE: c_int = 0x80;
const LZW_RESERVED: c_int = 0x60; /* reserved bits */
const CLEAR: c_int = 256; /* table clear output code*/
const FIRST: c_int = CLEAR + 1; /* first free entry */

type char_type = c_uchar;
type code_int = c_long;
type count_int = c_ulong;
type count_short = c_ushort;
type cmp_code_int = c_ulong;

macro_rules! MAXCODE {
    ($n:expr) => {
        (1 << $n)
    };
}

struct LZW_Compress<'a> {
    inbuf: [c_uchar; INBUFSIZ + INBUF_EXTRA],
    outbuf: [c_uchar; OUTBUFSIZ + OUTBUF_EXTRA],
    d_buf: [c_ushort; DIST_BUFSIZE],
    window: [c_uchar; 2 * WSIZE],
    tab_prefix: [c_ushort; 1 << BITS], /* prefix code for each entry */
    maxbits: c_int,                    /* max bits per code for LZW */
    method: c_int,                     /* compression method */
    exit_code: c_int,                  /* program exit code */
    last_member: c_int,                /* set for .zip and .Z files */
    bytes_in: usize,                   /* number of input bytes */
    bytes_out: usize,                  /* number of output bytes */
    ifname: [c_char; 128],             /* input file name */
    ifd: *mut FILE,                    /* input file descriptor */
    ofd: *mut FILE,                    /* output file descriptor */
    memptr: *mut *mut c_void,          /* memory location for uncompressed file */
    memsize: &'a mut usize,            /* size (bytes) of memory allocated for file */
    realloc_fn: Option<unsafe extern "C" fn(*mut c_void, usize) -> *mut c_void>, /* reallocation function */
    insize: usize,     /* valid bytes in inbuf */
    inptr: usize,      /* index of next byte to be processed in inbuf */
    block_mode: c_int, /* block compress mode -C compatible with 2.0 */
}

impl<'a> LZW_Compress<'a> {
    fn get_byte(&mut self) -> u8 {
        if self.inptr < self.insize {
            let v = self.inbuf[self.inptr];
            self.inptr += 1;
            v
        } else {
            self.fill_inbuf(0) as u8
        }
    }

    /* =========================================================================== */
    /// set if EOF acceptable as a result
    fn fill_inbuf(&mut self, eof_ok: c_int) -> c_int {
        let mut len: usize = 0;

        /* Read as much as possible from file */
        self.insize = 0;
        loop {
            len = unsafe {
                fread(
                    self.inbuf[self.insize..].as_mut_ptr() as *mut c_void,
                    1,
                    INBUFSIZ - self.insize,
                    self.ifd,
                )
            };
            if len == 0 || len == (EOF as usize) {
                break;
            }
            self.insize += len;

            if self.insize >= INBUFSIZ {
                break;
            }
        }

        if self.insize == 0 {
            if eof_ok != 0 {
                return EOF;
            }
            self.error("unexpected end of file");
            self.exit_code = ERROR;
            return ERROR;
        }

        self.bytes_in += self.insize;
        self.inptr = 1;
        self.inbuf[0] as c_int
    }

    /* =========================================================================== */
    ///  copy buffer into memory; allocate more memory if required
    unsafe fn write_buf(&mut self, cnt: usize) {
        unsafe {
            let buf = &mut self.outbuf[..];
            if self.realloc_fn.is_none() {
                /* append buffer to file */
                /* added 'unsigned' to get rid of compiler warning (WDP 1/1/99) */
                if fwrite(buf.as_ptr() as *mut c_void, 1, cnt, self.ofd) != cnt {
                    self.error("failed to write buffer to uncompressed output file (write_buf)");
                    self.exit_code = ERROR;
                }
            } else {
                /* get more memory if current buffer is too small */
                if self.bytes_out + cnt > *self.memsize {
                    *self.memptr = self.realloc_fn.unwrap()(*self.memptr, self.bytes_out + cnt);
                    *self.memsize = self.bytes_out + cnt; /* new memory buffer size */

                    if self.memptr.is_null() {
                        self.error("malloc failed while uncompressing (write_buf)");
                        self.exit_code = ERROR;
                        return;
                    }
                }
                /* copy  into memory buffer */
                memcpy(
                    (*self.memptr).add(self.bytes_out),
                    buf.as_ptr() as *mut c_void,
                    cnt,
                );
            }
        }
    }

    /* ======================================================================== */
    /// Error handler
    fn error(&self, m: &str) {
        ffpmsg_slice(&self.ifname);
        ffpmsg_str(m);
    }
}

/*--------------------------------------------------------------------------*/
/// Uncompress the file into memory.  Fill whatever amount of memory has
/// already been allocated, then realloc more memory, using the supplied
/// input function, if necessary.
pub(crate) unsafe fn zuncompress2mem(
    filename: &[c_char],   /* name of input file                 */
    indiskfile: *mut FILE, /* I - file pointer                        */
    buffptr: *mut *mut u8, /* IO - memory pointer                     */
    buffsize: &mut usize,  /* IO - size of buffer, in bytes           */
    mem_realloc: Option<unsafe extern "C" fn(p: *mut c_void, newsize: usize) -> *mut c_void>, /* function     */
    filesize: &mut usize, /* O - size of file, in bytes              */
    status: &mut c_int,   /* IO - error status                       */
) -> c_int {
    let mut magic: [c_char; 2] = [0; 2]; /* magic header */

    if *status != 0 {
        return *status;
    }

    let mut fn_buffer: [c_char; 128] = [0; 128]; /* buffer for file name */
    strncat_safe(&mut fn_buffer[..128], filename, 127);

    let mut lzw: LZW_Compress = LZW_Compress {
        inbuf: [0; INBUFSIZ + INBUF_EXTRA],
        outbuf: [0; OUTBUFSIZ + OUTBUF_EXTRA],
        d_buf: [0; DIST_BUFSIZE],
        window: [0; 2 * WSIZE],
        tab_prefix: [0; 1 << BITS],
        maxbits: 0,
        method: COMPRESSED,
        exit_code: 0,
        last_member: 1,
        bytes_in: 0,
        bytes_out: 0,
        ifname: fn_buffer,
        ifd: indiskfile,
        ofd: std::ptr::null_mut(),
        memptr: buffptr as *mut *mut c_void,
        memsize: buffsize,
        realloc_fn: mem_realloc,
        insize: 0,
        inptr: 0,
        block_mode: BLOCK_MODE,
    };

    magic[0] = lzw.get_byte() as c_char;
    magic[1] = lzw.get_byte() as c_char;

    if magic[0] != LZW_MAGIC[0] as c_char || magic[1] != LZW_MAGIC[1] as c_char {
        lzw.error("ERROR: input .Z file is in unrecognized compression format.\n");
        *status = ERROR;
        return ERROR;
    }

    if unlzw(&mut lzw, indiskfile, std::ptr::null_mut()) != OK {
        *status = DATA_DECOMPRESSION_ERR;
    }

    *filesize = lzw.bytes_out;

    *status
}

/*--------------------------------------------------------------------------*/
/// Decompress in to out.  This routine adapts to the codes in the
/// file building the "string" table on-the-fly; requiring no table to
/// be stored in the compressed file.
/// IN assertions: the buffer inbuf contains already the beginning of
/// the compressed data, from offsets iptr to insize-1 included.
/// The magic header has already been checked and skipped.
/// bytes_in and bytes_out have been initialized.
fn unlzw(lzw: &mut LZW_Compress, in_file: *mut FILE, out_file: *mut FILE) -> c_int {
    let mut inbits: c_long = 0;
    let mut code: code_int;
    let mut incode: code_int;
    let mut stackp: usize;

    lzw.ofd = out_file;
    lzw.maxbits = lzw.get_byte() as c_int;
    lzw.block_mode = lzw.maxbits & BLOCK_MODE;

    if lzw.maxbits & LZW_RESERVED != 0 {
        lzw.error("warning, unknown flags in unlzw decompression");
    }

    lzw.maxbits &= BIT_MASK;

    let maxmaxcode: code_int = MAXCODE!(lzw.maxbits);

    if lzw.maxbits > BITS {
        lzw.error("compressed with too many bits; cannot handle file");
        lzw.exit_code = ERROR;
        return ERROR;
    }

    let mut rsize: c_int = lzw.insize as c_int;
    let mut n_bits: c_int = INIT_BITS;
    let mut maxcode: code_int = MAXCODE!(n_bits) - 1;
    let mut bitmask: c_uint = (1 << n_bits) - 1;
    let mut oldcode: code_int = -1;
    let mut finchar: c_int = 0;
    let mut outpos: c_int = 0;
    let mut posbits: c_long = (lzw.inptr as c_long) << 3;

    let mut free_ent: code_int = if lzw.block_mode != 0 {
        FIRST as code_int
    } else {
        256
    };

    /* Initialize the first 256 entries in the table. */
    lzw.tab_prefix[..256].fill(0);

    for code in (0..256).rev() {
        lzw.window[code] = code as c_uchar;
    }

    let mut resetbuf = false;

    let mut i: c_int = 0;
    let mut e: c_int = 0;
    let mut o: c_int = 0;

    'outer_loop: loop {
        if !resetbuf {
            i = 0;
            e = 0;
            o = 0;
        }

        o = (posbits >> 3) as c_int;
        e = (lzw.insize - o as usize) as c_int;

        for i in 0..(e) {
            lzw.inbuf[i as usize] = lzw.inbuf[(o + i) as usize];
        }

        lzw.insize = e as usize;
        posbits = 0;

        if lzw.insize < INBUF_EXTRA {
            rsize = unsafe {
                fread(
                    lzw.inbuf[lzw.insize..].as_mut_ptr() as *mut c_void,
                    1,
                    INBUFSIZ,
                    in_file,
                ) as c_int
            };

            if rsize == EOF {
                lzw.error("unexpected end of file");
                lzw.exit_code = ERROR;
                return ERROR;
            }
            lzw.insize += rsize as usize;
            lzw.bytes_in += rsize as usize;
        }

        inbits = if rsize != 0 {
            ((lzw.insize as c_long) - (lzw.insize as c_long % n_bits as c_long)) << 3
        } else {
            ((lzw.insize as c_long) << 3) - (n_bits as c_long - 1)
        };

        while inbits > posbits {
            if free_ent > maxcode {
                posbits = (posbits - 1)
                    + (((n_bits as c_long) << 3)
                        - (posbits - 1 + (n_bits << 3) as c_long) % (n_bits << 3) as c_long);
                n_bits += 1;
                if n_bits == lzw.maxbits {
                    maxcode = maxmaxcode;
                } else {
                    maxcode = MAXCODE!(n_bits) - 1;
                }
                bitmask = (1 << n_bits) - 1;

                resetbuf = true;
                continue 'outer_loop;
            }

            // 'input' macro
            let p = &(lzw.inbuf)[((posbits >> 3) as usize)..];
            (code) = (((p[0] as c_long) | ((p[1] as c_long) << 8) | ((p[2] as c_long) << 16))
                >> ((posbits) & 0x7))
                & (bitmask as c_long);
            (posbits) += n_bits as c_long;

            if oldcode == -1 {
                if code >= 256 {
                    lzw.error("corrupt input.");
                    lzw.exit_code = ERROR;
                    return ERROR;
                }

                oldcode = code;
                finchar = code as c_int;
                lzw.outbuf[outpos as usize] = finchar as c_uchar;
                outpos += 1;

                continue;
            }

            if code == CLEAR as c_long && lzw.block_mode != 0 {
                lzw.tab_prefix[..256].fill(0);
                free_ent = (FIRST - 1) as c_long;
                posbits = (posbits - 1)
                    + (((n_bits as c_long) << 3)
                        - (posbits - 1 + (n_bits << 3) as c_long) % (n_bits << 3) as c_long);
                n_bits = INIT_BITS;
                maxcode = MAXCODE!(n_bits) - 1;
                bitmask = (1 << n_bits) - 1;

                resetbuf = true;
                continue 'outer_loop;
            }

            incode = code;

            let stack: &mut [char_type] = cast_slice_mut(&mut lzw.d_buf);
            stackp = (DIST_BUFSIZE - 1) * std::mem::size_of::<c_ushort>()
                / std::mem::size_of::<char_type>();

            if code >= free_ent {
                /* Special case for KwKwK string. */
                if code > free_ent {
                    if outpos > 0 {
                        unsafe {
                            lzw.write_buf(outpos as usize);
                        }
                        lzw.bytes_out += outpos as usize;
                    }
                    lzw.error("corrupt input.");
                    lzw.exit_code = ERROR;
                    return ERROR;
                }

                stackp -= 1;
                stack[stackp] = finchar as char_type;
                code = oldcode;
            }

            while (code as cmp_code_int) >= 256 {
                /* Generate output characters in reverse order */
                stackp -= 1;
                stack[stackp] = lzw.window[code as usize] as char_type;
                code = lzw.tab_prefix[code as usize] as code_int;
            }

            stackp -= 1;
            finchar = lzw.window[code as usize] as c_int;
            stack[stackp] = finchar as char_type;

            /* And put them out in forward order */
            {
                /*	REG1 int	i;   already defined above (WDP) */

                i = ((DIST_BUFSIZE - 1) * std::mem::size_of::<c_ushort>()
                    / std::mem::size_of::<char_type>()
                    - stackp) as c_int;

                if outpos + (i) >= OUTBUFSIZ as c_int {
                    loop {
                        if i > OUTBUFSIZ as c_int - outpos {
                            i = OUTBUFSIZ as c_int - outpos;
                        }

                        if i > 0 {
                            let s: &mut [c_char] = cast_slice_mut(&mut lzw.d_buf);

                            (&mut lzw.outbuf)[(outpos as usize)..((outpos + i) as usize)]
                                .copy_from_slice(cast_slice(&s[stackp..(stackp + i as usize)]));
                            outpos += i;
                        }

                        if outpos >= OUTBUFSIZ as c_int {
                            unsafe {
                                lzw.write_buf(outpos as usize);
                            }
                            lzw.bytes_out += outpos as usize;
                            outpos = 0;
                        }
                        stackp += i as usize;

                        i = ((DIST_BUFSIZE - 1) * std::mem::size_of::<c_ushort>()
                            / std::mem::size_of::<char_type>()
                            - stackp) as c_int;

                        if i <= 0 {
                            break;
                        }
                    }
                } else {
                    (&mut lzw.outbuf)[(outpos as usize)..((outpos + i) as usize)]
                        .copy_from_slice(cast_slice(&stack[stackp..(stackp + i as usize)]));
                    outpos += i;
                }
            }

            code = free_ent;
            if code < maxmaxcode {
                /* Generate the new entry. */

                lzw.tab_prefix[code as usize] = oldcode as c_ushort;
                lzw.window[code as usize] = finchar as char_type;
                free_ent = code + 1;
            }
            oldcode = incode; /* Remember previous code.	*/
        }

        if rsize == 0 {
            break;
        }
    }

    if outpos > 0 {
        unsafe { lzw.write_buf(outpos as usize) };
        lzw.bytes_out += outpos as usize;
    }

    OK
}

#[cfg(test)]
mod tests {

    use std::{ffi::CString, fs::File, io::Read};

    use crate::c_types::{c_char, c_int};
    use bytemuck::cast_slice;
    use libc::{fdopen, realloc};

    use crate::zuncompress::zuncompress2mem;

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_zuncompress2mem() {
        // Paths to the compressed and expected decompressed files
        let compressed_file_path = "test_resources/sample.txt.Z";
        let expected_file_path = "test_resources/sample.txt";

        // Open the compressed file
        let compressed_file =
            File::open(compressed_file_path).expect("Failed to open compressed file");
        let compressed_file_ptr = unsafe {
            #[cfg(windows)]
            {
                use libc::open_osfhandle;
                use std::os::windows::io::IntoRawHandle;

                let handle = compressed_file.into_raw_handle();
                let fd = open_osfhandle(handle as isize, 0);
                fdopen(fd, c"rb".as_ptr() as *const c_char)
            }

            #[cfg(unix)]
            {
                use std::os::fd::AsRawFd;

                let fd = compressed_file.as_raw_fd();

                fdopen(fd, c"rb".as_ptr() as *const c_char)
            }
        };

        // Prepare buffers and variables for decompression
        let mut decompressed_buffer: Vec<u8> = vec![0; 1024 * 1024 * 1024];
        let mut buffer_size: usize = 1024 * 1024 * 1024;
        let mut decompressed_size: usize = 0;
        let mut status: c_int = 0;

        // Call the zuncompress2mem function
        let result = unsafe {
            zuncompress2mem(
                cast_slice(
                    CString::new(compressed_file_path)
                        .unwrap()
                        .as_bytes_with_nul(),
                ),
                compressed_file_ptr,
                &mut decompressed_buffer.as_mut_ptr(),
                &mut buffer_size,
                Some(realloc),
                &mut decompressed_size,
                &mut status,
            )
        };

        // Ensure the decompression was successful
        assert_eq!(result, 0, "Decompression failed with status: {status}");

        // Read the expected decompressed content
        let mut expected_content = Vec::new();
        File::open(expected_file_path)
            .expect("Failed to open expected file")
            .read_to_end(&mut expected_content)
            .expect("Failed to read expected file");

        // Compare the decompressed content with the expected content
        assert_eq!(
            &decompressed_buffer[..decompressed_size],
            &expected_content,
            "Decompressed content does not match expected content"
        );

        // Clean up
        // unsafe {
        //libc::fclose(compressed_file_ptr);
        // }
    }
}
