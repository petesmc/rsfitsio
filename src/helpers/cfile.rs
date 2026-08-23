//! Adapters between C stream I/O and the Rust `Read` / `Write` traits.
#![warn(missing_docs)]

use std::{
    fs::File,
    io::{self, Read, Write},
};

use libc::{FILE, c_int, ferror, fflush, fread};

pub(crate) struct CFile {
    file: *mut FILE,
}

// Implement the trait to convert from a raw pointer to CFile
impl From<*mut FILE> for CFile {
    fn from(file: *mut FILE) -> Self {
        CFile { file }
    }
}

impl CFile {
    fn has_error(&self) -> bool {
        unsafe { ferror(self.file) != 0 }
    }

    fn flush_contents(&self) -> c_int {
        unsafe { fflush(self.file) }
    }
}

impl Read for CFile {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let nread = unsafe {
            fread(
                buf.as_mut_ptr().cast::<libc::c_void>(),
                1,
                buf.len(),
                self.file,
            ) as usize
        };
        if self.has_error() {
            Err(io::Error::last_os_error())
        } else {
            Ok(nread)
        }
    }
}

impl Write for CFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let nwritten = unsafe {
            libc::fwrite(buf.as_ptr().cast::<libc::c_void>(), 1, buf.len(), self.file) as usize
        };
        if self.has_error() {
            Err(io::Error::last_os_error())
        } else {
            Ok(nwritten)
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_contents();
        if self.has_error() {
            Err(io::Error::last_os_error())
        } else {
            Ok(())
        }
    }
}

// `std::io::ErrorKind` is retained deliberately: the `core::io` equivalent is still unstable.
#[allow(clippy::std_instead_of_core)]
pub(crate) fn fgets(buf: &mut [u8], size: usize, file: &mut File) -> Result<(), std::io::Error> {
    let mut byte = [0; 1];

    buf[buf.len() - 1] = 0; // Null-terminate the string

    for i in 0..(size - 1) {
        let n = file.read(&mut byte)?;

        if n == 0 {
            /* EOF. C's fgets returns NULL only when end-of-file is reached and
            no characters were read; if some chars were already read it returns
            them. Mirror that so `while fgets(...).is_ok()` loops terminate. */
            buf[i] = 0;
            if i == 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "fgets: end of file",
                ));
            }
            return Ok(());
        }

        buf[i] = byte[0];

        if byte[0] == 0 || byte[0] == (b'\n') {
            return Ok(());
        }
    }

    Ok(())
}
