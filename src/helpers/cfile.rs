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

pub(crate) fn fgets(buf: &mut [u8], size: usize, file: &mut File) -> Result<(), std::io::Error> {
    let mut byte = [0; 1];

    buf[buf.len() - 1] = 0; // Null-terminate the string

    for i in 0..(size - 1) {
        let read_result = file.read(&mut byte);

        if read_result.is_err() {
            return Err(read_result.unwrap_err());
        }

        let read_result = read_result.unwrap();

        buf[i] = byte[0];

        if byte[0] == 0 || byte[0] == (b'\n') {
            return Ok(());
        }
    }

    Ok(())
}
