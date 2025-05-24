use core::fmt;
use std::{cell::Cell, ptr};

use crate::c_types::c_int;

use super::io::{self, Write};

pub mod types;

thread_local! {
    /// The global `errno` variable used internally in relibc.
    pub(crate) static ERRNO: Cell<c_int> = const { Cell::new(0) };
}

pub(crate) trait WriteByte: fmt::Write {
    fn write_u8(&mut self, byte: u8) -> fmt::Result;
}

impl<W: WriteByte> WriteByte for &mut W {
    fn write_u8(&mut self, byte: u8) -> fmt::Result {
        (**self).write_u8(byte)
    }
}

pub(crate) struct SliceWriter<'a>(pub &'a mut [u8], pub usize);
impl<'a> Write for SliceWriter<'a> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let copy_size = buf.len().min(self.1);
        self.0[..copy_size].copy_from_slice(&buf[..copy_size]);
        self.1 -= copy_size;
        self.0[copy_size] = 0;
        Ok(copy_size)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl fmt::Write for SliceWriter<'_> {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        // can't fail
        self.write(s.as_bytes()).unwrap();
        Ok(())
    }
}
impl WriteByte for SliceWriter<'_> {
    fn write_u8(&mut self, byte: u8) -> fmt::Result {
        // can't fail
        self.write(&[byte]).unwrap();
        Ok(())
    }
}

pub(crate) struct StringWriter(pub *mut u8, pub usize);
impl Write for StringWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if self.1 > 1 {
            let copy_size = buf.len().min(self.1 - 1);
            unsafe {
                ptr::copy_nonoverlapping(buf.as_ptr(), self.0, copy_size);
                self.1 -= copy_size;

                self.0 = self.0.add(copy_size);
                *self.0 = 0;
            }
        }

        // Pretend the entire slice was written. This is because many functions
        // (like snprintf) expects a return value that reflects how many bytes
        // *would have* been written. So keeping track of this information is
        // good, and then if we want the *actual* written size we can just go
        // `cmp::min(written, maxlen)`.
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
impl fmt::Write for StringWriter {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        // can't fail
        self.write(s.as_bytes()).unwrap();
        Ok(())
    }
}
impl WriteByte for StringWriter {
    fn write_u8(&mut self, byte: u8) -> fmt::Result {
        // can't fail
        self.write(&[byte]).unwrap();
        Ok(())
    }
}

pub(crate) struct UnsafeStringWriter(pub *mut u8);
impl Write for UnsafeStringWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        unsafe {
            ptr::copy_nonoverlapping(buf.as_ptr(), self.0, buf.len());
            self.0 = self.0.add(buf.len());
            *self.0 = b'\0';
        }
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
impl fmt::Write for UnsafeStringWriter {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        // can't fail
        self.write(s.as_bytes()).unwrap();
        Ok(())
    }
}
impl WriteByte for UnsafeStringWriter {
    fn write_u8(&mut self, byte: u8) -> fmt::Result {
        // can't fail
        self.write(&[byte]).unwrap();
        Ok(())
    }
}

pub(crate) struct CountingWriter<T> {
    pub inner: T,
    pub written: usize,
}
impl<T> CountingWriter<T> {
    pub fn new(writer: T) -> Self {
        Self {
            inner: writer,
            written: 0,
        }
    }
}
impl<T: fmt::Write> fmt::Write for CountingWriter<T> {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.written += s.len();
        self.inner.write_str(s)
    }
}
impl<T: WriteByte> WriteByte for CountingWriter<T> {
    fn write_u8(&mut self, byte: u8) -> fmt::Result {
        self.written += 1;
        self.inner.write_u8(byte)
    }
}
impl<T: Write> Write for CountingWriter<T> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let res = self.inner.write(buf);
        if let Ok(written) = res {
            self.written += written;
        }
        res
    }
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self.inner.write_all(buf) {
            Ok(()) => (),
            Err(ref err) if err.kind() == io::ErrorKind::WriteZero => (),
            Err(err) => return Err(err),
        }
        self.written += buf.len();
        Ok(())
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}
