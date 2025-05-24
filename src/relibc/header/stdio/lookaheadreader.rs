struct LookAheadBuffer {
    buf: *const u8,
    pos: isize,
    look_ahead: isize,
}
impl LookAheadBuffer {
    fn look_ahead(&mut self) -> Result<Option<u8>, i32> {
        let byte = unsafe { *self.buf.offset(self.look_ahead) };
        if byte == 0 {
            Ok(None)
        } else {
            self.look_ahead += 1;
            Ok(Some(byte))
        }
    }

    fn commit(&mut self) {
        self.pos = self.look_ahead;
    }
}

impl From<*const u8> for LookAheadBuffer {
    fn from(buff: *const u8) -> LookAheadBuffer {
        LookAheadBuffer {
            buf: buff,
            pos: 0,
            look_ahead: 0,
        }
    }
}

enum LookAheadReaderEnum {
    // FILE(LookAheadFile<'a>),
    // (buffer, location)
    #[allow(clippy::upper_case_acronyms)]
    BUFFER(LookAheadBuffer),
}

pub struct LookAheadReader(LookAheadReaderEnum);

impl LookAheadReader {
    pub fn lookahead1(&mut self) -> Result<Option<u8>, i32> {
        match &mut self.0 {
            //LookAheadReaderEnum::FILE(f) => f.look_ahead(),
            LookAheadReaderEnum::BUFFER(b) => b.look_ahead(),
        }
    }
    pub fn commit(&mut self) {
        match &mut self.0 {
            //LookAheadReaderEnum::FILE(f) => f.commit(),
            LookAheadReaderEnum::BUFFER(b) => b.commit(),
        }
    }
}

impl From<*const u8> for LookAheadReader {
    fn from(buff: *const u8) -> LookAheadReader {
        LookAheadReader(LookAheadReaderEnum::BUFFER(buff.into()))
    }
}
