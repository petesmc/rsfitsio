use core::slice;
use std::cmp;

use crate::{
    cfileio::{ffflushx, ffread, ffread_int, ffseek, ffwrite, ffwrite_int},
    editcol::ffirow_safe,
    fitscore::{
        ffchdu, ffgext, ffghdn_safe, ffgtcl_safe, ffmahd_safe, ffpmsg_slice, ffpmsg_str,
        ffrdef_safe,
    },
    fitsio::*,
    fitsio2::*,
    int_snprintf,
    swapproc::{ffswap2, ffswap4, ffswap8},
};
use bytemuck::{cast_slice, cast_slice_mut};

use crate::c_types::{c_char, c_int, c_long, c_short, c_uchar};

/*--------------------------------------------------------------------------*/
/// Move to the input byte location in the file.
///
/// When writing to a file, a move
/// may sometimes be made to a position beyond the current EOF.  The err_mode
/// parameter determines whether such conditions should be returned as an error
/// or simply ignored.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmbyt(
    fptr: *mut fitsfile, /* I - FITS file pointer                */
    bytepos: LONGLONG,   /* I - byte position in file to move to */
    err_mode: c_int,     /* I - 1=ignore error, 0 = return error */
    status: *mut c_int,  /* IO - error status                    */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffmbyt_safe(fptr, bytepos, err_mode, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Move to the input byte location in the file.  When writing to a file, a move
/// may sometimes be made to a position beyond the current EOF.  The err_mode
/// parameter determines whether such conditions should be returned as an error
/// or simply ignored.
pub fn ffmbyt_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                */
    bytepos: LONGLONG,   /* I - byte position in file to move to */
    err_mode: c_int,     /* I - 1=ignore error, 0 = return error */
    status: &mut c_int,  /* IO - error status                    */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    if bytepos < 0 {
        *status = NEG_FILE_POS;
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let record = (bytepos / IOBUFLEN) as c_long; /* zero-indexed record number */

    /* if this is not the current record, then load it */
    if (fptr.Fptr.curbuf < 0) || (record != fptr.Fptr.bufrecnum[fptr.Fptr.curbuf as usize]) {
        ffldrc(fptr, record, err_mode, status);
    }

    if *status <= 0 {
        fptr.Fptr.bytepos = bytepos; /* save new file position */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the buffer of bytes to the output FITS file, starting at
/// the current file position.  Write large blocks of data directly to disk;
/// write smaller segments to intermediate IO buffers to improve efficiency.
pub(crate) fn ffpbyt(
    fptr: &mut fitsfile, /* I - FITS file pointer                    */
    nbytes: LONGLONG,    /* I - number of bytes to write             */
    buffer: &[u8],       /* I - buffer containing the bytes to write */
    status: &mut c_int,  /* IO - error status                        */
) -> c_int {
    let nbuff;
    let mut filepos: LONGLONG;
    let recstart: c_long;
    let recend: c_long;
    let mut bufpos: LONGLONG;
    let mut nspace: LONGLONG;
    let mut nwrite: LONGLONG;

    let mut j: usize = 0;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if nbytes > LONG_MAX as LONGLONG {
        ffpmsg_str("Number of bytes to write is greater than LONG_MAX (ffpbyt).");
        *status = WRITE_ERROR;
        return *status;
    }

    let mut ntodo = nbytes as LONGLONG;

    if fptr.Fptr.curbuf < 0 {
        /* no current data buffer for this file */
        /* so reload the last one that was used */
        ffldrc(
            fptr,
            ((fptr.Fptr.bytepos) / IOBUFLEN) as c_long,
            REPORT_EOF,
            status,
        );
    }

    if nbytes >= MINDIRECT {
        /* write large blocks of data directly to disk instead of via buffers */
        /* first, fill up the current IO buffer before flushing it to disk */

        nbuff = fptr.Fptr.curbuf as usize; /* current IO buffer number */
        filepos = fptr.Fptr.bytepos; /* save the write starting position */
        recstart = fptr.Fptr.bufrecnum[nbuff]; /* starting record */
        recend = ((filepos + nbytes - 1i64) / IOBUFLEN) as c_long; /* ending record   */

        /* bufpos is the starting position within the IO buffer */
        bufpos = (filepos - ((recstart as LONGLONG) * IOBUFLEN)) as LONGLONG;
        nspace = IOBUFLEN - bufpos; /* amount of space left in the buffer */

        if nspace != 0 {
            /* fill up the IO buffer */

            let _from = nbuff * IOBUFLEN as usize + bufpos as usize;
            let _to = _from + nspace as usize;
            fptr.Fptr.iobuffer[_from.._to]
                .copy_from_slice(cast_slice::<u8, c_char>(&buffer[j..(j + nspace as usize)]));

            ntodo -= nspace; /* decrement remaining number of bytes */
            j += nspace as usize; /* increment user buffer pointer */
            filepos += nspace; /* increment file position pointer */
            fptr.Fptr.dirty[nbuff] = TRUE as c_int; /* mark record as having been modified */
        }

        for ii in 0..(NIOBUF as usize) {
            /* flush any affected buffers to disk */

            if fptr.Fptr.bufrecnum[ii] >= recstart && fptr.Fptr.bufrecnum[ii] <= recend {
                if fptr.Fptr.dirty[ii] != 0 {
                    /* flush modified buffer to disk */
                    ffbfwt(&mut fptr.Fptr, ii as c_int, status);
                }

                fptr.Fptr.bufrecnum[ii] = -1; /* disassociate buffer from the file */
            }
        }

        /* move to the correct write position */
        if fptr.Fptr.io_pos != filepos {
            ffseek(&mut fptr.Fptr, filepos);
        }

        nwrite = ((ntodo - 1) / IOBUFLEN) * IOBUFLEN; /* don't write last buff */

        ffwrite(&mut fptr.Fptr, nwrite as c_long, &buffer[j..], status); /* write the data */
        ntodo -= nwrite; /* decrement remaining number of bytes */
        j += nwrite as usize; /* increment user buffer pointer */
        fptr.Fptr.io_pos = filepos + nwrite; /* update the file position */

        if fptr.Fptr.io_pos >= fptr.Fptr.filesize {
            /* at the EOF? */

            fptr.Fptr.filesize = fptr.Fptr.io_pos; /* increment file size */

            /* initialize the current buffer with the correct fill value */
            let _from = nbuff * IOBUFLEN as usize;
            let _n = IOBUFLEN as usize;
            if fptr.Fptr.hdutype == ASCII_TBL {
                fptr.Fptr.iobuffer[_from..(_from + _n)].fill(32); /* blank fill */
            } else {
                fptr.Fptr.iobuffer[_from..(_from + _n)].fill(0); /* zero fill */
            };
        } else {
            /* read next record */

            ffread_int(&mut fptr.Fptr, IOBUFLEN as usize, nbuff, status);
            // ffread(fptr.Fptr.as_mut(), IOBUFLEN, tmp, status);
            fptr.Fptr.io_pos += IOBUFLEN;
        }

        /* copy remaining bytes from user buffer into current IO buffer */
        let _from = nbuff * IOBUFLEN as usize;
        let _to = _from + ntodo as usize;
        fptr.Fptr.iobuffer[_from.._to]
            .copy_from_slice(cast_slice::<u8, c_char>(&buffer[j..(j + ntodo as usize)]));

        fptr.Fptr.dirty[nbuff] = TRUE as c_int; /* mark record as having been modified */
        fptr.Fptr.bufrecnum[nbuff] = recend; /* record number */
        fptr.Fptr.logfilesize =
            cmp::max(fptr.Fptr.logfilesize, ((recend + 1) as LONGLONG) * IOBUFLEN);
        fptr.Fptr.bytepos = filepos + nwrite + ntodo;
    } else {
        /* bufpos is the starting position in IO buffer */
        bufpos = (fptr.Fptr.bytepos
            - ((fptr.Fptr.bufrecnum[fptr.Fptr.curbuf as usize] as LONGLONG) * IOBUFLEN))
            as LONGLONG;
        nspace = IOBUFLEN - bufpos; /* amount of space left in the buffer */
        while ntodo != 0 {
            nwrite = cmp::min(ntodo, nspace);
            /* copy bytes from user's buffer to the IO buffer */

            let _from = fptr.Fptr.curbuf as usize * IOBUFLEN as usize + bufpos as usize;
            let _to = _from + nwrite as usize;
            fptr.Fptr.iobuffer[_from.._to]
                .copy_from_slice(cast_slice::<u8, c_char>(&buffer[j..(j + nwrite as usize)]));

            ntodo -= nwrite; /* decrement remaining number of bytes */
            j += nwrite as usize;
            fptr.Fptr.bytepos += nwrite; /* increment file position pointer */
            fptr.Fptr.dirty[fptr.Fptr.curbuf as usize] = TRUE as c_int; /* mark record as modified */
            if ntodo != 0 {
                /* load next record into a buffer */
                ffldrc(
                    fptr,
                    (fptr.Fptr.bytepos / IOBUFLEN) as c_long,
                    IGNORE_EOF,
                    status,
                );
                bufpos = 0;
                nspace = IOBUFLEN;
            };
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the buffer of bytes to the output FITS file, with an offset
/// between each group of bytes.  This function combines ffmbyt and ffpbyt
/// for increased efficiency.
pub(crate) fn ffpbytoff(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    gsize: c_long,       /* I - size of each group of bytes         */
    ngroups: c_long,     /* I - number of groups to write           */
    offset: c_long,      /* I - size of gap between groups          */
    buffer: &[u8],       /* I - buffer to be written                */
    status: &mut c_int,  /* IO - error status                       */
) -> c_int {
    let mut bcurrent: c_int;
    let mut bufpos: LONGLONG;
    let mut nspace: LONGLONG;
    let mut nwrite: LONGLONG;
    let mut record: c_long;

    let gsize = gsize as LONGLONG;
    let offset = offset as LONGLONG;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if fptr.Fptr.curbuf < 0 {
        /* no current data buffer for this file */

        /* so reload the last one that was used */
        ffldrc(
            fptr,
            ((fptr.Fptr.bytepos) / IOBUFLEN) as c_long,
            REPORT_EOF,
            status,
        );
    }

    bcurrent = fptr.Fptr.curbuf; /* number of the current IO buffer */
    record = fptr.Fptr.bufrecnum[bcurrent as usize]; /* zero-indexed record number */
    bufpos = (fptr.Fptr.bytepos - (record as LONGLONG * IOBUFLEN)) as LONGLONG; /* start pos */
    nspace = IOBUFLEN - bufpos; /* amount of space left in buffer */

    let mut cptr_index = 0; //buffer;
    let mut ioptr_index = (bcurrent as usize * IOBUFLEN as usize) + bufpos as usize; //iobuffer;

    let mut ii = 1;
    while ii < ngroups {
        /* write all but the last group */

        /* copy bytes from user's buffer to the IO buffer */
        nwrite = cmp::min(gsize, nspace);
        fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nwrite as usize)].copy_from_slice(
            cast_slice::<u8, c_char>(&buffer[cptr_index..(cptr_index + nwrite as usize)]),
        );
        cptr_index += nwrite as usize; /* increment buffer pointer */

        if nwrite < gsize {
            /* entire group did not fit */

            fptr.Fptr.dirty[bcurrent as usize] = TRUE as c_int; /* mark record as having been modified */
            record += 1;
            ffldrc(fptr, record, IGNORE_EOF, status); /* load next record */
            bcurrent = fptr.Fptr.curbuf;
            ioptr_index = bcurrent as usize * IOBUFLEN as usize;

            nwrite = gsize - nwrite;
            fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nwrite as usize)].copy_from_slice(
                cast_slice::<u8, c_char>(&buffer[cptr_index..(cptr_index + nwrite as usize)]),
            );
            cptr_index += nwrite as usize; /* increment buffer pointer */
            ioptr_index += (offset + nwrite) as usize; /* increment IO buffer pointer */

            nspace = IOBUFLEN - offset - nwrite; /* amount of space left */
        } else {
            ioptr_index += (offset + nwrite) as usize; /* increment IO bufer pointer */
            nspace -= offset + nwrite;
        }

        if nspace <= 0
        /* beyond current record? */
        {
            fptr.Fptr.dirty[bcurrent as usize] = TRUE as c_int;
            record += ((IOBUFLEN - nspace as LONGLONG) / IOBUFLEN) as c_long; /* new record number */
            ffldrc(fptr, record, IGNORE_EOF, status);
            bcurrent = fptr.Fptr.curbuf;

            bufpos = (-nspace) % IOBUFLEN; /* starting buffer pos */
            nspace = IOBUFLEN - bufpos;
            ioptr_index = (bcurrent as usize * IOBUFLEN as usize) + bufpos as usize;
            //iobuffer;
        }
        ii += 1;
    }

    /* now write the last group */
    nwrite = cmp::min(gsize, nspace);
    fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nwrite as usize)].copy_from_slice(cast_slice::<
        u8,
        c_char,
    >(
        &buffer[cptr_index..(cptr_index + nwrite as usize)],
    ));
    cptr_index += nwrite as usize; /* increment buffer pointer */

    if nwrite < gsize {
        /* entire group did not fit */

        fptr.Fptr.dirty[bcurrent as usize] = TRUE as c_int; /* mark record as having been modified */
        record += 1;
        ffldrc(fptr, record, IGNORE_EOF, status); /* load next record */
        bcurrent = fptr.Fptr.curbuf;
        ioptr_index = bcurrent as usize * IOBUFLEN as usize;

        nwrite = gsize - nwrite;
        fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nwrite as usize)].copy_from_slice(
            cast_slice::<u8, c_char>(&buffer[cptr_index..(cptr_index + nwrite as usize)]),
        );
    }

    fptr.Fptr.dirty[bcurrent as usize] = TRUE as c_int; /* mark record as having been modified */
    fptr.Fptr.bytepos =
        fptr.Fptr.bytepos + (ngroups as LONGLONG * gsize) + (ngroups as LONGLONG - 1) * offset;
    *status
}

/*--------------------------------------------------------------------------*/
/// get (read) the requested number of bytes from the file, starting at
/// the current file position.  Read large blocks of data directly from disk;
/// read smaller segments via intermediate IO buffers to improve efficiency.
pub(crate) fn ffgbyt(
    fptr: &mut fitsfile, /* I - FITS file pointer             */
    nbytes: LONGLONG,    /* I - number of bytes to read       */
    buffer: &mut [u8],   /* O - buffer to read into           */
    status: &mut c_int,  /* IO - error status                 */
) -> c_int {
    let mut ii: usize;
    let filepos: LONGLONG;
    let recstart: c_long;
    let recend: c_long;
    let mut ntodo: c_long;
    let mut bufpos: c_long;
    let mut nspace: c_long;
    let mut nread: c_long;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let mut cptr_index = 0; //buffer

    if nbytes >= MINDIRECT {
        /* read large blocks of data directly from disk instead of via buffers */
        filepos = fptr.Fptr.bytepos;

        /*  note that in this case, ffmbyt has not been called, and so
        bufrecnum[fptr.Fptr.curbuf] does not point to the intended
        output buffer */

        recstart = (filepos / IOBUFLEN) as c_long; /* starting record */
        recend = ((filepos + nbytes - 1i64) / IOBUFLEN) as c_long; /* ending record   */
        ii = 0;

        /* flush any affected buffers to disk */
        while ii < NIOBUF as usize {
            if fptr.Fptr.dirty[ii] != 0
                && fptr.Fptr.bufrecnum[ii] >= recstart
                && fptr.Fptr.bufrecnum[ii] <= recend
            {
                ffbfwt(&mut fptr.Fptr, ii as c_int, status); /* flush modified buffer to disk */
            };
            ii += 1
        }

        /* move to the correct read position */
        if fptr.Fptr.io_pos != filepos {
            ffseek(&mut fptr.Fptr, filepos);
        }

        ffread(&fptr.Fptr, nbytes as c_long, buffer, status); /* read the data */
        fptr.Fptr.io_pos = filepos + nbytes; /* update the file position */
    } else {
        /* read small chucks of data using the IO buffers for efficiency */

        if fptr.Fptr.curbuf < 0 {
            /* no current data buffer for this file */
            /* so reload the last one that was used */
            ffldrc(
                fptr,
                ((fptr.Fptr.bytepos) / IOBUFLEN) as c_long,
                REPORT_EOF,
                status,
            );
        }

        /* bufpos is the starting position in IO buffer */
        bufpos = (fptr.Fptr.bytepos
            - ((fptr.Fptr.bufrecnum[fptr.Fptr.curbuf as usize] as LONGLONG) * IOBUFLEN))
            as c_long;
        nspace = IOBUFLEN as c_long - bufpos; /* amount of space left in the buffer */

        ntodo = nbytes as c_long;
        while ntodo != 0 {
            nread = cmp::min(ntodo, nspace);

            /* copy bytes from IO buffer to user's buffer */
            let _from =
                ((LONGLONG::from(fptr.Fptr.curbuf) * IOBUFLEN) + bufpos as LONGLONG) as usize;
            let _to = _from + nread as usize;
            buffer[cptr_index..(cptr_index + nread as usize)]
                .copy_from_slice(cast_slice_mut(&mut fptr.Fptr.iobuffer[_from.._to]));

            ntodo -= nread; /* decrement remaining number of bytes */
            cptr_index += nread as usize;
            fptr.Fptr.bytepos += nread as LONGLONG; /* increment file position pointer */

            if ntodo != 0 {
                /* load next record into a buffer */
                ffldrc(
                    fptr,
                    (fptr.Fptr.bytepos / IOBUFLEN) as c_long,
                    REPORT_EOF,
                    status,
                );
                bufpos = 0;
                nspace = IOBUFLEN as c_long;
            };
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Get (read) the requested number of bytes from the file, starting at
/// the current file position.  This function combines ffmbyt and ffgbyt
/// for increased efficiency.
pub(crate) fn ffgbytoff(
    fptr: &mut fitsfile, /* I - FITS file pointer                   */
    gsize: c_long,       /* I - size of each group of bytes         */
    ngroups: c_long,     /* I - number of groups to read            */
    offset: c_long,      /* I - size of gap between groups (may be < 0) */
    buffer: &mut [u8],   /* I - buffer to be filled                 */
    status: &mut c_int,  /* IO - error status                       */
) -> c_int {
    let mut bcurrent: c_int;
    let mut bufpos: c_long;
    let mut nspace: c_long;
    let mut nread: c_long;
    let mut record: c_long;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if fptr.Fptr.curbuf < 0 {
        /* no current data buffer for this file */
        /* so reload the last one that was used */
        ffldrc(
            fptr,
            ((fptr.Fptr.bytepos) / IOBUFLEN) as c_long,
            REPORT_EOF,
            status,
        );
    }

    let mut cptr_index = 0; //buffer
    bcurrent = fptr.Fptr.curbuf; /* number of the current IO buffer */
    record = fptr.Fptr.bufrecnum[bcurrent as usize]; /* zero-indexed record number */
    bufpos = (fptr.Fptr.bytepos - (record as LONGLONG * IOBUFLEN)) as c_long; /* start pos */
    nspace = IOBUFLEN as c_long - bufpos; /* amount of space left in buffer */
    let mut ioptr_index = (bcurrent as usize * IOBUFLEN as usize) + bufpos as usize;

    for _ii in 1..(ngroups as usize) {
        /* read all but the last group */

        /* copy bytes from IO buffer to the user's buffer */
        nread = cmp::min(gsize, nspace);

        buffer[cptr_index..(cptr_index + nread as usize)].copy_from_slice(cast_slice_mut(
            &mut fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nread as usize)],
        ));
        cptr_index += nread as usize; /* increment buffer pointer */

        if nread < gsize {
            /* entire group did not fit */
            record += 1;
            ffldrc(fptr, record, REPORT_EOF, status); /* load next record */
            bcurrent = fptr.Fptr.curbuf;
            ioptr_index = bcurrent as usize * IOBUFLEN as usize;

            nread = gsize - nread;
            buffer[cptr_index..(cptr_index + nread as usize)].copy_from_slice(cast_slice_mut(
                &mut fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nread as usize)],
            ));
            cptr_index += nread as usize; /* increment buffer pointer */
            ioptr_index = ioptr_index
                .checked_add_signed((offset + nread) as isize)
                .unwrap(); /* increment IO buffer pointer */
            nspace = IOBUFLEN as c_long - offset - nread; /* amount of space left */
        } else {
            ioptr_index = ioptr_index
                .checked_add_signed((offset + nread) as isize)
                .unwrap(); /* increment IO buffer pointer */
            nspace -= offset + nread;
        }

        if nspace <= 0 || nspace > IOBUFLEN as c_long {
            /* beyond current record? */

            if nspace <= 0 {
                record += ((IOBUFLEN - nspace as LONGLONG) / IOBUFLEN) as c_long; /* new record number */
                bufpos = (-nspace) % IOBUFLEN as c_long; /* starting buffer pos */
            } else {
                record -= (nspace - 1) / IOBUFLEN as c_long; /* new record number */
                bufpos = (IOBUFLEN - (nspace as LONGLONG % IOBUFLEN)) as c_long; /* starting buffer pos */
            }

            ffldrc(fptr, record, REPORT_EOF, status);
            bcurrent = fptr.Fptr.curbuf;

            nspace = IOBUFLEN as c_long - bufpos;
            ioptr_index = (bcurrent as usize * IOBUFLEN as usize) + bufpos as usize;
        }
    }

    /* now read the last group */
    nread = cmp::min(gsize, nspace);
    buffer[cptr_index..(cptr_index + nread as usize)].copy_from_slice(cast_slice_mut(
        &mut fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nread as usize)],
    ));
    cptr_index += nread as usize; /* increment buffer pointer */

    if nread < gsize {
        /* entire group did not fit */

        record += 1;
        ffldrc(fptr, record, REPORT_EOF, status); /* load next record */
        bcurrent = fptr.Fptr.curbuf;
        ioptr_index = bcurrent as usize * IOBUFLEN as usize;

        nread = gsize - nread;
        buffer[cptr_index..(cptr_index + nread as usize)].copy_from_slice(cast_slice_mut(
            &mut fptr.Fptr.iobuffer[ioptr_index..(ioptr_index + nread as usize)],
        ));
    }

    fptr.Fptr.bytepos = fptr.Fptr.bytepos
        + (ngroups * gsize) as LONGLONG
        + (ngroups as LONGLONG - 1) * offset as LONGLONG;
    *status
}

/*--------------------------------------------------------------------------*/
/// low-level routine to load a specified record from a file into
/// a physical buffer, if it is not already loaded.  Reset all
/// pointers to make this the new current record for that file.
/// Update ages of all the physical buffers.
pub(crate) fn ffldrc(
    fptr: &mut fitsfile, /* I - FITS file pointer             */
    record: c_long,      /* I - record number to be loaded    */
    err_mode: c_int,     /* I - 1=ignore EOF, 0 = return EOF error */
    status: &mut c_int,  /* IO - error status                 */
) -> c_int {
    let mut updatebuf = false;
    let mut ibuff: c_int;
    let mut nbuff = 0;

    /* check if record is already loaded in one of the buffers */
    /* search from youngest to oldest buffer for efficiency */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status); /* EOF? */
    }

    ibuff = (NIOBUF - 1) as c_int;
    while ibuff >= 0 {
        nbuff = fptr.Fptr.ageindex[ibuff as usize]; /* blank fill */
        if record == fptr.Fptr.bufrecnum[nbuff as usize] {
            updatebuf = true;
            break;
        };
        ibuff -= 1;
    }

    if !updatebuf {
        /* record is not already loaded */
        let rstart = (record as LONGLONG) * IOBUFLEN;

        if err_mode == 0 && (rstart >= fptr.Fptr.logfilesize) {
            /* EOF? */
            *status = END_OF_FILE;
            return *status;
        }

        if ffwhbf(fptr, &mut nbuff) < 0 {
            /* which buffer should we reuse? */
            *status = TOO_MANY_FILES;
            return *status;
        }

        if fptr.Fptr.dirty[nbuff as usize] != 0 {
            /* write dirty buffer to disk */
            ffbfwt(&mut fptr.Fptr, nbuff, status);
        }

        if rstart >= fptr.Fptr.filesize {
            /* EOF? */

            /* initialize an empty buffer with the correct fill value */
            let _from = nbuff as usize * IOBUFLEN as usize;
            let _n = IOBUFLEN as usize;
            if fptr.Fptr.hdutype == ASCII_TBL {
                /* blank fill */
                fptr.Fptr.iobuffer[_from..(_from + _n)].fill(32);
            } else {
                /* zero fill */
                fptr.Fptr.iobuffer[_from..(_from + _n)].fill(0);
            }

            fptr.Fptr.logfilesize = cmp::max(fptr.Fptr.logfilesize, rstart + IOBUFLEN);
            fptr.Fptr.dirty[nbuff as usize] = TRUE as c_int; /* mark record as having been modified */
        } else {
            /* not EOF, so read record from disk */

            if fptr.Fptr.io_pos != rstart {
                ffseek(&mut fptr.Fptr, rstart);
            }

            ffread_int(&mut fptr.Fptr, IOBUFLEN as usize, nbuff as usize, status);
            //ffread(&fptr.Fptr, IOBUFLEN as c_long, tmp, status);
            fptr.Fptr.io_pos = rstart + IOBUFLEN; /* set new IO position */
        }
        fptr.Fptr.bufrecnum[nbuff as usize] = record; /* record number contained in buffer */
    }

    fptr.Fptr.curbuf = nbuff; /* this is the current buffer for this file */

    if ibuff < 0 {
        /* find the current position of the buffer in the age index */
        ibuff = 0;
        while ibuff < NIOBUF as c_int {
            if fptr.Fptr.ageindex[ibuff as usize] == nbuff {
                break;
            }
            ibuff += 1
        }
    }

    /* increment the age of all the buffers that were younger than it */
    ibuff += 1;
    while ibuff < NIOBUF as c_int {
        fptr.Fptr.ageindex[ibuff as usize - 1] = fptr.Fptr.ageindex[ibuff as usize];
        ibuff += 1
    }
    fptr.Fptr.ageindex[NIOBUF as usize - 1] = nbuff; /* this is now the youngest buffer */
    *status
}

/*--------------------------------------------------------------------------*/
/// decide which buffer to (re)use to hold a new file record
pub(crate) fn ffwhbf(
    fptr: &mut fitsfile, /* I - FITS file pointer             */
    nbuff: &mut c_int,   /* O - which buffer to use           */
) -> c_int {
    *nbuff = fptr.Fptr.ageindex[0]; /* return oldest buffer */
    *nbuff
}

/*--------------------------------------------------------------------------*/
/// Flush all the data in the current FITS file to disk. This ensures that if
/// the program subsequently dies, the disk FITS file will be closed correctly.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffflus(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffflus_safe(fptr, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Flush all the data in the current FITS file to disk. This ensures that if
/// the program subsequently dies, the disk FITS file will be closed correctly.
pub fn ffflus_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut hdunum: c_int = 0;
    let mut hdutype: c_int = 0;

    if *status > 0 {
        return *status;
    }

    ffghdn_safe(fptr, &mut hdunum); /* get the current HDU number */

    if ffchdu(fptr, status) > 0 {
        /* close out the current HDU */
        ffpmsg_str("ffflus could not close the current HDU.");
    }

    ffflsh_safe(fptr, false, status); /* flush any modified IO buffers to disk */

    if ffgext(fptr, hdunum - 1, Some(&mut hdutype), status) > 0 {
        /* reopen HDU */
        ffpmsg_str("ffflus could not reopen the current HDU.");
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// flush all dirty IO buffers associated with the file to disk
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffflsh(
    fptr: *mut fitsfile, /* I - FITS file pointer           */
    clearbuf: c_int,     /* I - also clear buffer contents? */
    status: *mut c_int,  /* IO - error status               */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffflsh_safe(fptr, clearbuf != 0, status)
    }
}

/*--------------------------------------------------------------------------*/
/// flush all dirty IO buffers associated with the file to disk
pub fn ffflsh_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    clearbuf: bool,      /* I - also clear buffer contents? */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    /*
       no need to move to a different HDU

        if (fptr->HDUposition != (fptr->Fptr)->curhdu)
            ffmahd(fptr, (fptr->HDUposition) + 1, NULL, status);
    */

    for ii in 0..(NIOBUF as usize) {
        /* flush system buffers to disk */
        if fptr.Fptr.bufrecnum[ii] >= 0 && fptr.Fptr.dirty[ii] != 0 {
            ffbfwt(&mut fptr.Fptr, ii as c_int, status);
        }
        if clearbuf {
            fptr.Fptr.bufrecnum[ii] = -1; /* set contents of buffer as undefined */
        };
    }

    if *status != READONLY_FILE {
        ffflushx(&mut fptr.Fptr); /* flush system buffers to disk */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// clear any buffers beyond the end of file
pub(crate) fn ffbfeof(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    for ii in 0..(NIOBUF as usize) {
        if (fptr.Fptr.bufrecnum[ii] as LONGLONG) * IOBUFLEN >= fptr.Fptr.filesize {
            fptr.Fptr.bufrecnum[ii] = -1; /* set contents of buffer as undefined */
        };
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// write contents of buffer to file;  If the position of the buffer
/// is beyond the current EOF, then the file may need to be extended
/// with fill values, and/or with the contents of some of the other
/// i/o buffers.
pub(crate) fn ffbfwt(
    Fptr: &mut FITSfile, /* I - FITS file pointer           */
    nbuff: c_int,        /* I - which buffer to write          */
    status: &mut c_int,  /* IO - error status                  */
) -> c_int {
    let mut ibuff: usize;
    let mut jj: c_long;
    let mut irec: c_long;
    let mut minrec: c_long;
    let mut nloop: c_long;
    let mut zeros: [i8; IOBUFLEN as usize] = [0; IOBUFLEN as usize]; /*  initialized to zero by default */

    if (Fptr.writemode) == 0 {
        ffpmsg_str("Error: trying to write to READONLY file.");
        if Fptr.driver == 8 {
            /* gzip compressed file */
            ffpmsg_str("Cannot write to a GZIP or COMPRESS compressed file.");
        }

        Fptr.dirty[nbuff as usize] = FALSE as c_int; /* reset buffer status to prevent later probs */
        *status = READONLY_FILE;
        return *status;
    }

    let mut filepos = (Fptr.bufrecnum[nbuff as usize] as LONGLONG) * IOBUFLEN;

    if filepos <= Fptr.filesize {
        /* record is located within current file, so just write it */

        /* move to the correct write position */
        if Fptr.io_pos != filepos {
            ffseek(Fptr, filepos);
        }

        ffwrite_int(Fptr, IOBUFLEN as usize, nbuff as usize, status);
        // ffwrite(Fptr, IOBUFLEN, tmp, status);

        Fptr.io_pos = filepos + IOBUFLEN;

        /* appended new record? */
        if filepos == Fptr.filesize {
            Fptr.filesize += IOBUFLEN; /* increment the file size */
        }

        Fptr.dirty[nbuff as usize] = FALSE as c_int;
    } else {
        /* if record is beyond the EOF, append any other records */
        /* and/or insert fill values if necessary */

        /* move to EOF */
        if Fptr.io_pos != Fptr.filesize {
            ffseek(Fptr, Fptr.filesize);
        }

        ibuff = NIOBUF as usize; /* initialize to impossible value */

        /* repeat until requested buffer is written */
        while ibuff != nbuff as usize {
            minrec = (Fptr.filesize / IOBUFLEN) as c_long;

            /* write lowest record beyond the EOF first */
            irec = Fptr.bufrecnum[nbuff as usize]; /* initially point to the requested buffer */
            ibuff = nbuff as usize;

            let mut ii: usize = 0;
            while ii < NIOBUF as usize {
                if Fptr.bufrecnum[ii] >= minrec && Fptr.bufrecnum[ii] < irec {
                    irec = Fptr.bufrecnum[ii]; /* found a lower record */
                    ibuff = ii;
                };
                ii += 1
            }
            filepos = (irec as LONGLONG) * IOBUFLEN; /* byte offset of record in file */

            /* append 1 or more fill records if necessary */
            if filepos > Fptr.filesize {
                nloop = ((filepos - (Fptr.filesize)) / IOBUFLEN) as c_long;
                jj = 0;
                while jj < nloop && (*status) == 0 {
                    ffwrite(Fptr, IOBUFLEN as c_long, cast_slice_mut(&mut zeros), status);
                    jj += 1
                }

                Fptr.filesize = filepos; /* increment the file size */
            }

            /* write the buffer itself */
            ffwrite_int(Fptr, IOBUFLEN as usize, ibuff, status);
            // ffwrite(Fptr, IOBUFLEN, tmp, status);

            Fptr.dirty[ibuff] = FALSE as c_int;
            Fptr.filesize += IOBUFLEN; /* increment the file size */
        } /* loop back if more buffers need to be written */

        Fptr.io_pos = Fptr.filesize; /* currently positioned at EOF */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Returns an optimal value for the number of rows in a binary table
/// or the number of pixels in an image that should be read or written
/// at one time for maximum efficiency. Accessing more data than this
/// may cause excessive flushing and rereading of buffers to/from disk.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgrsz(
    fptr: *mut fitsfile, /* I - FITS file pionter                        */
    ndata: *mut c_long,  /* O - optimal amount of data to access         */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let ndata = ndata.as_mut().expect(NULL_MSG);

        ffgrsz_safe(fptr, ndata, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Returns an optimal value for the number of rows in a binary table
/// or the number of pixels in an image that should be read or written
/// at one time for maximum efficiency. Accessing more data than this
/// may cause excessive flushing and rereading of buffers to/from disk.
pub fn ffgrsz_safe(
    fptr: &mut fitsfile, /* I - FITS file pionter                        */
    ndata: &mut c_long,  /* O - optimal amount of data to access         */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut typecode = 0;
    let bytesperpixel;

    /* There are NIOBUF internal buffers available each IOBUFLEN bytes long. */

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        /* rescan header to get hdu struct */

        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        /* calc pixels per buffer size */

        /* image pixels are in column 2 of the 'table' */
        ffgtcl_safe(fptr, 2, Some(&mut typecode), None, None, status);
        bytesperpixel = typecode / 10;
        *ndata = ((NIOBUF as c_long - 1) * IOBUFLEN as c_long) / c_long::from(bytesperpixel);
    } else {
        /* calc number of rows that fit in buffers */

        *ndata = ((NIOBUF as c_long - 1) * IOBUFLEN as c_long)
            / cmp::max(1, fptr.Fptr.rowlength) as c_long;
        *ndata = cmp::max(1, *ndata);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// read a consecutive string of bytes from an ascii or binary table.
/// This will span multiple rows of the table if nchars + firstchar is
/// greater than the length of a row.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtbb(
    fptr: *mut fitsfile,  /* I - FITS file pointer                 */
    firstrow: LONGLONG,   /* I - starting row (1 = first row)      */
    firstchar: LONGLONG,  /* I - starting byte in row (1=first)    */
    nchars: LONGLONG,     /* I - number of bytes to read           */
    values: *mut c_uchar, /* I - array of bytes to read            */
    status: *mut c_int,   /* IO - error status                     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let values = slice::from_raw_parts_mut(values, nchars as usize);

        ffgtbb_safe(fptr, firstrow, firstchar, nchars, values, status)
    }
}

/*--------------------------------------------------------------------------*/
/// read a consecutive string of bytes from an ascii or binary table.
/// This will span multiple rows of the table if nchars + firstchar is
/// greater than the length of a row.
pub fn ffgtbb_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                 */
    firstrow: LONGLONG,     /* I - starting row (1 = first row)      */
    firstchar: LONGLONG,    /* I - starting byte in row (1=first)    */
    nchars: LONGLONG,       /* I - number of bytes to read           */
    values: &mut [c_uchar], /* I - array of bytes to read            */
    status: &mut c_int,     /* IO - error status                     */
) -> c_int {
    if *status > 0 || nchars <= 0 {
        return *status;
    } else if firstrow < 1 {
        *status = BAD_ROW_NUM;
        return *status;
    } else if firstchar < 1 {
        *status = BAD_ELEM_NUM;
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    /* check that we do not exceed number of rows in the table */
    let endrow: LONGLONG = ((firstchar + nchars - 2) / fptr.Fptr.rowlength) + firstrow;
    if endrow > fptr.Fptr.numrows {
        ffpmsg_str("attempt to read past end of table (ffgtbb)");
        *status = BAD_ROW_NUM;
        return *status;
    }

    /* move the i/o pointer to the start of the sequence of characters */
    let bytepos: LONGLONG =
        fptr.Fptr.datastart + (fptr.Fptr.rowlength * (firstrow - 1)) + firstchar - 1;

    ffmbyt_safe(fptr, bytepos, REPORT_EOF, status);
    ffgbyt(fptr, nchars, values, status); /* read the bytes */

    *status
}

/*--------------------------------------------------------------------------*/
/// get (read) the array of values from the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffgi1b(
    fptr: &mut fitsfile,    /* I - FITS file pointer                         */
    byteloc: LONGLONG,      /* I - position within file to start reading     */
    nvals: c_long,          /* I - number of pixels to read                  */
    incre: c_long,          /* I - byte increment between pixels             */
    values: &mut [c_uchar], /* O - returned array of values           */
    status: &mut c_int,     /* IO - error status                             */
) -> c_int {
    let postemp: LONGLONG;

    if incre == 1 {
        /* read all the values at once (contiguous bytes) */

        if (nvals as LONGLONG) < MINDIRECT {
            /* read normally via IO buffers */

            ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
            ffgbyt(fptr, nvals as LONGLONG, values, status);
        } else {
            /* read directly from disk, bypassing IO buffers */

            postemp = fptr.Fptr.bytepos; /* store current file position */
            fptr.Fptr.bytepos = byteloc; /* set to the desired position */
            ffgbyt(fptr, nvals as LONGLONG, values, status);
            fptr.Fptr.bytepos = postemp; /* reset to original position */
        }
    } else {
        /* have to read each value individually (not contiguous ) */

        ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
        ffgbytoff(fptr, 1, nvals, incre - 1, values, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// get (read) the array of values from the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffgi2b(
    fptr: &mut fitsfile,    /* I - FITS file pointer                        */
    byteloc: LONGLONG,      /* I - position within file to start reading    */
    nvals: c_long,          /* I - number of pixels to read                 */
    incre: c_long,          /* I - byte increment between pixels            */
    values: &mut [c_short], /* O - returned array of values                 */
    status: &mut c_int,     /* IO - error status                            */
) -> c_int {
    if incre == 2 {
        /* read all the values at once (contiguous bytes) */

        if nvals as LONGLONG * 2 < MINDIRECT {
            /* read normally via IO buffers */

            ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
            ffgbyt(fptr, nvals as LONGLONG * 2, cast_slice_mut(values), status);
        } else {
            /* read directly from disk, bypassing IO buffers */

            let postemp = fptr.Fptr.bytepos; /* store current file position */
            fptr.Fptr.bytepos = byteloc; /* set to the desired position */
            ffgbyt(fptr, nvals as LONGLONG * 2, cast_slice_mut(values), status);
            fptr.Fptr.bytepos = postemp; /* reset to original position */
        }
    } else {
        /* have to read each value individually (not contiguous ) */

        ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
        ffgbytoff(fptr, 2, nvals, incre - 2, cast_slice_mut(values), status);
    }

    if BYTESWAPPED {
        ffswap2(values, nvals); /* reverse order of bytes in each value */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// get (read) the array of values from the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffgi4b(
    fptr: &mut fitsfile,     /* I - FITS file pointer                        */
    byteloc: LONGLONG,       /* I - position within file to start reading    */
    nvals: c_long,           /* I - number of pixels to read                 */
    incre: c_long,           /* I - byte increment between pixels            */
    values: &mut [INT32BIT], /* O - returned array of values                */
    status: &mut c_int,      /* IO - error status                            */
) -> c_int {
    if incre == 4 {
        /* read all the values at once (contiguous bytes) */

        if nvals as LONGLONG * 4 < MINDIRECT {
            /* read normally via IO buffers */

            ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
            ffgbyt(fptr, nvals as LONGLONG * 4, cast_slice_mut(values), status);
        } else {
            /* read directly from disk, bypassing IO buffers */

            let postemp = fptr.Fptr.bytepos; /* store current file position */
            fptr.Fptr.bytepos = byteloc; /* set to the desired position */
            ffgbyt(fptr, nvals as LONGLONG * 4, cast_slice_mut(values), status);
            fptr.Fptr.bytepos = postemp; /* reset to original position */
        }
    } else {
        /* have to read each value individually (not contiguous ) */
        ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
        ffgbytoff(fptr, 4, nvals, incre - 4, cast_slice_mut(values), status);
    }

    if BYTESWAPPED {
        ffswap4(values, nvals); /* reverse order of bytes in each value */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// get (read) the array of values from the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
///
/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// This routine reads 'nvals' 8-byte integers into 'values'.
/// This works both on platforms that have sizeof(long) = 64, and 32,
/// as long as 'values' has been allocated to large enough to hold
/// 8 * nvals bytes of data.
/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pub(crate) fn ffgi8b(
    fptr: &mut fitsfile,   /* I - FITS file pointer                        */
    byteloc: LONGLONG,     /* I - position within file to start reading    */
    nvals: c_long,         /* I - number of pixels to read                 */
    incre: c_long,         /* I - byte increment between pixels            */
    values: &mut [c_long], /* O - returned array of values                 */
    status: &mut c_int,    /* IO - error status                            */
) -> c_int {
    if incre == 8 {
        /* read all the values at once (contiguous bytes) */

        if nvals as LONGLONG * 8 < MINDIRECT {
            /* read normally via IO buffers */

            ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
            ffgbyt(fptr, nvals as LONGLONG * 8, cast_slice_mut(values), status);
        } else {
            /* read directly from disk, bypassing IO buffers */

            let postemp = fptr.Fptr.bytepos; /* store current file position */
            fptr.Fptr.bytepos = byteloc; /* set to the desired position */
            ffgbyt(fptr, nvals as LONGLONG * 8, cast_slice_mut(values), status);
            fptr.Fptr.bytepos = postemp; /* reset to original position */
        }
    } else {
        /* have to read each value individually (not contiguous ) */

        ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
        ffgbytoff(fptr, 8, nvals, incre - 8, cast_slice_mut(values), status);
    }

    if BYTESWAPPED {
        ffswap8(cast_slice_mut(values), nvals); /* reverse bytes in each value */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// get (read) the array of values from the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffgr4b(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    byteloc: LONGLONG,   /* I - position within file to start reading    */
    nvals: c_long,       /* I - number of pixels to read                 */
    incre: c_long,       /* I - byte increment between pixels            */
    values: &mut [f32],  /* O - returned array of values                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    if incre == 4 {
        /* read all the values at once (contiguous bytes) */

        if nvals as LONGLONG * 4 < MINDIRECT {
            /* read normally via IO buffers */

            ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
            ffgbyt(fptr, nvals as LONGLONG * 4, cast_slice_mut(values), status);
        } else {
            /* read directly from disk, bypassing IO buffers */

            let postemp = fptr.Fptr.bytepos; /* store current file position */
            fptr.Fptr.bytepos = byteloc; /* set to the desired position */
            ffgbyt(fptr, nvals as LONGLONG * 4, cast_slice_mut(values), status);
            fptr.Fptr.bytepos = postemp; /* reset to original position */
        }
    } else {
        /* have to read each value individually (not contiguous ) */

        ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
        ffgbytoff(fptr, 4, nvals, incre - 4, cast_slice_mut(values), status);
    }

    if CFITSIO_MACHINE == VAXVMS {
        todo!();
        //let ii = nvals; /* call VAX macro routine to convert */
        //ieevur(values, values, &ii); /* from  IEEE float -> F float       */
    } else if (CFITSIO_MACHINE == ALPHAVMS) && (FLOATTYPE == GFLOAT) {
        todo!("ALPHAVMS && GFLOAT not implemented");

        /*
        ffswap2(cast_slice_mut(values), nvals * 2); /* swap pairs of bytes */

        /* convert from IEEE float format to VMS GFLOAT float format */
        let mut sptr = 0;

        // Can't do cast_slice because purposely need to go around the borrow checker here...
        // TODO UNSAFE HACK
        let shortBuffer: &[c_short] = slice::from_raw_parts(values.as_ptr() as *const c_short, nvals as usize * 4);

        for ii in 0..(nvals as usize) {
            if fnan(shortBuffer[sptr]) == 0 {
                /* test for NaN or underflow */
                values[ii] *= 4.0;
            }
            sptr += 2;
        }
        */
    } else if BYTESWAPPED {
        ffswap4(cast_slice_mut(values), nvals); /* reverse order of bytes in values */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// get (read) the array of values from the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffgr8b(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    byteloc: LONGLONG,   /* I - position within file to start reading    */
    nvals: c_long,       /* I - number of pixels to read                 */
    incre: c_long,       /* I - byte increment between pixels            */
    values: &mut [f64],  /* O - returned array of values                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    if incre == 8 {
        /* read all the values at once (contiguous bytes) */

        if nvals as LONGLONG * 8 < MINDIRECT {
            /* read normally via IO buffers */

            ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
            ffgbyt(fptr, nvals as LONGLONG * 8, cast_slice_mut(values), status);
        } else {
            /* read directly from disk, bypassing IO buffers */

            let postemp = fptr.Fptr.bytepos; /* store current file position */
            fptr.Fptr.bytepos = byteloc; /* set to the desired position */
            ffgbyt(fptr, nvals as LONGLONG * 8, cast_slice_mut(values), status);
            fptr.Fptr.bytepos = postemp; /* reset to original position */
        }
    } else {
        /* have to read each value individually (not contiguous ) */

        ffmbyt_safe(fptr, byteloc, REPORT_EOF, status);
        ffgbytoff(fptr, 8, nvals, incre - 8, cast_slice_mut(values), status);
    }

    if CFITSIO_MACHINE == VAXVMS {
        /* call VAX macro routine to convert */
        todo!();
        //ieevud(values, values, &ii); /* from  IEEE float -> D float       */
    } else if (CFITSIO_MACHINE == ALPHAVMS) && (FLOATTYPE == GFLOAT) {
        todo!("ALPHAVMS && GFLOAT not implemented");

        /*
        ffswap2(cast_slice_mut(values), nvals * 4); /* swap pairs of bytes */

        /* convert from IEEE float format to VMS GFLOAT float format */
        let mut sptr = 0;

        // Can't do cast_slice because purposely need to go around the borrow checker here...
        // TODO UNSAFE HACK
        let shortBuffer: &[c_short] = slice::from_raw_parts(values.as_ptr() as *const c_short, nvals as usize * 4);

        for ii in 0..(nvals as usize) {
            if dnan(shortBuffer[sptr]) == 0 {
                /* test for NaN or underflow */
                values[ii] *= 4.0;
            }
            sptr += 4;
        }
        */
    } else if BYTESWAPPED {
        ffswap8(cast_slice_mut(values), nvals); /* reverse order of bytes in each value */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// write a consecutive string of bytes to an ascii or binary table.
/// This will span multiple rows of the table if nchars + firstchar is
/// greater than the length of a row.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffptbb(
    fptr: *mut fitsfile,    /* I - FITS file pointer                 */
    firstrow: LONGLONG,     /* I - starting row (1 = first row)      */
    firstchar: LONGLONG,    /* I - starting byte in row (1=first)    */
    nchars: LONGLONG,       /* I - number of bytes to write          */
    values: *const c_uchar, /* I - array of bytes to write           */
    status: *mut c_int,     /* IO - error status                     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let values: &[c_uchar] = cast_slice(slice::from_raw_parts(values, nchars as usize));

        ffptbb_safe(fptr, firstrow, firstchar, nchars, values, status)
    }
}

/*--------------------------------------------------------------------------*/
/// write a consecutive string of bytes to an ascii or binary table.
/// This will span multiple rows of the table if nchars + firstchar is
/// greater than the length of a row.
pub fn ffptbb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                 */
    firstrow: LONGLONG,  /* I - starting row (1 = first row)      */
    firstchar: LONGLONG, /* I - starting byte in row (1=first)    */
    nchars: LONGLONG,    /* I - number of bytes to write          */
    values: &[c_uchar],  /* I - array of bytes to write           */
    status: &mut c_int,  /* IO - error status                     */
) -> c_int {
    let nrows: LONGLONG;
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 || nchars <= 0 {
        return *status;
    } else if firstrow < 1 {
        *status = BAD_ROW_NUM;
        return *status;
    } else if firstchar < 1 {
        *status = BAD_ELEM_NUM;
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    } else if fptr.Fptr.datastart < 0 {
        /* rescan header if data undefined */
        ffrdef_safe(fptr, status);
    }

    let endrow: LONGLONG = ((firstchar + nchars - 2) / fptr.Fptr.rowlength) + firstrow;

    /* check if we are writing beyond the current end of table */
    if endrow > fptr.Fptr.numrows {
        /* if there are more HDUs following the current one, or */
        /* if there is a data heap, then we must insert space */
        /* for the new rows.  */
        if (fptr.Fptr.lasthdu) == 0 || fptr.Fptr.heapsize > 0 {
            nrows = endrow - (fptr.Fptr.numrows);

            /* ffirow also updates the heap address and numrows */
            if ffirow_safe(fptr, fptr.Fptr.numrows, nrows, status) > 0 {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "ffptbb failed to add space for {:.0} new rows in table.",
                    nrows as f64,
                );
                ffpmsg_slice(&message);
                return *status;
            }
        } else {
            /* manally update heap starting address */
            fptr.Fptr.heapstart += (endrow - fptr.Fptr.numrows) as LONGLONG * fptr.Fptr.rowlength;

            fptr.Fptr.numrows = endrow; /* update number of rows */
        }
    }

    /* move the i/o pointer to the start of the sequence of characters */
    let bytepos: LONGLONG =
        fptr.Fptr.datastart + (fptr.Fptr.rowlength * (firstrow - 1)) + firstchar - 1;

    ffmbyt_safe(fptr, bytepos, IGNORE_EOF, status);
    ffpbyt(fptr, nchars, values, status); /* write the bytes */

    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the array of values to the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffpi1b(
    fptr: &mut fitsfile, /* I - FITS file pointer                         */
    nvals: c_long,       /* I - number of pixels in the values array      */
    incre: c_long,       /* I - byte increment between pixels             */
    values: &[u8],       /* I - array of values to write           */
    status: &mut c_int,  /* IO - error status                             */
) -> c_int {
    if incre == 1 {
        /* write all the values at once (contiguous bytes) */

        ffpbyt(fptr, nvals as LONGLONG, values, status);
    } else {
        /* have to write each value individually (not contiguous ) */
        ffpbytoff(fptr, 1, nvals, incre - 1, values, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the array of values to the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffpi2b(
    fptr: &mut fitsfile, /* I - FITS file pointer                         */
    nvals: c_long,       /* I - number of pixels in the values array      */
    incre: c_long,       /* I - byte increment between pixels             */
    values: &[c_short],  /* I - array of values to write                  */
    status: &mut c_int,  /* IO - error status                             */
) -> c_int {
    // Copy slice so that we don't change the original values
    let mut v: Vec<c_short> = values.to_vec();

    if BYTESWAPPED {
        ffswap2(&mut v, nvals); /* reverse order of bytes in each value */
    }

    if incre == 2 {
        /* write all the values at once (contiguous bytes) */
        ffpbyt(fptr, nvals as LONGLONG * 2, cast_slice(&v), status);
    } else {
        /* have to write each value individually (not contiguous ) */
        ffpbytoff(fptr, 2, nvals, incre - 2, cast_slice(&v), status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the array of values to the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffpi4b(
    fptr: &mut fitsfile, /* I - FITS file pointer                         */
    nvals: c_long,       /* I - number of pixels in the values array      */
    incre: c_long,       /* I - byte increment between pixels             */
    values: &[INT32BIT], /* I - array of values to write                */
    status: &mut c_int,  /* IO - error status                             */
) -> c_int {
    // Copy slice so that we don't change the original values
    let mut v: Vec<INT32BIT> = values.to_vec();

    if BYTESWAPPED {
        ffswap4(&mut v, nvals); /* reverse order of bytes in each value */
    }

    if incre == 4 {
        /* write all the values at once (contiguous bytes) */

        ffpbyt(fptr, nvals as LONGLONG * 4, cast_slice(&v), status);
    } else {
        /* have to write each value individually (not contiguous ) */

        ffpbytoff(fptr, 4, nvals, incre - 4, cast_slice(&v), status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the array of values to the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
///
/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// This routine writes 'nvals' 8-byte integers from 'values'.
/// This works both on platforms that have sizeof(long) = 64, and 32,
/// as long as 'values' has been allocated to large enough to hold
/// 8 * nvals bytes of data.
/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pub(crate) fn ffpi8b(
    fptr: &mut fitsfile, /* I - FITS file pointer                         */
    nvals: c_long,       /* I - number of pixels in the values array      */
    incre: c_long,       /* I - byte increment between pixels             */
    values: &[c_long],   /* I - array of values to write                */
    status: &mut c_int,  /* IO - error status                             */
) -> c_int {
    // Copy slice so that we don't change the original values
    let mut v: Vec<c_long> = values.to_vec();

    if BYTESWAPPED {
        ffswap8(cast_slice_mut(&mut v), nvals); /* reverse bytes in each value */
    }

    if incre == 8 {
        /* write all the values at once (contiguous bytes) */
        ffpbyt(fptr, nvals as LONGLONG * 8, cast_slice(&v), status);
    } else {
        /* have to write each value individually (not contiguous ) */
        ffpbytoff(fptr, 8, nvals, incre - 8, cast_slice(&v), status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the array of values to the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffpr4b(
    fptr: &mut fitsfile, /* I - FITS file pointer                         */
    nvals: c_long,       /* I - number of pixels in the values array      */
    incre: c_long,       /* I - byte increment between pixels             */
    values: &[f32],      /* I - array of values to write                  */
    status: &mut c_int,  /* IO - error status                             */
) -> c_int {
    // Copy slice so that we don't change the original values
    let mut v: Vec<f32> = values.to_vec();

    if CFITSIO_MACHINE == VAXVMS {
        todo!(); /* call VAX macro routine to convert */
    //ieevpr(values, values, &ii);     /* from F float -> IEEE float        */
    } else if (CFITSIO_MACHINE == ALPHAVMS) && (FLOATTYPE == GFLOAT) {
        /* convert from VMS FFLOAT float format to IEEE float format */
        for ii in 0..(nvals as usize) {
            v[ii] *= 0.25;
        }

        ffswap2(cast_slice_mut(&mut v), nvals * 2); /* swap pairs of bytes */
    } else if BYTESWAPPED {
        ffswap4(cast_slice_mut(&mut v), nvals); /* reverse order of bytes in values */
    }

    if incre == 4 {
        /* write all the values at once (contiguous bytes) */
        ffpbyt(fptr, nvals as LONGLONG * 4, cast_slice(&v), status);
    } else {
        /* have to write each value individually (not contiguous ) */
        ffpbytoff(fptr, 4, nvals, incre - 4, cast_slice(&v), status);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// put (write) the array of values to the FITS file, doing machine dependent
/// format conversion (e.g. byte-swapping) if necessary.
pub(crate) fn ffpr8b(
    fptr: &mut fitsfile, /* I - FITS file pointer                         */
    nvals: c_long,       /* I - number of pixels in the values array      */
    incre: c_long,       /* I - byte increment between pixels             */
    values: &[f64],      /* I - array of values to write                  */
    status: &mut c_int,  /* IO - error status                             */
) -> c_int {
    // Copy slice so that we don't change the original values
    let mut v: Vec<f64> = values.to_vec();

    if CFITSIO_MACHINE == VAXVMS {
        /* call VAX macro routine to convert */
        todo!();
    //ieevpd(values, values, &ii);     /* from D float -> IEEE float        */
    } else if (CFITSIO_MACHINE == ALPHAVMS) && (FLOATTYPE == GFLOAT) {
        /* convert from VMS GFLOAT float format to IEEE float format */
        for ii in 0..(nvals as usize) {
            v[ii] *= 0.25;
        }

        ffswap2(cast_slice_mut(&mut v), nvals * 4); /* swap pairs of bytes */
    } else if BYTESWAPPED {
        ffswap8(cast_slice_mut(&mut v), nvals); /* reverse order of bytes in each value */
    }

    if incre == 8 {
        /* write all the values at once (contiguous bytes) */
        ffpbyt(fptr, nvals as LONGLONG * 8, cast_slice(&v), status);
    } else {
        /* have to write each value individually (not contiguous ) */
        ffpbytoff(fptr, 8, nvals, incre - 8, cast_slice(&v), status);
    }
    *status
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BAD_ELEM_NUM, BAD_ROW_NUM, BINARY_TBL, BYTE_IMG, DOUBLE_IMG, FLOAT_IMG,
        LONG_IMG, LONGLONG, LONGLONG_IMG, NEG_FILE_POS, READONLY, READWRITE, SHORT_IMG, TBYTE,
        TDOUBLE, TFLOAT, TINT, TLONG, TLONGLONG, TSBYTE, TSHORT, fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{c_char, c_int, c_long};

    const MINDIRECT_SIZE: usize = 8640;

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Build the (ttype, tform) reference vectors that ffcrtb expects.
    fn make_cols(ttype: &[&str], tform: &[&str]) -> (Vec<Option<Vec<c_char>>>, Vec<Vec<c_char>>) {
        let tt: Vec<Option<Vec<c_char>>> = ttype.iter().map(|s| Some(cc(s))).collect();
        let tf: Vec<Vec<c_char>> = tform.iter().map(|s| cc(s)).collect();
        (tt, tf)
    }

    fn crtb(
        f: &mut fitsfile,
        tbltype: c_int,
        nrows: LONGLONG,
        ttype: &[&str],
        tform: &[&str],
        status: &mut c_int,
    ) {
        let (tt, tf) = make_cols(ttype, tform);
        let tt_ref: Vec<Option<&[c_char]>> = tt.iter().map(|o| o.as_deref()).collect();
        let tf_ref: Vec<&[c_char]> = tf.iter().map(|v| v.as_slice()).collect();
        fits_create_tbl(
            f,
            tbltype,
            nrows,
            ttype.len() as c_int,
            &tt_ref,
            &tf_ref,
            None,
            None,
            status,
        );
    }

    #[test]
    fn test_ffgrsz_image() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [100, 100];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 2, &naxes, &mut status);
            assert_eq!(status, 0);
            let mut ndata: c_long = 0;
            fits_get_rowsize(f.as_deref_mut().unwrap(), &mut ndata, &mut status);
            assert_eq!(status, 0);
            assert!(ndata > 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_ffgrsz_binary_table() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                100,
                &["ACOL", "BCOL"],
                &["1E", "1D"],
                &mut status,
            );
            assert_eq!(status, 0);
            let mut ndata: c_long = 0;
            fits_get_rowsize(f.as_deref_mut().unwrap(), &mut ndata, &mut status);
            assert_eq!(status, 0);
            assert!(ndata > 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_ffflus() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [f32; 10] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 10, &data, &mut status);
            assert_eq!(status, 0);
            fits_flush_file(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_large_write() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10000];
            let data: Vec<f32> = (0..10000).map(|i| i as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 10000, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_cross_record_write() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: Vec<f32> = (0..1000).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                &["DATA"],
                &["1000E"],
                &mut status,
            );

            for i in 0..10 {
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    1,
                    i + 1,
                    1,
                    1000,
                    &data,
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result: Vec<f32> = vec![0.0; 1000];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                5,
                1,
                1000,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[999], 1000.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_multiple_hdu_access() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [f32; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 5, &data, &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                5,
                &["COL1"],
                &["1E"],
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 5];
            let mut anynull: c_int = 0;
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 1, None, &mut status);
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_write_interleaved() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [100];
            let data: Vec<f32> = (0..100).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                50,
                &data[..50],
                &mut status,
            );
            fits_flush_file(f.as_deref_mut().unwrap(), &mut status);
            fits_write_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                51,
                50,
                &data[50..],
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 100];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                100,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[99], 100.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_table_buffering() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [1.0, 2.0, 3.0];
            let nullval = 0.0f32;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                ASCII_TBL,
                100,
                &["A", "B", "C"],
                &["F10.2", "F10.2", "F10.2"],
                &mut status,
            );

            for i in 0..100 {
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    1,
                    i + 1,
                    1,
                    1,
                    &data[0..1],
                    &mut status,
                );
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    2,
                    i + 1,
                    1,
                    1,
                    &data[1..2],
                    &mut status,
                );
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    3,
                    i + 1,
                    1,
                    1,
                    &data[2..3],
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut readval = [0f32; 1];
            let mut anynull: c_int = 0;
            let mut check = |col: c_int, row: LONGLONG, lo: f32, hi: f32, s: &mut c_int| {
                fits_read_col_flt(
                    f.as_deref_mut().unwrap(),
                    col,
                    row,
                    1,
                    1,
                    nullval,
                    &mut readval,
                    Some(&mut anynull),
                    s,
                );
                assert!(!(readval[0] < lo || readval[0] > hi));
            };
            check(1, 1, 0.99, 1.01, &mut status);
            check(2, 1, 1.99, 2.01, &mut status);
            check(3, 1, 2.99, 3.01, &mut status);
            check(1, 50, 0.99, 1.01, &mut status);
            check(2, 50, 1.99, 2.01, &mut status);
            check(1, 100, 0.99, 1.01, &mut status);
            check(3, 100, 2.99, 3.01, &mut status);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_direct_write() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE / std::mem::size_of::<f32>() + 100) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let data: Vec<f32> = (0..nvals).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f32> = vec![0.0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[nvals as usize - 1], nvals as f32);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_direct_read() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE / std::mem::size_of::<f64>() + 100) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let data: Vec<f64> = (0..nvals).map(|i| (i + 1) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                DOUBLE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_img_dbl(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f64> = vec![0.0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[nvals as usize - 1], nvals as f64);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_table_byte_io() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: [u8; 100] = std::array::from_fn(|i| i as u8);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                &["DATA"],
                &["100B"],
                &mut status,
            );

            fits_write_tblbytes(f.as_deref_mut().unwrap(), 1, 1, 100, &data, &mut status);
            fits_write_tblbytes(f.as_deref_mut().unwrap(), 5, 1, 100, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut result = [0u8; 100];
            fits_read_tblbytes(
                f.as_deref_mut().unwrap(),
                1,
                1,
                100,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);

            fits_read_tblbytes(
                f.as_deref_mut().unwrap(),
                5,
                1,
                100,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffflsh() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [100];
            let data: Vec<f32> = (0..100).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 100, &data, &mut status);

            fits_flush_buffer(f.as_deref_mut().unwrap(), false, &mut status);
            fits_flush_buffer(f.as_deref_mut().unwrap(), true, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_buffer_management() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [100];
            let data: Vec<f32> = (0..100).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 100, &data, &mut status);

            fits_flush_buffer(f.as_deref_mut().unwrap(), false, &mut status);
            fits_flush_buffer(f.as_deref_mut().unwrap(), true, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 100];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                100,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[99], 100.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_beyond_eof() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [1000];
            let data: Vec<f32> = (0..1000).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 1000, &data, &mut status);

            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                100,
                &["COL"],
                &["1E"],
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 100, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f32> = vec![0.0; 1000];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1000,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[999], 1000.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_integer_types_io() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            // Recreated three times in one test: use the "!" overwrite prefix.
            let oname = to_buf(&format!("!{filename}"));
            let naxes: [c_long; 1] = [100];
            let sdata: Vec<i16> = (0..100).map(|i| (i + 1) as i16).collect();
            let jdata: Vec<c_long> = (0..100).map(|i| (i + 1) as c_long).collect();
            let lldata: Vec<LONGLONG> = (0..100).map(|i| (i + 1) as LONGLONG).collect();

            // SHORT_IMG
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &oname, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            fits_write_img_sht(f.as_deref_mut().unwrap(), 1, 1, 100, &sdata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut sresult = [0i16; 100];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                100,
                None,
                cast_slice_mut(&mut sresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(sresult[0], 1);
            assert_eq!(sresult[99], 100);
            fits_close_file(f.take().unwrap(), &mut status);

            // LONG_IMG (32-bit)
            f = None;
            fits_create_file(&mut f, &oname, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_write_img_lng(f.as_deref_mut().unwrap(), 1, 1, 100, &jdata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut jresult = [0 as c_long; 100];
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TLONG,
                1,
                100,
                None,
                cast_slice_mut(&mut jresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(jresult[0], 1);
            assert_eq!(jresult[99], 100);
            fits_close_file(f.take().unwrap(), &mut status);

            // LONGLONG_IMG (64-bit)
            f = None;
            fits_create_file(&mut f, &oname, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                LONGLONG_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_img_lnglng(f.as_deref_mut().unwrap(), 1, 1, 100, &lldata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut llresult = [0 as LONGLONG; 100];
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                1,
                100,
                None,
                cast_slice_mut(&mut llresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(llresult[0], 1);
            assert_eq!(llresult[99], 100);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_strided_io() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let fdata: Vec<f32> = (0..10).map(|i| (i + 1) as f32).collect();
            let ddata: Vec<f64> = (0..10).map(|i| (i + 1) as f64).collect();
            let jdata: Vec<c_int> = (0..10).map(|i| i + 1).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                &["A", "B", "C"],
                &["1E", "1D", "1J"],
                &mut status,
            );

            for i in 0..10 {
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    1,
                    i + 1,
                    1,
                    1,
                    &fdata[i as usize..i as usize + 1],
                    &mut status,
                );
                fits_write_col_dbl(
                    f.as_deref_mut().unwrap(),
                    2,
                    i + 1,
                    1,
                    1,
                    &ddata[i as usize..i as usize + 1],
                    &mut status,
                );
                fits_write_col_int(
                    f.as_deref_mut().unwrap(),
                    3,
                    i + 1,
                    1,
                    1,
                    &jdata[i as usize..i as usize + 1],
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut fresult = [0f32; 10];
            let mut dresult = [0f64; 10];
            let mut jresult = [0 as c_int; 10];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1,
                1,
                10,
                None,
                cast_slice_mut(&mut fresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                2,
                1,
                1,
                10,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TINT,
                3,
                1,
                1,
                10,
                None,
                cast_slice_mut(&mut jresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);

            for i in 0..10 {
                assert_eq!(fresult[i], fdata[i]);
                assert_eq!(dresult[i], ddata[i]);
                assert_eq!(jresult[i], jdata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_io() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [200];
            let data: [u8; 200] = std::array::from_fn(|i| (i % 256) as u8);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img_byt(f.as_deref_mut().unwrap(), 1, 1, 200, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 200];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                200,
                None,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_multi_record_spanning() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: Vec<f32> = (0..3000).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                5,
                &["BIGCOL"],
                &["3000E"],
                &mut status,
            );

            for i in 0..5 {
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    1,
                    i + 1,
                    1,
                    3000,
                    &data,
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut result: Vec<f32> = vec![0.0; 3000];
            let mut anynull: c_int = 0;
            for i in 0..5 {
                fits_read_col(
                    f.as_deref_mut().unwrap(),
                    TFLOAT,
                    1,
                    i + 1,
                    1,
                    3000,
                    None,
                    cast_slice_mut(&mut result),
                    Some(&mut anynull),
                    &mut status,
                );
                assert_eq!(result[0], 1.0);
                assert_eq!(result[2999], 3000.0);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_very_large_image() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5000];
            let data: Vec<f64> = (0..5000).map(|i| (i + 1) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                DOUBLE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_img_dbl(f.as_deref_mut().unwrap(), 1, 1, 5000, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f64> = vec![0.0; 5000];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                5000,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[4999], 5000.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffflus_table() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: [f32; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                5,
                &["COL1"],
                &["1E"],
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_flush_file(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_table_large() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                ASCII_TBL,
                200,
                &["VALUE"],
                &["E15.7"],
                &mut status,
            );

            for i in 0..200 {
                let data = [(i + 1) as f32];
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    1,
                    i + 1,
                    1,
                    1,
                    &data,
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut result = [0f32; 1];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(result[0], 1.0);
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                200,
                1,
                1,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(result[0], 200.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_auto_expand_table() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data = [42.0f32];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                0,
                &["DATA"],
                &["1E"],
                &mut status,
            );

            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 100, 1, 1, &data, &mut status);

            let mut nrows: LONGLONG = 0;
            fits_get_num_rowsll(f.as_deref_mut().unwrap(), &mut nrows, &mut status);
            assert_eq!(nrows, 100);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 1];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                100,
                1,
                1,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(result[0], 42.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_noncontiguous_column_io() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let fdata: Vec<f32> = (0..50).map(|i| (i + 1) as f32).collect();
            let jdata: Vec<c_int> = (0..50).map(|i| i + 1).collect();
            let ddata: Vec<f64> = (0..50).map(|i| (i + 1) as f64).collect();
            let idata: Vec<i16> = (0..50).map(|i| (i + 1) as i16).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                50,
                &["A", "B", "C", "D"],
                &["1E", "1J", "1D", "1I"],
                &mut status,
            );

            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 50, &fdata, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 2, 1, 1, 50, &jdata, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 3, 1, 1, 50, &ddata, &mut status);
            fits_write_col_sht(f.as_deref_mut().unwrap(), 4, 1, 1, 50, &idata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut fresult = [0f32; 50];
            let mut jresult = [0 as c_int; 50];
            let mut dresult = [0f64; 50];
            let mut iresult = [0i16; 50];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1,
                1,
                50,
                None,
                cast_slice_mut(&mut fresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TINT,
                2,
                1,
                1,
                50,
                None,
                cast_slice_mut(&mut jresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                3,
                1,
                1,
                50,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TSHORT,
                4,
                1,
                1,
                50,
                None,
                cast_slice_mut(&mut iresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);

            for i in 0..50 {
                assert_eq!(fresult[i], fdata[i]);
                assert_eq!(jresult[i], jdata[i]);
                assert_eq!(dresult[i], ddata[i]);
                assert_eq!(iresult[i], idata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_noncontiguous_spanning_records() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let ddata: Vec<f64> = (0..200).map(|i| (i + 1) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                200,
                &["PAD", "DATA"],
                &["492B", "1D"],
                &mut status,
            );

            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 200, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut dresult: Vec<f64> = vec![0.0; 200];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                2,
                1,
                1,
                200,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            for i in 0..200 {
                assert_eq!(dresult[i], ddata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_column_noncontiguous() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let bdata: [u8; 100] = std::array::from_fn(|i| (i % 256) as u8);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                100,
                &["PAD", "BYTES"],
                &["10B", "1B"],
                &mut status,
            );

            fits_write_col_byt(f.as_deref_mut().unwrap(), 2, 1, 1, 100, &bdata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut bresult = [0u8; 100];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TBYTE,
                2,
                1,
                1,
                100,
                None,
                &mut bresult,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(bresult, bdata);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_element_spanning_record() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let ddata: Vec<f64> = (0..50).map(|i| (i + 1) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                50,
                &["PAD", "DATA"],
                &["2862B", "1D"],
                &mut status,
            );

            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 50, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut dresult: Vec<f64> = vec![0.0; 50];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                2,
                1,
                1,
                50,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            for i in 0..50 {
                assert_eq!(dresult[i], ddata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_large_direct_write() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                ASCII_TBL,
                500,
                &["VALUE"],
                &["E20.10"],
                &mut status,
            );

            for i in 0..500 {
                let data = [(i + 1) as f32];
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    1,
                    i + 1,
                    1,
                    1,
                    &data,
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut result = [0f32; 1];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(result[0], 1.0);
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                500,
                1,
                1,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(result[0], 500.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_massive_direct_io() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals: LONGLONG = 20000;
            let naxes: [c_long; 1] = [20000];
            let data: Vec<f32> = (0..nvals).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f32> = vec![0.0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[nvals as usize - 1], nvals as f32);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_element_at_boundary() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let ddata: Vec<f64> = (0..30).map(|i| (i + 1) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                30,
                &["PAD", "VAL"],
                &["2876B", "1D"],
                &mut status,
            );

            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 30, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut dresult: Vec<f64> = vec![0.0; 30];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                2,
                1,
                1,
                30,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            for i in 0..30 {
                assert_eq!(dresult[i], ddata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_at_boundary() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let fdata: Vec<f32> = (0..40).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                40,
                &["PAD", "VAL"],
                &["2878B", "1E"],
                &mut status,
            );

            fits_write_col_flt(f.as_deref_mut().unwrap(), 2, 1, 1, 40, &fdata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut fresult: Vec<f32> = vec![0.0; 40];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                2,
                1,
                1,
                40,
                None,
                cast_slice_mut(&mut fresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            for i in 0..40 {
                assert_eq!(fresult[i], fdata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_at_boundary() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let sdata: Vec<i16> = (0..50).map(|i| (i + 1) as i16).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                50,
                &["PAD", "VAL"],
                &["2879B", "1I"],
                &mut status,
            );

            fits_write_col_sht(f.as_deref_mut().unwrap(), 2, 1, 1, 50, &sdata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut sresult: Vec<i16> = vec![0; 50];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TSHORT,
                2,
                1,
                1,
                50,
                None,
                cast_slice_mut(&mut sresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            for i in 0..50 {
                assert_eq!(sresult[i], sdata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_single_element_spanning() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let ddata = [42.5f64];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["PAD", "VAL"],
                &["2874B", "1D"],
                &mut status,
            );

            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 1, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut dresult = [0f64; 1];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                2,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(dresult[0], ddata[0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_existing_after_partial_write() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [1000];
            let mut data: Vec<f32> = (0..1000).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 1000, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READWRITE, &mut status);
            data[500] = 999.0;
            fits_write_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                501,
                1,
                &data[500..501],
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f32> = vec![0.0; 1000];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1000,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[500], 999.0);
            assert_eq!(result[999], 1000.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_image_read() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE / std::mem::size_of::<f32>() + 1000) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let data: Vec<f32> = (0..nvals).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f32> = vec![0.0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[nvals as usize - 1], nvals as f32);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_short_image() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE / std::mem::size_of::<i16>() + 1000) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let data: Vec<i16> = (0..nvals).map(|i| ((i + 1) % 30000) as i16).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            fits_write_img_sht(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<i16> = vec![0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[nvals as usize - 1], (nvals % 30000) as i16);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_int_image() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE / std::mem::size_of::<c_int>() + 1000) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let data: Vec<c_int> = (0..nvals).map(|i| (i + 1) as c_int).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_write_img_int(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<c_int> = vec![0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TINT,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[nvals as usize - 1], nvals as c_int);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_double_image() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE / std::mem::size_of::<f64>() + 1000) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let data: Vec<f64> = (0..nvals).map(|i| (i + 1) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                DOUBLE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_img_dbl(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f64> = vec![0.0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[nvals as usize - 1], nvals as f64);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_overwrite() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE / std::mem::size_of::<f32>() + 1000) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let mut data: Vec<f32> = (0..nvals).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            for (i, v) in data.iter_mut().enumerate() {
                *v = (i as LONGLONG + 10001) as f32;
            }

            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<f32> = vec![0.0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                nvals,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10001.0);
            assert_eq!(result[nvals as usize - 1], (nvals + 10000) as f32);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_byte_image() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let nvals = (MINDIRECT_SIZE + 1000) as LONGLONG;
            let naxes: [c_long; 1] = [nvals as c_long];
            let data: Vec<u8> = (0..nvals).map(|i| (i % 256) as u8).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img_byt(f.as_deref_mut().unwrap(), 1, 1, nvals, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: Vec<u8> = vec![0; nvals as usize];
            let mut anynull: c_int = 0;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                nvals,
                None,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 0);
            assert_eq!(result[255], 255);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_mixed_column_types() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let bdata = [255u8];
            let sdata = [32000i16];
            let idata = [2000000 as c_int];
            let ldata = [4000000000 as c_long];
            let fdata = [3.14159f32];
            let ddata = [2.71828182845f64];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["BYTE", "SHORT", "INT", "LONG", "FLOAT", "DOUBLE"],
                &["1B", "1I", "1J", "1K", "1E", "1D"],
                &mut status,
            );

            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &bdata, &mut status);
            fits_write_col_sht(f.as_deref_mut().unwrap(), 2, 1, 1, 1, &sdata, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 3, 1, 1, 1, &idata, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 4, 1, 1, 1, &ldata, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 5, 1, 1, 1, &fdata, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 6, 1, 1, 1, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut bresult = [0u8; 1];
            let mut sresult = [0i16; 1];
            let mut iresult = [0 as c_int; 1];
            let mut lresult = [0 as c_long; 1];
            let mut fresult = [0f32; 1];
            let mut dresult = [0f64; 1];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                1,
                1,
                1,
                None,
                &mut bresult,
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TSHORT,
                2,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut sresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TINT,
                3,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut iresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TLONG,
                4,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut lresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                5,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut fresult),
                Some(&mut anynull),
                &mut status,
            );
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                6,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);

            assert_eq!(bresult[0], bdata[0]);
            assert_eq!(sresult[0], sdata[0]);
            assert_eq!(iresult[0], idata[0]);
            assert_eq!(lresult[0], ldata[0]);
            assert_eq!(fresult[0], fdata[0]);
            assert_eq!(dresult[0], ddata[0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_delete_rows_truncate() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: Vec<f32> = (0..1000).map(|i| (i + 1) as f32).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                100,
                &["DATA"],
                &["1000E"],
                &mut status,
            );

            for i in 0..100 {
                fits_write_col_flt(
                    f.as_deref_mut().unwrap(),
                    1,
                    i + 1,
                    1,
                    1000,
                    &data,
                    &mut status,
                );
            }
            fits_flush_file(f.as_deref_mut().unwrap(), &mut status);

            // Read some data to load buffers.
            let mut anynull: c_int = 0;
            for i in 0..50 {
                let mut result: Vec<f32> = vec![0.0; 1000];
                fits_read_col(
                    f.as_deref_mut().unwrap(),
                    TFLOAT,
                    1,
                    i + 50,
                    1,
                    1000,
                    None,
                    cast_slice_mut(&mut result),
                    Some(&mut anynull),
                    &mut status,
                );
            }

            // Now delete many rows to trigger truncation and ffbfeof.
            fits_delete_rows(f.as_deref_mut().unwrap(), 1, 80, &mut status);

            let mut nrows: LONGLONG = 0;
            fits_get_num_rowsll(f.as_deref_mut().unwrap(), &mut nrows, &mut status);
            assert_eq!(nrows, 20);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffptbb_auto_extend_with_following_hdu() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 100] = std::array::from_fn(|i| (i % 256) as u8);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                &["DATA"],
                &["100B"],
                &mut status,
            );

            fits_write_tblbytes(f.as_deref_mut().unwrap(), 1, 1, 100, &data, &mut status);

            // Add another HDU after the table so the table is NOT the last HDU.
            fits_create_img(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);

            // Move back to the table.
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            // Now write to row 20 to trigger auto-extend with ffirow.
            fits_write_tblbytes(f.as_deref_mut().unwrap(), 20, 1, 100, &data, &mut status);

            let mut nrows: LONGLONG = 0;
            fits_get_num_rowsll(f.as_deref_mut().unwrap(), &mut nrows, &mut status);
            assert_eq!(nrows, 20);

            let mut result = [0u8; 100];
            fits_read_tblbytes(
                f.as_deref_mut().unwrap(),
                20,
                1,
                100,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffptbb_auto_extend_last_hdu() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: [u8; 100] = std::array::from_fn(|i| (i % 256) as u8);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                5,
                &["DATA"],
                &["100B"],
                &mut status,
            );

            fits_write_tblbytes(f.as_deref_mut().unwrap(), 1, 1, 100, &data, &mut status);
            fits_write_tblbytes(f.as_deref_mut().unwrap(), 50, 1, 100, &data, &mut status);

            let mut nrows: LONGLONG = 0;
            fits_get_num_rowsll(f.as_deref_mut().unwrap(), &mut nrows, &mut status);
            assert_eq!(nrows, 50);

            let mut result = [0u8; 100];
            fits_read_tblbytes(
                f.as_deref_mut().unwrap(),
                50,
                1,
                100,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_table_byte_io_errors() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut data = [0u8; 100];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                &["DATA"],
                &["100B"],
                &mut status,
            );

            // BAD_ROW_NUM for ffgtbb (row 0).
            status = 0;
            fits_read_tblbytes(f.as_deref_mut().unwrap(), 0, 1, 100, &mut data, &mut status);
            assert_eq!(status, BAD_ROW_NUM);

            // BAD_ELEM_NUM for ffgtbb (element 0).
            status = 0;
            fits_read_tblbytes(f.as_deref_mut().unwrap(), 1, 0, 100, &mut data, &mut status);
            assert_eq!(status, BAD_ELEM_NUM);

            // BAD_ROW_NUM for ffptbb (row 0).
            status = 0;
            fits_write_tblbytes(f.as_deref_mut().unwrap(), 0, 1, 100, &data, &mut status);
            assert_eq!(status, BAD_ROW_NUM);

            // BAD_ELEM_NUM for ffptbb (element 0).
            status = 0;
            fits_write_tblbytes(f.as_deref_mut().unwrap(), 1, 0, 100, &data, &mut status);
            assert_eq!(status, BAD_ELEM_NUM);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_past_table_end() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut data = [0u8; 100];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                &["DATA"],
                &["100B"],
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            // Try to read row 100 of a 10-row table.
            status = 0;
            fits_read_tblbytes(
                f.as_deref_mut().unwrap(),
                100,
                1,
                100,
                &mut data,
                &mut status,
            );
            assert_eq!(status, BAD_ROW_NUM);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_last_element_span_boundary() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let ddata = [123.456f64];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["PAD1", "PAD2", "VAL"],
                &["2800B", "76B", "1D"],
                &mut status,
            );

            fits_write_col_dbl(f.as_deref_mut().unwrap(), 3, 1, 1, 1, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut dresult = [0f64; 1];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                3,
                1,
                1,
                1,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(dresult[0], ddata[0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_multi_element_last_spans() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let ddata: Vec<f64> = (0..3).map(|i| (i + 1) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                3,
                &["PAD", "VAL"],
                &["2876B", "1D"],
                &mut status,
            );

            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 3, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut dresult = [0f64; 3];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                2,
                1,
                1,
                3,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            for i in 0..3 {
                assert_eq!(dresult[i], ddata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_buffered_write_wrap_last_chunk() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: Vec<u8> = (0..6000).map(|i| (i % 256) as u8).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                2,
                &["DATA"],
                &["6000B"],
                &mut status,
            );

            fits_write_tblbytes(f.as_deref_mut().unwrap(), 1, 1, 6000, &data, &mut status);
            fits_write_tblbytes(f.as_deref_mut().unwrap(), 2, 1, 6000, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result: Vec<u8> = vec![0; 6000];
            fits_read_tblbytes(
                f.as_deref_mut().unwrap(),
                1,
                1,
                6000,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_neg_file_pos() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [i8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                1,
                10,
                cast_slice(&data),
                &mut status,
            );
            assert_eq!(status, 0);

            // Try to move to negative position - should fail.
            // IGNORE_EOF = 1, REPORT_EOF = 0
            ffmbyt_safe(f.as_deref_mut().unwrap(), -1, 1, &mut status);
            assert_eq!(status, NEG_FILE_POS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_spanning_records() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: Vec<u8> = (0..10000).map(|i| (i % 256) as u8).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["BIGCOL"],
                &["10000B"],
                &mut status,
            );

            fits_write_tblbytes(f.as_deref_mut().unwrap(), 1, 1, 10000, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result: Vec<u8> = vec![0; 10000];
            fits_read_tblbytes(
                f.as_deref_mut().unwrap(),
                1,
                1,
                10000,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_last_element_exact_span() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let ddata: Vec<f64> = (0..4).map(|i| (i + 100) as f64).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            crtb(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                4,
                &["PAD", "VAL"],
                &["713B", "1D"],
                &mut status,
            );

            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 4, &ddata, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut dresult = [0f64; 4];
            let mut anynull: c_int = 0;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                2,
                1,
                1,
                4,
                None,
                cast_slice_mut(&mut dresult),
                Some(&mut anynull),
                &mut status,
            );
            for i in 0..4 {
                assert_eq!(dresult[i], ddata[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
