/*  This file, checksum.c, contains the checksum-related routines in the   */
/*  FITSIO library.                                                        */
/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
/*------------------------------------------------------------------------*/

use core::slice;
use std::cmp;

use crate::c_types::{c_char, c_int, c_long, c_uint, c_ulong};

use bytemuck::{cast_slice, cast_slice_mut};

use crate::cs;
use crate::fitscore::ffpmsg_str;
use crate::fitscore::{ffghadll_safe, ffpdfl, ffrdef_safe, ffuptf, ffwend};
use crate::fitsio::*;
use crate::fitsio2::*;
use crate::getkey::ffgkys_safe;
use crate::modkey::ffmkys_safe;
use crate::putkey::{ffgstm_safe, ffpkys_safe};
use crate::swapproc::ffswap2;
use crate::wrappers::*;
use crate::{buffers::*, int_snprintf};

/*--------------------------------------------------------------------------*/
/// Calculate a 32-bit 1's complement checksum of the FITS 2880-byte blocks.
///
/// This routine is based on the C algorithm developed by Rob
/// Seaman at NOAO that was presented at the 1994 ADASS conference,  
/// published in the Astronomical Society of the Pacific Conference Series.
/// This uses a 32-bit 1's complement checksum in which the overflow bits
/// are permuted back into the sum and therefore all bit positions are
/// sampled evenly.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcsum(
    fptr: *mut fitsfile, /* I - FITS file pointer                  */
    nrec: c_long,        /* I - number of 2880-byte blocks to sum  */
    sum: *mut c_ulong,   /* IO - accumulated checksum              */
    status: *mut c_int,  /* IO - error status                      */
) -> c_int {
    unsafe {
        if fptr.is_null() {
            *status = NULL_INPUT_PTR;
            return *status;
        }

        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let sum = sum.as_mut().expect(NULL_MSG);

        ffcsum_safe(fptr, nrec, sum, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Calculate a 32-bit 1's complement checksum of the FITS 2880-byte blocks.
/// This routine is based on the C algorithm developed by Rob
/// Seaman at NOAO that was presented at the 1994 ADASS conference,  
/// published in the Astronomical Society of the Pacific Conference Series.
/// This uses a 32-bit 1's complement checksum in which the overflow bits
/// are permuted back into the sum and therefore all bit positions are
/// sampled evenly.
pub fn ffcsum_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                  */
    nrec: c_long,        /* I - number of 2880-byte blocks to sum  */
    sum: &mut c_ulong,   /* IO - accumulated checksum              */
    status: &mut c_int,  /* IO - error status                      */
) -> c_int {
    let mut ii: usize;
    let mut sbuf: [u16; 1440] = [0; 1440];
    let mut hi: c_ulong;
    let mut lo: c_ulong;
    let mut hicarry: c_ulong;
    let mut locarry: c_ulong;

    let nrec = nrec as usize;

    if *status > 0 {
        return *status;
    }
    /*
      Sum the specified number of FITS 2880-byte records.  This assumes that
      the FITSIO file pointer points to the start of the records to be summed.
      Read each FITS block as 1440 short values (do byte swapping if needed).
    */
    for jj in 0..nrec {
        ffgbyt(fptr, 2880, cast_slice_mut(&mut sbuf), status);
        if BYTESWAPPED {
            /* reverse order of bytes in each value */
            ffswap2(cast_slice_mut(&mut sbuf), 1440);
        }

        hi = *sum >> 16;
        lo = *sum & 0xFFFF;
        ii = 0;
        while ii < 1440 {
            hi += sbuf[ii] as c_ulong;
            lo += sbuf[ii + 1] as c_ulong;
            ii += 2
        }

        hicarry = hi >> 16; /* fold carry bits in */
        locarry = lo >> 16;

        while hicarry | locarry != 0 {
            hi = (hi & 0xFFFF) + locarry;
            lo = (lo & 0xFFFF) + hicarry;
            hicarry = hi >> 16;
            locarry = lo >> 16;
        }
        *sum = (hi << 16) + lo;
    }
    *status
}

/*-------------------------------------------------------------------------*/
/// encode the 32 bit checksum by converting every
/// 2 bits of each byte into an ASCII character (32 bit word encoded
/// as 16 character string).   Only ASCII letters and digits are used
/// to encode the values (no ASCII punctuation characters).
///
/// If complm=TRUE, then the complement of the sum will be encoded.
///
/// This routine is based on the C algorithm developed by Rob
/// Seaman at NOAO that was presented at the 1994 ADASS conference,
/// published in the Astronomical Society of the Pacific Conference Series.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffesum(
    sum: c_ulong,             /* I - accumulated checksum                */
    complm: c_int,            /* I - = 1 to encode complement of the sum */
    ascii: *mut [c_char; 17], /* O - 16-char ASCII encoded checksum      */
) {
    unsafe {
        let ascii = ascii.as_mut().expect(NULL_MSG);
        ffesum_safe(sum, complm == 1, ascii)
    }
}

/*-------------------------------------------------------------------------*/
/// encode the 32 bit checksum by converting every
/// 2 bits of each byte into an ASCII character (32 bit word encoded
/// as 16 character string).   Only ASCII letters and digits are used
/// to encode the values (no ASCII punctuation characters).
///
/// If complm=TRUE, then the complement of the sum will be encoded.
///
/// This routine is based on the C algorithm developed by Rob
/// Seaman at NOAO that was presented at the 1994 ADASS conference,
/// published in the Astronomical Society of the Pacific Conference Series.
pub fn ffesum_safe(
    sum: c_ulong,             /* I - accumulated checksum                */
    complm: bool,             /* I - = 1 to encode complement of the sum */
    ascii: &mut [c_char; 17], /* O - 16-char ASCII encoded checksum      */
) {
    let exclude: [c_uint; 13] = [
        0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 0x40, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60,
    ];
    let mask: [c_ulong; 4] = [0xff000000, 0xff0000, 0xff00, 0xff];

    let offset = 0x30; /* ASCII 0 (zero) */

    let mut byte: c_int;
    let mut quotient: c_int;
    let mut remainder: c_int;
    let mut ch: [c_int; 4] = [0; 4];
    let mut check: c_int;
    let mut asc: [c_char; 32] = [0; 32];

    let value: c_ulong = if complm { 0xFFFFFFFF - sum } else { sum };

    for ii in 0..4 {
        byte = ((value & mask[ii]) >> (24 - (8 * ii))) as c_int;
        quotient = byte / 4 + offset;
        remainder = byte % 4;

        for jj in 0..4 {
            ch[jj] = quotient;
        }
        ch[0] += remainder;

        /* avoid ASCII  punctuation */
        check = 1;
        while check != 0 {
            check = 0;
            for kk in 0..13 {
                for jj in (0..4).step_by(2) {
                    if (ch[jj] as u8) == exclude[kk] as u8
                        || (ch[jj + 1] as u8) == exclude[kk] as u8
                    {
                        ch[jj] += 1;
                        ch[jj + 1] -= 1;
                        check += 1;
                    }
                }
            }
        }

        /* assign the bytes */
        for jj in 0..4 {
            asc[4 * jj + ii] = ch[jj] as c_char;
        }
    }

    /* shift the bytes 1 to the right */
    for ii in 0..16 {
        ascii[ii] = asc[(ii + 15) % 16];
    }

    ascii[16] = 0;
}

/*-------------------------------------------------------------------------*/
/// decode the 16-char ASCII encoded checksum into an unsigned 32-bit long.
/// If complm=TRUE, then the complement of the sum will be decoded.
///
/// This routine is based on the C algorithm developed by Rob
/// Seaman at NOAO that was presented at the 1994 ADASS conference,
/// published in the Astronomical Society of the Pacific Conference Series.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdsum(
    ascii: *const c_char, /* I - 16-char ASCII encoded checksum   */
    complm: c_int,        /* I - =1 to decode complement of the   */
    sum: *mut c_ulong,    /* O - 32-bit checksum           */
) -> c_ulong {
    unsafe {
        let sum = sum.as_mut().expect(NULL_MSG);
        let ascii = slice::from_raw_parts(ascii, 17).try_into().unwrap();

        ffdsum_safe(ascii, complm == 1, sum)
    }
}

/*-------------------------------------------------------------------------*/
/// decode the 16-char ASCII encoded checksum into an unsigned 32-bit long.
/// If complm=TRUE, then the complement of the sum will be decoded.
/// This routine is based on the C algorithm developed by Rob
/// Seaman at NOAO that was presented at the 1994 ADASS conference,
/// published in the Astronomical Society of the Pacific Conference Series.
pub fn ffdsum_safe(
    ascii: &[c_char; 17], /* I - 16-char ASCII encoded checksum   */
    complm: bool,         /* I - =1 to decode complement of the   */
    sum: &mut c_ulong,    /* O - 32-bit checksum           */
) -> c_ulong {
    let mut cbuf: [c_char; 16] = [0; 16];
    let mut hi: c_ulong = 0;
    let mut lo: c_ulong = 0;

    /* remove the permuted FITS byte alignment and the ASCII 0 offset */
    for ii in 0..16 {
        cbuf[ii] = ascii[(ii + 1) % 16];
        cbuf[ii] -= 48;
    }

    for ii in (0..16).step_by(4) {
        hi += ((cbuf[ii] as c_ulong) << 8) + (cbuf[ii + 1] as c_ulong);
        lo += ((cbuf[ii + 2] as c_ulong) << 8) + (cbuf[ii + 3] as c_ulong);
    }

    let mut hicarry = hi >> 16;
    let mut locarry = lo >> 16;
    while hicarry != 0 || locarry != 0 {
        hi = (hi & 0xFFFF) + locarry;
        lo = (lo & 0xFFFF) + hicarry;
        hicarry = hi >> 16;
        locarry = lo >> 16;
    }
    *sum = (hi << 16) + lo;
    if complm {
        *sum = 0xFFFFFFFF - *sum; /* complement each bit of the value */
    }
    *sum
}

/*------------------------------------------------------------------------*/
/// Create or update the checksum keywords in the CHDU.  These keywords
/// provide a checksum verification of the FITS HDU based on the ASCII
/// coded 1's complement checksum algorithm developed by Rob Seaman at NOAO.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcks(
    fptr: *mut fitsfile, /* I - FITS file pointer                  */
    status: *mut c_int,  /* IO - error status                      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffpcks_safe(fptr, status)
    }
}

/*------------------------------------------------------------------------*/
/// Create or update the checksum keywords in the CHDU.  These keywords
/// provide a checksum verification of the FITS HDU based on the ASCII
/// coded 1's complement checksum algorithm developed by Rob Seaman at NOAO.
pub fn ffpcks_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                  */
    status: &mut c_int,  /* IO - error status                      */
) -> c_int {
    let mut datestr: [c_char; 20] = [0; 20];
    let mut checksum: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut datasum: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut chkcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut datacomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut headstart: LONGLONG = 0;
    let mut datastart: LONGLONG = 0;
    let mut dataend: LONGLONG = 0;
    let mut dsum: c_ulong = 0;
    let mut olddsum: c_ulong = 0;
    let mut sum: c_ulong = 0;
    let mut tdouble: f64 = 0.0;

    /* inherit input status value if > 0 */
    if *status > 0 {
        return *status;
    }

    /* generate current date string and construct the keyword comments */
    ffgstm_safe(&mut datestr, None, status);
    strcpy_safe(&mut chkcomm, cs!(c"HDU checksum updated "));
    strcat_safe(&mut chkcomm, &datestr);
    strcpy_safe(&mut datacomm, cs!(c"data unit checksum updated "));
    strcat_safe(&mut datacomm, &datestr);

    /* write the CHECKSUM keyword if it does not exist */
    let mut tstatus = *status;
    if ffgkys_safe(
        fptr,
        cs!(c"CHECKSUM"),
        &mut checksum,
        Some(&mut comm),
        status,
    ) == KEY_NO_EXIST
    {
        *status = tstatus;
        strcpy_safe(&mut checksum, cs!(c"0000000000000000"));
        ffpkys_safe(fptr, cs!(c"CHECKSUM"), &checksum, Some(&chkcomm), status);
    }

    /* write the DATASUM keyword if it does not exist */
    tstatus = *status;
    if ffgkys_safe(fptr, cs!(c"DATASUM"), &mut datasum, Some(&mut comm), status) == KEY_NO_EXIST {
        *status = tstatus;
        olddsum = 0;
        ffpkys_safe(
            fptr,
            cs!(c"DATASUM"),
            cs!(c"         0"),
            Some(&datacomm),
            status,
        );

        /* set the CHECKSUM keyword as undefined, if it isn't already */
        if strcmp_safe(&checksum, cs!(c"0000000000000000")) != 0 {
            strcpy_safe(&mut checksum, cs!(c"0000000000000000"));
            ffmkys_safe(fptr, cs!(c"CHECKSUM"), &checksum, Some(&chkcomm), status);
        };
    } else {
        /* decode the datasum into an unsigned long variable */
        /* olddsum = strtoul(datasum, 0, 10); doesn't work on SUN OS */
        tdouble = atof_safe(&datasum);
        olddsum = tdouble as c_ulong;
    }

    /* close header: rewrite END keyword and following blank fill */
    /* and re-read the required keywords to determine the structure */
    if ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    /* update the variable length TFORM values */
    if fptr.Fptr.heapsize > 0 {
        ffuptf(fptr, status);
    }

    /* write the correct data fill values, if they are not already correct */
    if ffpdfl(fptr, status) > 0 {
        return *status;
    }

    /* calc size of data unit, in FITS 2880-byte blocks */
    if ffghadll_safe(
        fptr,
        Some(&mut headstart),
        Some(&mut datastart),
        Some(&mut dataend),
        status,
    ) > 0
    {
        return *status;
    }

    let mut nrec = ((dataend - datastart) / 2880i64) as c_long;

    if nrec > 0 {
        /* accumulate the 32-bit 1's complement checksum */
        ffmbyt_safe(fptr, datastart, REPORT_EOF, status);
        if ffcsum_safe(fptr, nrec, &mut dsum, status) > 0 {
            return *status;
        };
    }

    if dsum != olddsum {
        /* update the DATASUM keyword with the correct value */
        int_snprintf!(&mut datasum, FLEN_VALUE, "{}", dsum,);
        ffmkys_safe(fptr, cs!(c"DATASUM"), &datasum, Some(&datacomm), status);

        /* set the CHECKSUM keyword as undefined, if it isn't already */
        if strcmp_safe(&checksum, cs!(c"0000000000000000")) != 0 {
            strcpy_safe(&mut checksum, cs!(c"0000000000000000"));
            ffmkys_safe(fptr, cs!(c"CHECKSUM"), &checksum, Some(&chkcomm), status);
        };
    }

    if strcmp_safe(&checksum, cs!(c"0000000000000000")) != 0 {
        /* check if CHECKSUM is still OK; move to the start of the header */
        ffmbyt_safe(fptr, headstart, REPORT_EOF, status);

        /* accumulate the header checksum into the previous data checksum */
        nrec = ((datastart - headstart) / 2880) as c_long;
        sum = dsum;
        if ffcsum_safe(fptr, nrec, &mut sum, status) > 0 {
            return *status;
        }
        if sum == 0 || sum == 0xFFFFFFFF {
            /* CHECKSUM is correct */
            return *status;
        }

        /* Zero the CHECKSUM and recompute the new value */
        ffmkys_safe(
            fptr,
            cs!(c"CHECKSUM"),
            cs!(c"0000000000000000"),
            Some(&chkcomm),
            status,
        );
    }

    /* move to the start of the header */
    ffmbyt_safe(fptr, headstart, REPORT_EOF, status);

    /* accumulate the header checksum into the previous data checksum */
    nrec = ((datastart - headstart) / 2880i64) as c_long;
    sum = dsum;
    if ffcsum_safe(fptr, nrec, &mut sum, status) > 0 {
        return *status;
    }

    /* encode the COMPLEMENT of the checksum into a 16-character string */
    ffesum_safe(sum, true, (&mut checksum[..17]).try_into().unwrap());

    /* update the CHECKSUM keyword value with the new string */
    ffmkys_safe(fptr, cs!(c"CHECKSUM"), &checksum, Some(cs!(c"&")), status);

    *status
}

/*------------------------------------------------------------------------*/
/// Update the CHECKSUM keyword value.  This assumes that the DATASUM
/// keyword exists and has the correct value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffupck(
    fptr: *mut fitsfile, /* I - FITS file pointer                  */
    status: *mut c_int,  /* IO - error status                      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffupck_safe(fptr, status)
    }
}

/*------------------------------------------------------------------------*/
/// Update the CHECKSUM keyword value.  This assumes that the DATASUM
/// keyword exists and has the correct value.
pub fn ffupck_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                  */
    status: &mut c_int,  /* IO - error status                      */
) -> c_int {
    let mut datestr: [c_char; 20] = [0; 20];
    let mut chkcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut checksum: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut datasum: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut nrec: c_long = 0;
    let mut headstart: LONGLONG = 0;
    let mut datastart: LONGLONG = 0;
    let mut dataend: LONGLONG = 0;
    let mut sum: c_ulong = 0;

    /* inherit input status value if > 0 */
    if *status > 0 {
        return *status;
    }

    /* generate current date string and construct the keyword comments */
    ffgstm_safe(&mut datestr, None, status);
    strcpy_safe(&mut chkcomm, cs!(c"HDU checksum updated "));
    strcat_safe(&mut chkcomm, &datestr);

    /* get the DATASUM keyword and convert it to a unsigned long */
    if ffgkys_safe(fptr, cs!(c"DATASUM"), &mut datasum, Some(&mut comm), status) == KEY_NO_EXIST {
        ffpmsg_str("DATASUM keyword not found (ffupck");
        return *status;
    }

    /* read as a double as a workaround */
    let tdouble = atof_safe(&datasum);
    let dsum = tdouble as c_ulong;

    /* get size of the HDU */
    if ffghadll_safe(
        fptr,
        Some(&mut headstart),
        Some(&mut datastart),
        Some(&mut dataend),
        status,
    ) > 0
    {
        return *status;
    }

    /* get the checksum keyword, if it exists */
    let tstatus = *status;
    if ffgkys_safe(
        fptr,
        cs!(c"CHECKSUM"),
        &mut checksum,
        Some(&mut comm),
        status,
    ) == KEY_NO_EXIST
    {
        *status = tstatus;
        strcpy_safe(&mut checksum, cs!(c"0000000000000000"));
        ffpkys_safe(fptr, cs!(c"CHECKSUM"), &checksum, Some(&chkcomm), status);
    } else {
        /* check if CHECKSUM is still OK */
        /* rewrite END keyword and following blank fill */
        if ffwend(fptr, status) > 0 {
            return *status;
        }

        /* move to the start of the header */
        ffmbyt_safe(fptr, headstart, REPORT_EOF, status);

        /* accumulate the header checksum into the previous data checksum */
        nrec = ((datastart - headstart) / 2880) as c_long;
        sum = dsum;
        if ffcsum_safe(fptr, nrec, &mut sum, status) > 0 {
            return *status;
        }

        /* CHECKSUM is already correct */
        if sum == 0 || sum == 0xFFFFFFFF {
            return *status;
        }

        /* Zero the CHECKSUM and recompute the new value */
        ffmkys_safe(
            fptr,
            cs!(c"CHECKSUM"),
            cs!(c"0000000000000000"),
            Some(&chkcomm),
            status,
        );
    }

    /* move to the start of the header */
    ffmbyt_safe(fptr, headstart, REPORT_EOF, status);

    /* accumulate the header checksum into the previous data checksum */
    nrec = ((datastart - headstart) / 2880) as c_long;
    sum = dsum;
    if ffcsum_safe(fptr, nrec, &mut sum, status) > 0 {
        return *status;
    }

    /* encode the COMPLEMENT of the checksum into a 16-character string */
    ffesum_safe(sum, true, (&mut checksum[..17]).try_into().unwrap());

    /* update the CHECKSUM keyword value with the new string */
    ffmkys_safe(fptr, cs!(c"CHECKSUM"), &checksum, Some(cs!(c"&")), status);

    *status
}

/*------------------------------------------------------------------------*/
/// Verify the HDU by comparing the value of the computed checksums against
/// the values of the DATASUM and CHECKSUM keywords if they are present.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffvcks(
    fptr: *mut fitsfile,    /* I - FITS file pointer                  */
    datastatus: *mut c_int, /* O - data checksum status               */
    hdustatus: *mut c_int,  /* O - hdu checksum status                */
    /*     1  verification is correct         */
    /*     0  checksum keyword is not present */
    /*    -1 verification not correct         */
    status: *mut c_int, /* IO - error status                      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let datastatus = datastatus.as_mut().expect(NULL_MSG);
        let hdustatus = hdustatus.as_mut().expect(NULL_MSG);

        ffvcks_safe(fptr, datastatus, hdustatus, status)
    }
}

/*------------------------------------------------------------------------*/
/// Verify the HDU by comparing the value of the computed checksums against
/// the values of the DATASUM and CHECKSUM keywords if they are present.
pub fn ffvcks_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                  */
    datastatus: &mut c_int, /* O - data checksum status               */
    hdustatus: &mut c_int,  /* O - hdu checksum status                */
    /*     1  verification is correct         */
    /*     0  checksum keyword is not present */
    /*    -1 verification not correct         */
    status: &mut c_int, /* IO - error status                      */
) -> c_int {
    let mut datasum: c_ulong = 0;
    let mut hdusum: c_ulong = 0;
    let mut chksum: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT]; /* inherit input status value if > 0 */

    if *status > 0 {
        return *status;
    }

    *datastatus = -1;
    *hdustatus = -1;

    let tstatus = *status;
    if ffgkys_safe(fptr, cs!(c"CHECKSUM"), &mut chksum, Some(&mut comm), status) == KEY_NO_EXIST {
        /* CHECKSUM keyword does not exist */
        *hdustatus = 0;
        *status = tstatus;
    }

    if chksum[0] == 0 {
        /* all blank checksum means it is undefined */
        *hdustatus = 0;
    }

    if ffgkys_safe(fptr, cs!(c"DATASUM"), &mut chksum, Some(&mut comm), status) == KEY_NO_EXIST {
        /* DATASUM keyword does not exist */
        *datastatus = 0;
        *status = tstatus;
    }

    if chksum[0] == 0 {
        /* all blank checksum means it is undefined */
        *datastatus = 0;
    }

    if *status > 0 || ((*hdustatus) == 0 && (*datastatus) == 0) {
        /* return if neither keywords exist */
        return *status;
    }

    /* convert string to unsigned long */

    /* olddatasum = strtoul(chksum, 0, 10);  doesn't work w/ gcc on SUN OS */
    /* sscanf(chksum, "%u", &olddatasum);   doesn't work w/ cc on VAX/VMS */

    let tdouble = atof_safe(&chksum); /* read as a double as a workaround */
    let olddatasum: c_ulong = tdouble as c_ulong;

    /*  calculate the data checksum and the HDU checksum */
    if ffgcks_safe(fptr, &mut datasum, &mut hdusum, status) > 0 {
        return *status;
    }

    if *datastatus != 0 && datasum == olddatasum {
        *datastatus = 1;
    }

    if *hdustatus != 0 && (hdusum == 0 || hdusum == 0xFFFFFFFF) {
        *hdustatus = 1;
    }
    *status
}

/*------------------------------------------------------------------------*/
/// calculate the checksums of the data unit and the total HDU
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcks(
    fptr: *mut fitsfile,   /* I - FITS file pointer             */
    datasum: *mut c_ulong, /* O - data checksum                 */
    hdusum: *mut c_ulong,  /* O - hdu checksum                  */
    status: *mut c_int,    /* IO - error status                 */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let datasum = datasum.as_mut().expect(NULL_MSG);
        let hdusum = hdusum.as_mut().expect(NULL_MSG);

        ffgcks_safe(fptr, datasum, hdusum, status)
    }
}

/*------------------------------------------------------------------------*/
/// calculate the checksums of the data unit and the total HDU
pub fn ffgcks_safe(
    fptr: &mut fitsfile,   /* I - FITS file pointer             */
    datasum: &mut c_ulong, /* O - data checksum                 */
    hdusum: &mut c_ulong,  /* O - hdu checksum                  */
    status: &mut c_int,    /* IO - error status                 */
) -> c_int {
    let mut headstart: LONGLONG = 0;
    let mut datastart: LONGLONG = 0;
    let mut dataend: LONGLONG = 0;

    /* inherit input status value if > 0 */
    if *status > 0 {
        return *status;
    }

    /* get size of the HDU */
    if ffghadll_safe(
        fptr,
        Some(&mut headstart),
        Some(&mut datastart),
        Some(&mut dataend),
        status,
    ) > 0
    {
        return *status;
    }

    let mut nrec = ((dataend - datastart) / 2880i64) as c_long;

    *datasum = 0;

    if nrec > 0 {
        /* accumulate the 32-bit 1's complement checksum */
        ffmbyt_safe(fptr, datastart, REPORT_EOF, status);
        if ffcsum_safe(fptr, nrec, datasum, status) > 0 {
            return *status;
        };
    }

    /* move to the start of the header and calc. size of header */
    ffmbyt_safe(fptr, headstart, REPORT_EOF, status);
    nrec = ((datastart - headstart) / 2880) as c_long;

    /* accumulate the header checksum into the previous data checksum */
    *hdusum = *datasum;
    ffcsum_safe(fptr, nrec, hdusum, status);

    *status
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ffesum_safe_complm_false() {
        let sum: c_ulong = 0x12345678;
        let complm = false;
        let mut ascii: [c_char; 17] = [0; 17];

        ffesum_safe(sum, complm, &mut ascii);

        let result = std::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result, "N6AGN49EN4AEN49E"); // Replace with the expected encoded value
    }

    #[test]
    fn test_ffesum_safe_complm_true() {
        let sum: c_ulong = 0x12345678;
        let complm = true;
        let mut ascii: [c_char; 17] = [0; 17];

        ffesum_safe(sum, complm, &mut ascii);

        let result = std::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result, "QleaTkbTQkbZQkbZ"); // Replace with the expected encoded value
    }

    #[test]
    fn test_ffesum_safe_zero_sum() {
        let sum: c_ulong = 0;
        let complm = false;
        let mut ascii: [c_char; 17] = [0; 17];

        ffesum_safe(sum, complm, &mut ascii);

        let result = std::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result, "0000000000000000");
    }

    #[test]
    fn test_ffesum_safe_max_sum() {
        let sum: c_ulong = 0xFFFFFFFF;
        let complm = true;
        let mut ascii: [c_char; 17] = [0; 17];

        ffesum_safe(sum, complm, &mut ascii);

        let result = std::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result, "0000000000000000");
    }

    #[test]
    fn test_ffdsum_safe_complm_false() {
        let ascii: [u8; 17] = *b"N6AGN49EN4AEN49E\0";
        let ascii_cchar = ascii.map(|x| x as c_char);

        let complm = false;
        let mut sum: c_ulong = 0;

        ffdsum_safe(&ascii_cchar, complm, &mut sum);

        assert_eq!(sum, 0x12345678); // Replace with the expected decoded value
    }

    #[test]
    fn test_ffdsum_safe_complm_true() {
        let ascii: &[u8; 17] = b"QleaTkbTQkbZQkbZ\0";
        let ascii_cchar = ascii.map(|x| x as c_char);

        let complm = true;
        let mut sum: c_ulong = 0;

        ffdsum_safe(&ascii_cchar, complm, &mut sum);

        assert_eq!(sum, 0x12345678); // Replace with the expected decoded value
    }
}
