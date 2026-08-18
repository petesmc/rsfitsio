/*  This file, checksum.rs, contains the checksum-related routines in the   */
/*  FITSIO library.                                                        */
/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
/*------------------------------------------------------------------------*/

use core::slice;

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
    // FFI WRAPPER
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
    for _jj in 0..nrec {
        ffgbyt(fptr, 2880, cast_slice_mut(&mut sbuf), status);
        if BYTESWAPPED {
            /* reverse order of bytes in each value */
            ffswap2(cast_slice_mut(&mut sbuf), 1440);
        }

        hi = *sum >> 16;
        lo = *sum & 0xFFFF;
        ii = 0;
        while ii < 1440 {
            hi += c_ulong::from(sbuf[ii]);
            lo += c_ulong::from(sbuf[ii + 1]);
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    /* sscanf_u(&chksum, cs!(c"%u"), &mut olddatasum);   doesn't work w/ cc on VAX/VMS */

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
    // FFI WRAPPER
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

        let result = core::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
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

        let result = core::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
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

        let result = core::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
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

        let result = core::ffi::CStr::from_bytes_until_nul(cast_slice(&ascii))
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

    // --- Ported from test_checksum.c ---

    use crate::aliases::rust_api::*;
    use crate::fitsio::{BINARY_TBL, BYTE_IMG, KEY_NO_EXIST, READONLY, READWRITE, TBYTE, fitsfile};
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use libc::{c_char, c_int, c_long, c_ulong};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    #[test]
    fn test_write_checksum() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            assert_eq!(status, 0);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            assert_eq!(status, 0);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_verify_checksum() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            // setup: create + write checksum
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 1));
            assert!(!(hdustatus != 1));
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_update_checksum() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_update_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_get_checksums() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            let mut datasum: c_ulong = 0;
            let mut hdusum: c_ulong = 0;
            fits_get_chksum(
                f.as_deref_mut().unwrap(),
                &mut datasum,
                &mut hdusum,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_encode_decode_checksum() {
        let sum: c_ulong = 0x12345678;
        let mut decoded: c_ulong = 0;
        let mut ascii: [c_char; 17] = [0; 17];
        ffesum_safe(sum, false, &mut ascii);
        ffdsum_safe(&ascii, false, &mut decoded);
        assert!(!(decoded != sum));
    }

    #[test]
    fn test_encode_complement() {
        let sum: c_ulong = 0xABCDEF00;
        let mut ascii: [c_char; 17] = [0; 17];
        ffesum_safe(sum, true, &mut ascii);
        // strlen(ascii) == 16
        let len = ascii.iter().position(|&c| c == 0).unwrap();
        assert!(!(len != 16));
    }

    #[test]
    fn test_verify_no_keywords() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 0));
            assert!(!(hdustatus != 0));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_checksum_binary_table() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data: [c_long; 3] = [100, 200, 300];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                3,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, 0);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_update_existing_checksum() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
            let newdata: [u8; 10] = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READWRITE, &mut status);
            assert_eq!(status, 0);
            fits_write_img(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                10,
                &newdata,
                &mut status,
            );
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 1));
            assert!(!(hdustatus != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_data_checksum() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10000];
            let data: Vec<u8> = (0..10000).map(|i| (i % 256) as u8).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                10000,
                &data,
                &mut status,
            );
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 1));
            assert!(!(hdustatus != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_decode_with_carry() {
        // Test ffdsum with values that cause carry in the decode loop
        let mut sum: c_ulong = 0;
        let mut ascii: [c_char; 17] = [0; 17];
        // Use a value that will cause carries during decode
        ffesum_safe(0xFFFFFFFF, false, &mut ascii);
        ffdsum_safe(&ascii, false, &mut sum);
        assert!(!(sum != 0xFFFFFFFF));
        // Test with complement
        ffesum_safe(0xFFFFFFFF, true, &mut ascii);
        ffdsum_safe(&ascii, true, &mut sum);
        assert!(!(sum != 0xFFFFFFFF));
    }

    #[test]
    fn test_upck_no_checksum_keyword() {
        // Test ffupck when CHECKSUM keyword doesn't exist
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            // Write DATASUM but not CHECKSUM
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("DATASUM"),
                &cc("         0"),
                None,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_update_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_upck_no_datasum_keyword() {
        // Test ffupck when DATASUM keyword doesn't exist
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READWRITE, &mut status);
            status = 0;
            fits_update_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert!(!(status != KEY_NO_EXIST));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_verify_incorrect_checksum() {
        // Test ffvcks with incorrect checksum values
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let mut data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Corrupt the file by changing data after checksum was written
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            data[0] = 99;
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Verify should fail now
            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != -1));
            assert!(!(hdustatus != -1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_pcks_update_existing_datasum() {
        // Test ffpcks when DATASUM exists but data changes
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
            let newdata: [u8; 10] = [99, 98, 97, 96, 95, 94, 93, 92, 91, 90];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Open and change data, then recalculate checksum
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_write_img(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                10,
                &newdata,
                &mut status,
            );
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Verify new checksum is correct
            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 1));
            assert!(!(hdustatus != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_upck_existing_correct_checksum() {
        // Test ffupck when CHECKSUM already exists and is correct
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Call ffupck on file with correct checksum
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_update_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Verify checksum is still correct
            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 1));
            assert!(!(hdustatus != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_pcks_bad_status() {
        // Test that ffpcks returns immediately with bad status
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0);
            status = 1;
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert!(!(status != 1));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_vcks_bad_status() {
        // Test that ffvcks returns immediately with bad status
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0);
            status = 1;
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert!(!(status != 1));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_gcks_bad_status() {
        // Test that ffgcks returns immediately with bad status
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0);
            status = 1;
            let mut datasum: c_ulong = 0;
            let mut hdusum: c_ulong = 0;
            fits_get_chksum(
                f.as_deref_mut().unwrap(),
                &mut datasum,
                &mut hdusum,
                &mut status,
            );
            assert!(!(status != 1));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_upck_bad_status() {
        // Test that ffupck returns immediately with bad status
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0);
            status = 1;
            fits_update_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert!(!(status != 1));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_csum_bad_status() {
        // Test that ffcsum returns immediately with bad status
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0);
            status = 1;
            let mut sum: c_ulong = 0;
            ffcsum_safe(f.as_deref_mut().unwrap(), 1, &mut sum, &mut status);
            assert!(!(status != 1));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_empty_hdu_checksum() {
        // Test checksum on HDU with no data
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 1));
            assert!(!(hdustatus != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_varlen_column_checksum() {
        // Test checksum with variable-length column (heapsize > 0)
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let data1: [c_long; 3] = [1, 2, 3];
            let data2: [c_long; 5] = [10, 20, 30, 40, 50];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            let ttype = [Some(cc("VARDATA"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1PJ(100)")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                2,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data1, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 2, 1, 5, &data2, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcks_nonzero_checksum() {
        // Test ffpcks when CHECKSUM is not "0000..." and DATASUM exists
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img(f.as_deref_mut().unwrap(), TBYTE, 1, 10, &data, &mut status);
            // Write checksum first
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Reopen and call ffpcks again - checksum should be rechecked
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_write_chksum(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut datastatus: c_int = 0;
            let mut hdustatus: c_int = 0;
            fits_verify_chksum(
                f.as_deref_mut().unwrap(),
                &mut datastatus,
                &mut hdustatus,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(datastatus != 1));
            assert!(!(hdustatus != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_decode_large_values() {
        // Test ffdsum with values that cause carry in decode loop
        let mut sum: c_ulong = 0;
        let mut ascii: [c_char; 17] = [0; 17];
        let test_values: [c_ulong; 8] = [
            0xFFFFFFFF, 0xAAAAAAAA, 0x55555555, 0x80808080, 0xF0F0F0F0, 0x0F0F0F0F, 0x12345678,
            0x87654321,
        ];
        for &v in test_values.iter() {
            ffesum_safe(v, false, &mut ascii);
            ffdsum_safe(&ascii, false, &mut sum);
            assert!(!(sum != v));
        }
    }
}
