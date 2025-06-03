/*  This file, getcol.c, contains routines that read data elements from    */
/*  a FITS image or table.  There are generic datatype routines.           */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::{cmp, mem};

use crate::c_types::{c_char, c_int, c_long, c_short, c_uint, c_ulong, c_ushort, c_void};

use bytemuck::{cast_slice, cast_slice_mut};

use crate::NullValue;
use crate::bytes_per_datatype;
use crate::fitscore::{
    ffgidm_safe, ffgisz_safe, ffgiszll_safe, ffgnrwll_safe, ffgtclll_safe,
    fits_is_compressed_image_safe,
};
use crate::fitscore::{ffpmsg_slice, ffpmsg_str};
use crate::getcolb::{ffgclb, ffgpfb_safe, ffgpvb_safe, ffgsvb_safe};
use crate::getcold::{ffgcfm_safe, ffgcld, ffgpfd_safe, ffgpvd_safe, ffgsvd_safe};
use crate::getcole::{ffgcfc_safe, ffgcle, ffgpfe_safe, ffgpve_safe, ffgsve_safe};
use crate::getcoli::{ffgcli, ffgpfi_safe, ffgpvi_safe, ffgsvi_safe};
use crate::getcolj::{
    ffgclj, ffgcljj, ffgpfj_safe, ffgpfjj_safe, ffgpvj_safe, ffgpvjj_safe, ffgsvj_safe,
    ffgsvjj_safe,
};
use crate::getcolk::{ffgclk, ffgpfk_safe, ffgpvk_safe, ffgsvk_safe};
use crate::getcoll::{ffgcll, ffgcx_safe};
use crate::getcols::ffgcls;
use crate::getcolsb::{ffgclsb, ffgpfsb_safe, ffgpvsb_safe, ffgsvsb_safe};
use crate::getcolui::{ffgclui, ffgpfui_safe, ffgpvui_safe, ffgsvui_safe};
use crate::getcoluj::{
    ffgcluj, ffgclujj, ffgpfuj_safe, ffgpfujj_safe, ffgpvuj_safe, ffgpvujj_safe, ffgsvuj_safe,
    ffgsvujj_safe,
};
use crate::getcoluk::{ffgcluk, ffgpfuk_safe, ffgpvuk_safe, ffgsvuk_safe};
use crate::int_snprintf;
use crate::{NullCheckType, fitsio::*};
use crate::{buffers::*, calculate_subsection_length};

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpxv(
    fptr: *mut fitsfile,     /* I - FITS file pointer                       */
    datatype: c_int,         /* I - datatype of the value                   */
    firstpix: *const c_long, /* I - coord of first pixel to read (1s based) */
    nelem: LONGLONG,         /* I - number of values to read                */
    nulval: *const c_void,   /* I - value for undefined pixels              */
    array: *mut c_void,      /* O - array of values that are returned       */
    anynul: *mut c_int,      /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,      /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);
        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        /* get the size of the image */
        let mut naxis = 0;
        ffgidm_safe(fptr, &mut naxis, status);

        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        ffgpxv_safe(
            fptr, datatype, firstpix, nelem, nulval, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
pub fn ffgpxv_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    firstpix: &[c_long],        /* I - coord of first pixel to read (1s based) */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: Option<NullValue>,  /* I - value for undefined pixels              */
    array: &mut [u8],           /* O - array of values that are returned       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut tfirstpix: [LONGLONG; 99] = [0; 99];
    let mut naxis = 0;

    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* get the size of the image */
    ffgidm_safe(fptr, &mut naxis, status);

    for i in 0..(naxis as usize) {
        tfirstpix[i] = firstpix[i] as LONGLONG;
    }

    ffgpxvll_safe(
        fptr, datatype, &tfirstpix, nelem, nulval, array, anynul, status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpxvll(
    fptr: *mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstpix: *const LONGLONG, /* I - coord of first pixel to read (1s based) */
    nelem: LONGLONG,           /* I - number of values to read                */
    nulval: *const c_void,     /* I - value for undefined pixels              */
    array: *mut c_void,        /* O - array of values that are returned       */
    anynul: *mut c_int,        /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,        /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);
        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        if *status > 0 || nelem == 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* get the size of the image */
        let mut naxis = 0;
        ffgidm_safe(fptr, &mut naxis, status);

        if naxis == 0 {
            *status = BAD_DIMEN;
            return *status;
        }

        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        ffgpxvll_safe(
            fptr, datatype, firstpix, nelem, nulval, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffgpxvll_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    firstpix: &[LONGLONG],      /* I - coord of first pixel to read (1s based) */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: Option<NullValue>,  /* I - value for undefined pixels              */
    array: &mut [u8],           /* O - array of values that are returned       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut naxis = 0;

    let nullcheck = NullCheckType::SetPixel;
    let mut naxes: [LONGLONG; 9] = [0; 9];
    let mut trc: [LONGLONG; 9] = [1; 9];
    let inc: [c_long; 9] = [1; 9];
    let mut dimsize: LONGLONG = 1;
    let mut firstelem: LONGLONG = 0;

    let mut dummy_nularray = vec![0; nelem as usize];

    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* get the size of the image */
    ffgidm_safe(fptr, &mut naxis, status);

    ffgiszll_safe(fptr, 9, &mut naxes, status);

    if naxis == 0 || naxes[0] == 0 {
        *status = BAD_DIMEN;
        return *status;
    }

    /* calculate the position of the first element in the array */
    firstelem = 0;
    for ii in 0..(naxis as usize) {
        firstelem += (firstpix[ii] - 1) * dimsize;
        dimsize *= naxes[ii];
        trc[ii] = firstpix[ii];
    }
    firstelem += 1;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        /* test for special case of reading an integral number of */
        /* rows in a 2D or 3D image (which includes reading the whole image */

        if naxis > 1 && naxis < 4 && firstpix[0] == 1 && (nelem / naxes[0]) * naxes[0] == nelem {
            /* calculate coordinate of last pixel */
            trc[0] = naxes[0]; /* reading whole rows */
            trc[1] = firstpix[1] + (nelem / naxes[0] - 1);
            while trc[1] > naxes[1] {
                trc[1] -= naxes[1];
                trc[2] += 1; /* increment to next plane of cube */
            }
            todo!();
            //fits_read_compressed_img(fptr, datatype, firstpix, trc, inc, 1, nulval, array, None, anynul, status);
        } else {
            todo!();
            // fits_read_compressed_pixels(fptr, datatype, firstelem,           nelem, nullcheck, nulval, array, NULL, anynul, status);
        }

        return *status;
    }

    if datatype == TBYTE {
        match nulval {
            None => {
                ffgclb(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclb(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::UByte(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSBYTE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclsb(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclsb(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Byte(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclui(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclui(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::UShort(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcli(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcli(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Short(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcluk(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcluk(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::UInt(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclk(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclk(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Int(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcluj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcluj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::ULong(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Long(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONGLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclujj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclujj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::ULONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONGLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcljj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcljj(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::LONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TFLOAT {
        let array = cast_slice_mut(array);
        match nulval {
            None => {
                ffgcle(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0.0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcle(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Float(val) => val,
                        _ => 0.0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TDOUBLE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcld(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0.0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcld(
                    fptr,
                    2,
                    1,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Double(val) => val,
                        _ => 0.0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// The nullarray values will = 1 if the corresponding array value is null.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpxf(
    fptr: *mut fitsfile,     /* I - FITS file pointer                       */
    datatype: c_int,         /* I - datatype of the value                   */
    firstpix: *const c_long, /* I - coord of first pixel to read (1s based) */
    nelem: LONGLONG,         /* I - number of values to read            */
    array: *mut c_void,      /* O - array of values that are returned       */
    nullarray: *mut c_char,  /* O - returned array of null value flags      */
    anynul: *mut c_int,      /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,      /* IO - error status                           */
) -> c_int {
    unsafe {
        let mut tfirstpix: [LONGLONG; 99] = [0; 99];
        let mut naxis = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nullarray = slice::from_raw_parts_mut(nullarray, nelem as usize);

        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);

        if *status > 0 || nelem == 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* get the size of the image */
        ffgidm_safe(fptr, &mut naxis, status);

        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        for i in 0..(naxis as usize) {
            tfirstpix[i] = firstpix[i] as LONGLONG;
        }

        ffgpxfll_safe(
            fptr,
            datatype,
            &tfirstpix,
            nelem,
            array,
            nullarray,
            anynul.as_mut(),
            status,
        );

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// The nullarray values will = 1 if the corresponding array value is null.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpxfll(
    fptr: *mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstpix: *const LONGLONG, /* I - coord of first pixel to read (1s based) */
    nelem: LONGLONG,           /* I - number of values to read              */
    array: *mut c_void,        /* O - array of values that are returned       */
    nullarray: *mut c_char,    /* O - returned array of null value flags      */
    anynul: *mut c_int,        /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,        /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let nullarray = slice::from_raw_parts_mut(nullarray, nelem as usize);
        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);

        /* get the size of the image */
        let mut naxis = 0;
        ffgidm_safe(fptr, &mut naxis, status);

        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        ffgpxfll_safe(
            fptr,
            datatype,
            firstpix,
            nelem,
            array,
            nullarray,
            anynul.as_mut(),
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// The nullarray values will = 1 if the corresponding array value is null.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffgpxfll_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    firstpix: &[LONGLONG],      /* I - coord of first pixel to read (1s based) */
    nelem: LONGLONG,            /* I - number of values to read              */
    array: &mut [u8],           /* O - array of values that are returned       */
    nullarray: &mut [c_char],   /* O - returned array of null value flags      */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut naxis = 0;
    let nullcheck = NullCheckType::SetNullArray;
    let mut naxes: [LONGLONG; 9] = [0; 9];
    let mut dimsize: LONGLONG = 1;
    let mut firstelem: LONGLONG = 0;

    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* get the size of the image */
    ffgidm_safe(fptr, &mut naxis, status);
    ffgiszll_safe(fptr, 9, &mut naxes, status);

    /* calculate the position of the first element in the array */
    firstelem = 0;
    for ii in 0..(naxis as usize) {
        firstelem += (firstpix[ii] - 1) * dimsize;
        dimsize *= naxes[ii];
    }
    firstelem += 1;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        todo!();
        //fits_read_compressed_pixels(fptr, datatype, firstelem, nelem, nullcheck, None, array, nullarray, anynul, status);
        return *status;
    }

    if datatype == TBYTE {
        let array = cast_slice_mut(array);

        ffgclb(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TSBYTE {
        let array = cast_slice_mut(array);

        ffgclsb(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TUSHORT {
        let array = cast_slice_mut(array);

        ffgclui(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TSHORT {
        let array = cast_slice_mut(array);

        ffgcli(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TUINT {
        let array = cast_slice_mut(array);

        ffgcluk(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TINT {
        let array = cast_slice_mut(array);

        ffgclk(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TULONG {
        let array = cast_slice_mut(array);

        ffgcluj(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TLONG {
        let array = cast_slice_mut(array);

        ffgclj(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TULONGLONG {
        let array = cast_slice_mut(array);

        ffgclujj(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TLONGLONG {
        let array = cast_slice_mut(array);

        ffgcljj(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TFLOAT {
        let array = cast_slice_mut(array);

        ffgcle(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0.0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TDOUBLE {
        let array = cast_slice_mut(array);

        ffgcld(
            fptr,
            2,
            1,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            0.0,
            array,
            nullarray,
            anynul,
            status,
        );
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read an section of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsv(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    datatype: c_int,       /* I - datatype of the value                   */
    blc: *const c_long,    /* I - 'bottom left corner' of the subsection  */
    trc: *const c_long,    /* I - 'top right corner' of the subsection    */
    inc: *const c_long,    /* I - increment to be applied in each dim.    */
    nulval: *const c_void, /* I - value for undefined pixels              */
    array: *mut c_void,    /* O - array of values that are returned       */
    anynul: *mut c_int,    /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        let mut naxis = 0;
        /* get the size of the image */
        ffgidm_safe(fptr, &mut naxis, status);

        let blc = slice::from_raw_parts(blc, naxis as usize);
        let trc = slice::from_raw_parts(trc, naxis as usize);
        let inc = slice::from_raw_parts(inc, naxis as usize);

        let nelem = calculate_subsection_length(blc, trc, inc);
        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);

        ffgsv_safe(fptr, datatype, blc, trc, inc, nulval, array, anynul, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read an section of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
pub fn ffgsv_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    blc: &[c_long],             /* I - 'bottom left corner' of the subsection  */
    trc: &[c_long],             /* I - 'top right corner' of the subsection    */
    inc: &[c_long],             /* I - increment to be applied in each dim.    */
    nulval: Option<NullValue>,  /* I - value for undefined pixels              */
    array: &mut [u8],           /* O - array of values that are returned       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut naxis = 0;
    let mut naxes: [c_long; 9] = [0; 9];
    let mut nelem: LONGLONG = 1;
    let mut ii: usize = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* get the size of the image */
    ffgidm_safe(fptr, &mut naxis, status);
    ffgisz_safe(fptr, 9, &mut naxes, status);

    /* test for the important special case where we are reading the whole image */
    /* this is only useful for images that are not tile-compressed */
    if fits_is_compressed_image_safe(fptr, status) == 0 {
        while ii < naxis as usize {
            if inc[ii] != 1 || blc[ii] != 1 || trc[ii] != naxes[ii] {
                break;
            }

            nelem *= naxes[ii] as LONGLONG;
            ii += 1;
        }

        if ii == naxis as usize {
            /* read the whole image more efficiently */
            ffgpxv_safe(fptr, datatype, blc, nelem, nulval, array, anynul, status);
            return *status;
        }
    }

    if datatype == TBYTE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvb_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvb_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::UByte(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSBYTE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvsb_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvsb_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::Byte(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvui_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvui_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::UShort(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvi_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvi_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::Short(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvuk_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvuk_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::UInt(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvk_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvk_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::Int(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvuj_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvuj_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::ULong(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvj_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvj_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::Long(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONGLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvujj_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvujj_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::ULONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONGLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvjj_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvjj_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::LONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TFLOAT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsve_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0.0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsve_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::Float(val) => val,
                        _ => 0.0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TDOUBLE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgsvd_safe(
                    fptr, 1, naxis, &naxes, blc, trc, inc, 0.0, array, anynul, status,
                );
            }
            Some(val) => {
                ffgsvd_safe(
                    fptr,
                    1,
                    naxis,
                    &naxes,
                    blc,
                    trc,
                    inc,
                    match val {
                        NullValue::Double(val) => val,
                        _ => 0.0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
///
/// anynul can be a null pointer
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpv(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    datatype: c_int,       /* I - datatype of the value                   */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    nulval: *const c_void, /* I - value for undefined pixels              */
    array: *mut c_void,    /* O - array of values that are returned       */
    anynul: *mut c_int,    /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);

        let anynul = anynul.as_mut();
        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        ffgpv_safe(
            fptr, datatype, firstelem, nelem, nulval, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
///
/// anynul can be a null pointer
pub fn ffgpv_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: Option<NullValue>,  /* I - value for undefined pixels              */
    array: &mut [u8],           /* O - array of values that are returned       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if datatype == TBYTE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvb_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvb_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::UByte(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSBYTE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvsb_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvsb_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::Byte(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvui_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvui_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::UShort(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvi_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvi_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::Short(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvuk_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvuk_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::UInt(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvk_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvk_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::Int(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvuj_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvuj_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::ULong(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvj_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvj_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::Long(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONGLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvujj_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvujj_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::ULONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONGLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvjj_safe(fptr, 1, firstelem, nelem, 0, array, anynul, status);
            }
            Some(val) => {
                ffgpvjj_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::LONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TFLOAT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpve_safe(fptr, 1, firstelem, nelem, 0.0, array, anynul, status);
            }
            Some(val) => {
                ffgpve_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::Float(val) => val,
                        _ => 0.0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TDOUBLE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgpvd_safe(fptr, 1, firstelem, nelem, 0.0, array, anynul, status);
            }
            Some(val) => {
                ffgpvd_safe(
                    fptr,
                    1,
                    firstelem,
                    nelem,
                    match val {
                        NullValue::Double(val) => val,
                        _ => 0.0,
                    },
                    array,
                    anynul,
                    status,
                );
            }
        }
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// The nullarray values will = 1 if the corresponding array value is null.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpf(
    fptr: *mut fitsfile,    /* I - FITS file pointer                       */
    datatype: c_int,        /* I - datatype of the value                   */
    firstelem: LONGLONG,    /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,        /* I - number of values to read                */
    array: *mut c_void,     /* O - array of values that are returned       */
    nullarray: *mut c_char, /* O - array of null value flags               */
    anynul: *mut c_int,     /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,     /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let nullarray = slice::from_raw_parts_mut(nullarray, nelem as usize);
        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);

        let anynul = anynul.as_mut();

        ffgpf_safe(
            fptr, datatype, firstelem, nelem, array, nullarray, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// The nullarray values will = 1 if the corresponding array value is null.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffgpf_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [u8],           /* O - array of values that are returned       */
    nullarray: &mut [c_char],   /* O - array of null value flags               */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if datatype == TBYTE {
        let array = cast_slice_mut(array);

        ffgpfb_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TSBYTE {
        let array = cast_slice_mut(array);

        ffgpfsb_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TUSHORT {
        let array = cast_slice_mut(array);

        ffgpfui_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TSHORT {
        let array = cast_slice_mut(array);

        ffgpfi_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TUINT {
        let array = cast_slice_mut(array);

        ffgpfuk_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TINT {
        let array = cast_slice_mut(array);

        ffgpfk_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TULONG {
        let array = cast_slice_mut(array);

        ffgpfuj_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TLONG {
        let array = cast_slice_mut(array);

        ffgpfj_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TULONGLONG {
        let array = cast_slice_mut(array);

        ffgpfujj_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TLONGLONG {
        let array = cast_slice_mut(array);

        ffgpfjj_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TFLOAT {
        let array = cast_slice_mut(array);

        ffgpfe_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else if datatype == TDOUBLE {
        let array = cast_slice_mut(array);

        ffgpfd_safe(fptr, 1, firstelem, nelem, array, nullarray, anynul, status);
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a table column. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of true if any pixels are undefined.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcv(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    datatype: c_int,       /* I - datatype of the value                   */
    colnum: c_int,         /* I - number of column to read (1 = 1st col) */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)        */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    nulval: *const c_void, /* I - value for undefined pixels              */
    array: *mut c_void,    /* O - array of values that are returned       */
    anynul: *mut c_int,    /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let nulval = NullValue::from_raw_ptr(datatype, nulval);
        let bytes = bytes_per_datatype(datatype).unwrap();
        let array = slice::from_raw_parts_mut(array as *mut _, bytes * nelem as usize);

        ffgcv_safe(
            fptr,
            datatype,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            nulval,
            array,
            anynul,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a table column. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of true if any pixels are undefined.
pub fn ffgcv_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    colnum: c_int,              /* I - number of column to read (1 = 1st col) */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)        */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: Option<NullValue>,  /* I - value for undefined pixels              */
    array: &mut [u8],           /* O - array of values that are returned       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut cdummy: [c_char; 2] = [0; 2];

    let mut dummy_nularray = vec![0; nelem as usize];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if datatype == TBIT {
        let array = cast_slice_mut(array);

        ffgcx_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            status,
        );
    } else if datatype == TBYTE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclb(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclb(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::UByte(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSBYTE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclsb(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclsb(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Byte(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclui(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclui(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::UShort(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSHORT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcli(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcli(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Short(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TUINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcluk(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcluk(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::UInt(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TINT {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclk(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclk(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Int(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcluj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcluj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::ULong(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Long(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TULONGLONG {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgclujj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgclujj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::ULONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLONGLONG {
        let array = cast_slice_mut(array);
        match nulval {
            None => {
                ffgcljj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcljj(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::LONGLONG(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TFLOAT {
        let array = cast_slice_mut(array);
        match nulval {
            None => {
                ffgcle(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0.,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcle(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Float(val) => val,
                        _ => 0.,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TDOUBLE {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcld(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    0.,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcld(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Double(val) => val,
                        _ => 0.,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TCOMPLEX {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcle(
                    fptr,
                    colnum,
                    firstrow,
                    (firstelem - 1) * 2 + 1,
                    nelem * 2,
                    1,
                    NullCheckType::SetPixel,
                    0.,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcle(
                    fptr,
                    colnum,
                    firstrow,
                    (firstelem - 1) * 2 + 1,
                    nelem * 2,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Float(val) => val,
                        _ => 0.,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TDBLCOMPLEX {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcld(
                    fptr,
                    colnum,
                    firstrow,
                    (firstelem - 1) * 2 + 1,
                    nelem * 2,
                    1,
                    NullCheckType::SetPixel,
                    0.,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcld(
                    fptr,
                    colnum,
                    firstrow,
                    (firstelem - 1) * 2 + 1,
                    nelem * 2,
                    1,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Double(val) => val,
                        _ => 0.,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TLOGICAL {
        let array = cast_slice_mut(array);

        match nulval {
            None => {
                ffgcll(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    NullCheckType::SetPixel,
                    0,
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
            Some(val) => {
                ffgcll(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem,
                    nelem,
                    NullCheckType::SetPixel,
                    match val {
                        NullValue::Logical(val) => val,
                        _ => 0,
                    },
                    array,
                    &mut dummy_nularray,
                    anynul,
                    status,
                );
            }
        }
    } else if datatype == TSTRING {
        // SAFETY: TODO
        unsafe {
            let array =
                slice::from_raw_parts_mut(array.as_ptr() as *mut *mut c_char, nelem as usize);
            let mut v_array = Vec::new();
            for item in array {
                let array_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                v_array.push(array_item);
            }

            match nulval {
                None => {
                    cdummy[0] = 0;
                    let cdummy2 = cdummy;
                    ffgcls(
                        fptr,
                        colnum,
                        firstrow,
                        firstelem,
                        nelem,
                        NullCheckType::SetPixel,
                        Some(&cdummy),
                        &mut v_array,
                        &mut dummy_nularray,
                        anynul,
                        status,
                    );
                }
                Some(val) => match val {
                    NullValue::String(x) => {
                        ffgcls(
                            fptr,
                            colnum,
                            firstrow,
                            firstelem,
                            nelem,
                            NullCheckType::SetPixel,
                            Some(cast_slice(x.as_bytes_with_nul())),
                            &mut v_array,
                            &mut dummy_nularray,
                            anynul,
                            status,
                        );
                    }
                    _ => {
                        ffgcls(
                            fptr,
                            colnum,
                            firstrow,
                            firstelem,
                            nelem,
                            NullCheckType::SetPixel,
                            Some(&[0]),
                            &mut v_array,
                            &mut dummy_nularray,
                            anynul,
                            status,
                        );
                    }
                },
            }
        }
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read arrays of values from NCOLS table columns. This is an optimization
/// to read all columns in one pass through the table.  The datatypes of the
/// input arrays are defined by the 3rd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements for column i will be set equal to *(nulval[i]), unless nulval[i]=0
/// in which case no checking for undefined values will be performed.
/// anynul[i] is returned with a value of true if any pixels in column i are undefined.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcvn(
    fptr: *mut fitsfile,          /* I - FITS file pointer                       */
    ncols: c_int,                 /* I - number of columns to read               */
    datatype: *const c_int,       /* I - datatypes of the values                 */
    colnum: *const c_int,         /* I - columns numbers to read (1 = 1st col)   */
    firstrow: LONGLONG,           /* I - first row to read (1 = 1st row)     */
    nrows: LONGLONG,              /* I - number of rows to read              */
    nulval: *const *const c_void, /* I - array of pointers to values for undefined pixels */
    array: *mut *mut c_void,      /* O - array of pointers to values that are returned    */
    anynul: *mut c_int, /* O - anynul[i] set to 1 if any values in column i are null; else 0 */
    status: *mut c_int, /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts_mut(array, nrows as usize);

        let datatype = slice::from_raw_parts(datatype, nrows as usize);
        let colnum = slice::from_raw_parts(colnum, nrows as usize);

        let nulval = slice::from_raw_parts(nulval, nrows as usize);
        let nulval = nulval
            .iter()
            .zip(datatype.iter())
            .map(|(&val, &dtype)| NullValue::from_raw_ptr(dtype, val))
            .collect::<Vec<_>>();

        let anynul = match anynul.is_null() {
            false => Some(slice::from_raw_parts_mut(anynul, ncols as usize)),
            true => None,
        };

        ffgcvn_safer(
            fptr, ncols, datatype, colnum, firstrow, nrows, &nulval, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read arrays of values from NCOLS table columns. This is an optimization
/// to read all columns in one pass through the table.  The datatypes of the
/// input arrays are defined by the 3rd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements for column i will be set equal to *(nulval[i]), unless nulval[i]=0
/// in which case no checking for undefined values will be performed.
/// anynul[i] is returned with a value of true if any pixels in column i are undefined.
pub unsafe fn ffgcvn_safer(
    fptr: &mut fitsfile,          /* I - FITS file pointer                       */
    ncols: c_int,                 /* I - number of columns to read               */
    datatype: &[c_int],           /* I - datatypes of the values                 */
    colnum: &[c_int],             /* I - columns numbers to read (1 = 1st col)   */
    firstrow: LONGLONG,           /* I - first row to read (1 = 1st row)     */
    nrows: LONGLONG,              /* I - number of rows to read              */
    nulval: &[Option<NullValue>], /* I - array of pointers to values for undefined pixels */
    array: &mut [*mut c_void],    /* O - array of pointers to values that are returned    */
    mut anynul: Option<&mut [c_int]>, /* O - anynul[i] set to 1 if any values in column i are null; else 0 */
    status: &mut c_int,               /* IO - error status                           */
) -> c_int {
    unsafe {
        let mut ntotrows: LONGLONG = 0;
        let mut ndone: LONGLONG = 0;
        let mut nread: LONGLONG = 0;
        let mut currow: LONGLONG = 0;
        let mut nrowbuf: c_long = 0;

        let mut sizes: [usize; 255] = [0; 255];

        sizes[TBYTE as usize] = mem::size_of::<c_char>();
        sizes[TSBYTE as usize] = mem::size_of::<c_char>();
        sizes[TLOGICAL as usize] = mem::size_of::<c_char>();
        sizes[TUSHORT as usize] = mem::size_of::<c_short>();
        sizes[TSHORT as usize] = mem::size_of::<c_short>();
        sizes[TINT as usize] = mem::size_of::<c_int>();
        sizes[TUINT as usize] = mem::size_of::<c_int>();
        sizes[TLONG as usize] = mem::size_of::<c_long>();
        sizes[TULONG as usize] = mem::size_of::<c_long>();
        sizes[TLONGLONG as usize] = mem::size_of::<LONGLONG>();
        sizes[TULONGLONG as usize] = mem::size_of::<LONGLONG>();
        sizes[TFLOAT as usize] = mem::size_of::<f32>();
        sizes[TDOUBLE as usize] = mem::size_of::<f64>();
        sizes[TDBLCOMPLEX as usize] = 2 * mem::size_of::<f64>();

        if *status > 0 {
            return *status;
        }

        if ncols <= 0 {
            *status = 0;
            return *status;
        }

        let mut repeats: Vec<LONGLONG> = vec![0; ncols as usize];

        ffgnrwll_safe(fptr, &mut ntotrows, status);
        ffgrsz_safe(fptr, &mut nrowbuf, status);

        /* Retrieve column repeats */
        let mut icol: usize = 0;
        while (icol < ncols as usize) && (icol < 1000) {
            let mut typecode = 0;
            let mut repeat: LONGLONG = 0;
            let mut width: LONGLONG = 0;
            ffgtclll_safe(
                fptr,
                colnum[icol],
                Some(&mut typecode),
                Some(&mut repeat),
                Some(&mut width),
                status,
            );
            repeats[icol] = repeat;

            if datatype[icol] == TBIT
                || datatype[icol] == TSTRING
                || sizes[datatype[icol] as usize] == 0
            {
                ffpmsg_str("Cannot read from TBIT or TSTRING datatypes (ffgcvn)");
                *status = BAD_DATATYPE;
            }
            if typecode < 0 {
                ffpmsg_str("Cannot read from variable-length data (ffgcvn)");
                *status = BAD_DIMEN;
            }

            if *status > 0 {
                break;
            }
            icol += 1;
        }

        if *status > 0 {
            return *status;
        }

        /* Optimize for 1 column */
        if ncols == 1 {
            let bytes = match bytes_per_datatype(datatype[0]) {
                Some(x) => x * nrows as usize * repeats[0] as usize,
                None => {
                    *status = BAD_DATATYPE;
                    return *status;
                }
            };

            let arr = slice::from_raw_parts_mut(
                array[0] as *mut u8,
                bytes * (nrows * repeats[0]) as usize,
            );

            ffgcv_safe(
                fptr,
                datatype[0],
                colnum[0],
                firstrow,
                1,
                nrows * repeats[0],
                nulval[0].clone(),
                arr,
                anynul.map(|x| &mut x[0]),
                status,
            );

            return *status;
        }

        /* Scan through file, in chunks of nrowbuf */
        currow = firstrow;
        ndone = 0;
        while ndone < nrows {
            nread = nrows - ndone; /* Number of rows to read (not elements) */
            if nread > nrowbuf as LONGLONG {
                nread = nrowbuf as LONGLONG;
            }

            for icol in 0..(ncols as usize) {
                let nelem1: LONGLONG = nread * repeats[icol];

                let bytes = match bytes_per_datatype(datatype[icol]) {
                    Some(x) => x * nrows as usize * repeats[icol] as usize,
                    None => {
                        *status = BAD_DATATYPE;
                        return *status;
                    }
                };

                let arr = slice::from_raw_parts_mut(
                    array[icol] as *mut u8,
                    bytes * (nrows * repeats[icol]) as usize,
                );
                let array1 =
                    &mut arr[(repeats[icol] * ndone) as usize * (sizes[datatype[icol] as usize])..];

                ffgcv_safe(
                    fptr,
                    datatype[icol],
                    colnum[icol],
                    currow,
                    1,
                    nelem1,
                    nulval[icol].clone(),
                    array1,
                    anynul.as_deref_mut().map(|x| &mut x[icol]),
                    status,
                );

                if *status > 0 {
                    let mut errmsg: [c_char; 100] = [0; 100];
                    int_snprintf!(
                        &mut errmsg,
                        errmsg.len(),
                        "Failed to read column {} data rows {}-{} (ffgcvn)",
                        colnum[icol],
                        currow,
                        currow + nread - 1,
                    );
                    ffpmsg_slice(&errmsg);
                    break;
                }
            }

            if *status > 0 {
                break;
            }
            currow += nread;
            ndone += nread;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a table column. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// ANYNUL is returned with a value of true if any pixels are undefined.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcf(
    fptr: *mut fitsfile,    /* I - FITS file pointer                       */
    datatype: c_int,        /* I - datatype of the value                   */
    colnum: c_int,          /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,     /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,    /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,        /* I - number of values to read                */
    array: *mut c_void,     /* O - array of values that are returned       */
    nullarray: *mut c_char, /* O - array of null value flags               */
    anynul: *mut c_int,     /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,     /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let anynul = anynul.as_mut();

        let nullarray = slice::from_raw_parts_mut(nullarray, nelem as usize);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array: &mut [u8] = slice::from_raw_parts_mut(array as *mut u8, bytes);

        ffgcf_safe(
            fptr,
            datatype,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            nullarray,
            anynul,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a table column. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// ANYNUL is returned with a value of true if any pixels are undefined.
pub fn ffgcf_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    datatype: c_int,            /* I - datatype of the value                   */
    colnum: c_int,              /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,         /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [u8],           /* O - array of values that are returned       */
    nullarray: &mut [c_char],   /* O - array of null value flags               */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let nulval: f64 = 0.0;
    let cnulval: [c_char; 2] = [0; 2];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if datatype == TBIT {
        let array = cast_slice_mut(array);

        ffgcx_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            status,
        );
    } else if datatype == TBYTE {
        let array = cast_slice_mut(array);

        ffgclb(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as u8,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TSBYTE {
        let array = cast_slice_mut(array);

        ffgclsb(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as i8,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TUSHORT {
        let array = cast_slice_mut(array);

        ffgclui(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as c_ushort,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TSHORT {
        let array = cast_slice_mut(array);

        ffgcli(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as c_short,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TUINT {
        let array = cast_slice_mut(array);

        ffgcluk(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as c_uint,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TINT {
        let array = cast_slice_mut(array);

        ffgclk(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as c_int,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TULONG {
        let array = cast_slice_mut(array);

        ffgcluj(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as c_ulong,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TLONG {
        let array = cast_slice_mut(array);
        ffgclj(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as c_long,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TULONGLONG {
        let array = cast_slice_mut(array);

        ffgclujj(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as ULONGLONG,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TLONGLONG {
        let array = cast_slice_mut(array);
        ffgcljj(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as LONGLONG,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TFLOAT {
        let array = cast_slice_mut(array);
        ffgcle(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval as f32,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TDOUBLE {
        let array = cast_slice_mut(array);

        ffgcld(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            1,
            NullCheckType::SetNullArray,
            nulval,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TCOMPLEX {
        let array = cast_slice_mut(array);

        ffgcfc_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TDBLCOMPLEX {
        let array = cast_slice_mut(array);

        ffgcfm_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TLOGICAL {
        let array = cast_slice_mut(array);

        ffgcll(
            fptr,
            colnum,
            firstrow,
            firstelem,
            nelem,
            NullCheckType::SetNullArray,
            nulval as c_char,
            array,
            nullarray,
            anynul,
            status,
        );
    } else if datatype == TSTRING {
        unsafe {
            let array =
                slice::from_raw_parts_mut(array.as_mut_ptr() as *mut *mut _, nelem as usize);
            let mut v_array = Vec::new();
            for item in array {
                let array_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                v_array.push(array_item);
            }

            ffgcls(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                NullCheckType::SetNullArray,
                Some(&cnulval),
                &mut v_array,
                nullarray,
                anynul,
                status,
            );
        }
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}
