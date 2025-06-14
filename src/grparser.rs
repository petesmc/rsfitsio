/*  This file, grparser.c, contains the group parser template routines.      */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use std::ffi::CStr;

use crate::c_types::{c_char, c_int};
use crate::fitsio::{NULL_MSG, fitsfile};
use crate::raw_to_slice;

use bytemuck::cast_slice;

/*--------------------------------------------------------------------------*/
/// Execute template to fill in header keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_execute_template(
    ff: *mut fitsfile,           /* I - FITS file pointer */
    ngp_template: *const c_char, /* I - template string */
    status: *mut c_int,          /* IO - error status */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let ff = ff.as_mut().expect(NULL_MSG);

        raw_to_slice!(ngp_template);

        fits_execute_template_safer(ff, ngp_template, status)
    }
}

/// Execute template to fill in header keywords (safe version)
pub fn fits_execute_template_safer(
    ff: &mut fitsfile,       /* I - FITS file pointer */
    ngp_template: &[c_char], /* I - template string */
    status: &mut c_int,      /* IO - error status */
) -> c_int {
    todo!("fits_execute_template")
}
