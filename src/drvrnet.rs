use crate::c_types::{c_int, c_uint};

static mut SHOW_FITS_DOWNLOAD_PROGRESS: c_int = 0;
static mut NET_TIMEOUT: c_uint = 360; /* in seconds */

pub(crate) unsafe fn fits_net_timeout(sec: c_int) -> c_int {
    unsafe {
        /* If sec is 0 or negative, treat this as a 'get' call. */
        if sec > 0 {
            NET_TIMEOUT = sec as c_uint;
        }

        NET_TIMEOUT as c_int
    }
}

pub(crate) unsafe fn fits_dwnld_prog_bar(flag: c_int) {
    unsafe {
        if flag == 0 {
            SHOW_FITS_DOWNLOAD_PROGRESS = 0;
        } else {
            SHOW_FITS_DOWNLOAD_PROGRESS = 1;
        }
    }
}
