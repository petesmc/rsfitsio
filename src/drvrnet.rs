use crate::c_types::{c_int, c_uint};

static mut show_fits_download_progress: c_int = 0;
static mut net_timeout: c_uint = 360; /* in seconds */

pub(crate) unsafe fn fits_net_timeout(sec: c_int) -> c_int {
    unsafe {
        /* If sec is 0 or negative, treat this as a 'get' call. */
        if sec > 0 {
            net_timeout = sec as c_uint;
        }

        return net_timeout as c_int;
    }
}

pub(crate) unsafe fn fits_dwnld_prog_bar(flag: c_int) {
    unsafe {
        if (flag == 0) {
            show_fits_download_progress = 0;
        } else {
            show_fits_download_progress = 1;
        }
    }
}
