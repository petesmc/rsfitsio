use std::sync::Mutex;

use crate::c_types::{c_int, c_uint};

static SHOW_FITS_DOWNLOAD_PROGRESS: Mutex<c_int> = Mutex::new(0);
static NET_TIMEOUT: Mutex<c_uint> = Mutex::new(360); /* in seconds */

pub(crate) fn fits_net_timeout(sec: c_int) -> c_int {
    let mut net_timeout = NET_TIMEOUT.lock().unwrap();
    /* If sec is 0 or negative, treat this as a 'get' call. */
    if sec > 0 {
        *net_timeout = sec as c_uint;
    }

    *net_timeout as c_int
}

pub(crate) fn fits_dwnld_prog_bar(flag: c_int) {
    let mut download_progress = SHOW_FITS_DOWNLOAD_PROGRESS.lock().unwrap();
    if flag == 0 {
        *download_progress = 0;
    } else {
        *download_progress = 1;
    }
}
