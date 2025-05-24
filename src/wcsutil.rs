use crate::c_types::*;

use crate::bb;
use crate::fitsio::*;

const D2R: f64 = 0.01745329252;

#[allow(clippy::approx_constant)]
const TWOPI: f64 = 6.28318530717959; // Could use f64::consts::TAU

/*--------------------------------------------------------------------------*/
///This routine is based on the classic AIPS WCS routine.
///
/// It converts from pixel location to RA,Dec for 9 projective geometries:
/// "-CAR", "-SIN", "-TAN", "-ARC", "-NCP", "-GLS", "-MER", "-AIT" and "-STG".
///
///
/// routine to determine accurate position for pixel coordinates          
/// returns 0 if successful otherwise:                                    
/// 501 = angle too large for projection;                                 
/// does: -CAR, -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT  -STG projections
/// Input:                                                                
///   f   xpix    x pixel number  (RA or long without rotation)           
///   f   ypiy    y pixel number  (dec or lat without rotation)           
///   d   xref    x reference coordinate value (deg)                      
///   d   yref    y reference coordinate value (deg)                      
///   f   xrefpix x reference pixel                                       
///   f   yrefpix y reference pixel                                       
///   f   xinc    x coordinate increment (deg)                            
///   f   yinc    y coordinate increment (deg)                            
///   f   rot     rotation (deg)  (from N through E)                      
///   c  *ptype    projection ptype code e.g. "-SIN";                       
/// Output:                                                               
///   d   *xpos   x (RA) coordinate (deg)                                 
///   d   *ypos   y (dec) coordinate (deg)                                
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffwldp(
    xpix: f64,
    ypix: f64,
    xref: f64,
    yref: f64,
    xrefpix: f64,
    yrefpix: f64,
    xinc: f64,
    yinc: f64,
    rot: f64,
    ptype: *const [c_char; 5],
    xpos: *mut f64,
    ypos: *mut f64,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let xpos = xpos.as_mut().expect(NULL_MSG);
        let ypos = ypos.as_mut().expect(NULL_MSG);
        let ptype = ptype.as_ref().expect(NULL_MSG);

        ffwldp_safe(
            xpix, ypix, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, ptype, xpos, ypos, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
///This routine is based on the classic AIPS WCS routine.
///
/// It converts from pixel location to RA,Dec for 9 projective geometries:
/// "-CAR", "-SIN", "-TAN", "-ARC", "-NCP", "-GLS", "-MER", "-AIT" and "-STG".
///
///
/// routine to determine accurate position for pixel coordinates          
/// returns 0 if successful otherwise:                                    
/// 501 = angle too large for projection;                                 
/// does: -CAR, -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT  -STG projections
/// Input:                                                                
///   f   xpix    x pixel number  (RA or long without rotation)           
///   f   ypiy    y pixel number  (dec or lat without rotation)           
///   d   xref    x reference coordinate value (deg)                      
///   d   yref    y reference coordinate value (deg)                      
///   f   xrefpix x reference pixel                                       
///   f   yrefpix y reference pixel                                       
///   f   xinc    x coordinate increment (deg)                            
///   f   yinc    y coordinate increment (deg)                            
///   f   rot     rotation (deg)  (from N through E)                      
///   c  *ptype    projection ptype code e.g. "-SIN";                       
/// Output:                                                               
///   d   *xpos   x (RA) coordinate (deg)                                 
///   d   *ypos   y (dec) coordinate (deg)    
pub(crate) fn ffwldp_safe(
    xpix: f64,
    ypix: f64,
    xref: f64,
    yref: f64,
    xrefpix: f64,
    yrefpix: f64,
    xinc: f64,
    yinc: f64,
    rot: f64,
    ptype: &[c_char; 5],
    xpos: &mut f64,
    ypos: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut cosr: f64 = 0.0;
    let mut sinr: f64 = 0.0;
    let mut dx: f64 = 0.0;
    let mut dy: f64 = 0.0;
    let mut dz: f64 = 0.0;
    let mut temp: f64 = 0.0;
    let mut x: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut z: f64 = 0.0;
    let mut sins: f64 = 0.0;
    let mut coss: f64 = 0.0;
    let mut dect: f64 = 0.0;
    let mut rat: f64 = 0.0;
    let mut dt: f64 = 0.0;
    let mut l: f64 = 0.0;
    let mut m: f64 = 0.0;
    let mut mg: f64 = 0.0;
    let mut da: f64 = 0.0;
    let mut dd: f64 = 0.0;
    let mut cos0: f64 = 0.0;
    let mut sin0: f64 = 0.0;
    let mut dec0: f64 = 0.0;
    let mut ra0: f64 = 0.0;
    let mut geo1: f64 = 0.0;
    let mut geo2: f64 = 0.0;
    let mut geo3: f64 = 0.0;
    let deps: f64 = 1.0e-5;

    if *status > 0 {
        return *status;
    }

    /*   Offset from ref pixel  */
    dx = (xpix - xrefpix) * xinc;
    dy = (ypix - yrefpix) * yinc;

    /*   Take out rotation  */
    cosr = f64::cos(rot * D2R);
    sinr = f64::sin(rot * D2R);
    if rot != 0.0 {
        temp = dx * cosr - dy * sinr;
        dy = dy * cosr + dx * sinr;
        dx = temp;
    }

    /* convert to radians  */
    ra0 = xref * D2R;
    dec0 = yref * D2R;

    l = dx * D2R;
    m = dy * D2R;
    sins = l * l + m * m;
    cos0 = f64::cos(dec0);
    sin0 = f64::sin(dec0);

    if ptype[0] != bb(b'-') {
        /* unrecognized projection code */
        *status = BAD_WCS_PROJ;
        return *status;
    }

    let cptr = 1;

    if ptype[cptr] == bb(b'C') {
        /* linear -CAR */
        if ptype[cptr + 1] != bb(b'A') || ptype[cptr + 2] != bb(b'R') {
            *status = BAD_WCS_PROJ;
            return *status;
        }
        rat = ra0 + l;
        dect = dec0 + m;
    } else if ptype[cptr] == bb(b'T') {
        /* -TAN */
        if ptype[cptr + 1] != bb(b'A') || ptype[cptr + 2] != bb(b'N') {
            *status = BAD_WCS_PROJ;
            return *status;
        }
        x = cos0 * f64::cos(ra0) - l * f64::sin(ra0) - m * f64::cos(ra0) * sin0;
        y = cos0 * f64::sin(ra0) + l * f64::cos(ra0) - m * f64::sin(ra0) * sin0;
        z = sin0 + m * cos0;
        rat = f64::atan2(y, x);
        dect = f64::atan(z / f64::sqrt(x * x + y * y));
    } else if ptype[cptr] == bb(b'S') {
        if ptype[cptr + 1] == bb(b'I') && ptype[cptr + 2] == bb(b'N') {
            /* -SIN */
            if sins > 1.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            coss = f64::sqrt(1.0 - sins);
            dt = sin0 * coss + cos0 * m;
            if (dt > 1.0) || (dt < -1.0) {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            dect = f64::asin(dt);
            rat = cos0 * coss - sin0 * m;
            if (rat == 0.0) && (l == 0.0) {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            rat = f64::atan2(l, rat) + ra0;
        } else if ptype[cptr + 1] == bb(b'T') && ptype[cptr + 2] == bb(b'G') {
            /* -STG Sterographic*/
            dz = (4.0 - sins) / (4.0 + sins);
            if (dz.abs()) > 1.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            dect = dz * sin0 + m * cos0 * (1.0 + dz) / 2.0;
            if (dect.abs()) > 1.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            dect = f64::asin(dect);
            rat = f64::cos(dect);
            if (rat.abs()) < deps {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            rat = l * (1.0 + dz) / (2.0 * rat);
            if (rat.abs()) > 1.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            rat = f64::asin(rat);
            mg = 1.0 + f64::sin(dect) * sin0 + f64::cos(dect) * cos0 * f64::cos(rat);
            if (mg.abs()) < deps {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            mg = 2.0 * (f64::sin(dect) * cos0 - f64::cos(dect) * sin0 * f64::cos(rat)) / mg;
            if (mg - m).abs() > deps {
                rat = TWOPI / 2.0 - rat;
            }
            rat += ra0;
        } else {
            *status = BAD_WCS_PROJ;
            return *status;
        }
    } else if ptype[cptr] == bb(b'A') {
        if ptype[cptr + 1] == bb(b'R') && ptype[cptr + 2] == bb(b'C') {
            /* ARC */
            if sins >= TWOPI * TWOPI / 4.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            sins = f64::sqrt(sins);
            coss = f64::cos(sins);
            if sins != 0.0 {
                sins = f64::sin(sins) / sins;
            } else {
                sins = 1.0;
            }
            dt = m * cos0 * sins + sin0 * coss;
            if (dt > 1.0) || (dt < -1.0) {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            dect = f64::asin(dt);
            da = coss - dt * sin0;
            dt = l * sins * cos0;
            if (da == 0.0) && (dt == 0.0) {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            rat = ra0 + f64::atan2(dt, da);
        } else if ptype[cptr + 1] == bb(b'I') && ptype[cptr + 2] == bb(b'T') {
            /* -AIT Aitoff */
            dt = yinc * cosr + xinc * sinr;
            if dt == 0.0 {
                dt = 1.0;
            }
            dt *= D2R;
            dy = yref * D2R;
            dx = f64::sin(dy + dt) / f64::sqrt((1.0 + f64::cos(dy + dt)) / 2.0)
                - f64::sin(dy) / f64::sqrt((1.0 + f64::cos(dy)) / 2.0);
            if dx == 0.0 {
                dx = 1.0;
            }
            geo2 = dt / dx;
            dt = xinc * cosr - yinc * sinr;
            if dt == 0.0 {
                dt = 1.0;
            }
            dt *= D2R;
            dx = 2.0 * f64::cos(dy) * f64::sin(dt / 2.0);
            if dx == 0.0 {
                dx = 1.0;
            }
            geo1 = dt * f64::sqrt((1.0 + f64::cos(dy) * f64::cos(dt / 2.0)) / 2.0) / dx;
            geo3 = geo2 * f64::sin(dy) / f64::sqrt((1.0 + f64::cos(dy)) / 2.0);
            rat = ra0;
            dect = dec0;
            if (l != 0.0) || (m != 0.0) {
                dz = 4.0 - l * l / (4.0 * geo1 * geo1) - ((m + geo3) / geo2) * ((m + geo3) / geo2);
                if (dz > 4.0) || (dz < 2.0) {
                    *status = ANGLE_TOO_BIG;
                    return *status;
                }
                dz = 0.5 * f64::sqrt(dz);
                dd = (m + geo3) * dz / geo2;
                if (dd.abs()) > 1.0 {
                    *status = ANGLE_TOO_BIG;
                    return *status;
                }
                dd = f64::asin(dd);
                if (f64::cos(dd)).abs() < deps {
                    *status = ANGLE_TOO_BIG;
                    return *status;
                }
                da = l * dz / (2.0 * geo1 * f64::cos(dd));
                if (da.abs()) > 1.0 {
                    *status = ANGLE_TOO_BIG;
                    return *status;
                }
                da = f64::asin(da);
                rat = ra0 + 2.0 * da;
                dect = dd;
            }
        } else {
            *status = BAD_WCS_PROJ;
            return *status;
        }
    } else if ptype[cptr] == bb(b'N') {
        /* -NCP North celestial pole*/
        if ptype[cptr + 1] != bb(b'C') || ptype[cptr + 2] != bb(b'P') {
            *status = BAD_WCS_PROJ;
            return *status;
        }
        dect = cos0 - m * sin0;
        if dect == 0.0 {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        rat = ra0 + f64::atan2(l, dect);
        dt = f64::cos(rat - ra0);
        if dt == 0.0 {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        dect /= dt;
        if (dect > 1.0) || (dect < -1.0) {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        dect = f64::acos(dect);
        if dec0 < 0.0 {
            dect = -dect;
        }
    } else if ptype[cptr] == bb(b'G') {
        /* -GLS global sinusoid */
        if ptype[cptr + 1] != bb(b'L') || ptype[cptr + 2] != bb(b'S') {
            *status = BAD_WCS_PROJ;
            return *status;
        }
        dect = dec0 + m;
        if (dect.abs()) > TWOPI / 4.0 {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        coss = f64::cos(dect);
        if (l.abs()) > TWOPI * coss / 2.0 {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        rat = ra0;
        if coss > deps {
            rat += l / coss;
        }
    } else if ptype[cptr] == bb(b'M') {
        /* -MER mercator*/
        if ptype[cptr + 1] != bb(b'E') || ptype[cptr + 2] != bb(b'R') {
            *status = BAD_WCS_PROJ;
            return *status;
        }
        dt = yinc * cosr + xinc * sinr;
        if dt == 0.0 {
            dt = 1.0;
        }
        dy = (yref / 2.0 + 45.0) * D2R;
        dx = dy + dt / 2.0 * D2R;
        dy = f64::ln(f64::tan(dy));
        dx = f64::ln(f64::tan(dx));
        geo2 = dt * D2R / (dx - dy);
        geo3 = geo2 * dy;
        geo1 = f64::cos(yref * D2R);
        if geo1 <= 0.0 {
            geo1 = 1.0;
        }
        rat = l / geo1 + ra0;
        if (rat - ra0).abs() > TWOPI {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        dt = 0.0;
        if geo2 != 0.0 {
            dt = (m + geo3) / geo2;
        }
        dt = dt.exp();
        dect = 2.0 * f64::atan(dt) - TWOPI / 4.0;
    } else {
        *status = BAD_WCS_PROJ;
        return *status;
    }

    /*  correct for RA rollover  */
    if rat - ra0 > TWOPI / 2.0 {
        rat -= TWOPI;
    }
    if rat - ra0 < -TWOPI / 2.0 {
        rat += TWOPI;
    }
    if rat < 0.0 {
        rat += TWOPI;
    }

    /*  convert to degrees  */
    *xpos = rat / D2R;
    *ypos = dect / D2R;
    *status
}

/*--------------------------------------------------------------------------*/
/// This routine is based on the classic AIPS WCS routine.
///
/// It converts from RA,Dec to pixel location to for 9 projective geometries:
/// "-CAR", "-SIN", "-TAN", "-ARC", "-NCP", "-GLS", "-MER", "-AIT" and "-STG".
///
/// routine to determine accurate pixel coordinates for an RA and Dec     
/// returns 0 if successful otherwise:                                    
/// 501 = angle too large for projection;                                 
/// 502 = bad values                                                      
/// does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections            
/// anything else is linear                                               
/// Input:                                                                
///   d   xpos    x (RA) coordinate (deg)                                 
///   d   ypos    y (dec) coordinate (deg)                                
///   d   xref    x reference coordinate value (deg)                      
///   d   yref    y reference coordinate value (deg)                      
///   f   xrefpix x reference pixel                                       
///   f   yrefpix y reference pixel                                       
///   f   xinc    x coordinate increment (deg)                            
///   f   yinc    y coordinate increment (deg)                            
///   f   rot     rotation (deg)  (from N through E)                      
///   c  *ptype    projection ptype code e.g. "-SIN";                       
/// Output:                                                               
///   f  *xpix    x pixel number  (RA or long without rotation)           
///   f  *ypiy    y pixel number  (dec or lat without rotation)           
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffxypx(
    xpos: f64,
    ypos: f64,
    xref: f64,
    yref: f64,
    xrefpix: f64,
    yrefpix: f64,
    xinc: f64,
    yinc: f64,
    rot: f64,
    ptype: *const [c_char; 5],
    xpix: *mut f64,
    ypix: *mut f64,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let xpix = xpix.as_mut().expect(NULL_MSG);
        let ypix = ypix.as_mut().expect(NULL_MSG);
        let ptype = ptype.as_ref().expect(NULL_MSG);

        ffxypx_safe(
            xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, ptype, xpix, ypix, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// This routine is based on the classic AIPS WCS routine.
///
/// It converts from RA,Dec to pixel location to for 9 projective geometries:
/// "-CAR", "-SIN", "-TAN", "-ARC", "-NCP", "-GLS", "-MER", "-AIT" and "-STG".
///
/// routine to determine accurate pixel coordinates for an RA and Dec     
/// returns 0 if successful otherwise:                                    
/// 501 = angle too large for projection;                                 
/// 502 = bad values                                                      
/// does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections            
/// anything else is linear                                               
/// Input:                                                                
///   d   xpos    x (RA) coordinate (deg)                                 
///   d   ypos    y (dec) coordinate (deg)                                
///   d   xref    x reference coordinate value (deg)                      
///   d   yref    y reference coordinate value (deg)                      
///   f   xrefpix x reference pixel                                       
///   f   yrefpix y reference pixel                                       
///   f   xinc    x coordinate increment (deg)                            
///   f   yinc    y coordinate increment (deg)                            
///   f   rot     rotation (deg)  (from N through E)                      
///   c  *ptype    projection ptype code e.g. "-SIN";                       
/// Output:                                                               
///   f  *xpix    x pixel number  (RA or long without rotation)           
///   f  *ypiy    y pixel number  (dec or lat without rotation)       
pub(crate) fn ffxypx_safe(
    xpos: f64,
    ypos: f64,
    xref: f64,
    yref: f64,
    xrefpix: f64,
    yrefpix: f64,
    xinc: f64,
    yinc: f64,
    rot: f64,
    ptype: &[c_char; 5],
    xpix: &mut f64,
    ypix: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut dx: f64 = 0.0;
    let mut dy: f64 = 0.0;
    let mut dz: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut ra0: f64 = 0.0;
    let mut dec0: f64 = 0.0;
    let mut ra: f64 = 0.0;
    let mut dec: f64 = 0.0;
    let mut coss: f64 = 0.0;
    let mut sins: f64 = 0.0;
    let mut dt: f64 = 0.0;
    let mut da: f64 = 0.0;
    let mut dd: f64 = 0.0;
    let mut sint: f64 = 0.0;
    let mut l: f64 = 0.0;
    let mut m: f64 = 0.0;
    let mut geo1: f64 = 0.0;
    let mut geo2: f64 = 0.0;
    let mut geo3: f64 = 0.0;
    let mut sinr: f64 = 0.0;
    let mut cosr: f64 = 0.0;
    let mut cos0: f64 = 0.0;
    let mut sin0: f64 = 0.0;
    let deps: f64 = 1.0e-5;

    let mut xpos = xpos;

    if ptype[0] != bb(b'-') {
        /* unrecognized projection code */
        *status = BAD_WCS_PROJ;
        return *status;
    }

    let cptr = 1;

    dt = xpos - xref;
    if dt > 180.0 {
        xpos -= 360.0;
    }
    if dt < -180.0 {
        xpos += 360.0;
    }
    /* NOTE: changing input argument xpos is OK (call-by-value in C!) */

    /* default values - linear */
    dx = xpos - xref;
    dy = ypos - yref;

    /*  Correct for rotation */
    r = rot * D2R;
    cosr = f64::cos(r);
    sinr = f64::sin(r);
    dz = dx * cosr + dy * sinr;
    dy = dy * cosr - dx * sinr;
    dx = dz;

    /*     check axis increments - bail out if either 0 */
    if (xinc == 0.0) || (yinc == 0.0) {
        *xpix = 0.0;
        *ypix = 0.0;
        *status = BAD_WCS_VAL;
        return *status;
    }

    /*     convert to pixels  */
    *xpix = dx / xinc + xrefpix;
    *ypix = dy / yinc + yrefpix;

    if ptype[cptr] == bb(b'C') {
        /* linear -CAR */
        if ptype[cptr + 1] != bb(b'A') || ptype[cptr + 2] != bb(b'R') {
            *status = BAD_WCS_PROJ;
            return *status;
        }

        return *status; /* done if linear */
    }

    /* Non linear position */
    ra0 = xref * D2R;
    dec0 = yref * D2R;
    ra = xpos * D2R;
    dec = ypos * D2R;

    /* compute direction cosine */
    coss = f64::cos(dec);
    sins = f64::sin(dec);
    cos0 = f64::cos(dec0);
    sin0 = f64::sin(dec0);
    l = f64::sin(ra - ra0) * coss;
    sint = sins * sin0 + coss * cos0 * f64::cos(ra - ra0);

    /* process by case  */
    if ptype[cptr] == bb(b'T') {
        /* -TAN f64::tan */
        if ptype[cptr + 1] != bb(b'A') || ptype[cptr + 2] != bb(b'N') {
            *status = BAD_WCS_PROJ;
            return *status;
        }

        if sint <= 0.0 {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        if cos0 < 0.001 {
            /* Do a first order expansion around pole */
            m = (coss * f64::cos(ra - ra0)) / (sins * sin0);
            m = (-m + cos0 * (1.0 + m * m)) / sin0;
        } else {
            m = (sins / sint - sin0) / cos0;
        }
        if (f64::sin(ra0)).abs() < 0.3 {
            l = coss * f64::sin(ra) / sint - cos0 * f64::sin(ra0) + m * f64::sin(ra0) * sin0;
            l /= f64::cos(ra0);
        } else {
            l = coss * f64::cos(ra) / sint - cos0 * f64::cos(ra0) + m * f64::cos(ra0) * sin0;
            l /= -f64::sin(ra0);
        }
    } else if ptype[cptr] == bb(b'S') {
        if ptype[cptr + 1] == bb(b'I') && ptype[cptr + 2] == bb(b'N') {
            /* -SIN */
            if sint < 0.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            m = sins * f64::cos(dec0) - coss * f64::sin(dec0) * f64::cos(ra - ra0);
        } else if ptype[cptr + 1] == bb(b'T') && ptype[cptr + 2] == bb(b'G') {
            /* -STG Sterographic*/
            da = ra - ra0;
            if (dec.abs()) > TWOPI / 4.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            dd = 1.0 + sins * f64::sin(dec0) + coss * f64::cos(dec0) * f64::cos(da);
            if (dd.abs()) < deps {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            dd = 2.0 / dd;
            l *= dd;
            m = dd * (sins * f64::cos(dec0) - coss * f64::sin(dec0) * f64::cos(da));
        } else {
            *status = BAD_WCS_PROJ;
            return *status;
        }
    } else if ptype[cptr] == bb(b'A') {
        if ptype[cptr + 1] == bb(b'R') && ptype[cptr + 2] == bb(b'C') {
            /* ARC */
            m = sins * f64::sin(dec0) + coss * f64::cos(dec0) * f64::cos(ra - ra0);

            m = m.clamp(-1.0, 1.0);

            m = f64::acos(m);
            if m != 0.0 {
                m = m / f64::sin(m);
            } else {
                m = 1.0;
            }
            l *= m;
            m *= sins * f64::cos(dec0) - coss * f64::sin(dec0) * f64::cos(ra - ra0);
        } else if ptype[cptr + 1] == bb(b'I') && ptype[cptr + 2] == bb(b'T') {
            /* -AIT Aitoff */
            da = (ra - ra0) / 2.0;
            if (da.abs()) > TWOPI / 4.0 {
                *status = ANGLE_TOO_BIG;
                return *status;
            }
            dt = yinc * cosr + xinc * sinr;
            if dt == 0.0 {
                dt = 1.0;
            }
            dt *= D2R;
            dy = yref * D2R;
            dx = f64::sin(dy + dt) / f64::sqrt((1.0 + f64::cos(dy + dt)) / 2.0)
                - f64::sin(dy) / f64::sqrt((1.0 + f64::cos(dy)) / 2.0);
            if dx == 0.0 {
                dx = 1.0;
            }
            geo2 = dt / dx;
            dt = xinc * cosr - yinc * sinr;
            if dt == 0.0 {
                dt = 1.0;
            }
            dt *= D2R;
            dx = 2.0 * f64::cos(dy) * f64::sin(dt / 2.0);
            if dx == 0.0 {
                dx = 1.0;
            }
            geo1 = dt * f64::sqrt((1.0 + f64::cos(dy) * f64::cos(dt / 2.0)) / 2.0) / dx;
            geo3 = geo2 * f64::sin(dy) / f64::sqrt((1.0 + f64::cos(dy)) / 2.0);
            dt = f64::sqrt((1.0 + f64::cos(dec) * f64::cos(da)) / 2.0);
            if (dt.abs()) < deps {
                *status = WCS_ERROR;
                return *status;
            }
            l = 2.0 * geo1 * f64::cos(dec) * f64::sin(da) / dt;
            m = geo2 * f64::sin(dec) / dt - geo3;
        } else {
            *status = BAD_WCS_PROJ;
            return *status;
        }
    } else if ptype[cptr] == bb(b'N') {
        /* -NCP North celestial pole*/
        if ptype[cptr + 1] != bb(b'C') || ptype[cptr + 2] != bb(b'P') {
            *status = BAD_WCS_PROJ;
            return *status;
        }

        if dec0 == 0.0 {
            *status = ANGLE_TOO_BIG;
            return *status; /* can't stand the equator */
        } else {
            m = (f64::cos(dec0) - coss * f64::cos(ra - ra0)) / f64::sin(dec0);
        }
    } else if ptype[cptr] == bb(b'G') {
        /* -GLS global sinusoid */
        if ptype[cptr + 1] != bb(b'L') || ptype[cptr + 2] != bb(b'S') {
            *status = BAD_WCS_PROJ;
            return *status;
        }

        dt = ra - ra0;
        if (dec.abs()) > TWOPI / 4.0 {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        if (dec0.abs()) > TWOPI / 4.0 {
            *status = ANGLE_TOO_BIG;
            return *status;
        }
        m = dec - dec0;
        l = dt * coss;
    } else if ptype[cptr] == bb(b'M') {
        /* -MER mercator*/
        if ptype[cptr + 1] != bb(b'E') || ptype[cptr + 2] != bb(b'R') {
            *status = BAD_WCS_PROJ;
            return *status;
        }

        dt = yinc * cosr + xinc * sinr;
        if dt == 0.0 {
            dt = 1.0;
        }
        dy = (yref / 2.0 + 45.0) * D2R;
        dx = dy + dt / 2.0 * D2R;
        dy = f64::ln(f64::tan(dy));
        dx = f64::ln(f64::tan(dx));
        geo2 = dt * D2R / (dx - dy);
        geo3 = geo2 * dy;
        geo1 = f64::cos(yref * D2R);
        if geo1 <= 0.0 {
            geo1 = 1.0;
        }
        dt = ra - ra0;
        l = geo1 * dt;
        dt = dec / 2.0 + TWOPI / 8.0;
        dt = f64::tan(dt);
        if dt < deps {
            *status = BAD_WCS_VAL;
            return *status;
        }
        m = geo2 * f64::ln(dt) - geo3;
    } else {
        *status = BAD_WCS_PROJ;
        return *status;
    }

    /*   convert to degrees  */
    dx = l / D2R;
    dy = m / D2R;

    /*  Correct for rotation */
    dz = dx * cosr + dy * sinr;
    dy = dy * cosr - dx * sinr;
    dx = dz;

    /*     convert to pixels  */
    *xpix = dx / xinc + xrefpix;
    *ypix = dy / yinc + yrefpix;
    *status
}
