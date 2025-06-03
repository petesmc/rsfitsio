use std::cmp;
use std::ffi::CStr;
use std::fs::File;
use std::io::Read;

use crate::c_types::*;
use crate::helpers::boxed::box_try_new;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::NullValue;
use crate::aliases::rust_api::{fits_clear_errmark, fits_write_errmark};
use crate::cfileio::{ffclos_safer, ffopen_safer};
use crate::fitscore::{
    ffgcno_safe, ffmnhd_safe, ffpmsg_slice, ffpmsg_str, fits_strcasecmp, fits_strncasecmp,
};
use crate::fitsio::*;
use crate::getcol::ffgcv_safe;
use crate::getcold::ffgcvd_safe;
use crate::getcols::ffgcvs_safe;
use crate::getkey::{ffgky_safe, ffgtdm_safe};
use crate::wcssub::ffgtcs_safer;
use crate::wcsutil::{ffwldp_safe, ffxypx_safe};
use crate::{KeywordDatatypeMut, bb, cs};
use crate::{int_snprintf, wrappers::*};

#[allow(clippy::approx_constant)]
const MY_PI: f64 = 3.141_592_653_589_793;
const RAD_TO_DEG: f64 = 180.0 / MY_PI;

#[derive(Default, Debug, Copy, Clone)]
pub(crate) struct WCSdata {
    exists: bool,
    xrefval: f64,
    yrefval: f64,
    xrefpix: f64,
    yrefpix: f64,
    xinc: f64,
    yinc: f64,
    rot: f64,
    dtype: [c_char; 5],
}

#[derive(Default, Debug, PartialEq, Clone)]
enum ShapeType {
    #[default]
    Point,
    Line,
    Circle,
    Annulus,
    Ellipse,
    ElliptAnnulus,
    Box,
    BoxAnnulus,
    Rectangle,
    Diamond,
    Sector,
    Poly,
    Panda,
    EPanda,
    BPanda,
}

#[derive(Debug, PartialEq, Clone)]
enum CoordFmt {
    Pixel,
    Degree,

    #[allow(clippy::upper_case_acronyms)]
    HHMMSS,
}

#[derive(Default, Debug, Clone)]
pub(crate) struct RgnShape {
    sign: c_char,     /*  Include or exclude?        */
    shape: ShapeType, /*  Shape of this region       */
    comp: c_int,      /*  Component number for this region */

    /*  bounding box    */
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,

    /*  Parameters - In pixels     */
    genericParams: RgnShapeGeneric,
    polyParams: RgnShapePolygon,
}

#[derive(Default, Debug, Clone, Copy)]
struct RgnShapeGeneric {
    p: [f64; 11], /*  Region parameters       */

    /*  For rotated shapes      */
    sinT: f64,
    cosT: f64,

    /*  Extra scratch area      */
    a: f64,
    b: f64,
}

#[derive(Default, Debug, Clone)]
struct RgnShapePolygon {
    nPts: c_int,   /*  Number of Polygon pts   */
    Pts: Vec<f64>, /*  Polygon points          */
}

#[derive(Default, Debug)]
pub(crate) struct SAORegion {
    nShapes: c_int,
    Shapes: Vec<RgnShape>,
    wcs: WCSdata,
}

/*---------------------------------------------------------------------------*/
/// Read regions from either a FITS or ASCII region file and return the information
/// in the "SAORegion" structure.  If it is nonNULL, use wcs to convert the  
/// region coordinates to pixels.  Return an error if region is in degrees   
/// but no WCS data is provided.                                             
pub(crate) fn fits_read_rgnfile(
    filename: &[c_char],
    wcs: &mut WCSdata,
    Rgn: &mut Option<Box<SAORegion>>,
    status: &mut c_int,
) -> c_int {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut tstatus = 0;

    if *status != 0 {
        return *status;
    }

    /* try to open as a FITS file - if that doesn't work treat as an ASCII file */

    fits_write_errmark();

    // SAFETY: TODO
    let retval = unsafe { ffopen_safer(&mut fptr, filename, READONLY, &mut tstatus) };
    if retval != 0 {
        fits_clear_errmark();
        fits_read_ascii_region(filename, wcs, Rgn, status);
    } else {
        // SAFETY: TODO
        unsafe {
            let fptr = fptr.unwrap();
            fits_read_fits_region(fptr, wcs, Rgn, status);
        }
    }

    *status
}

/*---------------------------------------------------------------------------*/
/// Read regions from a SAO-style region file and return the information     
/// in the "SAORegion" structure.  If it is nonNULL, use wcs to convert the  
/// region coordinates to pixels.  Return an error if region is in degrees   
/// but no WCS data is provided.                                             
pub(crate) fn fits_read_ascii_region(
    filename: &[c_char],
    wcs: &mut WCSdata,
    Rgn: &mut Option<Box<SAORegion>>,
    status: &mut c_int,
) -> c_int {
    let mut namePtr: usize;
    let paramPtrSet = false;
    let mut paramPtr: usize;
    let mut pX: usize;
    let mut pY: usize;
    let endp = 0;

    let mut lineLen: usize;
    let mut hh: c_long;
    let mut mm: c_long;
    let mut dd: c_long;
    let mut coords: &mut [f64];
    let X: f64 = 0.0;
    let Y: f64 = 0.0;
    let x: f64 = 0.0;
    let y: f64 = 0.0;
    let ss: f64 = 0.0;
    let div: f64 = 0.0;
    let xsave: f64 = 0.;
    let ysave: f64 = 0.;
    let mut nParams: c_int;
    let mut nCoords: c_int;
    let mut negdec: c_int;
    let i: c_int = 0;
    let mut done: bool;

    let mut newShape: &mut RgnShape;
    let mut tmpShape: Box<RgnShape>;

    if *status != 0 {
        return *status;
    }

    let aRgn = box_try_new(SAORegion::default());
    if aRgn.is_err() {
        ffpmsg_str("Couldn't allocate memory to hold Region file contents.");
        *status = MEMORY_ALLOCATION;
        return *status;
    }
    let mut aRgn: Box<SAORegion> = aRgn.unwrap();

    aRgn.nShapes = 0;
    // aRgn.Shapes = NULL;

    if wcs.exists {
        aRgn.wcs = *wcs;
    } else {
        aRgn.wcs.exists = false;
    }

    let cFmt: CoordFmt = CoordFmt::Pixel; /* set default format */

    /*  Allocate Line Buffer  */

    let allocLen: usize = 512;
    let mut currLine: Vec<c_char> = Vec::new();
    if currLine.try_reserve_exact(allocLen).is_err() {
        ffpmsg_str("Couldn't allocate memory to hold Region file contents.");
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        currLine.resize(allocLen, 0);
    }

    /*  Open Region File  */
    let filename = CStr::from_bytes_until_nul(cast_slice(filename))
        .unwrap()
        .to_str()
        .unwrap();

    let rgnFile = File::options().read(true).write(false).open(filename);

    if rgnFile.is_err() {
        int_snprintf!(
            &mut currLine,
            allocLen,
            "Could not open Region file {}.",
            filename
        );
        ffpmsg_slice(&currLine);
        (*status = FILE_NOT_OPENED);
        return *status;
    }

    let rgnFile = rgnFile.unwrap();

    /*  Read in file, line by line  */
    /*  First, set error status in case file is empty */
    *status = FILE_NOT_OPENED;

    todo!(
        "Below code was badly translated from C. It previously used 'fgets' to read a line from a file. Rusts 'read' does not do this."
    );

    while rgnFile
        .read(cast_slice_mut(&mut currLine))
        .is_ok_and(|x| x > 0)
    {
        /* reset status if we got here */
        *status = 0;

        /*  Make sure we have a full line of text  */

        lineLen = strlen_safe(&currLine);
        while lineLen == allocLen - 1 && currLine[lineLen - 1] != bb(b'\n') {
            // TODO currLoc = (char *)realloc(currLine, 2 * allocLen * sizeof(char));
            if currLine.try_reserve_exact(allocLen).is_err() {
                ffpmsg_str("Couldn't allocate memory to hold Region file contents.");
                *status = MEMORY_ALLOCATION;
                fits_free_region(aRgn);
                return *status;
            } else {
                currLine.resize(2 * allocLen, 0);
            }

            let _ = rgnFile.read(cast_slice_mut(&mut currLine[lineLen..]));
            allocLen += allocLen;
            lineLen += strlen_safe(&currLine[lineLen..]);
        }

        let mut currLoc = 0;
        if currLine[currLoc] == bb(b'#') {
            /*  Look to see if it is followed by a format statement...  */
            /*  if not skip line                                        */

            currLoc += 1;
            while isspace(currLine[currLoc]) {
                currLoc += 1;
            }
            if fits_strncasecmp(&currLine[currLoc..], cs!(c"format:"), 7) == 0 {
                if aRgn.nShapes != 0 {
                    ffpmsg_str("Format code encountered after reading 1 or more shapes.");
                    *status = PARSE_SYNTAX_ERR;
                    fits_free_region(aRgn);
                    return *status;
                }
                currLoc += 7;
                while isspace(currLine[currLoc]) {
                    currLoc += 1;
                }
                if fits_strncasecmp(&currLine[currLoc..], cs!(c"pixel"), 5) == 0 {
                    cFmt = CoordFmt::Pixel;
                } else if fits_strncasecmp(&currLine[currLoc..], cs!(c"degree"), 6) == 0 {
                    cFmt = CoordFmt::Degree;
                } else if fits_strncasecmp(&currLine[currLoc..], cs!(c"hhmmss"), 6) == 0 {
                    cFmt = CoordFmt::HHMMSS;
                } else if fits_strncasecmp(&currLine[currLoc..], cs!(c"hms"), 3) == 0 {
                    cFmt = CoordFmt::HHMMSS;
                } else {
                    ffpmsg_str("Unknown format code encountered in region file.");
                    *status = PARSE_SYNTAX_ERR;
                    fits_free_region(aRgn);
                    return *status;
                }
            }
        } else if fits_strncasecmp(&currLine[currLoc..], cs!(c"glob"), 4) == 0 {
            /* skip lines that begin with the word 'global' */
        } else {
            while currLine[currLoc] != 0 {
                namePtr = currLoc;
                paramPtr = 0;
                nParams = 1;

                /*  Search for closing parenthesis  */

                done = false;
                while !done && *status == 0 && currLine[currLoc] != 0 {
                    match currLine[currLoc] as u8 {
                        b'(' => {
                            currLine[currLoc] = 0;
                            currLoc += 1;
                            if paramPtrSet {
                                /* Can't have two '(' in a region! */
                                *status = 1;
                            } else {
                                paramPtrSet = true;
                                paramPtr = currLoc;
                            }
                            break;
                        }
                        b')' => {
                            currLine[currLoc] = 0;
                            currLoc += 1;
                            if !paramPtrSet {
                                /* Can't have a ')' without a '(' first */
                                *status = 1;
                            } else {
                                done = true;
                            }
                            break;
                        }
                        b'#' | b'\n' => {
                            currLine[currLoc] = 0;
                            if !paramPtrSet {
                                /* Allow for a blank line */
                                done = true;
                            }
                            break;
                        }
                        b':' => {
                            currLoc += 1;
                            if paramPtrSet {
                                cFmt = CoordFmt::HHMMSS; /* set format if parameter has : */
                            }
                            break;
                        }
                        b'd' => {
                            currLoc += 1;
                            if paramPtrSet {
                                cFmt = CoordFmt::Degree; /* set format if parameter has d */
                            }
                            break;
                        }
                        b',' => {
                            nParams += 1;
                            currLoc += 1;
                        }
                        _ => {
                            currLoc += 1;
                            break;
                        }
                    }
                }

                if *status != 0 || !done {
                    ffpmsg_str("Error reading Region file");
                    *status = PARSE_SYNTAX_ERR;
                    fits_free_region(aRgn);
                    return *status;
                }

                /*  Skip white space in region name  */

                while isspace(currLine[namePtr]) {
                    namePtr += 1;
                }

                /*  Was this a blank line? Or the end of the current one  */

                if currLine[namePtr] == 0 && !paramPtrSet {
                    continue;
                }

                /*  Check for format code at beginning of the line */

                if fits_strncasecmp(&currLine[namePtr..], cs!(c"image;"), 6) == 0 {
                    namePtr += 6;
                    cFmt = CoordFmt::Pixel;
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"physical;"), 9) == 0 {
                    namePtr += 9;
                    cFmt = CoordFmt::Pixel;
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"linear;"), 7) == 0 {
                    namePtr += 7;
                    cFmt = CoordFmt::Pixel;
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"fk4;"), 4) == 0 {
                    namePtr += 4;
                    cFmt = CoordFmt::Degree;
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"fk5;"), 4) == 0 {
                    namePtr += 4;
                    cFmt = CoordFmt::Degree;
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"icrs;"), 5) == 0 {
                    namePtr += 5;
                    cFmt = CoordFmt::Degree;

                    /* the following 5 cases support region files created by POW
                    (or ds9 Version 4.x) which
                    may have lines containing  only a format code, not followed
                    by a ';' (and with no region specifier on the line).  We use
                    the 'continue' statement to jump to the end of the loop and
                    then continue reading the next line of the region file. */
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"fk5"), 3) == 0 {
                    cFmt = CoordFmt::Degree;
                    continue; /* supports POW region file format */
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"fk4"), 3) == 0 {
                    cFmt = CoordFmt::Degree;
                    continue; /* supports POW region file format */
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"icrs"), 4) == 0 {
                    cFmt = CoordFmt::Degree;
                    continue; /* supports POW region file format */
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"image"), 5) == 0 {
                    cFmt = CoordFmt::Pixel;
                    continue; /* supports POW region file format */
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"physical"), 8) == 0 {
                    cFmt = CoordFmt::Pixel;
                    continue; /* supports POW region file format */
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"galactic;"), 9) == 0 {
                    ffpmsg_str("Galactic region coordinates not supported");
                    ffpmsg_slice(&currLine[namePtr..]);
                    *status = PARSE_SYNTAX_ERR;
                    fits_free_region(aRgn);
                    return *status;
                } else if fits_strncasecmp(&currLine[namePtr..], cs!(c"ecliptic;"), 9) == 0 {
                    ffpmsg_str("ecliptic region coordinates not supported");
                    ffpmsg_slice(&currLine[namePtr..]);
                    *status = PARSE_SYNTAX_ERR;
                    fits_free_region(aRgn);
                    return *status;
                }

                /**************************************************/
                /*  We've apparently found a region... Set it up  */
                /**************************************************/

                if (aRgn.nShapes % 10) == 0 {
                    let len = aRgn.Shapes.len();

                    if aRgn.Shapes.try_reserve_exact(10).is_err() {
                        ffpmsg_str("Failed to allocate memory for Region data");
                        *status = MEMORY_ALLOCATION;
                        fits_free_region(aRgn);
                        return *status;
                    } else {
                        aRgn.Shapes.resize(len + 10, RgnShape::default());
                    }
                }
                newShape = &mut aRgn.Shapes[aRgn.nShapes as usize];
                aRgn.nShapes += 1;
                newShape.sign = 1;
                newShape.shape = ShapeType::Point;
                for i in 0..8 {
                    newShape.genericParams.p[i] = 0.0;
                }
                newShape.genericParams.a = 0.0;
                newShape.genericParams.b = 0.0;
                newShape.genericParams.sinT = 0.0;
                newShape.genericParams.cosT = 0.0;

                while isspace(currLine[namePtr]) {
                    namePtr += 1;
                }

                /*  Check for the shape's sign  */

                if currLine[namePtr] == bb(b'+') {
                    namePtr += 1;
                } else if currLine[namePtr] == bb(b'-') {
                    namePtr += 1;
                    newShape.sign = 0;
                }

                /* Skip white space in region name */

                while isspace(currLine[namePtr]) {
                    namePtr += 1;
                }
                if currLine[namePtr] == 0 {
                    ffpmsg_str("Error reading Region file");
                    *status = PARSE_SYNTAX_ERR;
                    fits_free_region(aRgn);
                    return *status;
                }
                lineLen = strlen_safe(&currLine[namePtr..]) - 1;
                while isspace(currLine[namePtr + lineLen]) {
                    currLine[namePtr + lineLen] = 0;
                    lineLen -= 1;
                }

                /*  Now identify the region  */

                if fits_strcasecmp(&currLine[namePtr..], cs!(c"circle")) == 0 {
                    newShape.shape = ShapeType::Circle;
                    if nParams != 3 {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"annulus")) == 0 {
                    newShape.shape = ShapeType::Annulus;
                    if nParams != 4 {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"ellipse")) == 0 {
                    if !(4..=8).contains(&nParams) {
                        *status = PARSE_SYNTAX_ERR;
                    } else if nParams < 6 {
                        newShape.shape = ShapeType::Ellipse;
                        newShape.genericParams.p[4] = 0.0;
                    } else {
                        newShape.shape = ShapeType::ElliptAnnulus;
                        newShape.genericParams.p[6] = 0.0;
                        newShape.genericParams.p[7] = 0.0;
                    }
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"elliptannulus")) == 0 {
                    newShape.shape = ShapeType::ElliptAnnulus;
                    if !(nParams == 8 || nParams == 6) {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    newShape.genericParams.p[6] = 0.0;
                    newShape.genericParams.p[7] = 0.0;
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"box")) == 0
                    || fits_strcasecmp(&currLine[namePtr..], cs!(c"rotbox")) == 0
                {
                    if !(4..=8).contains(&nParams) {
                        *status = PARSE_SYNTAX_ERR;
                    } else if nParams < 6 {
                        newShape.shape = ShapeType::Box;
                        newShape.genericParams.p[4] = 0.0;
                    } else {
                        newShape.shape = ShapeType::BoxAnnulus;
                        newShape.genericParams.p[6] = 0.0;
                        newShape.genericParams.p[7] = 0.0;
                    }
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"rectangle")) == 0
                    || fits_strcasecmp(&currLine[namePtr..], cs!(c"rotrectangle")) == 0
                {
                    newShape.shape = ShapeType::Rectangle;
                    if !(4..=5).contains(&nParams) {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    newShape.genericParams.p[4] = 0.0;
                    nCoords = 4;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"diamond")) == 0
                    || fits_strcasecmp(&currLine[namePtr..], cs!(c"rotdiamond")) == 0
                    || fits_strcasecmp(&currLine[namePtr..], cs!(c"rhombus")) == 0
                    || fits_strcasecmp(&currLine[namePtr..], cs!(c"rotrhombus")) == 0
                {
                    newShape.shape = ShapeType::Diamond;
                    if !(4..=5).contains(&nParams) {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    newShape.genericParams.p[4] = 0.0;
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"sector")) == 0
                    || fits_strcasecmp(&currLine[namePtr..], cs!(c"pie")) == 0
                {
                    newShape.shape = ShapeType::Sector;
                    if nParams != 4 {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"point")) == 0 {
                    newShape.shape = ShapeType::Point;
                    if nParams != 2 {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"line")) == 0 {
                    newShape.shape = ShapeType::Line;
                    if nParams != 4 {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    nCoords = 4;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"polygon")) == 0 {
                    newShape.shape = ShapeType::Poly;
                    if nParams < 6 || (nParams & 1) != 0 {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    nCoords = nParams;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"panda")) == 0 {
                    newShape.shape = ShapeType::Panda;
                    if nParams != 8 {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"epanda")) == 0 {
                    newShape.shape = ShapeType::EPanda;
                    if !(10..=11).contains(&nParams) {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    newShape.genericParams.p[10] = 0.0;
                    nCoords = 2;
                } else if fits_strcasecmp(&currLine[namePtr..], cs!(c"bpanda")) == 0 {
                    newShape.shape = ShapeType::BPanda;
                    if !(10..=11).contains(&nParams) {
                        *status = PARSE_SYNTAX_ERR;
                    }
                    newShape.genericParams.p[10] = 0.0;
                    nCoords = 2;
                } else {
                    ffpmsg_str("Unrecognized region found in region file:");
                    ffpmsg_slice(&currLine[namePtr..]);
                    *status = PARSE_SYNTAX_ERR;
                    fits_free_region(aRgn);
                    return *status;
                }
                if *status != 0 {
                    ffpmsg_str("Wrong number of parameters found for region");
                    ffpmsg_slice(&currLine[namePtr..]);
                    fits_free_region(aRgn);
                    return *status;
                }

                /*  Parse Parameter string... convert to pixels if necessary  */

                if newShape.shape == ShapeType::Poly {
                    if newShape
                        .polyParams
                        .Pts
                        .try_reserve_exact(nParams as usize)
                        .is_err()
                    {
                        ffpmsg_str("Could not allocate memory to hold polygon parameters");
                        *status = MEMORY_ALLOCATION;
                        fits_free_region(aRgn);
                        return *status;
                    } else {
                        newShape.polyParams.Pts.resize(nParams as usize, 0.0);
                    }

                    newShape.polyParams.nPts = nParams;
                    coords = &mut newShape.polyParams.Pts;
                } else {
                    coords = &mut newShape.genericParams.p;
                }

                /*  Parse the initial "WCS?" coordinates  */
                for i in (0..nCoords).step_by(2) {
                    pX = paramPtr; // Into currLine
                    while currLine[paramPtr] != bb(b',') {
                        paramPtr += 1;
                    }
                    currLine[paramPtr] = 0;
                    paramPtr += 1;

                    pY = paramPtr; // Into currLine
                    while currLine[paramPtr] != bb(b',') && currLine[paramPtr] != 0 {
                        paramPtr += 1;
                    }

                    currLine[paramPtr] = 0;
                    paramPtr += 1;

                    if strchr_safe(&currLine[pX..], bb(b':')).is_some() {
                        /*  Read in special format & convert to decimal degrees  */
                        cFmt = CoordFmt::HHMMSS;
                        mm = 0;
                        ss = 0.0;
                        // hh = strtol_safe(&currLine[pX..], &mut endp, 10);
                        let (r, p) = strtol_safe(&currLine[pX..]).unwrap();
                        hh = r;
                        let mut endp = p;

                        if endp != 0 && currLine[pX + endp] == bb(b':') {
                            pX += endp + 1;
                            // mm = strtol_safe(&currLine[pX..], &mut endp, 10);
                            let (r, p) = strtol_safe(&currLine[pX..]).unwrap();
                            mm = r;
                            endp = p;

                            if endp != 0 && currLine[pX + endp] == bb(b':') {
                                pX += endp + 1;
                                ss = atof_safe(&currLine[pX..]);
                            }
                        }
                        X = 15.0 * (hh as f64 + (mm as f64 / 60.0) + ss / 3600.0); /* convert to degrees */

                        mm = 0;
                        ss = 0.0;
                        negdec = 0;

                        while isspace(currLine[pY]) {
                            pY += 1;
                        }
                        if currLine[pY] == bb(b'-') {
                            negdec = 1;
                            pY += 1;
                        }

                        // dd = strtol_safe(&currLine[pY..], &mut endp, 10);
                        let (r, p) = strtol_safe(&currLine[pY..]).unwrap();
                        dd = r;
                        endp = p;

                        if endp != 0 && currLine[pY + endp] == bb(b':') {
                            pY += endp + 1;
                            // mm = strtol_safe(&currLine[pY..], &mut endp, 10);
                            let (r, p) = strtol_safe(&currLine[pY..]).unwrap();
                            mm = r;
                            endp = p;

                            if endp != 0 && currLine[pY + endp] == bb(b':') {
                                pY += endp + 1;
                                ss = atof_safe(&currLine[pY..]);
                            }
                        }

                        if negdec != 0 {
                            Y = -dd as f64 - (mm as f64 / 60.0) - ss / 3600.0; /* convert to degrees */
                        } else {
                            Y = dd as f64 + (mm as f64 / 60.0) + ss / 3600.0;
                        }
                    } else {
                        X = atof_safe(&currLine[pX..]);
                        Y = atof_safe(&currLine[pY..]);
                    }

                    if i == 0 {
                        /* save 1st coord. in case needed later */
                        xsave = X;
                        ysave = Y;
                    }

                    if cFmt != CoordFmt::Pixel {
                        /*  Convert to pixels  */
                        if !wcs.exists {
                            ffpmsg_str("WCS information needed to convert region coordinates.");
                            *status = NO_WCS_KEY;
                            fits_free_region(aRgn);
                            return *status;
                        }

                        if ffxypx_safe(
                            X,
                            Y,
                            wcs.xrefval,
                            wcs.yrefval,
                            wcs.xrefpix,
                            wcs.yrefpix,
                            wcs.xinc,
                            wcs.yinc,
                            wcs.rot,
                            &wcs.dtype,
                            &mut x,
                            &mut y,
                            status,
                        ) != 0
                        {
                            ffpmsg_str("Error converting region to pixel coordinates.");
                            fits_free_region(aRgn);
                            return *status;
                        }
                        X = x;
                        Y = y;
                    }
                    coords[i as usize] = X;
                    coords[i as usize + 1] = Y;
                }

                /*  Read in remaining parameters...  */
                while i < nParams {
                    pX = paramPtr; // Into currLine
                    while currLine[paramPtr] != bb(b',') && currLine[paramPtr] != 0 {
                        paramPtr += 1;
                    }

                    currLine[paramPtr] = 0;
                    paramPtr += 1;

                    coords[i as usize] = strtod_safe(&currLine[pX..], &mut endp);

                    if endp != 0
                        && (currLine[pX + endp] == bb(b'"')
                            || currLine[pX + endp] == bb(b'\'')
                            || currLine[pX + endp] == bb(b'd'))
                    {
                        div = 1.0;
                        if currLine[pX + endp] == bb(b'"') {
                            div = 3600.0;
                        }
                        if currLine[pX + endp] == bb(b'\'') {
                            div = 60.0;
                        }
                        /* parameter given in arcsec so convert to pixels. */
                        /* Increment first Y coordinate by this amount then calc */
                        /* the distance in pixels from the original coordinate. */
                        /* NOTE: This assumes the pixels are square!! */
                        if ysave < 0.0 {
                            Y = ysave + coords[i as usize] / div; /* don't exceed -90 */
                        } else {
                            Y = ysave - coords[i as usize] / div; /* don't exceed +90 */
                        }

                        X = xsave;
                        if ffxypx_safe(
                            X,
                            Y,
                            wcs.xrefval,
                            wcs.yrefval,
                            wcs.xrefpix,
                            wcs.yrefpix,
                            wcs.xinc,
                            wcs.yinc,
                            wcs.rot,
                            &wcs.dtype,
                            &mut x,
                            &mut y,
                            status,
                        ) != 0
                        {
                            ffpmsg_str("Error converting region to pixel coordinates.");
                            fits_free_region(aRgn);
                            return *status;
                        }

                        coords[i as usize] =
                            f64::sqrt(f64::powi(x - coords[0], 2) + f64::powi(y - coords[1], 2));
                    }
                    i += 1;
                }

                /* special case for elliptannulus and boxannulus if only one angle
                was given */

                if (newShape.shape == ShapeType::ElliptAnnulus
                    || newShape.shape == ShapeType::BoxAnnulus)
                    && nParams == 7
                {
                    coords[7] = coords[6];
                }

                /* Also, correct the position angle for any WCS rotation:  */
                /*    If regions are specified in WCS coordintes, then the angles */
                /*    are relative to the WCS system, not the pixel X,Y system */

                if cFmt != CoordFmt::Pixel {
                    match newShape.shape {
                        ShapeType::Sector | ShapeType::Panda => {
                            coords[2] += wcs.rot;
                            coords[3] += wcs.rot;
                            break;
                        }

                        ShapeType::Box
                        | ShapeType::Rectangle
                        | ShapeType::Diamond
                        | ShapeType::Ellipse => {
                            coords[4] += wcs.rot;
                            break;
                        }
                        ShapeType::BoxAnnulus | ShapeType::ElliptAnnulus => {
                            coords[6] += wcs.rot;
                            coords[7] += wcs.rot;
                            break;
                        }

                        ShapeType::EPanda | ShapeType::BPanda => {
                            coords[2] += wcs.rot;
                            coords[3] += wcs.rot;
                            coords[10] += wcs.rot;
                        }
                        _ => {}
                    }
                }

                /* do some precalculations to speed up tests */

                fits_setup_shape(newShape);
            } /* End of while( currLine[currLoc] ) */
            /*
              if (coords != 0)printf("%.8f %.8f %.8f %.8f %.8f\n",
               coords[0],coords[1],coords[2],coords[3],coords[4]);
            */
        } /* End of if...else parse line */
    } /* End of while( fgets(rgnFile) ) */

    /* set up component numbers */

    fits_set_region_components(&mut aRgn);

    // error:

    if *status != 0 {
        fits_free_region(aRgn);
    } else {
        *Rgn = Some(aRgn);
    }

    *status
}

/*---------------------------------------------------------------------------*/
/// Test if the given point is within the region described by Rgn.  X and
/// Y are in pixel coordinates.                                          
pub(crate) fn fits_in_region(X: f64, Y: f64, Rgn: &mut SAORegion) -> c_int {
    let mut x: f64;
    let mut y: f64;
    let mut dx: f64;
    let mut dy: f64;
    let mut xprime: f64;
    let mut yprime: f64;
    let mut r: f64;
    let th: f64;
    let mut Shapes: &mut [RgnShape];
    let mut i: c_int;
    let mut cur_comp: c_int;
    let mut result: bool = false;
    let mut comp_result: bool = false;

    cur_comp = Rgn.Shapes[0].comp;

    for (i, Shapes) in Rgn.Shapes.iter_mut().enumerate() {
        /* if this region has a different component number to the last one  */
        /*	then replace the accumulated selection logical with the union of */
        /*	the current logical and the total logical. Reinitialize the      */
        /* temporary logical.                                               */

        if i == 0 || Shapes.comp != cur_comp {
            result = result || comp_result;
            cur_comp = Shapes.comp;
            /* if an excluded region is given first, then implicitly   */
            /* assume a previous shape that includes the entire image. */
            comp_result = Shapes.sign == 0;
        }

        /* only need to test if  */
        /*   the point is not already included and this is an include region, */
        /* or the point is included and this is an excluded region */

        if (!comp_result && Shapes.sign != 0) || (comp_result && Shapes.sign == 0) {
            comp_result = true;

            match Shapes.shape {
                ShapeType::Box => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to region's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    dx = 0.5 * Shapes.genericParams.p[2];
                    dy = 0.5 * Shapes.genericParams.p[3];
                    if (x < -dx) || (x > dx) || (y < -dy) || (y > dy) {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::BoxAnnulus => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to region's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    dx = 0.5 * Shapes.genericParams.p[4];
                    dy = 0.5 * Shapes.genericParams.p[5];
                    if (x < -dx) || (x > dx) || (y < -dy) || (y > dy) {
                        comp_result = false;
                    } else {
                        /* Repeat test for inner box */
                        x = xprime * Shapes.genericParams.b + yprime * Shapes.genericParams.a;
                        y = -xprime * Shapes.genericParams.a + yprime * Shapes.genericParams.b;

                        dx = 0.5 * Shapes.genericParams.p[2];
                        dy = 0.5 * Shapes.genericParams.p[3];
                        if (x >= -dx) && (x <= dx) && (y >= -dy) && (y <= dy) {
                            comp_result = false;
                        }
                    }
                    break;
                }
                ShapeType::Rectangle => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[5];
                    yprime = Y - Shapes.genericParams.p[6];

                    /*  Rotate point to region's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    dx = Shapes.genericParams.a;
                    dy = Shapes.genericParams.b;
                    if (x < -dx) || (x > dx) || (y < -dy) || (y > dy) {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::Diamond => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to region's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    dx = 0.5 * Shapes.genericParams.p[2];
                    dy = 0.5 * Shapes.genericParams.p[3];
                    r = f64::abs(x / dx) + f64::abs(y / dy);
                    if r > 1.0 {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::Circle => {
                    /*  Shift origin to center of region  */
                    x = X - Shapes.genericParams.p[0];
                    y = Y - Shapes.genericParams.p[1];

                    r = x * x + y * y;
                    if r > Shapes.genericParams.a {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::Annulus => {
                    /*  Shift origin to center of region  */
                    x = X - Shapes.genericParams.p[0];
                    y = Y - Shapes.genericParams.p[1];

                    r = x * x + y * y;
                    if r < Shapes.genericParams.a || r > Shapes.genericParams.b {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::Sector => {
                    /*  Shift origin to center of region  */
                    x = X - Shapes.genericParams.p[0];
                    y = Y - Shapes.genericParams.p[1];

                    if x != 0.0 || y != 0.0 {
                        r = f64::atan2(y, x) * RAD_TO_DEG;
                        if Shapes.genericParams.p[2] <= Shapes.genericParams.p[3] {
                            if r < Shapes.genericParams.p[2] || r > Shapes.genericParams.p[3] {
                                comp_result = false;
                            }
                        } else if r < Shapes.genericParams.p[2] && r > Shapes.genericParams.p[3] {
                            comp_result = false;
                        }
                    }
                    break;
                }
                ShapeType::Ellipse => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to region's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    x /= Shapes.genericParams.p[2];
                    y /= Shapes.genericParams.p[3];
                    r = x * x + y * y;
                    if r > 1.0 {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::ElliptAnnulus => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to outer ellipse's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    x /= Shapes.genericParams.p[4];
                    y /= Shapes.genericParams.p[5];
                    r = x * x + y * y;
                    if r > 1.0 {
                        comp_result = false;
                    } else {
                        /*  Repeat test for inner ellipse  */
                        x = xprime * Shapes.genericParams.b + yprime * Shapes.genericParams.a;
                        y = -xprime * Shapes.genericParams.a + yprime * Shapes.genericParams.b;

                        x /= Shapes.genericParams.p[2];
                        y /= Shapes.genericParams.p[3];
                        r = x * x + y * y;
                        if r < 1.0 {
                            comp_result = false;
                        }
                    }
                    break;
                }
                ShapeType::Line => {
                    /*  Shift origin to first point of line  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to line's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    if !(-0.5..0.5).contains(&y) || (x < -0.5) || (x >= Shapes.genericParams.a) {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::Point => {
                    /*  Shift origin to center of region  */
                    x = X - Shapes.genericParams.p[0];
                    y = Y - Shapes.genericParams.p[1];

                    if !(-0.5..0.5).contains(&x) || !(-0.5..0.5).contains(&y) {
                        comp_result = false;
                    }
                    break;
                }
                ShapeType::Poly => {
                    if X < Shapes.xmin || X > Shapes.xmax || Y < Shapes.ymin || Y > Shapes.ymax {
                        comp_result = false;
                    } else {
                        comp_result =
                            Pt_in_Poly(X, Y, Shapes.polyParams.nPts, &mut Shapes.polyParams.Pts)
                                != 0;
                    }
                    break;
                }
                ShapeType::Panda => {
                    /*  Shift origin to center of region  */
                    x = X - Shapes.genericParams.p[0];
                    y = Y - Shapes.genericParams.p[1];

                    r = x * x + y * y;
                    if r < Shapes.genericParams.a || r > Shapes.genericParams.b {
                        comp_result = false;
                    } else if x != 0.0 || y != 0.0 {
                        th = f64::atan2(y, x) * RAD_TO_DEG;
                        if Shapes.genericParams.p[2] <= Shapes.genericParams.p[3] {
                            if th < Shapes.genericParams.p[2] || th > Shapes.genericParams.p[3] {
                                comp_result = false;
                            }
                        } else if th < Shapes.genericParams.p[2] && th > Shapes.genericParams.p[3] {
                            comp_result = false;
                        }
                    }
                    break;
                }
                ShapeType::EPanda => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to region's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;
                    xprime = x;
                    yprime = y;

                    /* outer region test */
                    x = xprime / Shapes.genericParams.p[7];
                    y = yprime / Shapes.genericParams.p[8];
                    r = x * x + y * y;
                    if r > 1.0 {
                        comp_result = false;
                    } else {
                        /* inner region test */
                        x = xprime / Shapes.genericParams.p[5];
                        y = yprime / Shapes.genericParams.p[6];
                        r = x * x + y * y;
                        if r < 1.0 {
                            comp_result = false;
                        } else {
                            /* angle test */
                            if xprime != 0.0 || yprime != 0.0 {
                                th = f64::atan2(yprime, xprime) * RAD_TO_DEG;
                                if Shapes.genericParams.p[2] <= Shapes.genericParams.p[3] {
                                    if th < Shapes.genericParams.p[2]
                                        || th > Shapes.genericParams.p[3]
                                    {
                                        comp_result = false;
                                    }
                                } else if th < Shapes.genericParams.p[2]
                                    && th > Shapes.genericParams.p[3]
                                {
                                    comp_result = false;
                                }
                            }
                        }
                    }
                    break;
                }
                ShapeType::BPanda => {
                    /*  Shift origin to center of region  */
                    xprime = X - Shapes.genericParams.p[0];
                    yprime = Y - Shapes.genericParams.p[1];

                    /*  Rotate point to region's orientation  */
                    x = xprime * Shapes.genericParams.cosT + yprime * Shapes.genericParams.sinT;
                    y = -xprime * Shapes.genericParams.sinT + yprime * Shapes.genericParams.cosT;

                    /* outer box test */
                    dx = 0.5 * Shapes.genericParams.p[7];
                    dy = 0.5 * Shapes.genericParams.p[8];
                    if (x < -dx) || (x > dx) || (y < -dy) || (y > dy) {
                        comp_result = false;
                    } else {
                        /* inner box test */
                        dx = 0.5 * Shapes.genericParams.p[5];
                        dy = 0.5 * Shapes.genericParams.p[6];
                        if (x >= -dx) && (x <= dx) && (y >= -dy) && (y <= dy) {
                            comp_result = false;
                        } else {
                            /* angle test */
                            if x != 0.0 || y != 0.0 {
                                th = f64::atan2(y, x) * RAD_TO_DEG;
                                if Shapes.genericParams.p[2] <= Shapes.genericParams.p[3] {
                                    if th < Shapes.genericParams.p[2]
                                        || th > Shapes.genericParams.p[3]
                                    {
                                        comp_result = false;
                                    }
                                } else if th < Shapes.genericParams.p[2]
                                    && th > Shapes.genericParams.p[3]
                                {
                                    comp_result = false;
                                }
                            }
                        }
                    }
                    break;
                }
            }

            if Shapes.sign == 0 {
                comp_result = !comp_result;
            }
        }
    }

    result = result || comp_result;

    if result { 1 } else { 0 }
}

/*---------------------------------------------------------------------------*/
///  Free up memory allocated to hold the region data.
///  This is more complicated for the case of polygons, which may be sharing
///  points arrays due to shallow copying (in fits_set_region_components) of
///  'exluded' regions.  We must ensure that these arrays are only freed once.
pub(crate) fn fits_free_region(Rgn: Box<SAORegion>) {
    let mut i: c_int;
    let mut j: c_int;

    let nFreedPoly: c_int = 0;
    let nPolyArraySize: c_int = 10;
    let mut freedPolyPtrs: *mut *mut f64;
    let mut ptsToFree: *mut f64;
    let isAlreadyFreed: bool = false;

    // TODO freedPolyPtrs = (double **)malloc(nPolyArraySize * sizeof(double *));
    /*
    i = 0;
    while i < Rgn.nShapes {
        if (Rgn.Shapes[i].shape == ShapeType::PolyRgn) {
            /* No shared arrays for 'include' polygons */
            if (Rgn.Shapes[i].sign) {
                free(Rgn.Shapes[i].polyParams.Pts);
            } else {
                ptsToFree = Rgn.Shapes[i].polyParams.Pts;
                isAlreadyFreed = false;
                j = 0;
                while j < nFreedPoly && isAlreadyFreed == 0 {
                    if (freedPolyPtrs[j] == ptsToFree) {
                        isAlreadyFreed = true;
                        j += 1;
                    }
                }
                if (!isAlreadyFreed) {
                    free(ptsToFree);
                    /* Now add pointer to array of freed points */
                    if (nFreedPoly == nPolyArraySize) {
                        nPolyArraySize *= 2;
                        // TODO freedPolyPtrs = (double **)realloc(freedPolyPtrs, nPolyArraySize * sizeof(double *));
                    }
                    freedPolyPtrs[nFreedPoly] = ptsToFree;
                    nFreedPoly += 1;
                }
            }
        }
        i += 1;
    }
    */

    drop(Rgn);
}

/*---------------------------------------------------------------------------*/
/// Internal routine for testing whether the coordinate x,y is within the
/// polygon region traced out by the array Pts.
fn Pt_in_Poly(x: f64, y: f64, nPts: c_int, Pts: &mut [f64]) -> c_int {
    let mut j: usize;
    let mut flag: c_int = 0;
    let mut prevX: f64;
    let mut prevY: f64;
    let mut nextX: f64;
    let mut nextY: f64;
    let mut dx: f64;
    let mut dy: f64;
    let mut Dy: f64;

    let nPts = nPts as usize;

    nextX = Pts[nPts - 2];
    nextY = Pts[nPts - 1];

    for i in (0..nPts).step_by(2) {
        prevX = nextX;
        prevY = nextY;

        nextX = Pts[i];
        nextY = Pts[i + 1];

        if (y > prevY && y >= nextY) || (y < prevY && y <= nextY) || (x > prevX && x >= nextX) {
            continue;
        }

        /* Check to see if x,y lies right on the segment */

        if x >= prevX || x > nextX {
            dy = y - prevY;
            Dy = nextY - prevY;

            if (Dy.abs()) < 1e-10 {
                if (dy.abs()) < 1e-10 {
                    return 1;
                } else {
                    continue;
                }
            }

            dx = prevX + ((nextX - prevX) / (Dy)) * dy - x;
            if dx < -1e-10 {
                continue;
            }
            if dx < 1e-10 {
                return 1;
            }
        }

        /* There is an intersection! Make sure it isn't a V point.  */

        if y != prevY {
            flag = 1 - flag;
        } else {
            j = i + 1; /* Point to Y component */
            loop {
                if j > 1 {
                    j -= 2;
                } else {
                    j = nPts - 1;
                }

                if y != Pts[j] {
                    break;
                }
            }

            if (nextY - y) * (y - Pts[j]) > 0.0 {
                flag = 1 - flag;
            }
        }
    }
    flag
}

/*---------------------------------------------------------------------------*/
/// Internal routine to turn a collection of regions read from an ascii file into
/// the more complex structure that is allowed by the FITS REGION extension with
/// multiple components. Regions are anded within components and ored between them
/// ie for a pixel to be selected it must be selected by at least one component
/// and to be selected by a component it must be selected by all that component's
/// shapes///
/// The algorithm is to replicate every exclude region after every include
/// region before it in the list. eg reg1, reg2, -reg3, reg4, -reg5 becomes
/// (reg1, -reg3, -reg5), (reg2, -reg5, -reg3), (reg4, -reg5) where the
/// parentheses designate components.
pub(crate) fn fits_set_region_components(aRgn: &mut SAORegion) {
    let mut i: usize;
    let mut j: isize;
    let mut k: isize;
    let mut icomp: c_int;

    /* loop round shapes */

    i = 0;
    while i < aRgn.nShapes as usize {
        /* first do the case of an exclude region */

        if aRgn.Shapes[i].sign == 0 {
            /* we need to run back through the list copying the current shape as
            required. start by findin the first include shape before this exclude */

            j = (i - 1) as isize;
            while j > 0 && aRgn.Shapes[j as usize].sign == 0 {
                j -= 1;
            }

            /* then go back one more shape */

            j -= 1;

            /* and loop back through the regions */

            while j >= 0 {
                /* if this is an include region then insert a copy of the exclude
                region immediately after it */

                /* Note that this makes shallow copies of a polygon's dynamically
                allocated Pts array -- the memory is shared.  This must be checked
                when freeing in fits_free_region. */

                if aRgn.Shapes[j as usize].sign != 0 {
                    // TODO aRgn.Shapes = (RgnShape *)realloc(aRgn.Shapes, (1 + aRgn.nShapes) * sizeof(RgnShape));
                    aRgn.nShapes += 1;

                    k = (aRgn.nShapes - 1) as isize;
                    while k > j + 1 {
                        aRgn.Shapes[k as usize] = aRgn.Shapes[(k - 1) as usize].clone();
                        k -= 1;
                    }

                    i += 1;
                    aRgn.Shapes[(j + 1) as usize] = aRgn.Shapes[i].clone();
                }

                j -= 1;
            }
        }

        i += 1;
    }

    /* now set the component numbers */

    icomp = 0;
    i = 0;
    while i < aRgn.Shapes.len() {
        if aRgn.Shapes[i].sign != 0 {
            icomp += 1;
        }
        aRgn.Shapes[i].comp = icomp;

        /*
        println!("i = {}, shape = {}, sign = {}, comp = {}", i, aRgn.Shapes[i].shape, aRgn.Shapes[i].sign, aRgn.Shapes[i].comp);
        */

        i += 1;
    }
}

/*---------------------------------------------------------------------------*/
pub(crate) fn fits_setup_shape(newShape: &mut RgnShape) {
    /* Perform some useful calculations now to speed up filter later             */

    let X: f64;
    let Y: f64;
    let mut R: f64;

    let mut coords: &mut [f64];
    let mut i: c_int;

    match newShape.shape {
        ShapeType::Circle => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.a = coords[2] * coords[2];
        }
        ShapeType::Annulus => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.a = coords[2] * coords[2];
            newShape.genericParams.b = coords[3] * coords[3];
        }
        ShapeType::Sector => {
            coords = &mut newShape.genericParams.p;
            while coords[2] > 180.0 {
                coords[2] -= 360.0;
            }
            while coords[2] <= -180.0 {
                coords[2] += 360.0;
            }
            while coords[3] > 180.0 {
                coords[3] -= 360.0;
            }
            while coords[3] <= -180.0 {
                coords[3] += 360.0;
            }
        }
        ShapeType::Ellipse => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.sinT = f64::sin(MY_PI * (coords[4] / 180.0));
            newShape.genericParams.cosT = f64::cos(MY_PI * (coords[4] / 180.0));
        }
        ShapeType::ElliptAnnulus => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.a = f64::sin(MY_PI * (coords[6] / 180.0));
            newShape.genericParams.b = f64::cos(MY_PI * (coords[6] / 180.0));
            newShape.genericParams.sinT = f64::sin(MY_PI * (coords[7] / 180.0));
            newShape.genericParams.cosT = f64::cos(MY_PI * (coords[7] / 180.0));
        }
        ShapeType::Box => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.sinT = f64::sin(MY_PI * (coords[4] / 180.0));
            newShape.genericParams.cosT = f64::cos(MY_PI * (coords[4] / 180.0));
        }
        ShapeType::BoxAnnulus => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.a = f64::sin(MY_PI * (coords[6] / 180.0));
            newShape.genericParams.b = f64::cos(MY_PI * (coords[6] / 180.0));
            newShape.genericParams.sinT = f64::sin(MY_PI * (coords[7] / 180.0));
            newShape.genericParams.cosT = f64::cos(MY_PI * (coords[7] / 180.0));
        }
        ShapeType::Rectangle => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.sinT = f64::sin(MY_PI * (coords[4] / 180.0));
            newShape.genericParams.cosT = f64::cos(MY_PI * (coords[4] / 180.0));
            X = 0.5 * (coords[2] - coords[0]);
            Y = 0.5 * (coords[3] - coords[1]);
            newShape.genericParams.a =
                f64::abs(X * newShape.genericParams.cosT + Y * newShape.genericParams.sinT);
            newShape.genericParams.b =
                f64::abs(Y * newShape.genericParams.cosT - X * newShape.genericParams.sinT);
            coords[5] = 0.5 * (coords[2] + coords[0]);
            coords[6] = 0.5 * (coords[3] + coords[1]);
        }
        ShapeType::Diamond => {
            coords = &mut newShape.genericParams.p;
            newShape.genericParams.sinT = f64::sin(MY_PI * (coords[4] / 180.0));
            newShape.genericParams.cosT = f64::cos(MY_PI * (coords[4] / 180.0));
        }
        ShapeType::Line => {
            coords = &mut newShape.genericParams.p;
            X = coords[2] - coords[0];
            Y = coords[3] - coords[1];
            R = f64::sqrt(X * X + Y * Y);
            newShape.genericParams.sinT = if R != 0.0 { Y / R } else { 0.0 };
            newShape.genericParams.cosT = if R != 0.0 { X / R } else { 1.0 };
            newShape.genericParams.a = R + 0.5;
        }
        ShapeType::Panda => {
            coords = &mut newShape.genericParams.p;
            while coords[2] > 180.0 {
                coords[2] -= 360.0;
            }
            while coords[2] <= -180.0 {
                coords[2] += 360.0;
            }
            while coords[3] > 180.0 {
                coords[3] -= 360.0;
            }
            while coords[3] <= -180.0 {
                coords[3] += 360.0;
            }
            newShape.genericParams.a = newShape.genericParams.p[5] * newShape.genericParams.p[5];
            newShape.genericParams.b = newShape.genericParams.p[6] * newShape.genericParams.p[6];
        }
        ShapeType::EPanda | ShapeType::BPanda => {
            coords = &mut newShape.genericParams.p;
            while coords[2] > 180.0 {
                coords[2] -= 360.0;
            }
            while coords[2] <= -180.0 {
                coords[2] += 360.0;
            }
            while coords[3] > 180.0 {
                coords[3] -= 360.0;
            }
            while coords[3] <= -180.0 {
                coords[3] += 360.0;
            }
            newShape.genericParams.sinT = f64::sin(MY_PI * (coords[10] / 180.0));
            newShape.genericParams.cosT = f64::cos(MY_PI * (coords[10] / 180.0));
        }
        _ => {}
    }

    /*  Set the xmin, xmax, ymin, ymax elements of the RgnShape structure */

    /* For everything which has first two parameters as center position just */
    /* find a circle that encompasses the region and use it to set the       */
    /* bounding box                                                          */

    R = -1.0;

    if newShape.shape == ShapeType::Poly {
        coords = &mut newShape.polyParams.Pts;
    } else {
        coords = &mut newShape.genericParams.p;
    }

    match newShape.shape {
        ShapeType::Circle => {
            R = coords[2];
        }
        ShapeType::Annulus => {
            R = coords[3];
        }
        ShapeType::Ellipse => {
            if coords[2] > coords[3] {
                R = coords[2];
            } else {
                R = coords[3];
            }
        }
        ShapeType::ElliptAnnulus => {
            if coords[4] > coords[5] {
                R = coords[4];
            } else {
                R = coords[5];
            }
        }
        ShapeType::Box => {
            R = f64::sqrt(coords[2] * coords[2] + coords[3] * coords[3]) / 2.0;
        }
        ShapeType::BoxAnnulus => {
            R = f64::sqrt(coords[4] * coords[5] + coords[4] * coords[5]) / 2.0;
        }
        ShapeType::Diamond => {
            if coords[2] > coords[3] {
                R = coords[2] / 2.0;
            } else {
                R = coords[3] / 2.0;
            }
        }
        ShapeType::Point => {
            R = 1.0;
        }
        ShapeType::Panda => {
            R = coords[6];
        }
        ShapeType::EPanda => {
            if coords[7] > coords[8] {
                R = coords[7];
            } else {
                R = coords[8];
            }
        }
        ShapeType::BPanda => {
            R = f64::sqrt(coords[7] * coords[8] + coords[7] * coords[8]) / 2.0;
        }
        _ => {}
    }

    if R > 0.0 {
        newShape.xmin = coords[0] - R;
        newShape.xmax = coords[0] + R;
        newShape.ymin = coords[1] - R;
        newShape.ymax = coords[1] + R;

        return;
    }

    /* Now do the rest of the shapes that require individual methods */

    match newShape.shape {
        ShapeType::Rectangle => {
            R = f64::sqrt(
                (coords[5] - coords[0]) * (coords[5] - coords[0])
                    + (coords[6] - coords[1]) * (coords[6] - coords[1]),
            );
            newShape.xmin = coords[5] - R;
            newShape.xmax = coords[5] + R;
            newShape.ymin = coords[6] - R;
            newShape.ymax = coords[6] + R;
        }
        ShapeType::Poly => {
            newShape.xmin = coords[0];
            newShape.xmax = coords[0];
            newShape.ymin = coords[1];
            newShape.ymax = coords[1];
            i = 2;
            while i < newShape.polyParams.nPts {
                if newShape.xmin > coords[i as usize] {
                    /* Min X */
                    newShape.xmin = coords[i as usize];
                }
                if newShape.xmax < coords[i as usize] {
                    /* Max X */
                    newShape.xmax = coords[i as usize];
                }
                i += 1;
                if newShape.ymin > coords[i as usize] {
                    /* Min Y */
                    newShape.ymin = coords[i as usize];
                }
                if newShape.ymax < coords[i as usize] {
                    /* Max Y */
                    newShape.ymax = coords[i as usize];
                }
                i += 1;
            }
        }
        ShapeType::Line => {
            if coords[0] > coords[2] {
                newShape.xmin = coords[2];
                newShape.xmax = coords[0];
            } else {
                newShape.xmin = coords[0];
                newShape.xmax = coords[2];
            }
            if coords[1] > coords[3] {
                newShape.ymin = coords[3];
                newShape.ymax = coords[1];
            } else {
                newShape.ymin = coords[1];
                newShape.ymax = coords[3];
            }
            /* sector doesn't have min and max so indicate by setting max < min */
        }
        ShapeType::Sector => {
            newShape.xmin = 1.0;
            newShape.xmax = -1.0;
            newShape.ymin = 1.0;
            newShape.ymax = -1.0;
        }
        _ => {}
    }
}

/*---------------------------------------------------------------------------*/
/// Read regions from a FITS region extension and return the information
/// in the "SAORegion" structure.  If it is nonNULL, use wcs to convert the
/// region coordinates to pixels.  Return an error if region is in degrees
/// but no WCS data is provided.
pub(crate) unsafe fn fits_read_fits_region(
    mut fptr: Box<fitsfile>,
    wcs: &WCSdata,
    Rgn: &mut Option<Box<SAORegion>>,
    status: &mut c_int,
) -> c_int {
    unsafe {
        let mut i: c_int;
        let mut j: c_int;
        let mut icol: [c_int; 6] = [0; 6];
        let mut idum: c_int = 0;
        let mut anynul: c_int = 0;
        let mut npos: c_int;
        let mut dotransform: bool = false;
        let mut got_component: bool = true;
        let mut tstatus: c_int = 0;
        let mut icsize: [c_long; 6] = [0; 6];

        let mut X: f64;
        let mut Y: f64;
        let mut Theta: f64;
        let mut Xsave: f64 = 0.0;
        let mut Ysave: f64 = 0.0;
        let mut Xpos: f64 = 0.0;
        let mut Ypos: f64 = 0.0;
        let mut coords: &mut [f64];

        let mut comment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let colname: [&[u8; FLEN_VALUE]; 6] = [
            b"X                                                                     \0",
            b"Y                                                                     \0",
            b"SHAPE                                                                 \0",
            b"R                                                                     \0",
            b"ROTANG                                                                \0",
            b"COMPONENT                                                             \0",
        ];
        let shapename: [&[u8; FLEN_VALUE]; 17] = [
            b"POINT                                                                 \0",
            b"CIRCLE                                                                \0",
            b"ELLIPSE                                                               \0",
            b"ANNULUS                                                               \0",
            b"ELLIPTANNULUS                                                         \0",
            b"BOX                                                                   \0",
            b"ROTBOX                                                                \0",
            b"BOXANNULUS                                                            \0",
            b"RECTANGLE                                                             \0",
            b"ROTRECTANGLE                                                          \0",
            b"POLYGON                                                               \0",
            b"PIE                                                                   \0",
            b"SECTOR                                                                \0",
            b"DIAMOND                                                               \0",
            b"RHOMBUS                                                               \0",
            b"ROTDIAMOND                                                            \0",
            b"ROTRHOMBUS                                                            \0",
        ];

        let shapetype: [ShapeType; 17] = [
            ShapeType::Point,
            ShapeType::Circle,
            ShapeType::Ellipse,
            ShapeType::Annulus,
            ShapeType::ElliptAnnulus,
            ShapeType::Box,
            ShapeType::Box,
            ShapeType::BoxAnnulus,
            ShapeType::Rectangle,
            ShapeType::Rectangle,
            ShapeType::Poly,
            ShapeType::Sector,
            ShapeType::Sector,
            ShapeType::Diamond,
            ShapeType::Diamond,
            ShapeType::Diamond,
            ShapeType::Diamond,
        ];

        let mut aRgn: &mut SAORegion;
        let newShape: &mut RgnShape;

        // TODO this should be unneccessary but the compiler is complaining because it doesn't
        // know that the variable `dotransform` is always set to true before it is used
        let mut regwcs: Box<WCSdata> = Box::default();

        if *status != 0 {
            return *status;
        }

        let aRgn = box_try_new(SAORegion::default());

        if aRgn.is_err() {
            ffpmsg_str("Couldn't allocate memory to hold Region file contents.");
            *status = MEMORY_ALLOCATION;
            return *status;
        }

        let mut aRgn = aRgn.unwrap();

        aRgn.nShapes = 0;
        aRgn.Shapes = Vec::new();
        if wcs.exists {
            aRgn.wcs = *wcs;
        } else {
            aRgn.wcs.exists = false;
        }

        /* See if we are already positioned to a region extension, else */
        /* move to the REGION extension (file is already open). */

        tstatus = 0;

        for i in 0..5 {
            ffgcno_safe(
                &mut fptr,
                CASEINSEN as c_int,
                cast_slice(colname[i]),
                &mut icol[i],
                &mut tstatus,
            );
        }

        if tstatus != 0 {
            /* couldn't find the required columns, so search for "REGION" extension */
            if ffmnhd_safe(&mut fptr, BINARY_TBL, cs!(c"REGION"), 1, status) != 0 {
                ffpmsg_str("Could not move to REGION extension.");
                fits_free_region(aRgn);
                return *status;
            }
        }

        /* get the number of shapes and allocate memory */

        if ffgky_safe(
            &mut fptr,
            KeywordDatatypeMut::TINT(&mut aRgn.nShapes),
            cs!(c"NAXIS2"),
            Some(&mut comment),
            status,
        ) != 0
        {
            ffpmsg_str("Could not read NAXIS2 keyword.");
            fits_free_region(aRgn);
            return *status;
        }

        if aRgn
            .Shapes
            .try_reserve_exact(aRgn.nShapes as usize)
            .is_err()
        {
            ffpmsg_str("Failed to allocate memory for Region data");
            *status = MEMORY_ALLOCATION;
            fits_free_region(aRgn);
            return *status;
        } else {
            aRgn.Shapes
                .resize(aRgn.nShapes as usize, RgnShape::default());
        }

        /* get the required column numbers */

        for i in 0..5 {
            if ffgcno_safe(
                &mut fptr,
                CASEINSEN as c_int,
                cast_slice(colname[i]),
                &mut icol[i],
                status,
            ) != 0
            {
                ffpmsg_str("Could not find column.");
                fits_free_region(aRgn);
                return *status;
            }
        }

        /* try to get the optional column numbers */

        if ffgcno_safe(
            &mut fptr,
            CASEINSEN as c_int,
            cast_slice(colname[5]),
            &mut icol[5],
            status,
        ) != 0
        {
            got_component = false;
        }

        /* if there was input WCS then read the WCS info for the region in case they */
        /* are different and we have to transform */

        dotransform = false;
        if aRgn.wcs.exists {
            let tmp_regwcs = box_try_new(WCSdata::default());
            if tmp_regwcs.is_err() {
                ffpmsg_str("Failed to allocate memory for Region WCS data");
                *status = MEMORY_ALLOCATION;
                fits_free_region(aRgn);
                return *status;
            }
            regwcs = tmp_regwcs.unwrap();

            regwcs.exists = true;

            if ffgtcs_safer(
                &mut fptr,
                icol[0],
                icol[1],
                &mut regwcs.xrefval,
                &mut regwcs.yrefval,
                &mut regwcs.xrefpix,
                &mut regwcs.yrefpix,
                &mut regwcs.xinc,
                &mut regwcs.yinc,
                &mut regwcs.rot,
                &mut regwcs.dtype,
                status,
            ) != 0
            {
                regwcs.exists = false;
                *status = 0;
            }

            if regwcs.exists
                && wcs.exists
                && (f64::abs(regwcs.xrefval - wcs.xrefval) > 1.0e-6
                    || f64::abs(regwcs.yrefval - wcs.yrefval) > 1.0e-6
                    || f64::abs(regwcs.xrefpix - wcs.xrefpix) > 1.0e-6
                    || f64::abs(regwcs.yrefpix - wcs.yrefpix) > 1.0e-6
                    || f64::abs(regwcs.xinc - wcs.xinc) > 1.0e-6
                    || f64::abs(regwcs.yinc - wcs.yinc) > 1.0e-6
                    || f64::abs(regwcs.rot - wcs.rot) > 1.0e-6
                    || strcmp_safe(&regwcs.dtype, &wcs.dtype) == 0)
            {
                dotransform = true;
            }
        }

        /* get the sizes of the X, Y, R, and ROTANG vectors */

        for i in 0..6 {
            if ffgtdm_safe(&mut fptr, icol[i], 1, &mut idum, &mut icsize[i..], status) != 0 {
                ffpmsg_str("Could not find vector size of column.");
                fits_free_region(aRgn);
                return *status;
            }
        }

        let mut cvalue: Vec<c_char> = vec![0; FLEN_VALUE + 1];

        /* loop over the shapes - note 1-based counting for rows in FITS files */

        #[allow(clippy::never_loop)]
        // Clippy thinks this never loops because aRgn.nShapes could be 0.
        for i in 1..=(aRgn.nShapes) {
            newShape = &mut aRgn.Shapes[i as usize - 1];
            for j in 0..8 {
                newShape.genericParams.p[j] = 0.0;
            }
            newShape.genericParams.a = 0.0;
            newShape.genericParams.b = 0.0;
            newShape.genericParams.sinT = 0.0;
            newShape.genericParams.cosT = 0.0;

            /* get the shape */

            if ffgcvs_safe(
                &mut fptr,
                icol[2],
                i as LONGLONG,
                1,
                1,
                Some(cs!(c" ")),
                &mut [&mut cvalue],
                Some(&mut anynul),
                status,
            ) != 0
            {
                ffpmsg_str("Could not read shape.");
                fits_free_region(aRgn);
                return *status;
            }

            /* set include or exclude */

            newShape.sign = 1;
            let mut cvalue2 = 0; // Index into cvalue
            if strncmp_safe(&cvalue, cs!(c"!"), 1) == 0 {
                newShape.sign = 0;
                cvalue2 += 1;
            }

            /* set the shape type */
            for j in 0..17 {
                if strcmp_safe(&cvalue[cvalue2..], cast_slice(shapename[j])) == 0 {
                    newShape.shape = shapetype[j].clone();
                }
            }

            /* allocate memory for polygon case and set coords pointer */

            if newShape.shape == ShapeType::Poly {
                if newShape
                    .polyParams
                    .Pts
                    .try_reserve_exact(2 * icsize[0] as usize)
                    .is_err()
                {
                    ffpmsg_str("Could not allocate memory to hold polygon parameters");
                    *status = MEMORY_ALLOCATION;
                    fits_free_region(aRgn);
                    return *status;
                } else {
                    newShape.polyParams.Pts.resize(2 * icsize[0] as usize, 0.0);
                }

                newShape.polyParams.nPts = 2 * icsize[0] as c_int;
            }

            if newShape.shape == ShapeType::Poly {
                coords = &mut newShape.polyParams.Pts;
            } else {
                coords = &mut newShape.genericParams.p;
            }

            /* read X and Y. Polygon and Rectangle require special cases */

            npos = 1;
            if newShape.shape == ShapeType::Poly {
                npos = newShape.polyParams.nPts / 2;
            }
            if newShape.shape == ShapeType::Rectangle {
                npos = 2;
            }

            let mut ci = 0; // coords index
            let mut j = 0;
            while j < npos {
                if ffgcvd_safe(
                    &mut fptr,
                    icol[0],
                    i as LONGLONG,
                    j as LONGLONG + 1,
                    1,
                    DOUBLENULLVALUE,
                    &mut coords[ci..],
                    Some(&mut anynul),
                    status,
                ) != 0
                {
                    ffpmsg_str("Failed to read X column for polygon region");
                    fits_free_region(aRgn);
                    return *status;
                }

                if coords[ci] == DOUBLENULLVALUE {
                    /* check for null value end of array marker */
                    npos = j;
                    newShape.polyParams.nPts = npos * 2;
                    break;
                }
                ci += 1;

                if ffgcvd_safe(
                    &mut fptr,
                    icol[1],
                    i as LONGLONG,
                    j as LONGLONG + 1,
                    1,
                    DOUBLENULLVALUE,
                    &mut coords[ci..],
                    Some(&mut anynul),
                    status,
                ) != 0
                {
                    ffpmsg_str("Failed to read Y column for polygon region");
                    fits_free_region(aRgn);
                    return *status;
                }

                if coords[ci] == DOUBLENULLVALUE {
                    /* check for null value end of array marker */
                    npos = j;
                    newShape.polyParams.nPts = npos * 2;
                    ci -= 1;
                    break;
                }
                ci += 1;

                if j == 0 {
                    /* save the first X and Y coordinate */
                    Xsave = coords[ci - 2];
                    Ysave = coords[ci - 1];
                } else if (Xsave == coords[ci]) && (Ysave == coords[ci - 1]) {
                    /* if point has same coordinate as first point, this marks the end of the array */
                    npos = j + 1;
                    newShape.polyParams.nPts = npos * 2;
                    break;
                }

                j += 1;
            }

            /* transform positions if the region and input wcs differ */

            if dotransform {
                ci -= npos as usize * 2;
                Xsave = coords[ci];
                Ysave = coords[ci + 1];
                for j in 0..npos {
                    ffwldp_safe(
                        coords[ci + 2 * j as usize],
                        coords[ci + (2 * j as usize) + 1],
                        regwcs.xrefval,
                        regwcs.yrefval,
                        regwcs.xrefpix,
                        regwcs.yrefpix,
                        regwcs.xinc,
                        regwcs.yinc,
                        regwcs.rot,
                        &regwcs.dtype,
                        &mut Xpos,
                        &mut Ypos,
                        status,
                    );

                    let mut tmpxpix = 0.0;
                    let mut tmpypix = 0.0;

                    ffxypx_safe(
                        Xpos,
                        Ypos,
                        wcs.xrefval,
                        wcs.yrefval,
                        wcs.xrefpix,
                        wcs.yrefpix,
                        wcs.xinc,
                        wcs.yinc,
                        wcs.rot,
                        &wcs.dtype,
                        &mut tmpxpix,
                        &mut tmpypix,
                        status,
                    );

                    coords[ci + (2 * j as usize)] = tmpxpix;
                    coords[ci + (2 * j as usize) + 1] = tmpypix;

                    if *status != 0 {
                        ffpmsg_str("Failed to transform coordinates");
                        fits_free_region(aRgn);
                        return *status;
                    }
                }
                ci += npos as usize * 2;
            }

            /* read R. Circle requires one number; Box, Diamond, Ellipse, Annulus, Sector
            and Panda two; Boxannulus and Elliptannulus four; Point, Rectangle and
            Polygon none. */

            npos = 0;
            match newShape.shape {
                ShapeType::Circle => {
                    npos = 1;
                }
                ShapeType::Box
                | ShapeType::Diamond
                | ShapeType::Ellipse
                | ShapeType::Annulus
                | ShapeType::Sector => {
                    npos = 2;
                }
                ShapeType::BoxAnnulus | ShapeType::ElliptAnnulus => {
                    npos = 4;
                }
                _ => {}
            }

            if npos > 0 {
                if ffgcvd_safe(
                    &mut fptr,
                    icol[3],
                    i as LONGLONG,
                    1,
                    npos.into(),
                    0.0,
                    &mut coords[ci..],
                    Some(&mut anynul),
                    status,
                ) != 0
                {
                    ffpmsg_str("Failed to read R column for region");
                    fits_free_region(aRgn);
                    return *status;
                }

                /* transform lengths if the region and input wcs differ */

                if dotransform {
                    for j in 0..npos {
                        if newShape.shape == ShapeType::Poly {
                            coords = &mut newShape.polyParams.Pts;
                        } else {
                            coords = &mut newShape.genericParams.p;
                        }
                        Y = Ysave + (coords[ci]);
                        X = Xsave;
                        ffwldp_safe(
                            X,
                            Y,
                            regwcs.xrefval,
                            regwcs.yrefval,
                            regwcs.xrefpix,
                            regwcs.yrefpix,
                            regwcs.xinc,
                            regwcs.yinc,
                            regwcs.rot,
                            &regwcs.dtype,
                            &mut Xpos,
                            &mut Ypos,
                            status,
                        );
                        ffxypx_safe(
                            Xpos,
                            Ypos,
                            wcs.xrefval,
                            wcs.yrefval,
                            wcs.xrefpix,
                            wcs.yrefpix,
                            wcs.xinc,
                            wcs.yinc,
                            wcs.rot,
                            &wcs.dtype,
                            &mut X,
                            &mut Y,
                            status,
                        );
                        if *status != 0 {
                            ffpmsg_str("Failed to transform coordinates");
                            fits_free_region(aRgn);
                            return *status;
                        }

                        if newShape.shape == ShapeType::Poly {
                            newShape.polyParams.Pts[ci] = f64::sqrt(
                                f64::powi(X - newShape.genericParams.p[0], 2)
                                    + f64::powi(Y - newShape.genericParams.p[1], 2),
                            );
                        } else {
                            coords[ci] = f64::sqrt(
                                f64::powi(X - coords[0], 2) + f64::powi(Y - coords[1], 2),
                            );
                        }

                        ci += 1;
                    }
                } else {
                    ci += npos as usize;
                }
            }

            /* read ROTANG. Requires two values for Boxannulus, Elliptannulus, Sector,
            Panda; one for Box, Diamond, Ellipse; and none for Circle, Point, Annulus,
            Rectangle, Polygon */

            npos = 0;
            match newShape.shape {
                ShapeType::Box | ShapeType::Diamond | ShapeType::Ellipse => {
                    npos = 1;
                    break;
                }
                ShapeType::BoxAnnulus | ShapeType::ElliptAnnulus | ShapeType::Sector => {
                    npos = 2;
                    break;
                }
                _ => {
                    break;
                }
            }

            if npos > 0 {
                if ffgcvd_safe(
                    &mut fptr,
                    icol[4],
                    i as LONGLONG,
                    1,
                    npos as LONGLONG,
                    0.0,
                    &mut coords[ci..],
                    Some(&mut anynul),
                    status,
                ) != 0
                {
                    ffpmsg_str("Failed to read ROTANG column for region");
                    fits_free_region(aRgn);
                    return *status;
                }

                /* transform angles if the region and input wcs differ */

                if dotransform {
                    Theta = (wcs.rot) - (regwcs.rot);
                    for j in 0..npos {
                        coords[ci] += Theta;
                        ci += 1;
                    }
                } else {
                    ci += npos as usize;
                }
            }

            /* read the component number */

            if got_component {
                if ffgcv_safe(
                    &mut fptr,
                    TINT,
                    icol[5],
                    i as LONGLONG,
                    1,
                    1,
                    Some(NullValue::Int(0)),
                    cast_slice_mut(&mut [(newShape.comp as usize)]),
                    Some(&mut anynul),
                    status,
                ) != 0
                {
                    ffpmsg_str("Failed to read COMPONENT column for region");
                    fits_free_region(aRgn);
                    return *status;
                }
            } else {
                newShape.comp = 1;
            }

            /* do some precalculations to speed up tests */

            fits_setup_shape(newShape);

            /* end loop over shapes */
        }

        if *status != 0 {
            fits_free_region(aRgn);
        } else {
            *Rgn = Some(aRgn);
        }

        ffclos_safer(fptr, status);

        *status
    }
}
