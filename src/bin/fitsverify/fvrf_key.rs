/* Transpiled from cfitsio/utilities/fvrf_key.c

The C walks the card with `char *' cursors; here every cursor is an index
into the card slice, and the `char **pt' in/out cursor parameters become
`pt: &mut usize'.  */

use rsfitsio::c_types::{c_char, c_int, c_ulong};
use rsfitsio::fitsio::{FLEN_CARD, FLEN_COMMENT, FLEN_KEYWORD, FLEN_VALUE};

use crate::common::*;
use crate::fvrf_misc::*;
use crate::{scat, spf};

/* Ref: Defininition of the Flexible Image Transport System(FITS),
    Sec. 5.1 and 5.2.
*/
pub(crate) fn fits_parse_card(
    out: Out,                           /* output file pointer */
    kpos: c_int,                        /* keyposition starting from 1 */
    card: &mut [c_char],                /* key card */
    kname: &mut [c_char; FLEN_KEYWORD], /* key name */
    ktype: &mut kwdtyp,                 /* key type */
    kvalue: &mut [c_char; FLEN_VALUE],  /* key value */
    kcomm: &mut [c_char],               /* comment */
) -> c_int {
    let mut vind: [c_char; 3] = [0; 3];
    let mut p: usize;
    let mut i: isize;
    let mut temp1: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut stat: c_ulong = 0;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    kname[0] = 0;
    kvalue[0] = 0;
    kcomm[0] = 0;
    *ktype = kwdtyp::UNKNOWN;

    if cstrlen(card) > FLEN_CARD - 1 {
        strncpy_c(&mut temp1, card, 20);
        temp1[21] = 0;
        spf!(errmes; "card ", CS(card), " is > 80.");
        wrterr(out, &errmes, 1);
        return 1;
    }
    card[FLEN_CARD - 1] = 0;

    /* get the kname */
    strncpy_c(kname, card, 8);
    kname[8] = 0;

    /* take out the trailing space */
    i = 7;
    p = 7;
    while i >= 0 && isspace_c(kname[p]) {
        kname[p] = 0;
        if p == 0 {
            i -= 1;
            break;
        }
        p -= 1;
        i -= 1;
    }

    /* Whether the keyword name is left justified */
    i = 0;
    p = 0;
    while isspace_c(kname[p]) && kname[p] != 0 {
        p += 1;
        i += 1;
    }
    if i < 8 && i > 0 {
        spf!(errmes; "Keyword #", kpos, ": Name ", CS(kname), " is not left justified.");
        wrterr(out, &errmes, 1);
    }
    /* Whether the characters in keyword name are valid */
    while kname[p] != 0 {
        let c = kname[p] as u8;
        if !(b'A'..=b'Z').contains(&c) && !(b'0'..=b'9').contains(&c) && c != b'-' && c != b'_' {
            spf!(errmes;
                "Keyword #", kpos, ": Name \"", CS(kname), "\" contains char \"", CHR(kname[p]),
                "\" which is not upper case letter, digit, \"-\", or \"_\".");
            wrterr(out, &errmes, 1);
            break;
        }
        p += 1;
        i += 1;
    }

    /* COMMENT, HISTORY, HIERARCH and "" keywords */
    if ceq(kname, b"COMMENT")
        || ceq(kname, b"HISTORY")
        || ceq(kname, b"HIERARCH")
        || ceq(kname, b"CONTINUE")
        || ceq(kname, b"")
    {
        *ktype = kwdtyp::COM_KEY;

        p = 8;
        strcpy_c(kcomm, &card[p..]);
        kcomm[FLEN_COMMENT - 1] = 0;
        while card[p] != 0 {
            if !isprint_c(card[p]) {
                spf!(errmes;
                    "Keyword #", kpos, ", ", CS(kname), ": String contains non-text characters.");
                wrterr(out, &errmes, 1);
                return 1;
            }
            p += 1;
        }
        p = 0;
        while !isspace_c(kname[p]) && kname[p] != 0 {
            p += 1;
        }
        kname[p] = 0;
        return 0;
    }

    /* End Keyword: 9-80 shall be filled with ASCII blanks \x20 */
    if ceq(kname, b"END") {
        *ktype = kwdtyp::COM_KEY;
        if card[3] == 0 {
            return 0;
        }
        p = 8;
        while card[p] != 0 {
            if card[p] as u8 != 0x20 {
                wrterr_str(out, "END keyword contains non-blank characters.", 1);
                return 1;
            }
            p += 1;
        }
        kname[3] = 0;
        return 0;
    }

    /* check for value indicator */
    p = 8;
    strncpy_c(&mut vind, &card[p..], 2);
    vind[2] = 0;
    if !ceq(&vind, b"= ") && !ceq(&vind, b"=") {
        /* no value indicator, so this is a commentary keyword */
        *ktype = kwdtyp::COM_KEY;
        strcpy_c(kcomm, &card[p..]);
        kcomm[FLEN_COMMENT - 1] = 0;
        while card[p] != 0 {
            if !isprint_c(card[p]) {
                spf!(errmes;
                    "Keyword #", kpos, ", ", CS(kname), ": String contains non-text characters.");
                wrterr(out, &errmes, 1);
                return 1;
            }
            p += 1;
        }
        p = 0;
        while !isspace_c(kname[p]) && kname[p] != 0 {
            p += 1;
        }
        kname[p] = 0;
        return 0;
    }

    p = 10;
    while isspace_c(card[p]) && card[p] != 0 {
        p += 1;
    }
    match card[p] as u8 {
        b'\'' => {
            /* string */
            get_str(card, &mut p, kvalue, &mut stat);
            *ktype = kwdtyp::STR_KEY;
            if card[p] != 0 {
                get_comm(card, &mut p, kcomm, &mut stat);
            }
        }
        b'T' | b'F' => {
            /*logical */
            get_log(card, &mut p, kvalue, &mut stat);
            *ktype = kwdtyp::LOG_KEY;
            if card[p] != 0 {
                get_comm(card, &mut p, kcomm, &mut stat);
            }
        }
        /* number */
        b'+' | b'-' | b'.' | b'0'..=b'9' => {
            get_num(card, &mut p, kvalue, ktype, &mut stat);
            if card[p] != 0 {
                get_comm(card, &mut p, kcomm, &mut stat);
            }
        }
        b'(' => {
            /* complex number */
            get_cmp(card, &mut p, kvalue, ktype, &mut stat);
            if card[p] != 0 {
                get_comm(card, &mut p, kcomm, &mut stat);
            }
        }
        b'/' => {
            /* comment */
            if card[p] != 0 {
                get_comm(card, &mut p, kcomm, &mut stat);
            }
            *ktype = kwdtyp::UNKNOWN;
        }
        _ => {
            get_unknown(card, &mut p, kvalue, ktype, &mut stat);
            if card[p] != 0 {
                get_comm(card, &mut p, kcomm, &mut stat);
            }
        }
    }
    /* take out the trailing blanks for non-string keys */
    if *ktype != kwdtyp::STR_KEY {
        i = cstrlen(kvalue) as isize;
        let mut q = i - 1;
        while q >= 0 && isspace_c(kvalue[q as usize]) && i > 0 {
            kvalue[q as usize] = 0;
            q -= 1;
            i -= 1;
        }
        if i == 0 && q >= 0 && isspace_c(kvalue[q as usize]) {
            kvalue[q as usize] = 0;
        }
    }
    pr_kval_err(out, kpos, kname, kvalue, stat);
    if stat != 0 {
        return 1;
    }
    0
}

/* parse And test the string keys */
pub(crate) fn get_str(
    card: &[c_char],       /* card string from character 11*/
    pt: &mut usize,        /* cursor into card */
    kvalue: &mut [c_char], /* key value string */
    stat: &mut c_ulong,    /* error number */
) {
    let mut pi: usize;
    let mut prev: u8; /* previous char */
    let mut nchar: isize;
    let mut p: usize;

    p = *pt;
    pi = p;
    p += 1;
    prev = b'a';
    while card[p] != 0 {
        if !isprint_c(card[p]) {
            *stat |= BAD_STR;
        }
        if prev == b'\'' && card[p] as u8 != b'\'' {
            break;
        }
        if prev == b'\'' && card[p] as u8 == b'\'' {
            /* skip the '' */
            p += 1;
            prev = b'a';
        } else {
            prev = card[p] as u8;
            p += 1;
        }
    }
    p -= 1;
    if card[p] as u8 != b'\'' {
        *stat |= NO_TRAIL_QUOTE;
    }
    pi += 1;
    nchar = p as isize - pi as isize; /* excluding the ' */
    if nchar < 0 {
        nchar = 0;
    }
    let nchar = std::cmp::min(nchar as usize, kvalue.len() - 1);
    strncpy_c(kvalue, &card[pi..], nchar);
    kvalue[nchar] = 0;
    let mut q = nchar as isize - 1;
    while q >= 0 && isspace_c(kvalue[q as usize]) {
        kvalue[q as usize] = 0;
        q -= 1;
    } /* delete the trailing space */
    p += 1; /* skip the  ' */
    while isspace_c(card[p]) && card[p] != 0 {
        p += 1;
    }
    *pt = p;
}

/* parse and test the logical keys */
pub(crate) fn get_log(
    card: &[c_char],       /* card string */
    pt: &mut usize,        /* cursor into card */
    kvalue: &mut [c_char], /* key value string */
    stat: &mut c_ulong,    /* error number */
) {
    let mut p: usize;

    p = *pt;
    kvalue[0] = card[p];
    kvalue[1] = 0;
    p += 1;
    while isspace_c(card[p]) {
        p += 1;
    }
    if card[p] as u8 != b'/' && card[p] != 0 {
        *stat |= BAD_LOGICAL;
    }
    *pt = p;
}

/* parse and test the numerical keys */
pub(crate) fn get_num(
    card: &[c_char],       /* card string */
    pt: &mut usize,        /* cursor into card */
    kvalue: &mut [c_char], /* comment string */
    ktype: &mut kwdtyp,
    stat: &mut c_ulong, /* error number */
) {
    let pi: usize;
    let mut set_deci = 0;
    let mut set_expo = 0;
    let nchar: usize;
    let mut p: usize;

    p = *pt;
    pi = p;
    *ktype = kwdtyp::INT_KEY;

    if card[p] as u8 != b'+'
        && card[p] as u8 != b'-'
        && !isdigit_c(card[p])
        && card[p] as u8 != b'.'
    {
        *stat |= BAD_NUM;
        return;
    }
    if card[p] as u8 == b'.' {
        *ktype = kwdtyp::FLT_KEY;
        set_deci = 1;
    }

    p += 1;
    while !isspace_c(card[p]) && card[p] != 0 && card[p] as u8 != b'/' {
        if card[p] as u8 == b'.' && set_deci == 0 {
            set_deci = 1;
            *ktype = kwdtyp::FLT_KEY;
            p += 1;
            continue;
        }
        if (card[p] as u8 == b'd' || card[p] as u8 == b'e') && set_expo == 0 {
            set_expo = 1;
            *ktype = kwdtyp::FLT_KEY;
            p += 1;
            if card[p] as u8 == b'+' || card[p] as u8 == b'-' {
                p += 1;
            }
            *stat |= LOWCASE_EXPO;
            continue;
        }
        if (card[p] as u8 == b'D' || card[p] as u8 == b'E') && set_expo == 0 {
            set_expo = 1;
            *ktype = kwdtyp::FLT_KEY;
            p += 1;
            if card[p] as u8 == b'+' || card[p] as u8 == b'-' {
                p += 1;
            }
            continue;
        }
        if !isdigit_c(card[p]) {
            *stat |= BAD_NUM;
        }
        p += 1;
    }
    nchar = std::cmp::min(p - pi, kvalue.len() - 1);
    strncpy_c(kvalue, &card[pi..], nchar);
    kvalue[nchar] = 0;
    while isspace_c(card[p]) && card[p] != 0 {
        p += 1;
    }
    *pt = p;
}

/* parse and test the complex keys */
pub(crate) fn get_cmp(
    incard: &[c_char],     /* card string */
    pt: &mut usize,        /* cursor into card */
    kvalue: &mut [c_char], /* comment string */
    ktype: &mut kwdtyp,
    stat: &mut c_ulong, /* error number */
) {
    let mut p: usize;
    let mut pr_beg: usize;
    /* pr_end / pi_beg are left unset by the C when there is no ',' in the
    value; it then dereferences a NULL pr_end.  Here they stay None and the
    terminating writes are skipped. */
    let mut pr_end: Option<usize> = None;
    let mut pi_beg: Option<usize> = None;
    let mut pi_end: usize = 0;
    let nchar: usize;
    let mut set_comm = 0;
    let mut set_paren = 0;

    let mut tr: c_ulong = 0;
    let mut ti: c_ulong = 0;
    let mut rtype = kwdtyp::UNKNOWN;
    let mut itype = kwdtyp::UNKNOWN;
    let mut temp: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    strcpy_c(&mut card, &incard[*pt..]); /* save the original */
    card[FLEN_CARD - 1] = 0;

    *ktype = kwdtyp::CMI_KEY; /* default: integer complex */
    p = 1;
    pr_beg = p;

    temp[0] = 0;
    while card[p] != 0 && card[p] as u8 != b'/' {
        if card[p] as u8 == b')' {
            set_paren = 1;
            pi_end = p;
            p += 1;
            break;
        }
        if set_comm == 0 && card[p] as u8 == b',' {
            set_comm = 1;
            pr_end = Some(p);
            pi_beg = Some(p + 1);
        } else if card[p] as u8 == b',' {
            *stat |= TOO_MANY_COMMA;
        }
        p += 1;
    }
    if set_comm == 0 {
        *stat |= NO_COMMA;
    }
    if set_paren == 0 {
        *stat |= NO_TRAIL_PAREN;
        pi_end = p;
        if pi_end > 0 {
            pi_end -= 1;
            while pi_end > 0 && isspace_c(card[pi_end]) {
                pi_end -= 1;
            }
            pi_end += 1;
        }
    }

    nchar = std::cmp::min(pi_end, kvalue.len() - 1);
    strncpy_c(kvalue, &card, nchar);
    kvalue[nchar] = 0;
    while isspace_c(card[p]) && card[p] != 0 {
        p += 1;
    }
    *pt += p;

    /* analyse the real and imagine part */
    if let Some(e) = pr_end {
        card[e] = 0;
    }
    card[pi_end] = 0;
    while isspace_c(card[pr_beg]) && card[pr_beg] != 0 {
        pr_beg += 1;
    }
    if let Some(mut b) = pi_beg {
        while isspace_c(card[b]) && card[b] != 0 {
            b += 1;
        }
        pi_beg = Some(b);
    }
    temp[0] = 0;
    let mut pp = pr_beg;
    get_num(&card, &mut pp, &mut temp, &mut rtype, &mut tr);
    if tr != 0 {
        *stat |= BAD_REAL;
    }
    temp[0] = 0;
    if let Some(b) = pi_beg {
        let mut pp = b;
        get_num(&card, &mut pp, &mut temp, &mut itype, &mut ti);
    }
    if ti != 0 {
        *stat |= BAD_IMG;
    }
    if rtype == kwdtyp::FLT_KEY || itype == kwdtyp::FLT_KEY {
        *ktype = kwdtyp::CMF_KEY;
    }
}

/* parse and test the comment keys */
fn get_comm(
    card: &[c_char],      /* card string */
    pt: &mut usize,       /* cursor into card */
    kcomm: &mut [c_char], /* comment string */
    stat: &mut c_ulong,   /* error number */
) {
    let pi: usize;
    let nchar: usize;
    let mut p: usize;

    p = *pt;
    pi = p;
    if card[p] as u8 != b'/' {
        *stat |= NO_START_SLASH;
    }
    p += 1;
    while card[p] != 0 {
        if !isprint_c(card[p]) {
            *stat |= BAD_COMMENT;
        }
        p += 1;
    }
    nchar = std::cmp::min(p - pi, kcomm.len() - 1);
    strncpy_c(kcomm, &card[pi..], nchar);
    kcomm[nchar] = 0;
}

/* parsing the unknown keyword */
pub(crate) fn get_unknown(
    card: &[c_char],       /* card string */
    pt: &mut usize,        /* cursor into card */
    kvalue: &mut [c_char], /* comment string */
    ktype: &mut kwdtyp,
    stat: &mut c_ulong, /* error number */
) {
    let mut p: usize;
    let mut p1: usize;
    let mut temp: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    p = *pt;
    strcpy_c(&mut temp, &card[*pt..]);
    p1 = 0;
    while card[p] != 0 && card[p] as u8 != b'/' {
        p += 1;
        p1 += 1;
    }
    temp[p1] = 0;
    *pt = p;

    strcpy_c(kvalue, &temp);
    *ktype = kwdtyp::UNKNOWN;
    *stat |= UNKNOWN_TYPE;
}

/* routine to print out the error of keyword value/comment */
pub(crate) fn pr_kval_err(
    out: Out,         /* output  FILE */
    kpos: c_int,      /* keyposition starting from 1 */
    kname: &[c_char], /* keyword name */
    kval: &[c_char],  /* keyword value */
    errnum: c_ulong,  /* error number */
) {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    if errnum == 0 {
        return;
    }
    if errnum & BAD_STR != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": String \"", CS(kval),
            "\"  contains non-text characters.");
        wrterr(out, &errmes, 1);
    }
    if errnum & NO_TRAIL_QUOTE != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname),
            ": The closing \"'\" is missing in the string.");
        wrterr(out, &errmes, 1);
    }
    if errnum & BAD_LOGICAL != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": Bad logical value \"", CS(kval), "\".");
        wrterr(out, &errmes, 1);
    }
    if errnum & BAD_NUM != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": Bad numerical value \"", CS(kval), "\".");
        wrterr(out, &errmes, 1);
    }
    if errnum & LOWCASE_EXPO != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname),
            ": lower-case exponent d or e is illegal in value ", CS(kval), ".");
        wrterr(out, &errmes, 1);
    }
    if errnum & NO_TRAIL_PAREN != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": Complex value \"", CS(kval),
            "\" misses closing \")\".");
        wrterr(out, &errmes, 1);
    }
    if errnum & NO_COMMA != 0 {
        spf!(errmes;
            "keyword #", kpos, ", ", CS(kname), " : Complex value \"", CS(kval),
            "\" misses \",\".");
        wrterr(out, &errmes, 1);
    }
    if errnum & TOO_MANY_COMMA != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname),
            ": Too many \",\" are in the complex value \"", CS(kval), "\".");
        wrterr(out, &errmes, 1);
    }
    if errnum & BAD_REAL != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": Real part of complex value \"", CS(kval),
            "\" is  bad.");
        wrterr(out, &errmes, 1);
    }
    if errnum & BAD_IMG != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": Imagine part of complex value \"", CS(kval),
            "\" is bad.");
        wrterr(out, &errmes, 1);
    }
    if errnum & NO_START_SLASH != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname),
            ": Value and Comment not separated by a \"/\".");
        wrterr(out, &errmes, 1);
    }
    if errnum & BAD_COMMENT != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": Comment contains non-text characters.");
        wrterr(out, &errmes, 1);
    }

    /* don't report null keywords as an error */
    if errnum & UNKNOWN_TYPE != 0 && kval[0] != 0 {
        spf!(errmes;
            "Keyword #", kpos, ", ", CS(kname), ": Type of value \"", CS(kval),
            "\" is unknown.");
        wrterr(out, &errmes, 1);
    }
}

pub(crate) fn check_str(pkey: &FitsKey, out: Out) -> c_int {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    if pkey.ktype == kwdtyp::UNKNOWN && pkey.kvalue[0] == 0 {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname),
            " has a null value; expected a string.");
        wrterr(out, &errmes, 1);
        return 0;
    } else if pkey.ktype != kwdtyp::STR_KEY {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname), ": \"", CS(&pkey.kvalue),
            "\" is not a string.");
        wrterr(out, &errmes, 1);
        return 0;
    }
    1
}

pub(crate) fn check_int(pkey: &FitsKey, out: Out) -> c_int {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    if pkey.ktype == kwdtyp::UNKNOWN && pkey.kvalue[0] == 0 {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname),
            " has a null value; expected an integer.");
        wrterr(out, &errmes, 1);
        return 0;
    } else if pkey.ktype != kwdtyp::INT_KEY {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname), ": value = ", CS(&pkey.kvalue),
            " is not an integer.");
        if pkey.ktype == kwdtyp::STR_KEY {
            scat!(errmes; " The value is entered as a string. ");
        }
        wrterr(out, &errmes, 1);
        return 0;
    }
    1
}

pub(crate) fn check_flt(pkey: &FitsKey, out: Out) -> c_int {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    if pkey.ktype == kwdtyp::UNKNOWN && pkey.kvalue[0] == 0 {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname),
            " has a null value; expected a float.");
        wrterr(out, &errmes, 1);
        return 0;
    } else if pkey.ktype != kwdtyp::INT_KEY && pkey.ktype != kwdtyp::FLT_KEY {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname), ": value = ", CS(&pkey.kvalue),
            " is not a floating point number.");
        if pkey.ktype == kwdtyp::STR_KEY {
            scat!(errmes; " The value is entered as a string. ");
        }
        wrterr(out, &errmes, 1);
        return 0;
    }
    1
}

fn check_cmi(pkey: &FitsKey, out: Out) -> c_int {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    if pkey.ktype != kwdtyp::CMI_KEY {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname), ": value = ", CS(&pkey.kvalue),
            " is not a integer complex number.");
        if pkey.ktype == kwdtyp::STR_KEY {
            scat!(errmes; " The value is entered as a string. ");
        }
        wrterr(out, &errmes, 1);
        return 0;
    }
    1
}

fn check_cmf(pkey: &FitsKey, out: Out) -> c_int {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    if pkey.ktype != kwdtyp::CMI_KEY && pkey.ktype != kwdtyp::CMF_KEY {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname), ": value = ", CS(&pkey.kvalue),
            " is not a floating point complex number.");
        if pkey.ktype == kwdtyp::STR_KEY {
            scat!(errmes; " The value is entered as a string. ");
        }
        wrterr(out, &errmes, 1);
        return 0;
    }
    1
}

pub(crate) fn check_log(pkey: &FitsKey, out: Out) -> c_int {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    if pkey.ktype != kwdtyp::LOG_KEY {
        spf!(errmes;
            "Keyword #", pkey.kindex, ", ", CS(&pkey.kname), ": value = ", CS(&pkey.kvalue),
            " is not a logical constant.");
        if pkey.ktype == kwdtyp::STR_KEY {
            scat!(errmes; " The value is entered as a string. ");
        }
        wrterr(out, &errmes, 1);
        return 0;
    }
    1
}

pub(crate) fn check_fixed_int(card: &[c_char], out: Out) -> c_int {
    let mut cptr: usize;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* fixed format integer must be right justified in columns 11-30 */

    cptr = 10;

    while card[cptr] as u8 == b' ' {
        cptr += 1;
    } /* skip leading spaces */

    if card[cptr] as u8 == b'-' {
        cptr += 1; /* skip leading minus sign */
    } else if card[cptr] as u8 == b'+' {
        cptr += 1; /* skip leading plus sign */
    }

    while isdigit_c(card[cptr]) {
        cptr += 1;
    } /* skip digits */

    /* should be pointing to column 31 of the card */

    if cptr != 30 {
        spf!(errmes; CSW(card, 0, Some(8)), " mandatory keyword is not in integer fixed format:");
        wrterr(out, &errmes, 1);
        print_fmt(out, card, 13);
        print_fmt_str(out, "          -------------------^", 13);

        return 0;
    }
    1
}

pub(crate) fn check_fixed_log(card: &[c_char], out: Out) -> c_int {
    let mut cptr: usize;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* fixed format logical must have T or F in column 30 */

    cptr = 10;

    while card[cptr] as u8 == b' ' {
        cptr += 1;
    } /* skip leading spaces */

    if card[cptr] as u8 != b'T' && card[cptr] as u8 != b'F' {
        spf!(errmes;
            CSW(card, 0, Some(8)), " mandatory keyword does not have T or F logical value.");
        wrterr(out, &errmes, 1);
        return 0;
    }

    /* should be pointing to column 31 of the card */

    if cptr != 29 {
        spf!(errmes; CSW(card, 0, Some(8)), " mandatory keyword is not in logical fixed format:");
        wrterr(out, &errmes, 1);
        print_fmt(out, card, 13);
        print_fmt_str(out, "          -------------------^", 13);

        return 0;
    }
    1
}

pub(crate) fn check_fixed_str(card: &[c_char], out: Out) -> c_int {
    let mut cptr: usize;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* fixed format string must have quotes in columns 11 and >= 20 */
    /* This only applys to the XTENSION and TFORMn keywords. */

    cptr = 10;

    if card[cptr] as u8 != b'\'' {
        spf!(errmes;
            CSW(card, 0, Some(8)), " mandatory string keyword does not start in col 11.");
        wrterr(out, &errmes, 1);
        print_fmt(out, card, 13);
        print_fmt_str(out, "          ^--------^", 13);
        return 0;
    }

    cptr += 1;

    while card[cptr] as u8 != b'\'' {
        if card[cptr] == 0 {
            spf!(errmes;
                CSW(card, 0, Some(8)),
                " mandatory string keyword missing closing quote character:");
            wrterr(out, &errmes, 1);
            print_fmt(out, card, 13);
            return 0;
        }
        cptr += 1;
    }

    if cptr < 19 {
        spf!(errmes; CSW(card, 0, Some(8)), " mandatory string keyword ends before column 20.");
        wrterr(out, &errmes, 1);
        print_fmt(out, card, 13);
        print_fmt_str(out, "          ^--------^", 13);

        return 0;
    }

    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fvrf_misc::reset_err_wrn;

    /* An 80-column card image, NUL-terminated as CFITSIO hands it over. */
    fn card(s: &str) -> [c_char; FLEN_CARD] {
        let mut b = [0 as c_char; FLEN_CARD];
        for i in 0..80 {
            b[i] = b' ' as c_char;
        }
        for (i, &c) in s.as_bytes().iter().enumerate() {
            b[i] = c as c_char;
        }
        b[80] = 0;
        b
    }

    fn parse(s: &str) -> (c_int, kwdtyp, String, String) {
        reset_err_wrn();
        let mut c = card(s);
        let mut kname = [0 as c_char; FLEN_KEYWORD];
        let mut ktype = kwdtyp::UNKNOWN;
        let mut kvalue = [0 as c_char; FLEN_VALUE];
        let mut kcomm = [0 as c_char; COMM_LEN];
        /* Out::Null suppresses reporting so only the parse result is under test */
        let r = fits_parse_card(
            Out::Null,
            1,
            &mut c,
            &mut kname,
            &mut ktype,
            &mut kvalue,
            &mut kcomm,
        );
        (
            r,
            ktype,
            String::from_utf8_lossy(cbytes(&kname)).into_owned(),
            String::from_utf8_lossy(cbytes(&kvalue)).into_owned(),
        )
    }

    #[test]
    fn test_parse_card_logical() {
        let (r, t, n, v) = parse("SIMPLE  =                    T / conforms to FITS standard");
        assert_eq!(r, 0);
        assert_eq!(t, kwdtyp::LOG_KEY);
        assert_eq!(n, "SIMPLE");
        assert_eq!(v, "T");
    }

    #[test]
    fn test_parse_card_integer() {
        let (r, t, n, v) = parse("BITPIX  =                   16 / bits per data value");
        assert_eq!(r, 0);
        assert_eq!(t, kwdtyp::INT_KEY);
        assert_eq!(n, "BITPIX");
        assert_eq!(v, "16");

        let (_, t, _, v) = parse("TZERO1  =                -32768");
        assert_eq!(t, kwdtyp::INT_KEY);
        assert_eq!(v, "-32768");
    }

    #[test]
    fn test_parse_card_float() {
        let (r, t, _, v) = parse("BSCALE  =                  1.5");
        assert_eq!(r, 0);
        assert_eq!(t, kwdtyp::FLT_KEY);
        assert_eq!(v, "1.5");

        /* upper-case exponent is legal */
        let (r, t, _, v) = parse("CRVAL1  =              1.25E+02");
        assert_eq!(r, 0);
        assert_eq!(t, kwdtyp::FLT_KEY);
        assert_eq!(v, "1.25E+02");

        /* a lower-case exponent is flagged, so the card parse fails */
        let (r, t, _, _) = parse("CRVAL1  =              1.25e+02");
        assert_eq!(r, 1);
        assert_eq!(t, kwdtyp::FLT_KEY);
    }

    #[test]
    fn test_parse_card_string() {
        /* get_str strips the quotes and the trailing blanks inside them */
        let (r, t, n, v) = parse("EXTNAME = 'FOO     '           / name of this HDU");
        assert_eq!(r, 0);
        assert_eq!(t, kwdtyp::STR_KEY);
        assert_eq!(n, "EXTNAME");
        assert_eq!(v, "FOO");

        /* '' is the escape for an embedded quote */
        let (r, t, _, v) = parse("OBJECT  = 'a''b'");
        assert_eq!(r, 0);
        assert_eq!(t, kwdtyp::STR_KEY);
        assert_eq!(v, "a''b");

        /* missing closing quote is an error */
        let (r, _, _, _) = parse("OBJECT  = 'unterminated");
        assert_eq!(r, 1);
    }

    #[test]
    fn test_parse_card_commentary() {
        let (r, t, n, _) = parse("COMMENT   some free text");
        assert_eq!(r, 0);
        assert_eq!(t, kwdtyp::COM_KEY);
        assert_eq!(n, "COMMENT");

        let (r, t, n, _) = parse("HISTORY   processed");
        assert_eq!((r, t), (0, kwdtyp::COM_KEY));
        assert_eq!(n, "HISTORY");

        /* no value indicator in columns 9-10 makes it commentary */
        let (r, t, _, _) = parse("FOO       bar");
        assert_eq!((r, t), (0, kwdtyp::COM_KEY));

        let (r, t, n, _) = parse("END");
        assert_eq!((r, t), (0, kwdtyp::COM_KEY));
        assert_eq!(n, "END");
    }

    #[test]
    fn test_parse_card_end_must_be_blank_filled() {
        let mut c = card("END");
        c[40] = b'X' as c_char;
        reset_err_wrn();
        let mut kname = [0 as c_char; FLEN_KEYWORD];
        let mut ktype = kwdtyp::UNKNOWN;
        let mut kvalue = [0 as c_char; FLEN_VALUE];
        let mut kcomm = [0 as c_char; COMM_LEN];
        let r = fits_parse_card(
            Out::Null,
            1,
            &mut c,
            &mut kname,
            &mut ktype,
            &mut kvalue,
            &mut kcomm,
        );
        assert_eq!(r, 1);
    }

    /* Fixed-format checks: the value has to end in column 30 for integers and
    sit in column 30 for logicals. */
    #[test]
    fn test_check_fixed_int() {
        reset_err_wrn();
        let good = card(&format!("BITPIX  ={}16", " ".repeat(19)));
        assert_eq!(check_fixed_int(&good, Out::Null), 1);

        let bad = card("BITPIX  = 16");
        assert_eq!(check_fixed_int(&bad, Out::Null), 0);

        let signed = card(&format!("TZERO1  ={}-8", " ".repeat(19)));
        assert_eq!(check_fixed_int(&signed, Out::Null), 1);
    }

    #[test]
    fn test_check_fixed_log() {
        reset_err_wrn();
        let good = card(&format!("SIMPLE  ={}T", " ".repeat(20)));
        assert_eq!(check_fixed_log(&good, Out::Null), 1);

        /* right column, wrong value */
        let notlog = card(&format!("SIMPLE  ={}X", " ".repeat(20)));
        assert_eq!(check_fixed_log(&notlog, Out::Null), 0);

        /* right value, wrong column */
        let shifted = card("SIMPLE  = T");
        assert_eq!(check_fixed_log(&shifted, Out::Null), 0);
    }

    #[test]
    fn test_check_fixed_str() {
        reset_err_wrn();
        let good = card("XTENSION= 'IMAGE   '");
        assert_eq!(check_fixed_str(&good, Out::Null), 1);

        /* closing quote before column 20 */
        let short = card("XTENSION= 'IM'");
        assert_eq!(check_fixed_str(&short, Out::Null), 0);

        /* does not start in column 11 */
        let late = card("XTENSION=  'IMAGE  '");
        assert_eq!(check_fixed_str(&late, Out::Null), 0);
    }

    /* The check_* helpers gate on the parsed keyword type. */
    #[test]
    fn test_check_type_helpers() {
        reset_err_wrn();
        let mut k = FitsKey {
            ktype: kwdtyp::INT_KEY,
            kindex: 1,
            ..Default::default()
        };
        set_cstr(&mut k.kname, b"NAXIS");
        set_cstr(&mut k.kvalue, b"2");
        assert_eq!(check_int(&k, Out::Null), 1);
        assert_eq!(check_flt(&k, Out::Null), 1); /* an integer is a valid float */
        assert_eq!(check_str(&k, Out::Null), 0);
        assert_eq!(check_log(&k, Out::Null), 0);

        k.ktype = kwdtyp::STR_KEY;
        assert_eq!(check_str(&k, Out::Null), 1);
        assert_eq!(check_int(&k, Out::Null), 0);

        k.ktype = kwdtyp::LOG_KEY;
        assert_eq!(check_log(&k, Out::Null), 1);

        k.ktype = kwdtyp::CMI_KEY;
        assert_eq!(check_cmi(&k, Out::Null), 1);
        assert_eq!(check_cmf(&k, Out::Null), 1);
        k.ktype = kwdtyp::CMF_KEY;
        assert_eq!(check_cmi(&k, Out::Null), 0);
        assert_eq!(check_cmf(&k, Out::Null), 1);
    }
}
