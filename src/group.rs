/*  This file, group.c, contains the grouping convention support routines.  */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
/*                                                                         */
/*  The group.c module of CFITSIO was written by Donald G. Jennings of     */
/*  the INTEGRAL Science Data Centre (ISDC) under NASA contract task       */
/*  66002J6. The above copyright laws apply. Copyright guidelines of The   */
/*  University of Geneva might also apply.                                 */

/*  The following routines are designed to create, read, and manipulate    */
/*  FITS Grouping Tables as defined in the FITS Grouping Convention paper  */
/*  by Jennings, Pence, Folk and Schlesinger. The development of the       */
/*  grouping structure was partially funded under the NASA AISRP Program.  */

use std::collections::VecDeque;

use crate::c_types::{c_char, c_int, c_long, c_uchar};
use bytemuck::cast_slice;

use std::ffi::CStr;

use crate::{
    aliases::rust_api::{fits_read_key_lng, fits_read_keyword},
    bb, cs,
    fitscore::*,
    fitsio::*,
    raw_to_slice,
    wrappers::{
        strcat_safe, strchr_safe, strcmp_safe, strcpy_safe, strlen_safe, strncpy_safe, strstr_safe,
    },
};

pub const HEX_ESCAPE: u8 = b'%';

pub const MAX_HDU_TRACKER: usize = 1000;

pub(crate) struct HDUtracker {
    nHDU: c_int,

    filename: [c_char; MAX_HDU_TRACKER],
    position: [c_int; MAX_HDU_TRACKER],

    newFilename: [c_char; MAX_HDU_TRACKER],
    newPosition: [c_int; MAX_HDU_TRACKER],
}

/*---------------------------------------------------------------------------*/
/// Create a grouping table at the end of the current FITS file.
///
/// This function makes the last HDU in the file the CHDU, then calls the
/// fits_insert_group() function to actually create the new grouping table.
/// grouping table information:
///   GT_ID_ALL_URI  0 ==> defualt (all columns)
///   GT_ID_REF      1 ==> ID by reference
///   GT_ID_POS      2 ==> ID by position
///   GT_ID_ALL      3 ==> ID by ref. and position
///   GT_ID_REF_URI 11 ==> (1) + URI info
///   GT_ID_POS_URI 12 ==> (2) + URI info  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtcr(
    fptr: *mut fitsfile,    /* FITS file pointer                         */
    grpname: *const c_char, /* name of the grouping table                */
    grouptype: c_int,       /* code specifying the type of  */
    status: *mut c_int,     /* return status code                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(grpname);

        ffgtcr_safe(fptr, grpname, grouptype, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Create a grouping table at the end of the current FITS file.
///
/// This function makes the last HDU in the file the CHDU, then calls the
/// fits_insert_group() function to actually create the new grouping table.
/// grouping table information:
///   GT_ID_ALL_URI  0 ==> defualt (all columns)
///   GT_ID_REF      1 ==> ID by reference
///   GT_ID_POS      2 ==> ID by position
///   GT_ID_ALL      3 ==> ID by ref. and position
///   GT_ID_REF_URI 11 ==> (1) + URI info
///   GT_ID_POS_URI 12 ==> (2) + URI info  
pub fn ffgtcr_safe(
    fptr: &mut fitsfile, /* FITS file pointer                         */
    grpname: &[c_char],  /* name of the grouping table                */
    grouptype: c_int,    /* code specifying the type of  */
    status: &mut c_int,  /* return status code                        */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Insert a grouping table just after the current HDU of the current FITS file.
///
/// This is the same as fits_create_group() only it allows the user to select
/// the place within the FITS file to add the grouping table.
/// grouping table information:
///   GT_ID_ALL_URI  0 ==> defualt (all columns)
///   GT_ID_REF      1 ==> ID by reference
///   GT_ID_POS      2 ==> ID by position
///   GT_ID_ALL      3 ==> ID by ref. and position
///   GT_ID_REF_URI 11 ==> (1) + URI info
///   GT_ID_POS_URI 12 ==> (2) + URI info  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtis(
    fptr: *mut fitsfile,    /* FITS file pointer                         */
    grpname: *const c_char, /* name of the grouping table                */
    grouptype: c_int,       /* code specifying the type of  */
    status: *mut c_int,     /* return status code                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(grpname);

        ffgtis_safe(fptr, grpname, grouptype, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Insert a grouping table just after the current HDU of the current FITS file.
///
/// This is the same as fits_create_group() only it allows the user to select
/// the place within the FITS file to add the grouping table.
/// grouping table information:
///   GT_ID_ALL_URI  0 ==> defualt (all columns)
///   GT_ID_REF      1 ==> ID by reference
///   GT_ID_POS      2 ==> ID by position
///   GT_ID_ALL      3 ==> ID by ref. and position
///   GT_ID_REF_URI 11 ==> (1) + URI info
///   GT_ID_POS_URI 12 ==> (2) + URI info  
pub fn ffgtis_safe(
    fptr: &mut fitsfile, /* FITS file pointer                         */
    grpname: &[c_char],  /* name of the grouping table                */
    grouptype: c_int,    /* code specifying the type of  */
    status: &mut c_int,  /* return status code                        */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Change the grouping table structure of the grouping table pointed to by gfptr.
///
/// The grouptype code specifies the new structure of the table. This
/// operation only adds or removes grouping table columns, it does not add
/// or delete group members (i.e., table rows). If the grouping table already
/// has the desired structure then no operations are performed and function   
/// simply returns with a (0) success status code. If the requested structure
/// change creates new grouping table columns, then the column values for all
/// existing members will be filled with the appropriate null values.
///   GT_ID_ALL_URI  0 ==> defualt (all columns)
///   GT_ID_REF      1 ==> ID by reference
///   GT_ID_POS      2 ==> ID by position
///   GT_ID_ALL      3 ==> ID by ref. and position
///   GT_ID_REF_URI 11 ==> (1) + URI info
///   GT_ID_POS_URI 12 ==> (2) + URI info  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtch(
    gfptr: *mut fitsfile, /* FITS file pointer                         */
    grouptype: c_int,     /* code specifying the type of  */
    status: *mut c_int,   /* return status code                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);

        ffgtch_safe(gfptr, grouptype, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Change the grouping table structure of the grouping table pointed to by gfptr.
///
/// The grouptype code specifies the new structure of the table. This
/// operation only adds or removes grouping table columns, it does not add
/// or delete group members (i.e., table rows). If the grouping table already
/// has the desired structure then no operations are performed and function   
/// simply returns with a (0) success status code. If the requested structure
/// change creates new grouping table columns, then the column values for all
/// existing members will be filled with the appropriate null values.
///   GT_ID_ALL_URI  0 ==> defualt (all columns)
///   GT_ID_REF      1 ==> ID by reference
///   GT_ID_POS      2 ==> ID by position
///   GT_ID_ALL      3 ==> ID by ref. and position
///   GT_ID_REF_URI 11 ==> (1) + URI info
///   GT_ID_POS_URI 12 ==> (2) + URI info  
pub fn ffgtch_safe(
    gfptr: &mut fitsfile, /* FITS file pointer                         */
    grouptype: c_int,     /* code specifying the type of  */
    status: &mut c_int,   /* return status code                        */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Remove a grouping table, and optionally all its members.
///
/// Any groups containing the grouping table are updated, and all members (if not
/// deleted) have their GRPIDn and GRPLCn keywords updated accordingly.
/// If the (deleted) members are members of another grouping table then those
/// tables are also updated. The CHDU of the FITS file pointed to by gfptr must
/// be positioned to the grouping table to be deleted.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtrm(
    gfptr: *mut fitsfile, /* FITS file pointer to group                   */
    rmopt: c_int,         /* code specifying if member
                          elements are to be deleted:
                          OPT_RM_GPT ==> remove only group table
                          OPT_RM_ALL ==> recursively remove members
                          and their members (if groups)                */
    status: *mut c_int, /* return status code                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);

        ffgtrm_safe(gfptr, rmopt, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Remove a grouping table, and optionally all its members.
///
/// Any groups containing the grouping table are updated, and all members (if not
/// deleted) have their GRPIDn and GRPLCn keywords updated accordingly.
/// If the (deleted) members are members of another grouping table then those
/// tables are also updated. The CHDU of the FITS file pointed to by gfptr must
/// be positioned to the grouping table to be deleted.
pub fn ffgtrm_safe(
    gfptr: &mut fitsfile, /* FITS file pointer to group                   */
    rmopt: c_int,         /* code specifying if member
                          elements are to be deleted:
                          OPT_RM_GPT ==> remove only group table
                          OPT_RM_ALL ==> recursively remove members
                          and their members (if groups)                */
    status: &mut c_int, /* return status code                           */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Copy a grouping table, and optionally all its members, to a new FITS file.
///
/// If the cpopt is set to OPT_GCP_GPT (copy grouping table only) then the
/// existing members have their GRPIDn and GRPLCn keywords updated to reflect
/// the existance of the new group, since they now belong to another group. If
/// cpopt is set to OPT_GCP_ALL (copy grouping table and members recursively)
/// then the original members are not updated; the new grouping table is
/// modified to include only the copied member HDUs and not the original members.
///
/// Note that the recursive version of this function, ffgtcpr(), is called
/// to perform the group table copy. In the case of cpopt == OPT_GCP_GPT
/// ffgtcpr() does not actually use recursion.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtcp(
    infptr: *mut fitsfile,  /* input FITS file pointer                     */
    outfptr: *mut fitsfile, /* output FITS file pointer                    */
    cpopt: c_int,           /* code specifying copy options:
                            OPT_GCP_GPT (0) ==> copy only grouping table
                            OPT_GCP_ALL (2) ==> recusrively copy members
                            and their members (if  groups)                  */
    status: *mut c_int, /* return status code                          */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffgtcp_safe(infptr, outfptr, cpopt, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Copy a grouping table, and optionally all its members, to a new FITS file.
///
/// If the cpopt is set to OPT_GCP_GPT (copy grouping table only) then the
/// existing members have their GRPIDn and GRPLCn keywords updated to reflect
/// the existance of the new group, since they now belong to another group. If
/// cpopt is set to OPT_GCP_ALL (copy grouping table and members recursively)
/// then the original members are not updated; the new grouping table is
/// modified to include only the copied member HDUs and not the original members.
///
/// Note that the recursive version of this function, ffgtcpr(), is called
/// to perform the group table copy. In the case of cpopt == OPT_GCP_GPT
/// ffgtcpr() does not actually use recursion.
pub fn ffgtcp_safe(
    infptr: &mut fitsfile,  /* input FITS file pointer                     */
    outfptr: &mut fitsfile, /* output FITS file pointer                    */
    cpopt: c_int,           /* code specifying copy options:
                            OPT_GCP_GPT (0) ==> copy only grouping table
                            OPT_GCP_ALL (2) ==> recusrively copy members
                            and their members (if  groups)                  */
    status: &mut c_int, /* return status code                          */
) -> c_int {
    todo!()
}

/*---------------------------------------------------------------------------*/
/// Merge two grouping tables by combining their members into a single table.
///
/// The source grouping table must be the CHDU of the fitsfile pointed to by
/// infptr, and the target grouping table must be the CHDU of the fitsfile to by
/// outfptr. All members of the source grouping table shall be copied to the
/// target grouping table. If the mgopt parameter is OPT_MRG_COPY then the source
/// grouping table continues to exist after the merge. If the mgopt parameter
/// is OPT_MRG_MOV then the source grouping table is deleted after the merge,
/// and all member HDUs are updated accordingly.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtmg(
    infptr: *mut fitsfile,  /* FITS file ptr to source grouping table      */
    outfptr: *mut fitsfile, /* FITS file ptr to target grouping table      */
    mgopt: c_int,           /* code specifying merge options:
                            OPT_MRG_COPY (0) ==> copy members to target group, leaving source group in place
                            OPT_MRG_MOV  (1) ==> move members to target group, source group is deleted after merge    */
    status: *mut c_int, /* return status code                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffgtmg_safe(infptr, outfptr, mgopt, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Merge two grouping tables by combining their members into a single table.
///
/// The source grouping table must be the CHDU of the fitsfile pointed to by
/// infptr, and the target grouping table must be the CHDU of the fitsfile to by
/// outfptr. All members of the source grouping table shall be copied to the
/// target grouping table. If the mgopt parameter is OPT_MRG_COPY then the source
/// grouping table continues to exist after the merge. If the mgopt parameter
/// is OPT_MRG_MOV then the source grouping table is deleted after the merge,
/// and all member HDUs are updated accordingly.
pub fn ffgtmg_safe(
    infptr: &mut fitsfile,  /* FITS file ptr to source grouping table      */
    outfptr: &mut fitsfile, /* FITS file ptr to target grouping table      */
    mgopt: c_int,           /* code specifying merge options:
                            OPT_MRG_COPY (0) ==> copy members to target group, leaving source group in place
                            OPT_MRG_MOV  (1) ==> move members to target group, source group is deleted after merge    */
    status: &mut c_int, /* return status code                         */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// "Compact" a group pointed to by the FITS file pointer gfptr.
///
/// This is achieved by flattening the tree structure of a group and its
/// (grouping table) members. All members HDUs of a grouping table which is
/// itself a member of the grouping table gfptr are added to gfptr. Optionally,
/// the grouping tables which are "compacted" are deleted. If the grouping
/// table contains no members that are themselves grouping tables then this
/// function performs a NOOP.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtcm(
    gfptr: *mut fitsfile, /* FITS file pointer to grouping table          */
    cmopt: c_int,         /* code specifying compact options
                          OPT_CMT_MBR      (1) ==> compact only direct members (if groups)
                          OPT_CMT_MBR_DEL (11) ==> (1) + delete all compacted groups    */
    status: *mut c_int, /* return status code                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);

        ffgtcm_safe(gfptr, cmopt, status)
    }
}

/*---------------------------------------------------------------------------*/
/// "Compact" a group pointed to by the FITS file pointer gfptr.
///
/// This is achieved by flattening the tree structure of a group and its
/// (grouping table) members. All members HDUs of a grouping table which is
/// itself a member of the grouping table gfptr are added to gfptr. Optionally,
/// the grouping tables which are "compacted" are deleted. If the grouping
/// table contains no members that are themselves grouping tables then this
/// function performs a NOOP.
pub fn ffgtcm_safe(
    gfptr: &mut fitsfile, /* FITS file pointer to grouping table          */
    cmopt: c_int,         /* code specifying compact options
                          OPT_CMT_MBR      (1) ==> compact only direct members (if groups)
                          OPT_CMT_MBR_DEL (11) ==> (1) + delete all compacted groups    */
    status: &mut c_int, /* return status code                           */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Check the integrity of a grouping table to make sure that all group members
/// are accessible and all the links to other grouping tables are valid. The
/// firstfailed parameter returns the member ID of the first member HDU to fail
/// verification if positive or the first group link to fail if negative;
/// otherwise firstfailed contains a return value of 0.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtvf(
    gfptr: *mut fitsfile,     /* FITS file pointer to group             */
    firstfailed: *mut c_long, /* Member ID (if positive) of first failed member HDU verify check or GRPID index (if negitive) of first failed group link verify check.                     */
    status: *mut c_int,       /* return status code                     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);
        let firstfailed = firstfailed.as_mut().expect(NULL_MSG);

        ffgtvf_safe(gfptr, firstfailed, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Check the integrity of a grouping table to make sure that all group members
/// are accessible and all the links to other grouping tables are valid. The
/// firstfailed parameter returns the member ID of the first member HDU to fail
/// verification if positive or the first group link to fail if negative;
/// otherwise firstfailed contains a return value of 0.
pub fn ffgtvf_safe(
    gfptr: &mut fitsfile,     /* FITS file pointer to group             */
    firstfailed: &mut c_long, /* Member ID (if positive) of first failed member HDU verify check or GRPID index (if negitive) of first failed group link verify check.                     */
    status: &mut c_int,       /* return status code                     */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Open the grouping table that contains the member HDU.
///
/// The member HDU must be the CHDU of the FITS file pointed to by mfptr, and the
/// grouping table is identified by the Nth index number of the GRPIDn keywords specified in
/// the member HDU's header. The fitsfile gfptr pointer is positioned with the
/// appropriate FITS file with the grouping table as the CHDU. If the group
/// grouping table resides in a file other than the member then an attempt
/// is first made to open the file readwrite, and failing that readonly.
///
/// Note that it is possible for the GRPIDn/GRPLCn keywords in a member
/// header to be non-continuous, e.g., GRPID1, GRPID2, GRPID5, GRPID6. In
/// such cases, the grpid index value specified in the function call shall
/// identify the (grpid)th GRPID value. In the above example, if grpid == 3,
/// then the group specified by GRPID5 would be opened.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtop(
    mfptr: *mut fitsfile,      /* FITS file pointer to the member HDU          */
    grpid: c_int,              /* group ID (GRPIDn index) within member HDU    */
    gfptr: *mut *mut fitsfile, /* FITS file pointer to grouping table HDU      */
    status: *mut c_int,        /* return status code                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let mfptr = mfptr.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);

        ffgtop_safe(mfptr, grpid, gfptr, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Open the grouping table that contains the member HDU.
///
/// The member HDU must be the CHDU of the FITS file pointed to by mfptr, and the
/// grouping table is identified by the Nth index number of the GRPIDn keywords specified in
/// the member HDU's header. The fitsfile gfptr pointer is positioned with the
/// appropriate FITS file with the grouping table as the CHDU. If the group
/// grouping table resides in a file other than the member then an attempt
/// is first made to open the file readwrite, and failing that readonly.
///
/// Note that it is possible for the GRPIDn/GRPLCn keywords in a member
/// header to be non-continuous, e.g., GRPID1, GRPID2, GRPID5, GRPID6. In
/// such cases, the grpid index value specified in the function call shall
/// identify the (grpid)th GRPID value. In the above example, if grpid == 3,
/// then the group specified by GRPID5 would be opened.
pub fn ffgtop_safe(
    mfptr: &mut fitsfile,      /* FITS file pointer to the member HDU          */
    grpid: c_int,              /* group ID (GRPIDn index) within member HDU    */
    gfptr: &mut *mut fitsfile, /* FITS file pointer to grouping table HDU      */
    status: &mut c_int,        /* return status code                           */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Add a member HDU to an existing grouping table.
///
/// The fitsfile pointer gfptr
/// must be positioned with the grouping table as the CHDU. The member HDU
/// may either be identifed with the fitsfile *mfptr (which must be positioned
/// to the member HDU) or the hdupos parameter (the HDU number of the member
/// HDU) if both reside in the same FITS file. The hdupos value is only used
/// if the mfptr parameter has a value of NULL (0). The new member HDU shall
/// have the appropriate GRPIDn and GRPLCn keywords created in its header.
///
/// Note that if the member HDU to be added to the grouping table is already
/// a member of the group then it will not be added a sceond time.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtam(
    gfptr: *mut fitsfile, /* FITS file pointer to grouping table HDU     */
    mfptr: *mut fitsfile, /* FITS file pointer to member HDU             */
    hdupos: c_int, /* member HDU position IF in the same file as the grouping table AND mfptr == NULL        */
    status: *mut c_int, /* return status code                          */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);
        let mfptr = mfptr.as_mut().expect(NULL_MSG);

        ffgtam_safe(gfptr, mfptr, hdupos, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Add a member HDU to an existing grouping table.
///
/// The fitsfile pointer gfptr
/// must be positioned with the grouping table as the CHDU. The member HDU
/// may either be identifed with the fitsfile *mfptr (which must be positioned
/// to the member HDU) or the hdupos parameter (the HDU number of the member
/// HDU) if both reside in the same FITS file. The hdupos value is only used
/// if the mfptr parameter has a value of NULL (0). The new member HDU shall
/// have the appropriate GRPIDn and GRPLCn keywords created in its header.
///
/// Note that if the member HDU to be added to the grouping table is already
/// a member of the group then it will not be added a sceond time.
pub fn ffgtam_safe(
    gfptr: &mut fitsfile, /* FITS file pointer to grouping table HDU     */
    mfptr: &mut fitsfile, /* FITS file pointer to member HDU             */
    hdupos: c_int, /* member HDU position IF in the same file as the grouping table AND mfptr == NULL        */
    status: &mut c_int, /* return status code                          */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Return the number of member HDUs in a grouping table.
///
/// The fitsfile pointer gfptr must be positioned with the grouping table as the CHDU.
/// The number of grouping table member HDUs is just the NAXIS2 value of the grouping
/// table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtnm(
    gfptr: *mut fitsfile,  /* FITS file pointer to grouping table        */
    nmembers: *mut c_long, /* member count of the grouping table         */
    status: *mut c_int,    /* return status code                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);
        let nmembers = nmembers.as_mut().expect(NULL_MSG);

        ffgtnm_safe(gfptr, nmembers, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Return the number of member HDUs in a grouping table.
///
/// The fitsfile pointer gfptr must be positioned with the grouping table as the CHDU.
/// The number of grouping table member HDUs is just the NAXIS2 value of the grouping
/// table.
pub fn ffgtnm_safe(
    gfptr: &mut fitsfile,  /* FITS file pointer to grouping table        */
    nmembers: &mut c_long, /* member count of the grouping table         */
    status: &mut c_int,    /* return status code                         */
) -> c_int {
    let mut keyvalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status != 0 {
        return *status;
    }

    *status = fits_read_keyword(
        gfptr,
        cs!(c"EXTNAME"),
        &mut keyvalue,
        Some(&mut comment),
        status,
    );

    if *status == KEY_NO_EXIST {
        *status = NOT_GROUP_TABLE;
    } else {
        prepare_keyvalue(&keyvalue);

        if fits_strcasecmp(&keyvalue, cs!(c"GROUPING")) != 0 {
            *status = NOT_GROUP_TABLE;
            ffpmsg_str("Specified HDU is not a Grouping table (ffgtnm)");
        }

        *status = fits_read_key_lng(gfptr, cs!(c"NAXIS2"), nmembers, Some(&mut comment), status);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Return the number of groups to which a HDU belongs, as defined by the number
/// of GRPIDn/GRPLCn keyword records that appear in the HDU header. The
/// fitsfile pointer mfptr must be positioned with the member HDU as the CHDU.
/// Each time this function is called, the indicies of the GRPIDn/GRPLCn
/// keywords are checked to make sure they are continuous (ie no gaps) and
/// are re-enumerated to eliminate gaps if gaps are found to be present.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgmng(
    mfptr: *mut fitsfile, /* FITS file pointer to member HDU            */
    ngroups: *mut c_long, /* total number of groups linked to HDU       */
    status: *mut c_int,   /* return status code                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let mfptr = mfptr.as_mut().expect(NULL_MSG);
        let ngroups = ngroups.as_mut().expect(NULL_MSG);

        ffgmng_safe(mfptr, ngroups, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Return the number of groups to which a HDU belongs, as defined by the number
/// of GRPIDn/GRPLCn keyword records that appear in the HDU header. The
/// fitsfile pointer mfptr must be positioned with the member HDU as the CHDU.
/// Each time this function is called, the indicies of the GRPIDn/GRPLCn
/// keywords are checked to make sure they are continuous (ie no gaps) and
/// are re-enumerated to eliminate gaps if gaps are found to be present.
pub fn ffgmng_safe(
    mfptr: &mut fitsfile, /* FITS file pointer to member HDU            */
    ngroups: *mut c_long, /* total number of groups linked to HDU       */
    status: &mut c_int,   /* return status code                         */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// open a grouping table member, returning a pointer to the member's FITS file
/// with the CHDU set to the member HDU. The grouping table must be the CHDU of
/// the FITS file pointed to by gfptr. The member to open is identified by its
/// row number within the grouping table (first row/member == 1).
///
/// If the member resides in a FITS file different from the grouping
/// table the member file is first opened readwrite and if this fails then
/// it is opened readonly. For access type of FILE:// the member file is
/// searched for assuming (1) an absolute path is given, (2) a path relative
/// to the CWD is given, and (3) a path relative to the grouping table file
/// but not relative to the CWD is given. If all of these fail then the
/// error FILE_NOT_FOUND is returned.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgmop(
    gfptr: *mut fitsfile,      /* FITS file pointer to grouping table          */
    member: c_long,            /* member ID (row num) within grouping table    */
    mfptr: *mut *mut fitsfile, /* FITS file pointer to member HDU              */
    status: *mut c_int,        /* return status code                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);
        let mfptr = mfptr.as_mut().expect(NULL_MSG);

        ffgmop_safe(gfptr, member, mfptr, status)
    }
}

/*---------------------------------------------------------------------------*/
/// open a grouping table member, returning a pointer to the member's FITS file
/// with the CHDU set to the member HDU. The grouping table must be the CHDU of
/// the FITS file pointed to by gfptr. The member to open is identified by its
/// row number within the grouping table (first row/member == 1).
///
/// If the member resides in a FITS file different from the grouping
/// table the member file is first opened readwrite and if this fails then
/// it is opened readonly. For access type of FILE:// the member file is
/// searched for assuming (1) an absolute path is given, (2) a path relative
/// to the CWD is given, and (3) a path relative to the grouping table file
/// but not relative to the CWD is given. If all of these fail then the
/// error FILE_NOT_FOUND is returned.
pub fn ffgmop_safe(
    gfptr: *mut fitsfile,      /* FITS file pointer to grouping table          */
    member: c_long,            /* member ID (row num) within grouping table    */
    mfptr: *mut *mut fitsfile, /* FITS file pointer to member HDU              */
    status: *mut c_int,        /* return status code                           */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// copy a member HDU of a grouping table to a new FITS file. The grouping table
/// must be the CHDU of the FITS file pointed to by gfptr. The copy of the
/// group member shall be appended to the end of the FITS file pointed to by
/// mfptr. If the cpopt parameter is set to OPT_MCP_ADD then the copy of the
/// member is added to the grouping table as a new member, if OPT_MCP_NADD
/// then the copied member is not added to the grouping table, and if
/// OPT_MCP_REPL then the copied member is used to replace the original member.
/// The copied member HDU also has its EXTVER value updated so that its
/// combination of XTENSION, EXTNAME and EXVTER is unique within its new
/// FITS file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgmcp(
    gfptr: *mut fitsfile, /* FITS file pointer to group                   */
    mfptr: *mut fitsfile, /* FITS file pointer to new member FITS file                                    */
    member: c_long,       /* member ID (row num) within grouping table    */
    cpopt: c_int,         /* code specifying copy options:
                          OPT_MCP_ADD  (0) ==> add copied member to the
                                              grouping table
                          OPT_MCP_NADD (1) ==> do not add member copy to
                                              the grouping table
                          OPT_MCP_REPL (2) ==> replace current member
                                              entry with member copy  */
    status: *mut c_int, /* return status code                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let mfptr = mfptr.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);

        ffgmcp_safe(gfptr, mfptr, member, cpopt, status)
    }
}

/*---------------------------------------------------------------------------*/
/// copy a member HDU of a grouping table to a new FITS file. The grouping table
/// must be the CHDU of the FITS file pointed to by gfptr. The copy of the
/// group member shall be appended to the end of the FITS file pointed to by
/// mfptr. If the cpopt parameter is set to OPT_MCP_ADD then the copy of the
/// member is added to the grouping table as a new member, if OPT_MCP_NADD
/// then the copied member is not added to the grouping table, and if
/// OPT_MCP_REPL then the copied member is used to replace the original member.
/// The copied member HDU also has its EXTVER value updated so that its
/// combination of XTENSION, EXTNAME and EXVTER is unique within its new
/// FITS file.
pub fn ffgmcp_safe(
    gfptr: &mut fitsfile, /* FITS file pointer to group                   */
    mfptr: &mut fitsfile, /* FITS file pointer to new member FITS file                                    */
    member: c_long,       /* member ID (row num) within grouping table    */
    cpopt: c_int,         /* code specifying copy options:
                          OPT_MCP_ADD  (0) ==> add copied member to the
                                              grouping table
                          OPT_MCP_NADD (1) ==> do not add member copy to
                                              the grouping table
                          OPT_MCP_REPL (2) ==> replace current member
                                              entry with member copy  */
    status: &mut c_int, /* return status code                           */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// transfer a group member from one grouping table to another. The source
/// grouping table must be the CHDU of the fitsfile pointed to by infptr, and
/// the destination grouping table must be the CHDU of the fitsfile to by
/// outfptr. If the tfopt parameter is OPT_MCP_ADD then the member is made a
/// member of the target group and remains a member of the source group. If
/// the tfopt parameter is OPT_MCP_MOV then the member is deleted from the
/// source group after the transfer to the destination group. The member to be
/// transfered is identified by its row number within the source grouping table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgmtf(
    infptr: *mut fitsfile,  /* FITS file pointer to source grouping table */
    outfptr: *mut fitsfile, /* FITS file pointer to target grouping table */
    member: c_long,         /* member ID within source grouping table     */
    tfopt: c_int,           /* code specifying transfer opts:
                            OPT_MCP_ADD (0) ==> copy member to dest.
                            OPT_MCP_MOV (3) ==> move member to dest.   */
    status: *mut c_int, /* return status code                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffgmtf_safe(infptr, outfptr, member, tfopt, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Transfer a group member from one grouping table to another. The source
/// grouping table must be the CHDU of the fitsfile pointed to by infptr, and
/// the destination grouping table must be the CHDU of the fitsfile to by
/// outfptr. If the tfopt parameter is OPT_MCP_ADD then the member is made a
/// member of the target group and remains a member of the source group. If
/// the tfopt parameter is OPT_MCP_MOV then the member is deleted from the
/// source group after the transfer to the destination group. The member to be
/// transfered is identified by its row number within the source grouping table.
pub fn ffgmtf_safe(
    infptr: &mut fitsfile,  /* FITS file pointer to source grouping table */
    outfptr: &mut fitsfile, /* FITS file pointer to target grouping table */
    member: c_long,         /* member ID within source grouping table     */
    tfopt: c_int,           /* code specifying transfer opts:
                            OPT_MCP_ADD (0) ==> copy member to dest.
                            OPT_MCP_MOV (3) ==> move member to dest.   */
    status: &mut c_int, /* return status code                         */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// remove a member HDU from a grouping table. The fitsfile pointer gfptr must
/// be positioned with the grouping table as the CHDU, and the member to
/// delete is identified by its row number in the table (first member == 1).
/// The rmopt parameter determines if the member entry is deleted from the
/// grouping table (in which case GRPIDn and GRPLCn keywords in the member
/// HDU's header shall be updated accordingly) or if the member HDU shall
/// itself be removed from its FITS file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgmrm(
    gfptr: *mut fitsfile, /* FITS file pointer to group table             */
    member: c_long,       /* member ID (row num) in the group             */
    rmopt: c_int,         /* code specifying the delete option:
                          OPT_RM_ENTRY ==> delete the member entry
                          OPT_RM_MBR   ==> delete entry and member HDU */
    status: *mut c_int, /* return status code                          */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let gfptr = gfptr.as_mut().expect(NULL_MSG);

        ffgmrm_safe(gfptr, member, rmopt, status)
    }
}

/*---------------------------------------------------------------------------*/
/// remove a member HDU from a grouping table. The fitsfile pointer gfptr must
/// be positioned with the grouping table as the CHDU, and the member to
/// delete is identified by its row number in the table (first member == 1).
/// The rmopt parameter determines if the member entry is deleted from the
/// grouping table (in which case GRPIDn and GRPLCn keywords in the member
/// HDU's header shall be updated accordingly) or if the member HDU shall
/// itself be removed from its FITS file.
pub fn ffgmrm_safe(
    gfptr: &mut fitsfile, /* FITS file pointer to group table             */
    member: c_long,       /* member ID (row num) in the group             */
    rmopt: c_int,         /* code specifying the delete option:
                          OPT_RM_ENTRY ==> delete the member entry
                          OPT_RM_MBR   ==> delete entry and member HDU */
    status: &mut c_int, /* return status code                          */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------
                 Grouping Table support functions
---------------------------------------------------------------------------*/

/****************************************************************************/
/// Examine the grouping table pointed to by gfptr and determine the column
/// index ID of each possible grouping column. If a column is not found then
/// an index of 0 is returned. the grptype parameter returns the structure
/// of the grouping table ==> what columns are defined.
pub(crate) fn ffgtgc(
    gfptr: &mut fitsfile,    /* pointer to the grouping table                */
    xtensionCol: &mut c_int, /* column ID of the MEMBER_XTENSION column      */
    extnameCol: &mut c_int,  /* column ID of the MEMBER_NAME column          */
    extverCol: &mut c_int,   /* column ID of the MEMBER_VERSION column       */
    positionCol: &mut c_int, /* column ID of the MEMBER_POSITION column      */
    locationCol: &mut c_int, /* column ID of the MEMBER_LOCATION column      */
    uriCol: &mut c_int,      /* column ID of the MEMBER_URI_TYPE column      */
    grptype: &mut c_int,     /* group structure type code specifying the
                             grouping table columns that are defined:
                             GT_ID_ALL_URI  (0) ==> all columns defined
                             GT_ID_REF      (1) ==> reference cols only
                             GT_ID_POS      (2) ==> position col only
                             GT_ID_ALL      (3) ==> ref & pos cols
                             GT_ID_REF_URI (11) ==> ref & loc cols
                             GT_ID_POS_URI (12) ==> pos & loc cols        */
    status: &mut c_int, /* return status code                           */
) -> c_int {
    todo!();
}

/*****************************************************************************/
/// Perform validation on column formats to ensure this matches the grouping
/// format the get functions expect.  Particularly want to check widths of
/// string columns.
pub(crate) fn ffvcfm(
    gfptr: &mut fitsfile,
    xtensionCol: c_int,
    extnameCol: c_int,
    extverCol: c_int,
    positionCol: c_int,
    locationCol: c_int,
    uriCol: c_int,
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*****************************************************************************/
/// Create the TTYPE and TFORM values for the grouping table according to the
/// value of the grouptype parameter and the values of the *col flags. The
/// resulting TTYPE and TFORM are returned in ttype[] and tform[] respectively.
/// The number of TTYPE and TFORMs returned is given by ncols. Both the TTYPE[]
/// and TTFORM[] arrays must contain enough pre-allocated strings to hold
/// the returned information.
pub(crate) fn ffgtdc(
    grouptype: c_int,            /* code specifying the type of
                                 grouping table information:
                                 GT_ID_ALL_URI  0 ==> defualt (all columns)
                                 GT_ID_REF      1 ==> ID by reference
                                 GT_ID_POS      2 ==> ID by position
                                 GT_ID_ALL      3 ==> ID by ref. and position
                                 GT_ID_REF_URI 11 ==> (1) + URI info
                                 GT_ID_POS_URI 12 ==> (2) + URI info       */
    xtensioncol: c_int,          /* does MEMBER_XTENSION already exist?         */
    extnamecol: c_int,           /* does MEMBER_NAME aleady exist?              */
    extvercol: c_int,            /* does MEMBER_VERSION already exist?          */
    positioncol: c_int,          /* does MEMBER_POSITION already exist?         */
    locationcol: c_int,          /* does MEMBER_LOCATION already exist?         */
    uricol: c_int,               /* does MEMBER_URI_TYPE aleardy exist?         */
    ttype: &[Option<&[c_char]>], /* array of grouping table column TTYPE names to define (if *col var false)               */
    tform: &[Option<&[c_char]>], /* array of grouping table column TFORM values to define (if*col variable false)           */
    ncols: &mut c_int,           /* number of TTYPE and TFORM values returned   */
    status: &mut c_int,          /* return status code                          */
) -> c_int {
    todo!();
}

/*****************************************************************************/
/// Examine all the GRPIDn and GRPLCn keywords in the member HDUs header
/// and remove the member from the grouping tables referenced; This
/// effectively "unlinks" the member from all of its groups. The rmopt
/// specifies if the GRPIDn/GRPLCn keywords are to be removed from the
/// member HDUs header after the unlinking.
pub(crate) fn ffgmul(
    mfptr: &mut fitsfile, /* pointer to the grouping table member HDU    */
    rmopt: c_int,         /* 0 ==> leave GRPIDn/GRPLCn keywords,
                          1 ==> remove GRPIDn/GRPLCn keywords         */
    status: &mut c_int, /* return status code                          */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Try to find the entry for the member HDU defined by the xtension, extname,
/// extver, position, and location parameters within the grouping table
/// pointed to by gfptr. If the member HDU is found then its ID (row number)
/// within the grouping table is returned in the member variable; if not
/// found then member is returned with a value of 0 and the status return
/// code will be set to MEMBER_NOT_FOUND.
///
/// Note that the member HDU postion information is used to obtain a member
/// match only if the grouping table type is GT_ID_POS_URI or GT_ID_POS. This
/// is because the position information can become invalid much more
/// easily then the reference information for a group member.
pub(crate) fn ffgmf(
    gfptr: &mut fitsfile,        /* pointer to grouping table HDU to search       */
    xtension: &[c_char],         /* XTENSION value for member HDU                */
    extname: &[c_char],          /* EXTNAME value for member HDU                 */
    extver: c_int,               /* EXTVER value for member HDU                  */
    position: c_int,             /* HDU position value for member HDU            */
    location: &mut [c_char],     /* FITS file location value for member HDU      */
    member: Option<&mut c_long>, /* member HDU ID within group table (if found)  */
    status: &mut c_int,          /* return status code                           */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------
                      Recursive Group Functions
--------------------------------------------------------------------------*/

/****************************************************************************/
/// Recursively remove a grouping table and all its members. Each member of
/// the grouping table pointed to by gfptr it processed. If the member is itself
/// a grouping table then ffgtrmr() is recursively called to process all
/// of its members. The HDUtracker struct *HDU is used to make sure a member
/// is not processed twice, thus avoiding an infinite loop (e.g., a grouping
/// table contains itself as a member).
pub(crate) fn ffgtrmr(
    gfptr: &mut fitsfile, /* FITS file pointer to group               */
    HDU: &mut HDUtracker, /* list of processed HDUs                   */
    status: &mut c_int,   /* return status code                       */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Copy a Group to a new FITS file. If the cpopt parameter is set to
/// OPT_GCP_GPT (copy grouping table only) then the existing members have their
/// GRPIDn and GRPLCn keywords updated to reflect the existance of the new group,
/// since they now belong to another group. If cpopt is set to OPT_GCP_ALL
/// (copy grouping table and members recursively) then the original members are
/// not updated; the new grouping table is modified to include only the copied
/// member HDUs and not the original members.
///
/// Note that this function is recursive. When copt is OPT_GCP_ALL it will call
/// itself whenever a member HDU of the current grouping table is itself a
/// grouping table (i.e., EXTNAME = 'GROUPING').
pub(crate) fn ffgtcpr(
    infptr: &mut fitsfile,  /* input FITS file pointer                 */
    outfptr: &mut fitsfile, /* output FITS file pointer                */
    cpopt: c_int,           /* code specifying copy options:
                            OPT_GCP_GPT (0) ==> cp only grouping table
                            OPT_GCP_ALL (2) ==> recusrively copy
                            members and their members (if groups)   */
    HDU: &mut [HDUtracker], /* list of already copied HDUs             */
    status: &mut c_int,     /* return status code                      */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------
              HDUtracker struct manipulation functions
--------------------------------------------------------------------------*/

/****************************************************************************/
/// add an HDU to the HDUtracker struct pointed to by HDU. The HDU is only
/// added if it does not already reside in the HDUtracker. If it already
/// resides in the HDUtracker then the new HDU postion and file name are
/// returned in  newPosition and newFileName (if != NULL)
pub(crate) fn fftsad(
    mfptr: &mut fitsfile,       /* pointer to an member HDU             */
    HDU: &mut HDUtracker,       /* pointer to an HDU tracker struct     */
    newPosition: &mut c_int,    /* new HDU position of the member HDU   */
    newFileName: &mut [c_char], /* file containing member HDU           */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Update the HDU information in the HDUtracker struct pointed to by HDU. The
/// HDU to update is pointed to by mfptr. If non-zero, the value of newPosition
/// is used to update the HDU->newPosition[] value for the mfptr, and if
/// non-NULL the newFileName value is used to update the HDU->newFilename[]
/// value for mfptr.
pub(crate) fn fftsud(
    mfptr: &mut fitsfile,       /* pointer to an member HDU             */
    HDU: &mut HDUtracker,       /* pointer to an HDU tracker struct     */
    newPosition: c_int,         /* new HDU position of the member HDU   */
    newFileName: &mut [c_char], /* file containing member HDU           */
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Strip off all single quote characters "'" and blank spaces from a keyword
/// value retrieved via fits_read_key*() routines
///
/// this is necessary so that a standard comparision of keyword values may
/// be made
pub(crate) fn prepare_keyvalue(keyvalue: &[c_char]) /* string containing keyword value     */
{
    todo!();
}

/*---------------------------------------------------------------------------
      Host dependent directory path to/from URL functions
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/// Convert a file path into its Unix-style equivelent for URL
/// purposes. Note that this process is platform dependent. This
/// function supports Unix, MSDOS/WIN32, VMS and Macintosh platforms.
/// The plaform dependant code is conditionally compiled depending upon
/// the setting of the appropriate C preprocessor macros.
pub(crate) fn fits_path2url(
    inpath: &[c_char],      /* input file path string                  */
    maxlength: usize, /* I max number of chars that can be written to output, including terminating NULL */
    outpath: &mut [c_char], /* output file path string                 */
    status: &mut c_int,
) -> c_int {
    let mut buff: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];

    if cfg!(target_os = "winnt") {
        /*
        Microsoft Windows NT case. We assume input file paths of the form:

        //disk/path/filename

        All path segments may be null, so that a single file name is the
        simplist case.

        The leading "//" becomes a single "/" if present. If no "//" is present,
        then make sure the resulting URL path is relative, i.e., does not
        begin with a "/". In other words, the only way that an absolute URL
        file path may be generated is if the drive specification is given.
        */

        if *status > 0 {
            return *status;
        }

        if inpath[0] == bb(b'/') {
            strcpy_safe(&mut buff, &inpath[1..]);
        } else {
            strcpy_safe(&mut buff, inpath);
        }
    } else if cfg!(target_os = "windows") {
        /*
        MSDOS or Microsoft windows/NT case. The assumed form of the
        input path is:

        disk:\path\filename

        All path segments may be null, so that a single file name is the
        simplist case.

        All back-slashes '\' become slashes '/'; if the path starts with a
        string of the form "X:" then it is replaced with "/X/"
        */

        let mut i: usize = 0;
        let mut j: usize = 0;
        let mut k: isize = 0;
        let mut size = 0;

        if *status != 0 {
            return *status;
        }

        size = strlen_safe(inpath);
        buff[0] = 0;

        while i < size {
            match inpath[i] as u8 {
                b':' => {
                    /*
                    must be a disk desiginator; add a slash '/' at the start of
                    outpath to designate that the path is absolute, then change
                    the colon ':' to a slash '/'
                    */

                    k = j as isize;
                    while k >= 0 {
                        buff[(k + 1) as usize] = buff[k as usize];
                        k -= 1;
                    }

                    buff[0] = bb(b'/');
                    strcat_safe(&mut buff, cs!(c"/"));
                    i += 1;
                }

                b'\\' => {
                    /* just replace the '\' with a '/' IF its not the first character */

                    if i != 0 && buff[if j == 0 { 0 } else { j - 1 }] != bb(b'/') {
                        buff[j] = bb(b'/');
                        buff[j + 1] = 0;
                    }

                    i += 1;
                }
                _ => {
                    /* copy the character from inpath to buff as is */

                    buff[j] = inpath[i];
                    buff[j + 1] = 0;
                    i += 1;
                }
            }
            j = strlen_safe(&buff);
        }
    } else if cfg!(target_os = "vms") {
        /*
        VMS case. Assumed format of the input path is:

        node::disk:[path]filename.ext;version

        Any part of the file path may be missing, so that in the simplist
        case a single file name/extension is given.

        all brackets "[", "]" and dots "." become "/"; dashes "-" become "..",
        all single colons ":" become ":/", all double colons "::" become
        "FILE://"
        */
        todo!();
    /*
    int i,j,k;
    int done;
    int size;

    if(*status > 0) return(*status);

    /* see if inpath contains a directory specification */

    if(strchr(inpath,']') == NULL)
    done = 1;
    else
    done = 0;

    for(i = 0, j = 0, size = strlen(inpath), buff[0] = 0;
                 i < size && j < FLEN_FILENAME - 8; j = strlen(buff))
    {
    switch(inpath[i])
    {

    case ':':

    /*
    must be a logical/symbol separator or (in the case of a double
    colon "::") machine node separator
    */

    if(inpath[i+1] == ':')
    {
    /* insert a "FILE://" at the start of buff ==> machine given */

    for(k = j; k >= 0; --k) buff[k+7] = buff[k];
    strncpy(buff,"FILE://",7);
    i += 2;
    }
    else if(strstr(buff,"FILE://") == NULL)
    {
    /* insert a "/" at the start of buff ==> absolute path */

    for(k = j; k >= 0; --k) buff[k+1] = buff[k];
    buff[0] = '/';
    ++i;
    }
    else
    ++i;

    /* a colon always ==> path separator */

    strcat(buff,"/");

    break;

    case ']':

    /* end of directory spec, file name spec begins after this */

    done = 1;

    buff[j]   = '/';
    buff[j+1] = 0;
    ++i;

    break;

    case '[':

    /*
    begin directory specification; add a '/' only if the last char
    is not '/'
    */

    if(i != 0 && buff[(j == 0 ? 0 : j-1)] != '/')
    {
    buff[j]   = '/';
    buff[j+1] = 0;
    }

    ++i;

    break;

    case '.':

    /*
    directory segment separator or file name/extension separator;
    we decide which by looking at the value of done
    */

    if(!done)
    {
    /* must be a directory segment separator */
    if(inpath[i-1] == '[')
    {
    strcat(buff,"./");
    ++j;
    }
    else
    buff[j] = '/';
    }
    else
    /* must be a filename/extension separator */
    buff[j] = '.';

    buff[j+1] = 0;

    ++i;

    break;

    case '-':

    /*
    a dash is the same as ".." in Unix speak, but lets make sure
    that its not part of the file name first!
    */

    if(!done)
    /* must be part of the directory path specification */
    strcat(buff,"..");
    else
    {
    /* the dash is part of the filename, so just copy it as is */
    buff[j] = '-';
    buff[j+1] = 0;
    }

    ++i;

    break;

    default:

    /* nothing special, just copy the character as is */

    buff[j]   = inpath[i];
    buff[j+1] = 0;

    ++i;

    break;

    }
    }

    if(j > FLEN_FILENAME - 8)
    {
    *status = URL_PARSE_ERROR;
    ffpmsg("resulting path to URL conversion too big (fits_path2url)");
    }
    */
    } else if cfg!(target_os = "macos") {
        /*
        MacOS case. The assumed form of the input path is:

        disk:path:filename

        It is assumed that all paths are absolute with disk and path specified,
        unless no colons ":" are supplied with the string ==> a single file name
        only. All colons ":" become slashes "/", and if one or more colon is
        encountered then the path is specified as absolute.
        */
        todo!();
    /*
    int i,j,k;
    int firstColon;
    int size;

    if(*status > 0) return(*status);

    for(i = 0, j = 0, firstColon = 1, size = strlen(inpath), buff[0] = 0;
                                         i < size; j = strlen(buff))
    {
    switch(inpath[i])
    {

    case ':':

    /*
    colons imply path separators. If its the first colon encountered
    then assume that its the disk designator and add a slash to the
    beginning of the buff string
    */

    if(firstColon)
    {
    firstColon = 0;

    for(k = j; k >= 0; --k) buff[k+1] = buff[k];
    buff[0] = '/';
    }

    /* all colons become slashes */

    strcat(buff,"/");

    ++i;

    break;

    default:

    /* copy the character from inpath to buff as is */

    buff[j]   = inpath[i];
    buff[j+1] = 0;

    ++i;

    break;
    }
    }
    */
    } else {
        /*
        Default Unix case.

        Nothing special to do here except to remove the double or more // and
        replace them with single /
        */

        let mut ii = 0;
        let mut jj = 0;

        if *status > 0 {
            return *status;
        }

        while inpath[ii] > 0 {
            if inpath[ii] == bb(b'/') && inpath[ii + 1] == bb(b'/') {
                /* do nothing */
            } else {
                buff[jj] = inpath[ii];
                jj += 1;
            }
            ii += 1;
        }
        buff[jj] = 0;
        /* printf("buff is %s\ninpath is %s\n",buff,inpath); */
        /* strcpy(buff,inpath); */
    }

    /*
    encode all "unsafe" and "reserved" URL characters
    */

    *status = fits_encode_url(&buff, maxlength, outpath, status);

    *status
}

/*---------------------------------------------------------------------------*/
/// Convert a Unix-style URL into a platform dependent directory path.
/// Note that this process is platform dependent. This
/// function supports Unix, MSDOS/WIN32, VMS and Macintosh platforms. Each
/// platform dependent code segment is conditionally compiled depending
/// upon the setting of the appropriate C preprocesser macros.
pub(crate) fn fits_url2path(
    inpath: &[c_char],      /* input file path string  */
    outpath: &mut [c_char], /* output file path string */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/****************************************************************************/
/// Retrieve the string containing the current working directory absolute
/// path in Unix-like URL standard notation. It is assumed that the CWD
/// string has a size of at least FLEN_FILENAME.
///
/// Note that this process is platform dependent. This
/// function supports Unix, MSDOS/WIN32, VMS and Macintosh platforms. Each
/// platform dependent code segment is conditionally compiled depending
/// upon the setting of the appropriate C preprocesser macros.
pub(crate) fn fits_get_cwd(
    cwd: &mut [c_char], /* IO current working directory string */
    status: &mut c_int,
) -> c_int {
    let mut buff: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
    let current_dir;

    if cwd.len() < FLEN_FILENAME {
        let msg = format!(
            "fits_get_cwd expects the buffer to be at least FLEN_FILENAME={FLEN_FILENAME} bytes long"
        );
        ffpmsg_str(&msg);
        *status = URL_PARSE_ERROR;
        return *status;
    }

    if *status != 0 {
        return *status;
    }

    if cfg!(target_os = "macos") {
        /*
        MacOS case. Currently unknown !!!!
        */
        todo!();
        buff[0] = 0;
    } else {
        /*
        Good old getcwd() seems to work with all other platforms
        */

        current_dir = std::env::current_dir();

        if current_dir.is_err() {
            cwd[0] = 0;
            ffpmsg_str("Path and file name too long (fits_get_cwd)");
            *status = URL_PARSE_ERROR;
            return *status;
        }
    }

    /*
    convert the cwd string to a URL standard path string
    */
    let current_dir = current_dir.unwrap();
    let cwd_bytes = current_dir.to_str().unwrap().as_bytes();
    for i in 0..cwd_bytes.len() {
        buff[i] = cwd_bytes[i] as c_char;
    }
    buff[cwd_bytes.len()] = 0;

    fits_path2url(&buff, FLEN_FILENAME, cwd, status);

    *status
}

/*---------------------------------------------------------------------------*/
/// For grouping convention purposes, determine the URL of the FITS file
/// associated with the fitsfile pointer fptr. The true access type (file://,
/// mem://, shmem://, root://), starting "official" access type, and iostate
/// (0 ==> readonly, 1 ==> readwrite) are also returned.
///
/// It is assumed that the url string has enough room to hold the resulting
/// URL, and the the accessType string has enough room to hold the access type.
pub(crate) fn fits_get_url(
    fptr: &mut fitsfile,      /* I ptr to FITS file to evaluate    */
    realURL: &mut [c_char],   /* O URL of real FITS file           */
    startURL: &mut [c_char],  /* O URL of starting FITS file       */
    realAccess: &mut c_char,  /* O true access method of FITS file */
    startAccess: &mut c_char, /* O "official" access of FITS file  */
    iostate: &mut c_int,      /* O can this file be modified?      */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------
                         URL parse support functions
--------------------------------------------------------------------------*/

// HAVE NOT INCLUDED THIS SET OF FUNTIONS.

/*--------------------------------------------------------------------------*/
/// Clean the URL by eliminating any ".." or "." specifiers in the inURL
/// string, and write the output to the outURL string.
///
/// Note that this function must have a valid Unix-style URL as input; platform
/// dependent path strings are not allowed.
pub(crate) fn fits_clean_url(
    inURL: &[c_char],      /* I input URL string                      */
    outURL: &mut [c_char], /* O output URL string                     */
    status: &mut c_int,
) -> c_int {
    let mut mystack: VecDeque<&[c_char]> = VecDeque::new(); /* stack to hold pieces of URL */

    let mut inURL_ptr = inURL;

    if *status > 0 {
        return *status;
    }

    outURL[0] = 0;

    loop {
        /* handle URL scheme and domain if they exist */
        let tmp = strstr_safe(inURL_ptr, cs!(c"://"));

        if let Some(tmp) = tmp {
            // tmp is now the index into the current inURL_ptr

            /* there is a URL scheme, so look for the end of the domain too */
            let tmp_inner = strchr_safe(&inURL_ptr[(tmp + 3)..], b'/' as c_char);

            if let Some(tmp_inner) = tmp_inner {
                // tmp_inner is the index into tmp
                // adjust so that its the index into inURL_ptr

                let tmp_inner = tmp + tmp_inner + 3;

                /* tmp is now the end of the domain, so copy URL scheme and domain as is, and terminate by hand */
                let string_size = tmp_inner;

                strncpy_safe(outURL, &inURL_ptr[tmp_inner..], string_size);
                outURL[string_size] = 0;

                /* now advance the input pointer to just after the domain and go on */
                inURL_ptr = &inURL_ptr[tmp_inner..];
            } else {
                /* '/' was not found, which means there are no path-like
                 * portions, so copy whole inURL to outURL and we're done */
                strcpy_safe(outURL, inURL);
                continue; /* while(0) */
            }
        }

        /* explicitly copy a leading / (absolute path) */
        if bb(b'/') == inURL_ptr[0] {
            strcat_safe(outURL, cs!(c"/"));
        }

        /* now clean the remainder of the inURL. push URL segments onto stack, dealing with .. and . as we go */
        let first_null = inURL_ptr
            .iter()
            .position(|&x| x == 0)
            .unwrap_or(inURL_ptr.len());
        let tokens = inURL_ptr[..first_null].split(|f| (*f == b'/' as c_char) || (*f == 0));

        for tmp in tokens {
            if tmp.is_empty() {
                continue;
            }

            if strcmp_safe(tmp, cast_slice(b"..")) == 0 {
                // WARNING tokens are not null-terminated
                /*
                discard previous URL segment, if there was one. if not,
                add the .. to the stack if this is *not* an absolute path
                (for absolute paths, leading .. has no effect, so skip it)
                */
                if !mystack.is_empty() {
                    mystack.pop_back();
                } else if bb(b'/') != inURL_ptr[0] {
                    mystack.push_back(tmp);
                }
            } else {
                /* always just skip ., but otherwise add segment to stack */
                if strcmp_safe(tmp, cast_slice(b".")) > 0 {
                    mystack.push_back(tmp);
                }
            }
        }

        /*
        stack now has pieces of cleaned URL, so just catenate them
        onto output string until stack is empty
        */
        while !mystack.is_empty() {
            let tmp = mystack.pop_front().unwrap(); // This will be missing all the null terminators

            if strlen_safe(outURL) + tmp.len() + 1 > FLEN_FILENAME - 1 {
                outURL[0] = 0;
                ffpmsg_str("outURL is too long (fits_clean_url)");
                *status = URL_PARSE_ERROR;
                return *status;
            }

            // Note: We need to manually keep track of the len since we aren't using C null terminators here
            let len = strlen_safe(outURL) + tmp.len();

            let mut c = Vec::from(tmp);
            c.push(0); // Add a null terminator to the end of the tmp slice
            let tmp = c.as_slice(); // Convert back to a slice

            strcat_safe(outURL, tmp);

            outURL[len] = 0; // Manually add the null terminator

            strcat_safe(outURL, cs!(c"/"));
        }

        let outurl_len = strlen_safe(outURL);
        if outurl_len > 0 {
            outURL[strlen_safe(outURL) - 1] = 0; /* blank out trailing / */
        }
        break;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Create a relative URL to the file referenced by absURL with respect to the
/// reference URL refURL. The relative URL is returned in relURL.
///
/// Both refURL and absURL must be absolute URL strings; i.e. either begin
/// with an access method specification "XXX://" or with a '/' character
/// signifiying that they are absolute file paths.
///
/// Note that it is possible to make a relative URL from two input URLs
/// (absURL and refURL) that are not compatable. This function does not
/// check to see if the resulting relative URL makes any sence. For instance,
/// it is impossible to make a relative URL from the following two inputs:
///
/// absURL = ftp://a.b.c.com/x/y/z/foo.fits
/// refURL = /a/b/c/ttt.fits
///
/// The resulting relURL will be:
///
/// ../../../ftp://a.b.c.com/x/y/z/foo.fits
///
/// Which is syntically correct but meaningless. The problem is that a file
/// with an access method of ftp:// cannot be expressed a a relative URL to
/// a local disk file.
pub(crate) fn fits_url2relurl(
    refURL: &[c_char],     /* I reference URL string             */
    absURL: &[c_char],     /* I absoulute URL string to process  */
    relURL: &mut [c_char], /* O resulting relative URL string    */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Create an absolute URL from a relative url and a reference URL. The
/// reference URL is given by the FITS file pointed to by fptr.
///
/// The construction of the absolute URL from the partial and reference URl
/// is performed using the rules set forth in:
///
/// http://www.w3.org/Addressing/URL/URL_TOC.html
/// and
/// http://www.w3.org/Addressing/URL/4_3_Partial.html
///
/// Note that the relative URL string relURL must conform to the Unix-like
/// URL syntax; host dependent partial URL strings are not allowed.
pub fn fits_relurl2url(
    cwd: [c_char; FLEN_FILENAME], /* I reference URL string             */
    filename: *mut c_char,        /* I relative URL string to process   */
    abs_url: [c_char; FLEN_FILENAME], /* O absolute URL string              */
    status: &mut c_int,
) -> c_int {
    todo!()
}

/*--------------------------------------------------------------------------*/
/// Encode all URL "unsafe" and "reserved" characters using the "%XX"
/// convention, where XX stand for the two hexidecimal digits of the
/// encode character's ASCII code.
///
/// Note that the outpath length, as specified by the maxlength argument,
/// should be at least as large as inpath and preferably larger (to hold
/// any characters that need encoding).  If more than maxlength chars are
/// required for outpath, including the terminating NULL, outpath will
/// be set to size 0 and an error status will be returned.
///
/// This function was adopted from code in the libwww.a library available
/// via the W3 consortium <URL: http://www.w3.org>
pub(crate) fn fits_encode_url(
    inpath: &[c_char],      /* I URL  to be encoded                  */
    maxlength: usize, /* I max number of chars that may be copied to outpath, including terminating NULL. */
    outpath: &mut [c_char], /* O output encoded URL                  */
    status: &mut c_int,
) -> c_int {
    let mut a: c_uchar;

    let mut p = 0;
    let mut q = 0;
    let hex: &[c_char] = cast_slice(b"0123456789ABCDEF");
    let mut iout = 0;

    let isAcceptable: [c_uchar; 96] = [
        /* 0x0 0x1 0x2 0x3 0x4 0x5 0x6 0x7 0x8 0x9 0xA 0xB 0xC 0xD 0xE 0xF */
        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0xF, 0xE, 0x0, 0xF, 0xF, 0xC,
        /* 2x  !"#$%&'()*+,-./   */
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0x8, 0x0, 0x0, 0x0, 0x0, 0x0,
        /* 3x 0123456789:;<=>?   */
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        /* 4x @ABCDEFGHIJKLMNO   */
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0x0, 0x0, 0x0, 0x0, 0xF,
        /* 5X PQRSTUVWXYZ[\]^_   */
        0x0, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        /* 6x `abcdefghijklmno   */
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0x0, 0x0, 0x0, 0x0,
        0x0, /* 7X pqrstuvwxyz{\}~DEL */
    ];

    if *status != 0 {
        return *status;
    }

    /* loop over all characters in inpath until 0 is encountered */

    while inpath[p] > 0 && (iout < maxlength - 1) {
        a = inpath[p] as u8;

        /* if the charcter requires encoding then process it */

        if !(a >= 32 && a < 128 && (isAcceptable[a as usize - 32] > 0)) {
            if iout + 2 < maxlength - 1 {
                /* add a '%' character to the outpath */
                outpath[q] = HEX_ESCAPE as c_char;
                q += 1;

                /* add the most significant ASCII code hex value */
                outpath[q] = hex[(a >> 4) as usize];
                q += 1;

                /* add the least significant ASCII code hex value */
                outpath[q] = hex[(a & 15) as usize];
                q += 1;
                iout += 3;
            } else {
                ffpmsg_str("URL input is too long to encode (fits_encode_url)");
                *status = URL_PARSE_ERROR;
                outpath[0] = 0;
                return *status;
            }
        }
        /* else just copy the character as is */
        else {
            outpath[q] = inpath[p];
            q += 1;
            iout += 1;
        }
        p += 1;
    }

    /* null terminate the outpath string */

    if inpath[p] > 0 && (iout == maxlength - 1) {
        ffpmsg_str("URL input is too long to encode (fits_encode_url)");
        *status = URL_PARSE_ERROR;
        outpath[0] = 0;
        return *status;
    }

    outpath[q] = 0;
    q += 1;

    *status
}

/*---------------------------------------------------------------------------*/
/// unencode all URL "unsafe" and "reserved" characters to their actual
/// ASCII representation. All tokens of the form "%XX" where XX is the
/// hexidecimal code for an ASCII character, are searched for and
/// translated into the actuall ASCII character (so three chars become
/// 1 char).
///
/// It is assumed that OUTPATH has enough room to hold the unencoded
/// URL.
///
/// This function was adopted from code in the libwww.a library available
/// via the W3 consortium <URL: http://www.w3.org>
pub(crate) fn fits_unencode_url(
    inpath: &[c_char],      /* I input URL with encoding            */
    outpath: &mut [c_char], /* O unencoded URL                      */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*---------------------------------------------------------------------------*/
/// Return a True (1) or False (0) value indicating whether or not the passed
/// URL string contains an access method specifier or not. Note that this is
/// a boolean function and it neither reads nor returns the standard error
/// status parameter
/// Note that this doens't validate the URL, it just checks for the presence
pub fn fits_is_url_absolute(url: &[c_char]) -> c_int {
    let reserved: [c_char; 10] = [
        bb(b':'),
        bb(b';'),
        bb(b'/'),
        bb(b'?'),
        bb(b'@'),
        bb(b'&'),
        bb(b'='),
        bb(b'+'),
        bb(b'$'),
        bb(b','),
    ];

    /*
     The rule for determing if an URL is relative or absolute is that it (1)
     must have a colon ":" and (2) that the colon must appear before any other
     reserved URL character in the URL string. We first see if a colon exists,
     get its position in the string, and then check to see if any of the other
     reserved characters exists and if their position in the string is greater
     than that of the colons.
    */

    if let Some(colon_pos) = url.iter().position(|&c| c == reserved[0]) {
        if reserved[1..].iter().all(|&c| {
            url.iter()
                .position(|&x| x == c)
                .is_none_or(|pos| pos > colon_pos)
        }) {
            return 1;
        }
    }
    0
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fits_is_url_absolute() {
        assert_eq!(fits_is_url_absolute(cs!(c"http://example.com")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"ftp://example.com")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"file:///path/to/file")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"mem://memory")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"shmem://shared_memory")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"root://root_path")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"/absolute/path")), 0); // No colon
        assert_eq!(fits_is_url_absolute(cs!(c"relative/path")), 0);
        assert_eq!(fits_is_url_absolute(cs!(c"C:/absolute/path")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"")), 0);
        assert_eq!(fits_is_url_absolute(cs!(c"no_colon_here")), 0);
        assert_eq!(fits_is_url_absolute(cs!(c"http:/example.com")), 1);
        assert_eq!(fits_is_url_absolute(cs!(c"http:example.com")), 1);
    }

    #[test]
    fn test_fits_clean_url() {
        let mut status = 0;
        let mut out_url = [0 as c_char; FLEN_FILENAME];

        // Test cases
        let test_cases = vec![
            // ("/dir/../file.fits", "/file.fits", 0),
            // ("/home/user/../file.fits", "/home/file.fits", 0),
            // ("/home//user///file.fits", "/home/user/file.fits", 0),
            //("", "", 0),
            ("././.", "", 0),
            // ("dir/../file?.fits", "file?.fits", 0),
            // ("/dir/./file.fits", "/dir/file.fits", 0),
            // ("/dir/././file.fits", "/dir/file.fits", 0),
            // ("/dir/../dir2/../file.fits", "/file.fits", 0),
            // ("/dir/dir2/../../file.fits", "/file.fits", 0),
            // ("dir/./file.fits", "dir/file.fits", 0),
            // ("dir/../dir2/file.fits", "dir2/file.fits", 0),
            // ("dir/dir2/../../file.fits", "file.fits", 0),
            // ("dir/dir2/../../../file.fits", "../file.fits", 0),
            // ("/dir/dir2/../../../file.fits", "/file.fits", 0),
            // ("dir/dir2/../../dir3/./file.fits", "dir3/file.fits", 0),
            // ("dir/dir2/../../dir3/../file.fits", "file.fits", 0),
            // ("dir/dir2/../../dir3/.././file.fits", "file.fits", 0),
            // ("/dir/dir2/../../dir3/./file.fits", "/dir3/file.fits", 0),
            // ("/dir/dir2/../../dir3/../file.fits", "/file.fits", 0),
            // ("/dir/dir2/../../dir3/.././file.fits", "/file.fits", 0),
        ];

        for (input, expected, expected_status) in test_cases {
            let mut input_cstr = [0 as c_char; FLEN_FILENAME];
            for (i, &byte) in input.as_bytes().iter().enumerate() {
                input_cstr[i] = byte as c_char;
            }

            fits_clean_url(&input_cstr, &mut out_url, &mut status);

            let output = CStr::from_bytes_until_nul(cast_slice(&out_url))
                .unwrap()
                .to_str()
                .unwrap();
            assert_eq!(output, expected);
            assert_eq!(status, expected_status);
        }
    }
}
