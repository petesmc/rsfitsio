#![allow(deprecated)]
#![allow(unused_imports)]

use crate::c_types::{c_char, c_int};
use crate::fitsio::CFITSIO_SONAME;
use crate::fitsio::fitsfile;

pub(crate) use crate::cfileio::ffourl as fits_parse_output_url;

pub mod c_api {
    use super::*;

    pub use crate::cfileio::ffexist as fits_file_exists;
    pub use crate::cfileio::ffifile as fits_parse_input_filename;
    pub use crate::cfileio::ffiurl as fits_parse_input_url;
    pub use crate::cfileio::ffrtnm as fits_parse_rootname;

    pub use crate::cfileio::ffextn as fits_parse_extnum;
    pub use crate::cfileio::ffexts as fits_parse_extspec;
    pub use crate::cfileio::ffomem as fits_open_memfile;
    pub use crate::editcol::ffrwrg as fits_parse_range;
    pub use crate::editcol::ffrwrgll as fits_parse_rangell;
    pub use crate::histo::ffbinr as fits_parse_binrange;
    pub use crate::histo::ffbins as fits_parse_binspec;

    /*
    use the following special macro to test that the fitsio.h include file
    that was used to build the CFITSIO library is compatible with the version
    as included when compiling the application program
    */
    #[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
    pub unsafe extern "C" fn fits_open_file(
        A: *mut Option<Box<fitsfile>>,
        B: *const c_char,
        C: c_int,
        D: *mut c_int,
    ) -> c_int {
        unsafe {
            use crate::fitsio::NULL_MSG;
            use bytemuck::cast_slice;
            use std::ffi::CStr;

            let A = A.as_mut().expect(NULL_MSG);
            let D = D.as_mut().expect(NULL_MSG);
            let B = cast_slice(CStr::from_ptr(B).to_bytes_with_nul());

            crate::cfileio::ffopentest_safe(CFITSIO_SONAME as c_int, A, B, C, D)
        }
    }

    pub use crate::buffers::ffflsh as fits_flush_buffer;
    pub use crate::buffers::ffflus as fits_flush_file;
    pub use crate::cfileio::ffclos as fits_close_file;
    pub use crate::cfileio::ffdelt as fits_delete_file;
    pub use crate::cfileio::ffdkinit as fits_create_diskfile;
    pub use crate::cfileio::ffdkopn as fits_open_diskfile;
    pub use crate::cfileio::ffdopn as fits_open_data;
    pub use crate::cfileio::ffeopn as fits_open_extlist;
    pub use crate::cfileio::ffimem as fits_create_memfile;
    pub use crate::cfileio::ffinit as fits_create_file;
    pub use crate::cfileio::ffiopn as fits_open_image;
    pub use crate::cfileio::ffreopen as fits_reopen_file;
    pub use crate::cfileio::fftopn as fits_open_table;
    pub use crate::cfileio::fftplt as fits_create_template;
    pub use crate::cfileio::ffurlt as fits_url_type;
    pub use crate::fitscore::ffflmd as fits_file_mode;
    pub use crate::fitscore::ffflnm as fits_file_name;

    pub use crate::buffers::ffgrsz as fits_get_rowsize;
    pub use crate::cfileio::ffrprt as fits_report_error;
    pub use crate::fitscore::ffasfm as fits_ascii_tform;
    pub use crate::fitscore::ffbnfm as fits_binary_tform;
    pub use crate::fitscore::ffbnfmll as fits_binary_tformll;
    pub use crate::fitscore::ffcmps as fits_compare_str;
    pub use crate::fitscore::ffcmrk as fits_clear_errmark;
    pub use crate::fitscore::ffcmsg as fits_clear_errmsg;
    pub use crate::fitscore::ffdtyp as fits_get_keytype;
    pub use crate::fitscore::ffgabc as fits_get_tbcol;
    pub use crate::fitscore::ffgerr as fits_get_errstatus;
    pub use crate::fitscore::ffgkcl as fits_get_keyclass;
    pub use crate::fitscore::ffgmsg as fits_read_errmsg;
    pub use crate::fitscore::ffgthd as fits_parse_template;
    pub use crate::fitscore::ffinttyp as fits_get_inttype;
    pub use crate::fitscore::ffkeyn as fits_make_keyn;
    pub use crate::fitscore::ffmkky as fits_make_key;
    pub use crate::fitscore::ffnkey as fits_make_nkey;
    pub use crate::fitscore::ffpmrk as fits_write_errmark;
    pub use crate::fitscore::ffpmsg as fits_write_errmsg;
    pub use crate::fitscore::ffpsvc as fits_parse_value;
    pub use crate::fitscore::fftkey as fits_test_keyword;
    pub use crate::fitscore::fftrec as fits_test_record;
    pub use crate::fitscore::ffupch as fits_uppercase;
    pub use crate::fitscore::ffvers as fits_get_version;
    pub use crate::getcols::ffgcdw as fits_get_col_display_width;
    pub use crate::getkey::ffgknm as fits_get_keyname;
    pub use crate::getkey::ffnchk as fits_null_check;

    pub use crate::editcol::ffcpky as fits_copy_key;
    pub use crate::modkey::ffpunt as fits_write_key_unit;
    pub use crate::putkey::ffdt2s as fits_date2str;
    pub use crate::putkey::ffgsdt as fits_get_system_date;
    pub use crate::putkey::ffgstm as fits_get_system_time;
    pub use crate::putkey::ffpcom as fits_write_comment;
    pub use crate::putkey::ffpdat as fits_write_date;
    pub use crate::putkey::ffphbn as fits_write_btblhdr;
    pub use crate::putkey::ffphext as fits_write_exthdr;
    pub use crate::putkey::ffphis as fits_write_history;
    pub use crate::putkey::ffphpr as fits_write_grphdr;
    pub use crate::putkey::ffphprll as fits_write_grphdrll;
    pub use crate::putkey::ffphps as fits_write_imghdr;
    pub use crate::putkey::ffphpsll as fits_write_imghdrll;
    pub use crate::putkey::ffphtb as fits_write_atblhdr;
    pub use crate::putkey::ffpkfc as fits_write_key_fixcmp;
    pub use crate::putkey::ffpkfm as fits_write_key_fixdblcmp;
    pub use crate::putkey::ffpkls as fits_write_key_longstr;
    pub use crate::putkey::ffpknd as fits_write_keys_dbl;
    pub use crate::putkey::ffpkne as fits_write_keys_flt;
    pub use crate::putkey::ffpknf as fits_write_keys_fixflt;
    pub use crate::putkey::ffpkng as fits_write_keys_fixdbl;
    pub use crate::putkey::ffpknj as fits_write_keys_lng;
    pub use crate::putkey::ffpknl as fits_write_keys_log;
    pub use crate::putkey::ffpkns as fits_write_keys_str;
    pub use crate::putkey::ffpktp as fits_write_key_template;
    pub use crate::putkey::ffpky as fits_write_key;
    pub use crate::putkey::ffpkyc as fits_write_key_cmp;
    pub use crate::putkey::ffpkyd as fits_write_key_dbl;
    pub use crate::putkey::ffpkye as fits_write_key_flt;
    pub use crate::putkey::ffpkyf as fits_write_key_fixflt;
    pub use crate::putkey::ffpkyg as fits_write_key_fixdbl;
    pub use crate::putkey::ffpkyj as fits_write_key_lng;
    pub use crate::putkey::ffpkyl as fits_write_key_log;
    pub use crate::putkey::ffpkym as fits_write_key_dblcmp;
    pub use crate::putkey::ffpkys as fits_write_key_str;
    pub use crate::putkey::ffpkyt as fits_write_key_triple;
    pub use crate::putkey::ffpkyu as fits_write_key_null;
    pub use crate::putkey::ffpkyuj as fits_write_key_ulng;
    pub use crate::putkey::ffplsw as fits_write_key_longwarn;
    pub use crate::putkey::ffprec as fits_write_record;
    pub use crate::putkey::ffptdm as fits_write_tdim;
    pub use crate::putkey::ffptdmll as fits_write_tdimll;
    pub use crate::putkey::ffs2dt as fits_str2date;
    pub use crate::putkey::ffs2tm as fits_str2time;
    pub use crate::putkey::fftm2s as fits_time2str;

    pub use crate::getkey::ffghps as fits_get_hdrpos;
    pub use crate::getkey::ffghsp as fits_get_hdrspace;
    pub use crate::getkey::ffgnxk as fits_find_nextkey;
    pub use crate::getkey::ffmaky as fits_movabs_key;
    pub use crate::getkey::ffmrky as fits_movrel_key;

    pub use crate::getkey::ffcnvthdr2str as fits_convert_hdr2str;
    pub use crate::getkey::ffdtdm as fits_decode_tdim;
    pub use crate::getkey::ffdtdmll as fits_decode_tdimll;
    pub use crate::getkey::fffree as fits_free_memory;
    pub use crate::getkey::ffgcrd as fits_read_card;
    pub use crate::getkey::ffghbn as fits_read_btblhdr;
    pub use crate::getkey::ffghbnll as fits_read_btblhdrll;
    pub use crate::getkey::ffghpr as fits_read_imghdr;
    pub use crate::getkey::ffghprll as fits_read_imghdrll;
    pub use crate::getkey::ffghtb as fits_read_atblhdr;
    pub use crate::getkey::ffghtbll as fits_read_atblhdrll;
    pub use crate::getkey::ffgkcsl as fits_get_key_com_strlen;
    pub use crate::getkey::ffgkey as fits_read_keyword;
    pub use crate::getkey::ffgkls as fits_read_key_longstr;
    pub use crate::getkey::ffgknd as fits_read_keys_dbl;
    pub use crate::getkey::ffgkne as fits_read_keys_flt;
    pub use crate::getkey::ffgknj as fits_read_keys_lng;
    pub use crate::getkey::ffgknjj as fits_read_keys_lnglng;
    pub use crate::getkey::ffgknl as fits_read_keys_log;
    pub use crate::getkey::ffgkns as fits_read_keys_str;
    pub use crate::getkey::ffgksl as fits_get_key_strlen;
    pub use crate::getkey::ffgky as fits_read_key;
    pub use crate::getkey::ffgkyc as fits_read_key_cmp;
    pub use crate::getkey::ffgkyd as fits_read_key_dbl;
    pub use crate::getkey::ffgkye as fits_read_key_flt;
    pub use crate::getkey::ffgkyj as fits_read_key_lng;
    pub use crate::getkey::ffgkyjj as fits_read_key_lnglng;
    pub use crate::getkey::ffgkyl as fits_read_key_log;
    pub use crate::getkey::ffgkym as fits_read_key_dblcmp;
    pub use crate::getkey::ffgkyn as fits_read_keyn;
    pub use crate::getkey::ffgkys as fits_read_key_str;
    pub use crate::getkey::ffgkyt as fits_read_key_triple;
    pub use crate::getkey::ffgkyujj as fits_read_key_ulnglng;
    pub use crate::getkey::ffgrec as fits_read_record;
    pub use crate::getkey::ffgsky as fits_read_string_key;
    pub use crate::getkey::ffgskyc as fits_read_string_key_com;
    pub use crate::getkey::ffgstr as fits_read_str;
    pub use crate::getkey::ffgtdm as fits_read_tdim;
    pub use crate::getkey::ffgtdmll as fits_read_tdimll;
    pub use crate::getkey::ffgunt as fits_read_key_unit;
    pub use crate::getkey::ffhdr2str as fits_hdr2str;

    pub use crate::modkey::ffucrd as fits_update_card;
    pub use crate::modkey::ffukfc as fits_update_key_fixcmp;
    pub use crate::modkey::ffukfm as fits_update_key_fixdblcmp;
    pub use crate::modkey::ffukls as fits_update_key_longstr;
    pub use crate::modkey::ffuky as fits_update_key;
    pub use crate::modkey::ffukyc as fits_update_key_cmp;
    pub use crate::modkey::ffukyd as fits_update_key_dbl;
    pub use crate::modkey::ffukye as fits_update_key_flt;
    pub use crate::modkey::ffukyf as fits_update_key_fixflt;
    pub use crate::modkey::ffukyg as fits_update_key_fixdbl;
    pub use crate::modkey::ffukyj as fits_update_key_lng;
    pub use crate::modkey::ffukyl as fits_update_key_log;
    pub use crate::modkey::ffukym as fits_update_key_dblcmp;
    pub use crate::modkey::ffukys as fits_update_key_str;
    pub use crate::modkey::ffukyu as fits_update_key_null;
    pub use crate::modkey::ffukyuj as fits_update_key_ulng;

    pub use crate::modkey::ffmcom as fits_modify_comment;
    pub use crate::modkey::ffmcrd as fits_modify_card;
    pub use crate::modkey::ffmkfc as fits_modify_key_fixcmp;
    pub use crate::modkey::ffmkfm as fits_modify_key_fixdblcmp;
    pub use crate::modkey::ffmkls as fits_modify_key_longstr;
    pub use crate::modkey::ffmkyc as fits_modify_key_cmp;
    pub use crate::modkey::ffmkyd as fits_modify_key_dbl;
    pub use crate::modkey::ffmkye as fits_modify_key_flt;
    pub use crate::modkey::ffmkyf as fits_modify_key_fixflt;
    pub use crate::modkey::ffmkyg as fits_modify_key_fixdbl;
    pub use crate::modkey::ffmkyj as fits_modify_key_lng;
    pub use crate::modkey::ffmkyl as fits_modify_key_log;
    pub use crate::modkey::ffmkym as fits_modify_key_dblcmp;
    pub use crate::modkey::ffmkys as fits_modify_key_str;
    pub use crate::modkey::ffmkyu as fits_modify_key_null;
    pub use crate::modkey::ffmkyuj as fits_modify_key_ulng;
    pub use crate::modkey::ffmnam as fits_modify_name;
    pub use crate::modkey::ffmrec as fits_modify_record;

    pub use crate::modkey::ffikey as fits_insert_card;
    pub use crate::modkey::ffikfc as fits_insert_key_fixcmp;
    pub use crate::modkey::ffikfm as fits_insert_key_fixdblcmp;
    pub use crate::modkey::ffikls as fits_insert_key_longstr;
    pub use crate::modkey::ffikyc as fits_insert_key_cmp;
    pub use crate::modkey::ffikyd as fits_insert_key_dbl;
    pub use crate::modkey::ffikye as fits_insert_key_flt;
    pub use crate::modkey::ffikyf as fits_insert_key_fixflt;
    pub use crate::modkey::ffikyg as fits_insert_key_fixdbl;
    pub use crate::modkey::ffikyj as fits_insert_key_lng;
    pub use crate::modkey::ffikyl as fits_insert_key_log;
    pub use crate::modkey::ffikym as fits_insert_key_dblcmp;
    pub use crate::modkey::ffikys as fits_insert_key_str;
    pub use crate::modkey::ffikyu as fits_insert_key_null;
    pub use crate::modkey::ffirec as fits_insert_record;

    pub use crate::fitscore::ffghad as fits_get_hduaddr;
    pub use crate::fitscore::ffghadll as fits_get_hduaddrll;
    pub use crate::fitscore::ffghdn as fits_get_hdu_num;
    pub use crate::fitscore::ffghdt as fits_get_hdu_type;
    pub use crate::fitscore::ffghof as fits_get_hduoff;
    pub use crate::modkey::ffdkey as fits_delete_key;
    pub use crate::modkey::ffdrec as fits_delete_record;
    pub use crate::modkey::ffdstr as fits_delete_str;

    pub use crate::fitscore::ffgipr as fits_get_img_param;
    pub use crate::fitscore::ffgiprll as fits_get_img_paramll;

    pub use crate::fitscore::ffgidm as fits_get_img_dim;
    pub use crate::fitscore::ffgidt as fits_get_img_type;
    pub use crate::fitscore::ffgiet as fits_get_img_equivtype;
    pub use crate::fitscore::ffgisz as fits_get_img_size;
    pub use crate::fitscore::ffgiszll as fits_get_img_sizell;

    pub use crate::editcol::ffrsim as fits_resize_img;
    pub use crate::editcol::ffrsimll as fits_resize_imgll;
    pub use crate::edithdu::ffibin as fits_insert_btbl;
    pub use crate::edithdu::ffiimg as fits_insert_img;
    pub use crate::edithdu::ffiimgll as fits_insert_imgll;
    pub use crate::edithdu::ffitab as fits_insert_atbl;
    pub use crate::fitscore::ffcrhd as fits_create_hdu;
    pub use crate::fitscore::ffmahd as fits_movabs_hdu;
    pub use crate::fitscore::ffmnhd as fits_movnam_hdu;
    pub use crate::fitscore::ffmrhd as fits_movrel_hdu;
    pub use crate::fitscore::ffthdu as fits_get_num_hdus;
    pub use crate::putkey::ffcrim as fits_create_img;
    pub use crate::putkey::ffcrimll as fits_create_imgll;
    pub use crate::putkey::ffcrtb as fits_create_tbl;

    pub use crate::edithdu::ffcopy as fits_copy_hdu;
    pub use crate::edithdu::ffcpdt as fits_copy_data;
    pub use crate::edithdu::ffcpfl as fits_copy_file;
    pub use crate::edithdu::ffcphd as fits_copy_header;
    pub use crate::edithdu::ffcpht as fits_copy_hdutab;
    pub use crate::edithdu::ffdhdu as fits_delete_hdu;
    pub use crate::edithdu::ffwrhdu as fits_write_hdu;

    pub use crate::fitscore::ffhdef as fits_set_hdrsize;
    pub use crate::fitscore::ffrdef as fits_set_hdustruc;
    pub use crate::scalnull::ffpthp as fits_write_theap;

    pub use crate::checksum::ffdsum as fits_decode_chksum;
    pub use crate::checksum::ffesum as fits_encode_chksum;
    pub use crate::checksum::ffgcks as fits_get_chksum;
    pub use crate::checksum::ffpcks as fits_write_chksum;
    pub use crate::checksum::ffupck as fits_update_chksum;
    pub use crate::checksum::ffvcks as fits_verify_chksum;

    pub use crate::scalnull::ffpnul as fits_set_imgnull;
    pub use crate::scalnull::ffpscl as fits_set_bscale;
    pub use crate::scalnull::ffsnul as fits_set_atblnull;
    pub use crate::scalnull::fftnul as fits_set_btblnull;
    pub use crate::scalnull::fftscl as fits_set_tscale;

    pub use crate::fitscore::ffeqty as fits_get_eqcoltype;
    pub use crate::fitscore::ffeqtyll as fits_get_eqcoltypell;
    pub use crate::fitscore::ffgacl as fits_get_acolparms;
    pub use crate::fitscore::ffgbcl as fits_get_bcolparms;
    pub use crate::fitscore::ffgbclll as fits_get_bcolparmsll;
    pub use crate::fitscore::ffgcnn as fits_get_colname;
    pub use crate::fitscore::ffgcno as fits_get_colnum;
    pub use crate::fitscore::ffgncl as fits_get_num_cols;
    pub use crate::fitscore::ffgnrw as fits_get_num_rows;
    pub use crate::fitscore::ffgnrwll as fits_get_num_rowsll;
    pub use crate::fitscore::ffgtcl as fits_get_coltype;
    pub use crate::fitscore::ffgtclll as fits_get_coltypell;

    pub use crate::putcol::ffiter as fits_iterate_data;

    pub use crate::getcolb::ffggpb as fits_read_grppar_byt;
    pub use crate::getcold::ffggpd as fits_read_grppar_dbl;
    pub use crate::getcole::ffggpe as fits_read_grppar_flt;
    pub use crate::getcoli::ffggpi as fits_read_grppar_sht;
    pub use crate::getcolj::ffggpj as fits_read_grppar_lng;
    pub use crate::getcolj::ffggpjj as fits_read_grppar_lnglng;
    pub use crate::getcolk::ffggpk as fits_read_grppar_int;
    pub use crate::getcolsb::ffggpsb as fits_read_grppar_sbyt;
    pub use crate::getcolui::ffggpui as fits_read_grppar_usht;
    pub use crate::getcoluj::ffggpuj as fits_read_grppar_ulng;
    pub use crate::getcoluj::ffggpujj as fits_read_grppar_ulnglng;
    pub use crate::getcoluk::ffggpuk as fits_read_grppar_uint;

    pub use crate::getcol::ffgpf as fits_read_imgnull;
    pub use crate::getcol::ffgpv as fits_read_img;
    pub use crate::getcol::ffgpxf as fits_read_pixnull;
    pub use crate::getcol::ffgpxfll as fits_read_pixnullll;
    pub use crate::getcol::ffgpxv as fits_read_pix;
    pub use crate::getcol::ffgpxvll as fits_read_pixll;
    pub use crate::getcolb::ffgpvb as fits_read_img_byt;
    pub use crate::getcold::ffgpvd as fits_read_img_dbl;
    pub use crate::getcole::ffgpve as fits_read_img_flt;
    pub use crate::getcoli::ffgpvi as fits_read_img_sht;
    pub use crate::getcolj::ffgpvj as fits_read_img_lng;
    pub use crate::getcolj::ffgpvjj as fits_read_img_lnglng;
    pub use crate::getcolk::ffgpvk as fits_read_img_int;
    pub use crate::getcolsb::ffgpvsb as fits_read_img_sbyt;
    pub use crate::getcolui::ffgpvui as fits_read_img_usht;
    pub use crate::getcoluj::ffgpvuj as fits_read_img_ulng;
    pub use crate::getcoluj::ffgpvujj as fits_read_img_ulnglng;
    pub use crate::getcoluk::ffgpvuk as fits_read_img_uint;

    pub use crate::getcolb::ffgpfb as fits_read_imgnull_byt;
    pub use crate::getcold::ffgpfd as fits_read_imgnull_dbl;
    pub use crate::getcole::ffgpfe as fits_read_imgnull_flt;
    pub use crate::getcoli::ffgpfi as fits_read_imgnull_sht;
    pub use crate::getcolj::ffgpfj as fits_read_imgnull_lng;
    pub use crate::getcolj::ffgpfjj as fits_read_imgnull_lnglng;
    pub use crate::getcolk::ffgpfk as fits_read_imgnull_int;
    pub use crate::getcolsb::ffgpfsb as fits_read_imgnull_sbyt;
    pub use crate::getcolui::ffgpfui as fits_read_imgnull_usht;
    pub use crate::getcoluj::ffgpfuj as fits_read_imgnull_ulng;
    pub use crate::getcoluj::ffgpfujj as fits_read_imgnull_ulnglng;
    pub use crate::getcoluk::ffgpfuk as fits_read_imgnull_uint;

    pub use crate::getcolb::ffg2db as fits_read_2d_byt;
    pub use crate::getcold::ffg2dd as fits_read_2d_dbl;
    pub use crate::getcole::ffg2de as fits_read_2d_flt;
    pub use crate::getcoli::ffg2di as fits_read_2d_sht;
    pub use crate::getcolj::ffg2dj as fits_read_2d_lng;
    pub use crate::getcolj::ffg2djj as fits_read_2d_lnglng;
    pub use crate::getcolk::ffg2dk as fits_read_2d_int;
    pub use crate::getcolsb::ffg2dsb as fits_read_2d_sbyt;
    pub use crate::getcolui::ffg2dui as fits_read_2d_usht;
    pub use crate::getcoluj::ffg2duj as fits_read_2d_ulng;
    pub use crate::getcoluj::ffg2dujj as fits_read_2d_ulnglng;
    pub use crate::getcoluk::ffg2duk as fits_read_2d_uint;

    pub use crate::getcolb::ffg3db as fits_read_3d_byt;
    pub use crate::getcold::ffg3dd as fits_read_3d_dbl;
    pub use crate::getcole::ffg3de as fits_read_3d_flt;
    pub use crate::getcoli::ffg3di as fits_read_3d_sht;
    pub use crate::getcolj::ffg3dj as fits_read_3d_lng;
    pub use crate::getcolj::ffg3djj as fits_read_3d_lnglng;
    pub use crate::getcolk::ffg3dk as fits_read_3d_int;
    pub use crate::getcolsb::ffg3dsb as fits_read_3d_sbyt;
    pub use crate::getcolui::ffg3dui as fits_read_3d_usht;
    pub use crate::getcoluj::ffg3duj as fits_read_3d_ulng;
    pub use crate::getcoluj::ffg3dujj as fits_read_3d_ulnglng;
    pub use crate::getcoluk::ffg3duk as fits_read_3d_uint;

    pub use crate::getcol::ffgsv as fits_read_subset;
    pub use crate::getcolb::ffgsvb as fits_read_subset_byt;
    pub use crate::getcold::ffgsvd as fits_read_subset_dbl;
    pub use crate::getcole::ffgsve as fits_read_subset_flt;
    pub use crate::getcoli::ffgsvi as fits_read_subset_sht;
    pub use crate::getcolj::ffgsvj as fits_read_subset_lng;
    pub use crate::getcolj::ffgsvjj as fits_read_subset_lnglng;
    pub use crate::getcolk::ffgsvk as fits_read_subset_int;
    pub use crate::getcolsb::ffgsvsb as fits_read_subset_sbyt;
    pub use crate::getcolui::ffgsvui as fits_read_subset_usht;
    pub use crate::getcoluj::ffgsvuj as fits_read_subset_ulng;
    pub use crate::getcoluj::ffgsvujj as fits_read_subset_ulnglng;
    pub use crate::getcoluk::ffgsvuk as fits_read_subset_uint;

    pub use crate::getcolb::ffgsfb as fits_read_subsetnull_byt;
    pub use crate::getcold::ffgsfd as fits_read_subsetnull_dbl;
    pub use crate::getcole::ffgsfe as fits_read_subsetnull_flt;
    pub use crate::getcoli::ffgsfi as fits_read_subsetnull_sht;
    pub use crate::getcolj::ffgsfj as fits_read_subsetnull_lng;
    pub use crate::getcolj::ffgsfjj as fits_read_subsetnull_lnglng;
    pub use crate::getcolk::ffgsfk as fits_read_subsetnull_int;
    pub use crate::getcolsb::ffgsfsb as fits_read_subsetnull_sbyt;
    pub use crate::getcolui::ffgsfui as fits_read_subsetnull_usht;
    pub use crate::getcoluj::ffgsfuj as fits_read_subsetnull_ulng;
    pub use crate::getcoluj::ffgsfujj as fits_read_subsetnull_ulnglng;
    pub use crate::getcoluk::ffgsfuk as fits_read_subsetnull_uint;

    pub use crate::cfileio::fits_copy_image_section as ffcpimg;
    // Not part of external API
    // pub use 	fits_comp_img as fits_compress_img;
    // pub use 	fits_decomp_img as fits_decompress_img;

    pub use crate::getcol::ffgcf as fits_read_colnull;
    pub use crate::getcol::ffgcv as fits_read_col;
    pub use crate::getcol::ffgcvn as fits_read_cols;
    pub use crate::getcolb::ffgcvb as fits_read_col_byt;
    pub use crate::getcold::ffgcvd as fits_read_col_dbl;
    pub use crate::getcold::ffgcvm as fits_read_col_dblcmp;
    pub use crate::getcole::ffgcvc as fits_read_col_cmp;
    pub use crate::getcole::ffgcve as fits_read_col_flt;
    pub use crate::getcoli::ffgcvi as fits_read_col_sht;
    pub use crate::getcolj::ffgcvj as fits_read_col_lng;
    pub use crate::getcolj::ffgcvjj as fits_read_col_lnglng;
    pub use crate::getcolk::ffgcvk as fits_read_col_int;
    pub use crate::getcoll::ffgcvl as fits_read_col_log;
    pub use crate::getcoll::ffgcx as fits_read_col_bit;
    pub use crate::getcoll::ffgcxui as fits_read_col_bit_usht;
    pub use crate::getcoll::ffgcxuk as fits_read_col_bit_uint;
    pub use crate::getcols::ffgcvs as fits_read_col_str;
    pub use crate::getcolsb::ffgcvsb as fits_read_col_sbyt;
    pub use crate::getcolui::ffgcvui as fits_read_col_usht;
    pub use crate::getcoluj::ffgcvuj as fits_read_col_ulng;
    pub use crate::getcoluj::ffgcvujj as fits_read_col_ulnglng;
    pub use crate::getcoluk::ffgcvuk as fits_read_col_uint;

    pub use crate::getcolb::ffgcfb as fits_read_colnull_byt;
    pub use crate::getcold::ffgcfd as fits_read_colnull_dbl;
    pub use crate::getcold::ffgcfm as fits_read_colnull_dblcmp;
    pub use crate::getcole::ffgcfc as fits_read_colnull_cmp;
    pub use crate::getcole::ffgcfe as fits_read_colnull_flt;
    pub use crate::getcoli::ffgcfi as fits_read_colnull_sht;
    pub use crate::getcolj::ffgcfj as fits_read_colnull_lng;
    pub use crate::getcolj::ffgcfjj as fits_read_colnull_lnglng;
    pub use crate::getcolk::ffgcfk as fits_read_colnull_int;
    pub use crate::getcoll::ffgcfl as fits_read_colnull_log;
    pub use crate::getcols::ffgcfs as fits_read_colnull_str;
    pub use crate::getcolsb::ffgcfsb as fits_read_colnull_sbyt;
    pub use crate::getcolui::ffgcfui as fits_read_colnull_usht;
    pub use crate::getcoluj::ffgcfuj as fits_read_colnull_ulng;
    pub use crate::getcoluj::ffgcfujj as fits_read_colnull_ulnglng;
    pub use crate::getcoluk::ffgcfuk as fits_read_colnull_uint;

    pub use crate::buffers::ffgtbb as fits_read_tblbytes;
    pub use crate::fitscore::ffgdes as fits_read_descript;
    pub use crate::fitscore::ffgdesll as fits_read_descriptll;
    pub use crate::fitscore::ffgdess as fits_read_descripts;
    pub use crate::fitscore::ffgdessll as fits_read_descriptsll;

    pub use crate::putcolb::ffpgpb as fits_write_grppar_byt;
    pub use crate::putcold::ffpgpd as fits_write_grppar_dbl;
    pub use crate::putcole::ffpgpe as fits_write_grppar_flt;
    pub use crate::putcoli::ffpgpi as fits_write_grppar_sht;
    pub use crate::putcolj::ffpgpj as fits_write_grppar_lng;
    pub use crate::putcolj::ffpgpjj as fits_write_grppar_lnglng;
    pub use crate::putcolk::ffpgpk as fits_write_grppar_int;
    pub use crate::putcolsb::ffpgpsb as fits_write_grppar_sbyt;
    pub use crate::putcolui::ffpgpui as fits_write_grppar_usht;
    pub use crate::putcoluj::ffpgpuj as fits_write_grppar_ulng;
    pub use crate::putcoluj::ffpgpujj as fits_write_grppar_ulnglng;
    pub use crate::putcoluk::ffpgpuk as fits_write_grppar_uint;

    pub use crate::putcol::ffppr as fits_write_img;
    pub use crate::putcol::ffppx as fits_write_pix;
    pub use crate::putcol::ffppxll as fits_write_pixll;
    pub use crate::putcol::ffppxn as fits_write_pixnull;
    pub use crate::putcol::ffppxnll as fits_write_pixnullll;
    pub use crate::putcolb::ffpprb as fits_write_img_byt;
    pub use crate::putcold::ffpprd as fits_write_img_dbl;
    pub use crate::putcole::ffppre as fits_write_img_flt;
    pub use crate::putcoli::ffppri as fits_write_img_sht;
    pub use crate::putcolj::ffpprj as fits_write_img_lng;
    pub use crate::putcolj::ffpprjj as fits_write_img_lnglng;
    pub use crate::putcolk::ffpprk as fits_write_img_int;
    pub use crate::putcolsb::ffpprsb as fits_write_img_sbyt;
    pub use crate::putcolui::ffpprui as fits_write_img_usht;
    pub use crate::putcoluj::ffppruj as fits_write_img_ulng;
    pub use crate::putcoluj::ffpprujj as fits_write_img_ulnglng;
    pub use crate::putcoluk::ffppruk as fits_write_img_uint;

    pub use crate::putcol::ffppn as fits_write_imgnull;
    pub use crate::putcolb::ffppnb as fits_write_imgnull_byt;
    pub use crate::putcold::ffppnd as fits_write_imgnull_dbl;
    pub use crate::putcole::ffppne as fits_write_imgnull_flt;
    pub use crate::putcoli::ffppni as fits_write_imgnull_sht;
    pub use crate::putcolj::ffppnj as fits_write_imgnull_lng;
    pub use crate::putcolj::ffppnjj as fits_write_imgnull_lnglng;
    pub use crate::putcolk::ffppnk as fits_write_imgnull_int;
    pub use crate::putcolsb::ffppnsb as fits_write_imgnull_sbyt;
    pub use crate::putcolui::ffppnui as fits_write_imgnull_usht;
    pub use crate::putcoluj::ffppnuj as fits_write_imgnull_ulng;
    pub use crate::putcoluj::ffppnujj as fits_write_imgnull_ulnglng;
    pub use crate::putcoluk::ffppnuk as fits_write_imgnull_uint;

    pub use crate::putcolu::ffpprn as fits_write_null_img;
    pub use crate::putcolu::ffppru as fits_write_img_null;

    pub use crate::putcolb::ffp2db as fits_write_2d_byt;
    pub use crate::putcold::ffp2dd as fits_write_2d_dbl;
    pub use crate::putcole::ffp2de as fits_write_2d_flt;
    pub use crate::putcoli::ffp2di as fits_write_2d_sht;
    pub use crate::putcolj::ffp2dj as fits_write_2d_lng;
    pub use crate::putcolj::ffp2djj as fits_write_2d_lnglng;
    pub use crate::putcolk::ffp2dk as fits_write_2d_int;
    pub use crate::putcolsb::ffp2dsb as fits_write_2d_sbyt;
    pub use crate::putcolui::ffp2dui as fits_write_2d_usht;
    pub use crate::putcoluj::ffp2duj as fits_write_2d_ulng;
    pub use crate::putcoluj::ffp2dujj as fits_write_2d_ulnglng;
    pub use crate::putcoluk::ffp2duk as fits_write_2d_uint;

    pub use crate::putcolb::ffp3db as fits_write_3d_byt;
    pub use crate::putcold::ffp3dd as fits_write_3d_dbl;
    pub use crate::putcole::ffp3de as fits_write_3d_flt;
    pub use crate::putcoli::ffp3di as fits_write_3d_sht;
    pub use crate::putcolj::ffp3dj as fits_write_3d_lng;
    pub use crate::putcolj::ffp3djj as fits_write_3d_lnglng;
    pub use crate::putcolk::ffp3dk as fits_write_3d_int;
    pub use crate::putcolsb::ffp3dsb as fits_write_3d_sbyt;
    pub use crate::putcolui::ffp3dui as fits_write_3d_usht;
    pub use crate::putcoluj::ffp3duj as fits_write_3d_ulng;
    pub use crate::putcoluj::ffp3dujj as fits_write_3d_ulnglng;
    pub use crate::putcoluk::ffp3duk as fits_write_3d_uint;

    pub use crate::putcol::ffpss as fits_write_subset;
    pub use crate::putcolb::ffpssb as fits_write_subset_byt;
    pub use crate::putcold::ffpssd as fits_write_subset_dbl;
    pub use crate::putcole::ffpsse as fits_write_subset_flt;
    pub use crate::putcoli::ffpssi as fits_write_subset_sht;
    pub use crate::putcolj::ffpssj as fits_write_subset_lng;
    pub use crate::putcolj::ffpssjj as fits_write_subset_lnglng;
    pub use crate::putcolk::ffpssk as fits_write_subset_int;
    pub use crate::putcolsb::ffpsssb as fits_write_subset_sbyt;
    pub use crate::putcolui::ffpssui as fits_write_subset_usht;
    pub use crate::putcoluj::ffpssuj as fits_write_subset_ulng;
    pub use crate::putcoluj::ffpssujj as fits_write_subset_ulnglng;
    pub use crate::putcoluk::ffpssuk as fits_write_subset_uint;

    pub use crate::putcol::ffpcl as fits_write_col;
    pub use crate::putcol::ffpcln as fits_write_cols;
    pub use crate::putcolb::ffpclb as fits_write_col_byt;
    pub use crate::putcold::ffpcld as fits_write_col_dbl;
    pub use crate::putcold::ffpclm as fits_write_col_dblcmp;
    pub use crate::putcole::ffpclc as fits_write_col_cmp;
    pub use crate::putcole::ffpcle as fits_write_col_flt;
    pub use crate::putcoli::ffpcli as fits_write_col_sht;
    pub use crate::putcolj::ffpclj as fits_write_col_lng;
    pub use crate::putcolj::ffpcljj as fits_write_col_lnglng;
    pub use crate::putcolk::ffpclk as fits_write_col_int;
    pub use crate::putcoll::ffpcll as fits_write_col_log;
    pub use crate::putcoll::ffpclx as fits_write_col_bit;
    pub use crate::putcols::ffpcls as fits_write_col_str;
    pub use crate::putcolsb::ffpclsb as fits_write_col_sbyt;
    pub use crate::putcolu::ffpclu as fits_write_col_null;
    pub use crate::putcolu::ffprwu as fits_write_nulrows;
    pub use crate::putcolu::ffprwu as fits_write_nullrows;
    pub use crate::putcolui::ffpclui as fits_write_col_usht;
    pub use crate::putcoluj::ffpcluj as fits_write_col_ulng;
    pub use crate::putcoluj::ffpclujj as fits_write_col_ulnglng;
    pub use crate::putcoluk::ffpcluk as fits_write_col_uint;

    pub use crate::putcol::ffpcn as fits_write_colnull;
    pub use crate::putcolb::ffpcnb as fits_write_colnull_byt;
    pub use crate::putcold::ffpcnd as fits_write_colnull_dbl;
    pub use crate::putcole::ffpcne as fits_write_colnull_flt;
    pub use crate::putcoli::ffpcni as fits_write_colnull_sht;
    pub use crate::putcolj::ffpcnj as fits_write_colnull_lng;
    pub use crate::putcolj::ffpcnjj as fits_write_colnull_lnglng;
    pub use crate::putcolk::ffpcnk as fits_write_colnull_int;
    pub use crate::putcoll::ffpcnl as fits_write_colnull_log;
    pub use crate::putcols::ffpcns as fits_write_colnull_str;
    pub use crate::putcolsb::ffpcnsb as fits_write_colnull_sbyt;
    pub use crate::putcolui::ffpcnui as fits_write_colnull_usht;
    pub use crate::putcoluj::ffpcnuj as fits_write_colnull_ulng;
    pub use crate::putcoluj::ffpcnujj as fits_write_colnull_ulnglng;
    pub use crate::putcoluk::ffpcnuk as fits_write_colnull_uint;

    pub use crate::getcolb::ffgextn as fits_read_ext;
    pub use crate::putcolb::ffpextn as fits_write_ext;

    pub use crate::fitscore::ffcmph as fits_compress_heap;
    pub use crate::fitscore::ffpdes as fits_write_descript;
    pub use crate::fitscore::fftheap as fits_test_heap;

    pub use crate::buffers::ffptbb as fits_write_tblbytes;
    pub use crate::editcol::ffccls as fits_copy_cols;
    pub use crate::editcol::ffcpcl as fits_copy_col;
    pub use crate::editcol::ffcprw as fits_copy_rows;
    pub use crate::editcol::ffcpsr as fits_copy_selrows;
    pub use crate::editcol::ffdcol as fits_delete_col;
    pub use crate::editcol::ffdrow as fits_delete_rows;
    pub use crate::editcol::ffdrrg as fits_delete_rowrange;
    pub use crate::editcol::ffdrws as fits_delete_rowlist;
    pub use crate::editcol::ffdrwsll as fits_delete_rowlistll;
    pub use crate::editcol::fficls as fits_insert_cols;
    pub use crate::editcol::fficol as fits_insert_col;
    pub use crate::editcol::ffirow as fits_insert_rows;
    pub use crate::editcol::ffmvec as fits_modify_vector_len;

    pub use crate::wcssub::ffgics as fits_read_img_coord;
    pub use crate::wcssub::ffgicsa as fits_read_img_coord_version;
    pub use crate::wcssub::ffgtcs as fits_read_tbl_coord;
    pub use crate::wcsutil::ffwldp as fits_pix_to_world;
    pub use crate::wcsutil::ffxypx as fits_world_to_pix;

    pub use crate::wcssub::ffgiwcs as fits_get_image_wcs_keys;
    pub use crate::wcssub::ffgtwcs as fits_get_table_wcs_keys;

    pub use crate::eval_f::ffcalc as fits_calculator;
    pub use crate::eval_f::ffcalc_rng as fits_calculator_rng;
    pub use crate::eval_f::ffcrow as fits_calc_rows;
    pub use crate::eval_f::ffffrw as fits_find_first_row;
    pub use crate::eval_f::fffrow as fits_find_rows;
    pub use crate::eval_f::fffrwc as fits_find_rows_cmp;
    pub use crate::eval_f::ffsrow as fits_select_rows;
    pub use crate::eval_f::fftexp as fits_test_expr;

    pub use crate::group::ffgtam as fits_add_group_member;
    pub use crate::group::ffgtch as fits_change_group;
    pub use crate::group::ffgtcm as fits_compact_group;
    pub use crate::group::ffgtcp as fits_copy_group;
    pub use crate::group::ffgtcr as fits_create_group;
    pub use crate::group::ffgtis as fits_insert_group;
    pub use crate::group::ffgtmg as fits_merge_groups;
    pub use crate::group::ffgtnm as fits_get_num_members;
    pub use crate::group::ffgtop as fits_open_group;
    pub use crate::group::ffgtrm as fits_remove_group;
    pub use crate::group::ffgtvf as fits_verify_group;

    pub use crate::group::ffgmcp as fits_copy_member;
    pub use crate::group::ffgmng as fits_get_num_groups;
    pub use crate::group::ffgmop as fits_open_member;
    pub use crate::group::ffgmrm as fits_remove_member;
    pub use crate::group::ffgmtf as fits_transfer_member;

    pub use crate::cfileio::ffgtmo as fits_get_timeout;
    pub use crate::cfileio::ffshdwn as fits_show_download_progress;
    pub use crate::cfileio::ffstmo as fits_set_timeout;

    pub use crate::cfileio::ffchtps as fits_cleanup_https;
    pub use crate::cfileio::ffihtps as fits_init_https;
    pub use crate::cfileio::ffvhtps as fits_verbose_https;
}

#[deny(unsafe_code, unsafe_op_in_unsafe_fn, deprecated)]
pub mod rust_api {
    use super::*;

    pub use crate::cfileio::ffexist_safer as fits_file_exists;
    pub use crate::cfileio::ffifile_safer as fits_parse_input_filename;
    pub use crate::cfileio::ffiurl_safer as fits_parse_input_url;
    pub use crate::cfileio::ffrtnm_safe as fits_parse_rootname;

    pub use crate::cfileio::ffextn_safer as fits_parse_extnum;
    pub use crate::cfileio::ffexts_safer as fits_parse_extspec;
    pub use crate::cfileio::ffomem_safer as fits_open_memfile;
    pub use crate::editcol::ffrwrg_safe as fits_parse_range;
    pub use crate::editcol::ffrwrgll_safe as fits_parse_rangell;
    pub use crate::histo::ffbinr_safer as fits_parse_binrange;
    pub use crate::histo::ffbins_safe as fits_parse_binspec;

    /*
    use the following special macro to test that the fitsio.h include file
    that was used to build the CFITSIO library is compatible with the version
    as included when compiling the application program
    */
    pub fn fits_open_file(
        A: &mut Option<Box<fitsfile>>,
        B: &[c_char],
        C: c_int,
        D: &mut c_int,
    ) -> c_int {
        crate::cfileio::ffopentest_safe(CFITSIO_SONAME as c_int, A, B, C, D)
    }

    pub use crate::buffers::ffflsh_safe as fits_flush_buffer;
    pub use crate::buffers::ffflus_safer as fits_flush_file;
    pub use crate::cfileio::ffclos_safer as fits_close_file;
    pub use crate::cfileio::ffdelt_safer as fits_delete_file;
    pub use crate::cfileio::ffdkinit_safer as fits_create_diskfile;
    pub use crate::cfileio::ffdkopn_safer as fits_open_diskfile;
    pub use crate::cfileio::ffdopn_safer as fits_open_data;
    pub use crate::cfileio::ffeopn_safer as fits_open_extlist;
    pub use crate::cfileio::ffimem_safer as fits_create_memfile;
    pub use crate::cfileio::ffinit_safer as fits_create_file;
    pub use crate::cfileio::ffiopn_safer as fits_open_image;
    pub use crate::cfileio::ffreopen_safer as fits_reopen_file;
    pub use crate::cfileio::fftopn_safer as fits_open_table;
    pub use crate::cfileio::fftplt_safer as fits_create_template;
    pub use crate::cfileio::ffurlt_safe as fits_url_type;
    pub use crate::fitscore::ffflmd_safe as fits_file_mode;
    pub use crate::fitscore::ffflnm_safe as fits_file_name;

    pub use crate::buffers::ffgrsz_safe as fits_get_rowsize;
    pub use crate::cfileio::ffrprt_safer as fits_report_error;
    pub use crate::fitscore::ffasfm_safe as fits_ascii_tform;
    pub use crate::fitscore::ffbnfm_safe as fits_binary_tform;
    pub use crate::fitscore::ffbnfmll_safe as fits_binary_tformll;
    pub use crate::fitscore::ffcmps_safe as fits_compare_str;
    pub use crate::fitscore::ffcmrk_safe as fits_clear_errmark;
    pub use crate::fitscore::ffcmsg_safe as fits_clear_errmsg;
    pub use crate::fitscore::ffdtyp_safe as fits_get_keytype;
    pub use crate::fitscore::ffgabc_safe as fits_get_tbcol;
    pub use crate::fitscore::ffgerr_safe as fits_get_errstatus;
    pub use crate::fitscore::ffgkcl_safe as fits_get_keyclass;
    pub use crate::fitscore::ffgmsg_safe as fits_read_errmsg;
    pub use crate::fitscore::ffgthd_safe as fits_parse_template;
    pub use crate::fitscore::ffinttyp_safe as fits_get_inttype;
    pub use crate::fitscore::ffkeyn_safe as fits_make_keyn;
    pub use crate::fitscore::ffmkky_safe as fits_make_key;
    pub use crate::fitscore::ffnkey_safe as fits_make_nkey;
    pub use crate::fitscore::ffpmrk_safe as fits_write_errmark;
    pub use crate::fitscore::ffpmsg_safer as fits_write_errmsg;
    pub use crate::fitscore::ffpsvc_safe as fits_parse_value;
    pub use crate::fitscore::fftkey_safe as fits_test_keyword;
    pub use crate::fitscore::fftrec_safe as fits_test_record;
    pub use crate::fitscore::ffupch_safe as fits_uppercase;
    pub use crate::fitscore::ffvers_safe as fits_get_version;
    pub use crate::getcols::ffgcdw_safe as fits_get_col_display_width;
    pub use crate::getkey::ffgknm_safe as fits_get_keyname;
    pub use crate::getkey::ffnchk_safe as fits_null_check;

    pub use crate::editcol::ffcpky_safe as fits_copy_key;
    pub use crate::modkey::ffpunt_safe as fits_write_key_unit;
    pub use crate::putkey::ffdt2s_safe as fits_date2str;
    pub use crate::putkey::ffgsdt_safe as fits_get_system_date;
    pub use crate::putkey::ffgstm_safe as fits_get_system_time;
    pub use crate::putkey::ffpcom_safe as fits_write_comment;
    pub use crate::putkey::ffpdat_safe as fits_write_date;
    pub use crate::putkey::ffphbn_safe as fits_write_btblhdr;
    pub use crate::putkey::ffphext_safe as fits_write_exthdr;
    pub use crate::putkey::ffphis_safe as fits_write_history;
    pub use crate::putkey::ffphpr_safe as fits_write_grphdr;
    pub use crate::putkey::ffphprll_safe as fits_write_grphdrll;
    pub use crate::putkey::ffphps_safe as fits_write_imghdr;
    pub use crate::putkey::ffphpsll_safe as fits_write_imghdrll;
    pub use crate::putkey::ffphtb_safe as fits_write_atblhdr;
    pub use crate::putkey::ffpkfc_safe as fits_write_key_fixcmp;
    pub use crate::putkey::ffpkfm_safe as fits_write_key_fixdblcmp;
    pub use crate::putkey::ffpkls_safe as fits_write_key_longstr;
    pub use crate::putkey::ffpknd_safe as fits_write_keys_dbl;
    pub use crate::putkey::ffpkne_safe as fits_write_keys_flt;
    pub use crate::putkey::ffpknf_safe as fits_write_keys_fixflt;
    pub use crate::putkey::ffpkng_safe as fits_write_keys_fixdbl;
    pub use crate::putkey::ffpknj_safe as fits_write_keys_lng;
    pub use crate::putkey::ffpknl_safe as fits_write_keys_log;
    pub use crate::putkey::ffpkns_safe as fits_write_keys_str;
    pub use crate::putkey::ffpktp_safe as fits_write_key_template;
    pub use crate::putkey::ffpky_safe as fits_write_key;
    pub use crate::putkey::ffpkyc_safe as fits_write_key_cmp;
    pub use crate::putkey::ffpkyd_safe as fits_write_key_dbl;
    pub use crate::putkey::ffpkye_safe as fits_write_key_flt;
    pub use crate::putkey::ffpkyf_safe as fits_write_key_fixflt;
    pub use crate::putkey::ffpkyg_safe as fits_write_key_fixdbl;
    pub use crate::putkey::ffpkyj_safe as fits_write_key_lng;
    pub use crate::putkey::ffpkyl_safe as fits_write_key_log;
    pub use crate::putkey::ffpkym_safe as fits_write_key_dblcmp;
    pub use crate::putkey::ffpkys_safe as fits_write_key_str;
    pub use crate::putkey::ffpkyt_safe as fits_write_key_triple;
    pub use crate::putkey::ffpkyu_safe as fits_write_key_null;
    pub use crate::putkey::ffpkyuj_safe as fits_write_key_ulng;
    pub use crate::putkey::ffplsw_safe as fits_write_key_longwarn;
    pub use crate::putkey::ffprec_safe as fits_write_record;
    pub use crate::putkey::ffptdm_safe as fits_write_tdim;
    pub use crate::putkey::ffptdmll_safe as fits_write_tdimll;
    pub use crate::putkey::ffs2dt_safe as fits_str2date;
    pub use crate::putkey::ffs2tm_safe as fits_str2time;
    pub use crate::putkey::fftm2s_safe as fits_time2str;

    pub use crate::getkey::ffghps_safe as fits_get_hdrpos;
    pub use crate::getkey::ffghsp_safe as fits_get_hdrspace;
    pub use crate::getkey::ffgnxk_safe as fits_find_nextkey;
    pub use crate::getkey::ffmaky_safe as fits_movabs_key;
    pub use crate::getkey::ffmrky_safe as fits_movrel_key;

    pub use crate::getkey::ffcnvthdr2str_safer as fits_convert_hdr2str;
    pub use crate::getkey::ffdtdm_safe as fits_decode_tdim;
    pub use crate::getkey::ffdtdmll_safe as fits_decode_tdimll;
    pub use crate::getkey::fffree_safer as fits_free_memory;
    pub use crate::getkey::ffgcrd_safe as fits_read_card;
    pub use crate::getkey::ffghbn_safer as fits_read_btblhdr;
    pub use crate::getkey::ffghbnll_safer as fits_read_btblhdrll;
    pub use crate::getkey::ffghpr_safe as fits_read_imghdr;
    pub use crate::getkey::ffghprll_safe as fits_read_imghdrll;
    pub use crate::getkey::ffghtb_safer as fits_read_atblhdr;
    pub use crate::getkey::ffghtbll_safer as fits_read_atblhdrll;
    pub use crate::getkey::ffgkcsl_safe as fits_get_key_com_strlen;
    pub use crate::getkey::ffgkey_safe as fits_read_keyword;
    pub use crate::getkey::ffgkls_safe as fits_read_key_longstr;
    pub use crate::getkey::ffgknd_safe as fits_read_keys_dbl;
    pub use crate::getkey::ffgkne_safe as fits_read_keys_flt;
    pub use crate::getkey::ffgknj_safe as fits_read_keys_lng;
    pub use crate::getkey::ffgknjj_safe as fits_read_keys_lnglng;
    pub use crate::getkey::ffgknl_safe as fits_read_keys_log;
    pub use crate::getkey::ffgkns_safe as fits_read_keys_str;
    pub use crate::getkey::ffgksl_safe as fits_get_key_strlen;
    pub use crate::getkey::ffgky_safe as fits_read_key;
    pub use crate::getkey::ffgkyc_safe as fits_read_key_cmp;
    pub use crate::getkey::ffgkyd_safe as fits_read_key_dbl;
    pub use crate::getkey::ffgkye_safe as fits_read_key_flt;
    pub use crate::getkey::ffgkyj_safe as fits_read_key_lng;
    pub use crate::getkey::ffgkyjj_safe as fits_read_key_lnglng;
    pub use crate::getkey::ffgkyl_safe as fits_read_key_log;
    pub use crate::getkey::ffgkym_safe as fits_read_key_dblcmp;
    pub use crate::getkey::ffgkyn_safe as fits_read_keyn;
    pub use crate::getkey::ffgkys_safe as fits_read_key_str;
    pub use crate::getkey::ffgkyt_safe as fits_read_key_triple;
    pub use crate::getkey::ffgkyujj_safe as fits_read_key_ulnglng;
    pub use crate::getkey::ffgrec_safe as fits_read_record;
    pub use crate::getkey::ffgsky_safe as fits_read_string_key;
    pub use crate::getkey::ffgskyc_safe as fits_read_string_key_com;
    pub use crate::getkey::ffgstr_safe as fits_read_str;
    pub use crate::getkey::ffgtdm_safe as fits_read_tdim;
    pub use crate::getkey::ffgtdmll_safer as fits_read_tdimll;
    pub use crate::getkey::ffgunt_safe as fits_read_key_unit;
    pub use crate::getkey::ffhdr2str_safe as fits_hdr2str;

    pub use crate::modkey::ffucrd_safe as fits_update_card;
    pub use crate::modkey::ffukfc_safe as fits_update_key_fixcmp;
    pub use crate::modkey::ffukfm_safe as fits_update_key_fixdblcmp;
    pub use crate::modkey::ffukls_safe as fits_update_key_longstr;
    pub use crate::modkey::ffuky_safe as fits_update_key;
    pub use crate::modkey::ffukyc_safe as fits_update_key_cmp;
    pub use crate::modkey::ffukyd_safe as fits_update_key_dbl;
    pub use crate::modkey::ffukye_safe as fits_update_key_flt;
    pub use crate::modkey::ffukyf_safe as fits_update_key_fixflt;
    pub use crate::modkey::ffukyg_safe as fits_update_key_fixdbl;
    pub use crate::modkey::ffukyj_safe as fits_update_key_lng;
    pub use crate::modkey::ffukyl_safe as fits_update_key_log;
    pub use crate::modkey::ffukym_safe as fits_update_key_dblcmp;
    pub use crate::modkey::ffukys_safe as fits_update_key_str;
    pub use crate::modkey::ffukyu_safe as fits_update_key_null;
    pub use crate::modkey::ffukyuj_safe as fits_update_key_ulng;

    pub use crate::modkey::ffmcom_safe as fits_modify_comment;
    pub use crate::modkey::ffmcrd_safe as fits_modify_card;
    pub use crate::modkey::ffmkfc_safe as fits_modify_key_fixcmp;
    pub use crate::modkey::ffmkfm_safe as fits_modify_key_fixdblcmp;
    pub use crate::modkey::ffmkls_safe as fits_modify_key_longstr;
    pub use crate::modkey::ffmkyc_safe as fits_modify_key_cmp;
    pub use crate::modkey::ffmkyd_safe as fits_modify_key_dbl;
    pub use crate::modkey::ffmkye_safe as fits_modify_key_flt;
    pub use crate::modkey::ffmkyf_safe as fits_modify_key_fixflt;
    pub use crate::modkey::ffmkyg_safe as fits_modify_key_fixdbl;
    pub use crate::modkey::ffmkyj_safe as fits_modify_key_lng;
    pub use crate::modkey::ffmkyl_safe as fits_modify_key_log;
    pub use crate::modkey::ffmkym_safe as fits_modify_key_dblcmp;
    pub use crate::modkey::ffmkys_safe as fits_modify_key_str;
    pub use crate::modkey::ffmkyu_safe as fits_modify_key_null;
    pub use crate::modkey::ffmkyuj_safe as fits_modify_key_ulng;
    pub use crate::modkey::ffmnam_safe as fits_modify_name;
    pub use crate::modkey::ffmrec_safe as fits_modify_record;

    pub use crate::modkey::ffikey_safe as fits_insert_card;
    pub use crate::modkey::ffikfc_safer as fits_insert_key_fixcmp;
    pub use crate::modkey::ffikfm_safer as fits_insert_key_fixdblcmp;
    pub use crate::modkey::ffikls_safe as fits_insert_key_longstr;
    pub use crate::modkey::ffikyc_safer as fits_insert_key_cmp;
    pub use crate::modkey::ffikyd_safer as fits_insert_key_dbl;
    pub use crate::modkey::ffikye_safer as fits_insert_key_flt;
    pub use crate::modkey::ffikyf_safe as fits_insert_key_fixflt;
    pub use crate::modkey::ffikyg_safer as fits_insert_key_fixdbl;
    pub use crate::modkey::ffikyj_safe as fits_insert_key_lng;
    pub use crate::modkey::ffikyl_safer as fits_insert_key_log;
    pub use crate::modkey::ffikym_safe as fits_insert_key_dblcmp;
    pub use crate::modkey::ffikys_safer as fits_insert_key_str;
    pub use crate::modkey::ffikyu_safer as fits_insert_key_null;
    pub use crate::modkey::ffirec_safe as fits_insert_record;

    pub use crate::fitscore::ffghad_safe as fits_get_hduaddr;
    pub use crate::fitscore::ffghadll_safe as fits_get_hduaddrll;
    pub use crate::fitscore::ffghdn_safe as fits_get_hdu_num;
    pub use crate::fitscore::ffghdt_safe as fits_get_hdu_type;
    pub use crate::fitscore::ffghof_safe as fits_get_hduoff;
    pub use crate::modkey::ffdkey_safe as fits_delete_key;
    pub use crate::modkey::ffdrec_safe as fits_delete_record;
    pub use crate::modkey::ffdstr_safe as fits_delete_str;

    pub use crate::fitscore::ffgipr_safe as fits_get_img_param;
    pub use crate::fitscore::ffgiprll_safe as fits_get_img_paramll;

    pub use crate::fitscore::ffgidm_safe as fits_get_img_dim;
    pub use crate::fitscore::ffgidt_safe as fits_get_img_type;
    pub use crate::fitscore::ffgiet_safe as fits_get_img_equivtype;
    pub use crate::fitscore::ffgisz_safe as fits_get_img_size;
    pub use crate::fitscore::ffgiszll_safe as fits_get_img_sizell;

    pub use crate::editcol::ffrsim_safe as fits_resize_img;
    pub use crate::editcol::ffrsimll_safe as fits_resize_imgll;
    pub use crate::edithdu::ffibin_safe as fits_insert_btbl;
    pub use crate::edithdu::ffiimg_safer as fits_insert_img;
    pub use crate::edithdu::ffiimgll_safer as fits_insert_imgll;
    pub use crate::edithdu::ffitab_safe as fits_insert_atbl;
    pub use crate::fitscore::ffcrhd_safer as fits_create_hdu;
    pub use crate::fitscore::ffmahd_safe as fits_movabs_hdu;
    pub use crate::fitscore::ffmnhd_safe as fits_movnam_hdu;
    pub use crate::fitscore::ffmrhd_safe as fits_movrel_hdu;
    pub use crate::fitscore::ffthdu_safe as fits_get_num_hdus;
    pub use crate::putkey::ffcrim_safer as fits_create_img;
    pub use crate::putkey::ffcrimll_safer as fits_create_imgll;
    pub use crate::putkey::ffcrtb_safer as fits_create_tbl;

    pub use crate::edithdu::ffcopy_safer as fits_copy_hdu;
    pub use crate::edithdu::ffcpdt_safe as fits_copy_data;
    pub use crate::edithdu::ffcpfl_safer as fits_copy_file;
    pub use crate::edithdu::ffcphd_safer as fits_copy_header;
    pub use crate::edithdu::ffcpht_safer as fits_copy_hdutab;
    pub use crate::edithdu::ffdhdu_safer as fits_delete_hdu;
    pub use crate::edithdu::ffwrhdu_safer as fits_write_hdu;

    pub use crate::fitscore::ffhdef_safe as fits_set_hdrsize;
    pub use crate::fitscore::ffrdef_safe as fits_set_hdustruc;
    pub use crate::scalnull::ffpthp_safe as fits_write_theap;

    pub use crate::checksum::ffdsum_safe as fits_decode_chksum;
    pub use crate::checksum::ffesum_safe as fits_encode_chksum;
    pub use crate::checksum::ffgcks_safe as fits_get_chksum;
    pub use crate::checksum::ffpcks_safe as fits_write_chksum;
    pub use crate::checksum::ffupck_safe as fits_update_chksum;
    pub use crate::checksum::ffvcks_safe as fits_verify_chksum;

    pub use crate::scalnull::ffpnul_safe as fits_set_imgnull;
    pub use crate::scalnull::ffpscl_safe as fits_set_bscale;
    pub use crate::scalnull::ffsnul_safe as fits_set_atblnull;
    pub use crate::scalnull::fftnul_safe as fits_set_btblnull;
    pub use crate::scalnull::fftscl_safe as fits_set_tscale;

    pub use crate::fitscore::ffeqty_safe as fits_get_eqcoltype;
    pub use crate::fitscore::ffeqtyll_safe as fits_get_eqcoltypell;
    pub use crate::fitscore::ffgacl_safer as fits_get_acolparms;
    pub use crate::fitscore::ffgbcl_safe as fits_get_bcolparms;
    pub use crate::fitscore::ffgbclll_safe as fits_get_bcolparmsll;
    pub use crate::fitscore::ffgcnn_safe as fits_get_colname;
    pub use crate::fitscore::ffgcno_safe as fits_get_colnum;
    pub use crate::fitscore::ffgncl_safe as fits_get_num_cols;
    pub use crate::fitscore::ffgnrw_safe as fits_get_num_rows;
    pub use crate::fitscore::ffgnrwll_safe as fits_get_num_rowsll;
    pub use crate::fitscore::ffgtcl_safe as fits_get_coltype;
    pub use crate::fitscore::ffgtclll_safe as fits_get_coltypell;

    pub use crate::putcol::ffiter_safe as fits_iterate_data;

    pub use crate::getcolb::ffggpb_safe as fits_read_grppar_byt;
    pub use crate::getcold::ffggpd_safe as fits_read_grppar_dbl;
    pub use crate::getcole::ffggpe_safe as fits_read_grppar_flt;
    pub use crate::getcoli::ffggpi_safe as fits_read_grppar_sht;
    pub use crate::getcolj::ffggpj_safe as fits_read_grppar_lng;
    pub use crate::getcolj::ffggpjj_safe as fits_read_grppar_lnglng;
    pub use crate::getcolk::ffggpk_safe as fits_read_grppar_int;
    pub use crate::getcolsb::ffggpsb_safe as fits_read_grppar_sbyt;
    pub use crate::getcolui::ffggpui_safe as fits_read_grppar_usht;
    pub use crate::getcoluj::ffggpuj_safe as fits_read_grppar_ulng;
    pub use crate::getcoluj::ffggpujj_safe as fits_read_grppar_ulnglng;
    pub use crate::getcoluk::ffggpuk_safe as fits_read_grppar_uint;

    pub use crate::getcol::ffgpf_safe as fits_read_imgnull;
    pub use crate::getcol::ffgpv_safe as fits_read_img;
    pub use crate::getcol::ffgpxf_safe as fits_read_pixnull;
    pub use crate::getcol::ffgpxfll_safe as fits_read_pixnullll;
    pub use crate::getcol::ffgpxv_safe as fits_read_pix;
    pub use crate::getcol::ffgpxvll_safe as fits_read_pixll;
    pub use crate::getcolb::ffgpvb_safe as fits_read_img_byt;
    pub use crate::getcold::ffgpvd_safe as fits_read_img_dbl;
    pub use crate::getcole::ffgpve_safe as fits_read_img_flt;
    pub use crate::getcoli::ffgpvi_safe as fits_read_img_sht;
    pub use crate::getcolj::ffgpvj_safe as fits_read_img_lng;
    pub use crate::getcolj::ffgpvjj_safe as fits_read_img_lnglng;
    pub use crate::getcolk::ffgpvk_safe as fits_read_img_int;
    pub use crate::getcolsb::ffgpvsb_safe as fits_read_img_sbyt;
    pub use crate::getcolui::ffgpvui_safe as fits_read_img_usht;
    pub use crate::getcoluj::ffgpvuj_safe as fits_read_img_ulng;
    pub use crate::getcoluj::ffgpvujj_safe as fits_read_img_ulnglng;
    pub use crate::getcoluk::ffgpvuk_safe as fits_read_img_uint;

    pub use crate::getcolb::ffgpfb_safe as fits_read_imgnull_byt;
    pub use crate::getcold::ffgpfd_safe as fits_read_imgnull_dbl;
    pub use crate::getcole::ffgpfe_safe as fits_read_imgnull_flt;
    pub use crate::getcoli::ffgpfi_safe as fits_read_imgnull_sht;
    pub use crate::getcolj::ffgpfj_safe as fits_read_imgnull_lng;
    pub use crate::getcolj::ffgpfjj_safe as fits_read_imgnull_lnglng;
    pub use crate::getcolk::ffgpfk_safe as fits_read_imgnull_int;
    pub use crate::getcolsb::ffgpfsb_safe as fits_read_imgnull_sbyt;
    pub use crate::getcolui::ffgpfui_safe as fits_read_imgnull_usht;
    pub use crate::getcoluj::ffgpfuj_safe as fits_read_imgnull_ulng;
    pub use crate::getcoluj::ffgpfujj_safe as fits_read_imgnull_ulnglng;
    pub use crate::getcoluk::ffgpfuk_safe as fits_read_imgnull_uint;

    pub use crate::getcolb::ffg2db_safe as fits_read_2d_byt;
    pub use crate::getcold::ffg2dd_safe as fits_read_2d_dbl;
    pub use crate::getcole::ffg2de_safe as fits_read_2d_flt;
    pub use crate::getcoli::ffg2di_safe as fits_read_2d_sht;
    pub use crate::getcolj::ffg2dj_safe as fits_read_2d_lng;
    pub use crate::getcolj::ffg2djj_safe as fits_read_2d_lnglng;
    pub use crate::getcolk::ffg2dk_safe as fits_read_2d_int;
    pub use crate::getcolsb::ffg2dsb_safe as fits_read_2d_sbyt;
    pub use crate::getcolui::ffg2dui_safe as fits_read_2d_usht;
    pub use crate::getcoluj::ffg2duj_safe as fits_read_2d_ulng;
    pub use crate::getcoluj::ffg2dujj_safe as fits_read_2d_ulnglng;
    pub use crate::getcoluk::ffg2duk_safe as fits_read_2d_uint;

    pub use crate::getcolb::ffg3db_safe as fits_read_3d_byt;
    pub use crate::getcold::ffg3dd_safe as fits_read_3d_dbl;
    pub use crate::getcole::ffg3de_safe as fits_read_3d_flt;
    pub use crate::getcoli::ffg3di_safe as fits_read_3d_sht;
    pub use crate::getcolj::ffg3dj_safe as fits_read_3d_lng;
    pub use crate::getcolj::ffg3djj_safe as fits_read_3d_lnglng;
    pub use crate::getcolk::ffg3dk_safe as fits_read_3d_int;
    pub use crate::getcolsb::ffg3dsb_safe as fits_read_3d_sbyt;
    pub use crate::getcolui::ffg3dui_safe as fits_read_3d_usht;
    pub use crate::getcoluj::ffg3duj_safe as fits_read_3d_ulng;
    pub use crate::getcoluj::ffg3dujj_safe as fits_read_3d_ulnglng;
    pub use crate::getcoluk::ffg3duk_safe as fits_read_3d_uint;

    pub use crate::getcol::ffgsv_safe as fits_read_subset;
    pub use crate::getcolb::ffgsvb_safe as fits_read_subset_byt;
    pub use crate::getcold::ffgsvd_safe as fits_read_subset_dbl;
    pub use crate::getcole::ffgsve_safe as fits_read_subset_flt;
    pub use crate::getcoli::ffgsvi_safe as fits_read_subset_sht;
    pub use crate::getcolj::ffgsvj_safe as fits_read_subset_lng;
    pub use crate::getcolj::ffgsvjj_safe as fits_read_subset_lnglng;
    pub use crate::getcolk::ffgsvk_safe as fits_read_subset_int;
    pub use crate::getcolsb::ffgsvsb_safe as fits_read_subset_sbyt;
    pub use crate::getcolui::ffgsvui_safe as fits_read_subset_usht;
    pub use crate::getcoluj::ffgsvuj_safe as fits_read_subset_ulng;
    pub use crate::getcoluj::ffgsvujj_safe as fits_read_subset_ulnglng;
    pub use crate::getcoluk::ffgsvuk_safe as fits_read_subset_uint;

    pub use crate::getcolb::ffgsfb_safe as fits_read_subsetnull_byt;
    pub use crate::getcold::ffgsfd_safe as fits_read_subsetnull_dbl;
    pub use crate::getcole::ffgsfe_safe as fits_read_subsetnull_flt;
    pub use crate::getcoli::ffgsfi_safe as fits_read_subsetnull_sht;
    pub use crate::getcolj::ffgsfj_safe as fits_read_subsetnull_lng;
    pub use crate::getcolj::ffgsfjj_safe as fits_read_subsetnull_lnglng;
    pub use crate::getcolk::ffgsfk_safe as fits_read_subsetnull_int;
    pub use crate::getcolsb::ffgsfsb_safe as fits_read_subsetnull_sbyt;
    pub use crate::getcolui::ffgsfui_safe as fits_read_subsetnull_usht;
    pub use crate::getcoluj::ffgsfuj_safe as fits_read_subsetnull_ulng;
    pub use crate::getcoluj::ffgsfujj_safe as fits_read_subsetnull_ulnglng;
    pub use crate::getcoluk::ffgsfuk_safe as fits_read_subsetnull_uint;

    pub use crate::cfileio::fits_copy_image_section_safer as ffcpimg;

    pub use crate::getcol::ffgcf_safe as fits_read_colnull;
    pub use crate::getcol::ffgcv_safe as fits_read_col;
    pub use crate::getcol::ffgcvn_safer as fits_read_cols;
    pub use crate::getcolb::ffgcvb_safe as fits_read_col_byt;
    pub use crate::getcold::ffgcvd_safe as fits_read_col_dbl;
    pub use crate::getcold::ffgcvm_safe as fits_read_col_dblcmp;
    pub use crate::getcole::ffgcvc_safe as fits_read_col_cmp;
    pub use crate::getcole::ffgcve_safe as fits_read_col_flt;
    pub use crate::getcoli::ffgcvi_safe as fits_read_col_sht;
    pub use crate::getcolj::ffgcvj_safe as fits_read_col_lng;
    pub use crate::getcolj::ffgcvjj_safe as fits_read_col_lnglng;
    pub use crate::getcolk::ffgcvk_safe as fits_read_col_int;
    pub use crate::getcoll::ffgcvl_safe as fits_read_col_log;
    pub use crate::getcoll::ffgcx_safe as fits_read_col_bit;
    pub use crate::getcoll::ffgcxui_safe as fits_read_col_bit_usht;
    pub use crate::getcoll::ffgcxuk_safe as fits_read_col_bit_uint;
    pub use crate::getcols::ffgcvs_safe as fits_read_col_str;
    pub use crate::getcolsb::ffgcvsb_safe as fits_read_col_sbyt;
    pub use crate::getcolui::ffgcvui_safe as fits_read_col_usht;
    pub use crate::getcoluj::ffgcvuj_safe as fits_read_col_ulng;
    pub use crate::getcoluj::ffgcvujj_safe as fits_read_col_ulnglng;
    pub use crate::getcoluk::ffgcvuk_safe as fits_read_col_uint;

    pub use crate::getcolb::ffgcfb_safe as fits_read_colnull_byt;
    pub use crate::getcold::ffgcfd_safe as fits_read_colnull_dbl;
    pub use crate::getcole::ffgcfc_safe as fits_read_colnull_cmp;
    pub use crate::getcole::ffgcfe_safe as fits_read_colnull_flt;
    pub use crate::getcoli::ffgcfi_safe as fits_read_colnull_sht;
    pub use crate::getcolj::ffgcfj_safe as fits_read_colnull_lng;
    pub use crate::getcolj::ffgcfjj_safe as fits_read_colnull_lnglng;
    pub use crate::getcolk::ffgcfk_safe as fits_read_colnull_int;
    pub use crate::getcoll::ffgcfl_safe as fits_read_colnull_log;
    pub use crate::getcols::ffgcfs_safe as fits_read_colnull_str;
    pub use crate::getcolsb::ffgcfsb_safe as fits_read_colnull_sbyt;
    pub use crate::getcolui::ffgcfui_safe as fits_read_colnull_usht;
    pub use crate::getcoluj::ffgcfuj_safe as fits_read_colnull_ulng;
    pub use crate::getcoluj::ffgcfujj_safe as fits_read_colnull_ulnglng;
    pub use crate::getcoluk::ffgcfuk_safe as fits_read_colnull_uint;

    pub use crate::buffers::ffgtbb_safe as fits_read_tblbytes;
    pub use crate::fitscore::ffgdes_safe as fits_read_descript;
    pub use crate::fitscore::ffgdesll_safe as fits_read_descriptll;
    pub use crate::fitscore::ffgdess_safe as fits_read_descripts;
    pub use crate::fitscore::ffgdessll_safe as fits_read_descriptsll;

    pub use crate::putcolb::ffpgpb_safe as fits_write_grppar_byt;
    pub use crate::putcold::ffpgpd_safe as fits_write_grppar_dbl;
    pub use crate::putcole::ffpgpe_safe as fits_write_grppar_flt;
    pub use crate::putcoli::ffpgpi_safe as fits_write_grppar_sht;
    pub use crate::putcolj::ffpgpj_safe as fits_write_grppar_lng;
    pub use crate::putcolj::ffpgpjj_safe as fits_write_grppar_lnglng;
    pub use crate::putcolk::ffpgpk_safe as fits_write_grppar_int;
    pub use crate::putcolsb::ffpgpsb_safe as fits_write_grppar_sbyt;
    pub use crate::putcolui::ffpgpui_safe as fits_write_grppar_usht;
    pub use crate::putcoluj::ffpgpuj_safe as fits_write_grppar_ulng;
    pub use crate::putcoluj::ffpgpujj_safe as fits_write_grppar_ulnglng;
    pub use crate::putcoluk::ffpgpuk_safe as fits_write_grppar_uint;

    pub use crate::putcol::ffppr_safe as fits_write_img;
    pub use crate::putcol::ffppx_safe as fits_write_pix;
    pub use crate::putcol::ffppxll_safe as fits_write_pixll;
    pub use crate::putcol::ffppxn_safe as fits_write_pixnull;
    pub use crate::putcol::ffppxnll_safe as fits_write_pixnullll;
    pub use crate::putcolb::ffpprb_safe as fits_write_img_byt;
    pub use crate::putcold::ffpprd_safe as fits_write_img_dbl;
    pub use crate::putcole::ffppre_safe as fits_write_img_flt;
    pub use crate::putcoli::ffppri_safe as fits_write_img_sht;
    pub use crate::putcolj::ffpprj_safe as fits_write_img_lng;
    pub use crate::putcolj::ffpprjj_safe as fits_write_img_lnglng;
    pub use crate::putcolk::ffpprk_safe as fits_write_img_int;
    pub use crate::putcolsb::ffpprsb_safe as fits_write_img_sbyt;
    pub use crate::putcolui::ffpprui_safe as fits_write_img_usht;
    pub use crate::putcoluj::ffppruj_safe as fits_write_img_ulng;
    pub use crate::putcoluj::ffpprujj_safe as fits_write_img_ulnglng;
    pub use crate::putcoluk::ffppruk_safe as fits_write_img_uint;

    pub use crate::putcol::ffppn_safe as fits_write_imgnull;
    pub use crate::putcolb::ffppnb_safe as fits_write_imgnull_byt;
    pub use crate::putcold::ffppnd_safe as fits_write_imgnull_dbl;
    pub use crate::putcole::ffppne_safe as fits_write_imgnull_flt;
    pub use crate::putcoli::ffppni_safe as fits_write_imgnull_sht;
    pub use crate::putcolj::ffppnj_safe as fits_write_imgnull_lng;
    pub use crate::putcolj::ffppnjj_safe as fits_write_imgnull_lnglng;
    pub use crate::putcolk::ffppnk_safe as fits_write_imgnull_int;
    pub use crate::putcolsb::ffppnsb_safe as fits_write_imgnull_sbyt;
    pub use crate::putcolui::ffppnui_safe as fits_write_imgnull_usht;
    pub use crate::putcoluj::ffppnuj_safe as fits_write_imgnull_ulng;
    pub use crate::putcoluj::ffppnujj_safe as fits_write_imgnull_ulnglng;
    pub use crate::putcoluk::ffppnuk_safe as fits_write_imgnull_uint;

    pub use crate::putcolu::ffpprn_safe as fits_write_null_img;
    pub use crate::putcolu::ffppru_safe as fits_write_img_null;

    pub use crate::putcolb::ffp2db_safe as fits_write_2d_byt;
    pub use crate::putcold::ffp2dd_safe as fits_write_2d_dbl;
    pub use crate::putcole::ffp2de_safe as fits_write_2d_flt;
    pub use crate::putcoli::ffp2di_safe as fits_write_2d_sht;
    pub use crate::putcolj::ffp2dj_safe as fits_write_2d_lng;
    pub use crate::putcolj::ffp2djj_safe as fits_write_2d_lnglng;
    pub use crate::putcolk::ffp2dk_safe as fits_write_2d_int;
    pub use crate::putcolsb::ffp2dsb_safe as fits_write_2d_sbyt;
    pub use crate::putcolui::ffp2dui_safe as fits_write_2d_usht;
    pub use crate::putcoluj::ffp2duj_safe as fits_write_2d_ulng;
    pub use crate::putcoluj::ffp2dujj_safe as fits_write_2d_ulnglng;
    pub use crate::putcoluk::ffp2duk_safe as fits_write_2d_uint;

    pub use crate::putcolb::ffp3db_safe as fits_write_3d_byt;
    pub use crate::putcold::ffp3dd_safe as fits_write_3d_dbl;
    pub use crate::putcole::ffp3de_safe as fits_write_3d_flt;
    pub use crate::putcoli::ffp3di_safe as fits_write_3d_sht;
    pub use crate::putcolj::ffp3dj_safe as fits_write_3d_lng;
    pub use crate::putcolj::ffp3djj_safe as fits_write_3d_lnglng;
    pub use crate::putcolk::ffp3dk_safe as fits_write_3d_int;
    pub use crate::putcolsb::ffp3dsb_safe as fits_write_3d_sbyt;
    pub use crate::putcolui::ffp3dui_safe as fits_write_3d_usht;
    pub use crate::putcoluj::ffp3duj_safe as fits_write_3d_ulng;
    pub use crate::putcoluj::ffp3dujj_safe as fits_write_3d_ulnglng;
    pub use crate::putcoluk::ffp3duk_safe as fits_write_3d_uint;

    pub use crate::putcol::ffpss_safe as fits_write_subset;
    pub use crate::putcolb::ffpssb_safe as fits_write_subset_byt;
    pub use crate::putcold::ffpssd_safe as fits_write_subset_dbl;
    pub use crate::putcole::ffpsse_safe as fits_write_subset_flt;
    pub use crate::putcoli::ffpssi_safe as fits_write_subset_sht;
    pub use crate::putcolj::ffpssj_safe as fits_write_subset_lng;
    pub use crate::putcolj::ffpssjj_safe as fits_write_subset_lnglng;
    pub use crate::putcolk::ffpssk_safe as fits_write_subset_int;
    pub use crate::putcolsb::ffpsssb_safe as fits_write_subset_sbyt;
    pub use crate::putcolui::ffpssui_safe as fits_write_subset_usht;
    pub use crate::putcoluj::ffpssuj_safe as fits_write_subset_ulng;
    pub use crate::putcoluj::ffpssujj_safe as fits_write_subset_ulnglng;
    pub use crate::putcoluk::ffpssuk_safe as fits_write_subset_uint;

    pub use crate::putcol::ffpcl_safe as fits_write_col;
    pub use crate::putcol::ffpcln_safe as fits_write_cols;
    pub use crate::putcolb::ffpclb_safe as fits_write_col_byt;
    pub use crate::putcold::ffpcld_safe as fits_write_col_dbl;
    pub use crate::putcold::ffpclm_safe as fits_write_col_dblcmp;
    pub use crate::putcole::ffpclc_safe as fits_write_col_cmp;
    pub use crate::putcole::ffpcle_safe as fits_write_col_flt;
    pub use crate::putcoli::ffpcli_safe as fits_write_col_sht;
    pub use crate::putcolj::ffpclj_safe as fits_write_col_lng;
    pub use crate::putcolj::ffpcljj_safe as fits_write_col_lnglng;
    pub use crate::putcolk::ffpclk_safe as fits_write_col_int;
    pub use crate::putcoll::ffpcll_safe as fits_write_col_log;
    pub use crate::putcoll::ffpclx_safe as fits_write_col_bit;
    pub use crate::putcols::ffpcls_safe as fits_write_col_str;
    pub use crate::putcolsb::ffpclsb_safe as fits_write_col_sbyt;
    pub use crate::putcolu::ffpclu_safe as fits_write_col_null;
    pub use crate::putcolu::ffprwu_safe as fits_write_nulrows;
    pub use crate::putcolu::ffprwu_safe as fits_write_nullrows;
    pub use crate::putcolui::ffpclui_safe as fits_write_col_usht;
    pub use crate::putcoluj::ffpcluj_safe as fits_write_col_ulng;
    pub use crate::putcoluj::ffpclujj_safe as fits_write_col_ulnglng;
    pub use crate::putcoluk::ffpcluk_safe as fits_write_col_uint;

    pub use crate::putcol::ffpcn_safer as fits_write_colnull;
    pub use crate::putcolb::ffpcnb_safe as fits_write_colnull_byt;
    pub use crate::putcold::ffpcnd_safe as fits_write_colnull_dbl;
    pub use crate::putcole::ffpcne_safe as fits_write_colnull_flt;
    pub use crate::putcoli::ffpcni_safe as fits_write_colnull_sht;
    pub use crate::putcolj::ffpcnj_safe as fits_write_colnull_lng;
    pub use crate::putcolj::ffpcnjj_safe as fits_write_colnull_lnglng;
    pub use crate::putcolk::ffpcnk_safe as fits_write_colnull_int;
    pub use crate::putcoll::ffpcnl_safe as fits_write_colnull_log;
    pub use crate::putcols::ffpcns_safe as fits_write_colnull_str;
    pub use crate::putcolsb::ffpcnsb_safe as fits_write_colnull_sbyt;
    pub use crate::putcolui::ffpcnui_safe as fits_write_colnull_usht;
    pub use crate::putcoluj::ffpcnuj_safe as fits_write_colnull_ulng;
    pub use crate::putcoluj::ffpcnujj_safe as fits_write_colnull_ulnglng;
    pub use crate::putcoluk::ffpcnuk_safe as fits_write_colnull_uint;

    pub use crate::getcolb::ffgextn_safe as fits_read_ext;
    pub use crate::putcolb::ffpextn_safe as fits_write_ext;

    pub use crate::fitscore::ffcmph_safer as fits_compress_heap;
    pub use crate::fitscore::ffpdes_safe as fits_write_descript;
    pub use crate::fitscore::fftheap_safe as fits_test_heap;

    pub use crate::buffers::ffptbb_safe as fits_write_tblbytes;
    pub use crate::editcol::ffccls_safe as fits_copy_cols;
    pub use crate::editcol::ffcpcl_safe as fits_copy_col;
    pub use crate::editcol::ffcprw_safe as fits_copy_rows;
    pub use crate::editcol::ffcpsr_safe as fits_copy_selrows;
    pub use crate::editcol::ffdcol_safe as fits_delete_col;
    pub use crate::editcol::ffdrow_safe as fits_delete_rows;
    pub use crate::editcol::ffdrrg_safe as fits_delete_rowrange;
    pub use crate::editcol::ffdrws_safe as fits_delete_rowlist;
    pub use crate::editcol::ffdrwsll_safe as fits_delete_rowlistll;
    pub use crate::editcol::fficls_safe as fits_insert_cols;
    pub use crate::editcol::fficol_safe as fits_insert_col;
    pub use crate::editcol::ffirow_safe as fits_insert_rows;
    pub use crate::editcol::ffmvec_safe as fits_modify_vector_len;

    pub use crate::wcssub::ffgics_safe as fits_read_img_coord;
    pub use crate::wcssub::ffgicsa_safe as fits_read_img_coord_version;
    pub use crate::wcssub::ffgtcs_safer as fits_read_tbl_coord;
    pub use crate::wcsutil::ffwldp_safe as fits_pix_to_world;
    pub use crate::wcsutil::ffxypx_safe as fits_world_to_pix;

    pub use crate::wcssub::ffgiwcs_safe as fits_get_image_wcs_keys;
    pub use crate::wcssub::ffgtwcs_safe as fits_get_table_wcs_keys;

    pub use crate::eval_f::ffcalc_rng_safe as fits_calculator_rng;
    pub use crate::eval_f::ffcalc_safe as fits_calculator;
    pub use crate::eval_f::ffcrow_safe as fits_calc_rows;
    pub use crate::eval_f::ffffrw_safer as fits_find_first_row;
    pub use crate::eval_f::fffrow_safe as fits_find_rows;
    pub use crate::eval_f::fffrwc_safe as fits_find_rows_cmp;
    pub use crate::eval_f::ffsrow_safe as fits_select_rows;
    pub use crate::eval_f::fftexp_safe as fits_test_expr;

    pub use crate::group::ffgtam_safe as fits_add_group_member;
    pub use crate::group::ffgtch_safe as fits_change_group;
    pub use crate::group::ffgtcm_safe as fits_compact_group;
    pub use crate::group::ffgtcp_safe as fits_copy_group;
    pub use crate::group::ffgtcr_safe as fits_create_group;
    pub use crate::group::ffgtis_safe as fits_insert_group;
    pub use crate::group::ffgtmg_safe as fits_merge_groups;
    pub use crate::group::ffgtnm_safe as fits_get_num_members;
    pub use crate::group::ffgtop_safe as fits_open_group;
    pub use crate::group::ffgtrm_safe as fits_remove_group;
    pub use crate::group::ffgtvf_safe as fits_verify_group;

    pub use crate::group::ffgmcp_safe as fits_copy_member;
    pub use crate::group::ffgmng_safe as fits_get_num_groups;
    pub use crate::group::ffgmop_safe as fits_open_member;
    pub use crate::group::ffgmrm_safe as fits_remove_member;
    pub use crate::group::ffgmtf_safe as fits_transfer_member;

    pub use crate::cfileio::ffgtmo_safer as fits_get_timeout;
    pub use crate::cfileio::ffshdwn_safe as fits_show_download_progress;
    pub use crate::cfileio::ffstmo_safer as fits_set_timeout;

    pub use crate::cfileio::ffchtps_safer as fits_cleanup_https;
    pub use crate::cfileio::ffihtps_safer as fits_init_https;
    pub use crate::cfileio::ffvhtps_safer as fits_verbose_https;
}
