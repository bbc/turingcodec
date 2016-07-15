/*
Copyright (C) 2016 British Broadcasting Corporation, Parabola Research
and Queen Mary University of London.

This file is part of the Turing codec.

The Turing codec is free software; you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as
published by the Free Software Foundation.

The Turing codec is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Commercial support and intellectual property rights for
the Turing codec are also available under a proprietary license.
For more information, contact us at info @ turingcodec.org.
 */

#ifndef INCLUDED_RangeLimits_h
#define INCLUDED_RangeLimits_h

#include "Global.h"



// Used as baseclass for RangeLimits<> if no limits are defined
struct NoRangeLimits { };


// Defines minimum and maximum limits of legal values for value V.
template <class V> struct RangeLimits : NoRangeLimits { };


// Used as baseclass for RangeLimits<> specialisation if parsing should be aborted if value out of range
// In general, this is only used where an invalid value could cause a crash, for example if it was used as an
// index to a static array.
struct Fatal { };


// CondCheck 7.4.3.1-C
template <> struct RangeLimits<vps_max_sub_layers_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 6; }
    static const char *clause() { return "7.4.3.1"; }
};


// CondCheck 7.4.3.2-A
template <> struct RangeLimits<sps_max_sub_layers_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 6; }
    static const char *clause() { return "7.4.3.2"; }
};


// CondCheck 7.4.3.2-N
template <> struct RangeLimits<log2_max_pic_order_cnt_lsb_minus4>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 12; }
    static const char *clause() { return "7.4.3.2"; }
};


// CondCheck 7.4.3.1-J
template <> struct RangeLimits<vps_max_latency_increase_plus1>
{
    template <class H> static unsigned min(H &) { return 0; }
    template <class H> static unsigned max(H &) { return 0xfffffffe; }
    static const char *clause() { return "7.4.3.1"; }
};


// CondCheck 7.4.3.1-P
template <> struct RangeLimits<vps_num_ticks_poc_diff_one_minus1>
{
    template <class H> static unsigned min(H &) { return 0; }
    template <class H> static unsigned max(H &) { return 0xfffffffe; }
    static const char *clause() { return "7.4.3.1"; }
};


// CondCheck 7.4.3.2-Y
template <> struct RangeLimits<max_transform_hierarchy_depth_inter>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[CtbLog2SizeY()] - h[MinTbLog2SizeY()]; }
    static const char *clause() { return "7.4.3.2"; }
};


// CondCheck 7.4.3.2-Z
template <> struct RangeLimits<max_transform_hierarchy_depth_intra>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[CtbLog2SizeY()] - h[MinTbLog2SizeY()]; }
    static const char *clause() { return "7.4.3.2"; }
};


// CondCheck 7.4.3.2-AE
template <> struct RangeLimits<num_short_term_ref_pic_sets>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 64; }
    static const char *clause() { return "7.4.3.2"; }
};


// CondCheck 7.4.3.2-AF
template <> struct RangeLimits<num_long_term_ref_pics_sps>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 32; }
    static const char *clause() { return "7.4.3.2"; }
};


// CondCheck 7.4.3.3-D
template <> struct RangeLimits<num_ref_idx_l0_default_active_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 14; }
    static const char *clause() { return "7.4.3.3"; }
};


// CondCheck 7.4.3.3-E
template <> struct RangeLimits<num_ref_idx_l1_default_active_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 14; }
    static const char *clause() { return "7.4.3.3"; }
};


// CondCheck 7.4.3.3-H
template <> struct RangeLimits<pps_cb_qp_offset>
{
    template <class H> static int min(H &) { return -12; }
    template <class H> static int max(H &) { return 12; }
    static const char *clause() { return "7.4.3.3"; }
};


// CondCheck 7.4.3.3-I
template <> struct RangeLimits<pps_cr_qp_offset>
{
    template <class H> static int min(H &) { return -12; }
    template <class H> static int max(H &) { return 12; }
    static const char *clause() { return "7.4.3.3"; }
};


// CondCheck 7.4.3.3-S
template <> struct RangeLimits<pps_beta_offset_div2>
{
    template <class H> static int min(H &) { return -6; }
    template <class H> static int max(H &) { return 6; }
    static const char *clause() { return "7.4.3.3"; }
};


// CondCheck 7.4.3.3-T
template <> struct RangeLimits<pps_tc_offset_div2>
{
    template <class H> static int min(H &) { return -6; }
    template <class H> static int max(H &) { return 6; }
    static const char *clause() { return "7.4.3.3"; }
};


// CondCheck 7.4.5-B
template <> struct RangeLimits<scaling_list_dc_coef_minus8>
{
    template <class H> static int min(H &) { return -7; }
    template <class H> static int max(H &) { return 247; }
    static const char *clause() { return "7.4.5"; }
};


// CondCheck 7.4.5-C
template <> struct RangeLimits<scaling_list_delta_coef>
{
    template <class H> static int min(H &) { return -128; }
    template <class H> static int max(H &) { return 127; }
    static const char *clause() { return "7.4.5"; }
};


// CondCheck 7.4.7.1-N
template <> struct RangeLimits<slice_segment_address>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[PicSizeInCtbsY()] - 1; }
    static const char *clause() { return "7.4.7.1"; }
};


// CondCheck 7.4.7.1-S
template <> struct RangeLimits<colour_plane_id>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 2; }
    static const char *clause() { return "7.4.7.1"; }
};


// CondCheck 7.4.7.1-V
template <> struct RangeLimits<short_term_ref_pic_set_idx>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[num_short_term_ref_pic_sets()] - 1; }
    static const char *clause() { return "7.4.7.1"; }
};

// CondCheck 7.4.7.1-AC
template <> struct RangeLimits<num_ref_idx_l0_active_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 14; }
    static const char *clause() { return "7.4.7.1"; }
};


// CondCheck 7.4.7.1-AD
template <> struct RangeLimits<num_ref_idx_l1_active_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 14; }
    static const char *clause() { return "7.4.7.1"; }
};


// CondCheck 7.4.7.1-AN
template <> struct RangeLimits<slice_beta_offset_div2>
{
    template <class H> static int min(H &) { return -6; }
    template <class H> static int max(H &) { return 6; }
    static const char *clause() { return "7.4.7.1"; }
};


// CondCheck 7.4.7.1-AO
template <> struct RangeLimits<slice_tc_offset_div2>
{
    template <class H> static int min(H &) { return -6; }
    template <class H> static int max(H &) { return 6; }
    static const char *clause() { return "7.4.7.1"; }
};


template <> struct RangeLimits<offset_len_minus1> // CondCheck 7.4.7.1-AS
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 31; }
    static const char *clause() { return "7.4.7.1"; }
};

// CondCheck 7.4.7.1-AZ
template <> struct RangeLimits<slice_segment_header_extension_length>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 256; }
    static const char *clause() { return "7.4.7.1"; }
};

// CondCheck 7.4.7.2-A
template <> struct RangeLimits<list_entry_l0>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[NumPocTotalCurr()] - 1; }
    static const char *clause() { return "7.4.7.2"; }
};

// CondCheck 7.4.7.2-B
template <> struct RangeLimits<list_entry_l1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[NumPocTotalCurr()] - 1; }
    static const char *clause() { return "7.4.7.2"; }
};

// CondCheck 7.4.7.3-A
template <> struct RangeLimits<luma_log2_weight_denom>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 7; }
    static const char *clause() { return "7.4.7.3"; }
};

// CondCheck 7.4.7.3-D
template <> struct RangeLimits<luma_offset_l0>
{
    template <class H> static int min(H &) { return -128; }
    template <class H> static int max(H &) { return 127; }
    static const char *clause() { return "7.4.7.3"; }
};

// CondCheck 7.4.7.3-F
template <> struct RangeLimits<delta_chroma_offset_l0>
{
    template <class H> static int min(H &) { return -512; }
    template <class H> static int max(H &) { return 511; }
    static const char *clause() { return "7.4.7.3"; }
};

// CondCheck 7.4.7.3-H
template <> struct RangeLimits<luma_offset_l1>
{
    template <class H> static int min(H &) { return -128; }
    template <class H> static int max(H &) { return 127; }
    static const char *clause() { return "7.4.7.3"; }
};

// CondCheck 7.4.7.3-J
template <> struct RangeLimits<delta_chroma_offset_l1>
{
    template <class H> static int min(H &) { return -512; }
    template <class H> static int max(H &) { return 511; }
    static const char *clause() { return "7.4.7.3"; }
};

// CondCheck 7.4.7.2-Q
template <> struct RangeLimits<delta_idx_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return static_cast<short_term_ref_pic_set *>(h)->stRpsIdx - 1; }
    static const char *clause() { return "7.4.7.2"; }
};

// CondCheck 7.4.7.2-R
template <> struct RangeLimits<abs_delta_rps_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return (1 << 15) - 1; }
    static const char *clause() { return "7.4.7.2"; }
};

// CondCheck 7.4.7.2-S
template <> struct RangeLimits<num_negative_pics>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[sps_max_dec_pic_buffering_minus1(h[sps_max_sub_layers_minus1()]) ]; }
    static const char *clause() { return "7.4.7.2"; }
};

// CondCheck 7.4.7.2-T
template <> struct RangeLimits<num_positive_pics>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h) { return h[sps_max_dec_pic_buffering_minus1(h[sps_max_sub_layers_minus1()]) ] - h[num_negative_pics()]; }
    static const char *clause() { return "7.4.7.2"; }
};

// CondCheck 7.4.7.2-U
template <> struct RangeLimits<delta_poc_s0_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return (1 << 15) - 1; }
    static const char *clause() { return "7.4.7.2"; }
};

// CondCheck 7.4.7.2-V
template <> struct RangeLimits<delta_poc_s1_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return (1 << 15) - 1; }
    static const char *clause() { return "7.4.7.2"; }
};


// CondCheck E.2.2-A
template <> struct RangeLimits<elemental_duration_in_tc_minus1>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 2047; }
    static const char *clause() { return "E.2.2"; }
};

// CondCheck E.2.2-D
template <> struct RangeLimits<cpb_cnt_minus1> :
    Fatal
    {
        template <class H> static int min(H &) { return 0; }
        template <class H> static int max(H &) { return 31; }
        static const char *clause() { return "E.2.2, E.3.2"; }
    };

// CondCheck 7.4.9.11-A
template <> struct RangeLimits<last_sig_coeff_x_prefix>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h)
    {
        const residual_coding &rc = *static_cast<residual_coding *>(h);
        return (rc.log2TrafoSize  <<  1) - 1;
    }
    static const char *clause() { return "7.4.9.11"; }
};

// CondCheck 7.4.9.11-B
template <> struct RangeLimits<last_sig_coeff_y_prefix>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h)
    {
        const residual_coding &rc = *static_cast<residual_coding *>(h);
        return (rc.log2TrafoSize  <<  1) - 1;
    }
    static const char *clause() { return "7.4.9.11"; }
};

// CondCheck 7.4.9.11-C
template <> struct RangeLimits<last_sig_coeff_x_suffix>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h)
    {
        const residual_coding &rc = *static_cast<residual_coding *>(h);
        return (1 << ((h[last_sig_coeff_x_prefix()] >> 1) - 1)) - 1;
    }
    static const char *clause() { return "7.4.9.11"; }
};

// CondCheck 7.4.9.11-D
template <> struct RangeLimits<last_sig_coeff_y_suffix>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &h)
    {
        const residual_coding &rc = *static_cast<residual_coding *>(h);
        return ( 1  <<  ( ( h[last_sig_coeff_y_prefix()]  >>  1 ) - 1 ) ) - 1;
    }
    static const char *clause() { return "7.4.9.11"; }
};

// CondCheck D.3.2-E0
template <> struct RangeLimits<bp_seq_parameter_set_id>
{
    template <class H> static int min(H &) { return 0; }
    template <class H> static int max(H &) { return 15; }
    static const char *clause() { return "D.3.2"; }
};

// review: prevent parsing of unimplemented syntax
template <> struct RangeLimits<pps_multilayer_extension_flag> :
    Fatal
    {
        template <class H> static int min(H &) { return 0; }
        template <class H> static int max(H &) { return 0; }
        static const char *clause() { return "unimplemented"; }
    };

#endif
