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

#include "Cabac.h"


// sao()
template<>
int ContextInit<sao_merge_X_flag>::initValue[3] =
{
        153, 153, 153
};

template<>
int ContextInit<sao_type_idx_X>::initValue[3] =
{
        200, 185, 160
};

// coding_quadtree()
template<>
int ContextInit<split_cu_flag>::initValue[9] =
{
        139, 141, 157, 107, 139, 126, 107, 139, 126
};


// coding_unit()
template<>
int ContextInit<cu_transquant_bypass_flag>::initValue[3] =
{
        154, 154, 154
};

template<>
int ContextInit<cu_skip_flag>::initValue[6] =
{
        197, 185, 201, 197, 185, 201
};

template<>
int ContextInit<cu_qp_delta_abs>::initValue[6] =
{
        154, 154, 154, 154, 154, 154
};

template<>
int ContextInit<cu_chroma_qp_offset_flag>::initValue[3] =
{
        154, 154, 154
};

template<>
int ContextInit<cu_chroma_qp_offset_idx>::initValue[3] =
{
        154, 154, 154
};

template<>
int ContextInit<pred_mode_flag>::initValue[2] =
{
        149, 134
};

template<>
int ContextInit<part_mode>::initValue[9] =
{
        184, 154, 139, 154, 154, 154, 139, 154, 154
};


// prediction_unit()
template<>
int ContextInit<prev_intra_luma_pred_flag>::initValue[3] =
{
        184, 154, 183
};

template<>
int ContextInit<intra_chroma_pred_mode>::initValue[3] =
{
        63, 152, 152
};

template<>
int ContextInit<merge_flag>::initValue[2] =
{
        110, 154
};

template<>
int ContextInit<merge_idx>::initValue[2] =
{
        122, 137
};

template<>
int ContextInit<inter_pred_idc>::initValue[10] =
{
        95, 79, 63, 31, 31, 95, 79, 63, 31, 31
};

template<>
int ContextInit<ref_idx_lX>::initValue[4] =
{
        153, 153, 153, 153
};

template<>
int ContextInit<abs_mvd_greaterX_flag>::initValue[4] =
{
        140, 169, 198, 198
};

template<>
int ContextInit<mvp_lX_flag>::initValue[2] =
{
        168, 168
};


// transform_tree()
template<>
int ContextInit<rqt_root_cbf>::initValue[2] =
{
        79, 79
};

template<>
int ContextInit<split_transform_flag>::initValue[9] =
{
        153, 138, 138, 124, 138, 94, 224, 167, 122  // different to standard text?
};

template<>
int ContextInit<cbf_luma>::initValue[6] =
{
        111, 141, 153, 111, 153, 111
};

template<>
int ContextInit<cbf_cX>::initValue[12] =
{
        94, 138, 182, 154, 149, 107, 167, 154, 149, 92, 167, 154
};


// residual_coding()
template<>
int ContextInit<transform_skip_flag_X>::initValue[6] =
{
        139, 139, 139, 139, 139, 139
};

template<>
int ContextInit<last_sig_coeff_x_prefix>::initValue[54] =
{
        110, 110, 124, 125, 140, 153, 125, 127, 140, 109, 111, 143, 127, 111, 79, 108, 123, 63,
        125, 110, 94, 110, 95, 79, 125, 111, 110, 78, 110, 111, 111, 95, 94, 108, 123, 108,
        125, 110, 124, 110, 95, 94, 125, 111, 111, 79, 125, 126, 111, 111, 79, 108, 123, 93
};

template<>
int ContextInit<last_sig_coeff_y_prefix>::initValue[54] =
{
        110, 110, 124, 125, 140, 153, 125, 127, 140, 109, 111, 143, 127, 111, 79, 108, 123, 63,
        125, 110, 94, 110, 95, 79, 125, 111, 110, 78, 110, 111, 111, 95, 94, 108, 123, 108,
        125, 110, 124, 110, 95, 94, 125, 111, 111, 79, 125, 126, 111, 111, 79, 108, 123, 93
};

template<>
int ContextInit<coded_sub_block_flag>::initValue[12] =
{
        91, 171, 134, 141, 121, 140, 61, 154, 121, 140, 61, 154
};

template<>
int ContextInit<sig_coeff_flag>::initValue[132] =
{
        111, 111, 125, 110, 110, 94, 124, 108, 124, 107, 125, 141, 179, 153, 125, 107,
        125, 141, 179, 153, 125, 107, 125, 141, 179, 153, 125, 140, 139, 182, 182, 152,
        136, 152, 136, 153, 136, 139, 111, 136, 139, 111, 141, 111,

        155, 154, 139, 153, 139, 123,
        123, 63, 153, 166, 183, 140, 136, 153, 154, 166, 183, 140, 136, 153, 154, 166,
        183, 140, 136, 153, 154, 170, 153, 123, 123, 107, 121, 107, 121, 167, 151, 183,
        140, 151, 183, 140, 140, 140,

        170, 154, 139, 153, 139, 123, 123, 63, 124, 166, 183, 140,
        136, 153, 154, 166, 183, 140, 136, 153, 154, 166, 183, 140, 136, 153, 154, 170,
        153, 138, 138, 122, 121, 122, 121, 167, 151, 183, 140, 151, 183, 140, 140, 140
};

template<>
int ContextInit<coeff_abs_level_greater1_flag>::initValue[72] =
{
        140, 92, 137, 138, 140, 152, 138, 139, 153, 74, 149, 92, 139, 107, 122, 152,
        140, 179, 166, 182, 140, 227, 122, 197, 154, 196, 196, 167, 154, 152, 167, 182,
        182, 134, 149, 136, 153, 121, 136, 137, 169, 194, 166, 167, 154, 167, 137, 182,
        154, 196, 167, 167, 154, 152, 167, 182, 182, 134, 149, 136, 153, 121, 136, 122,
        169, 208, 166, 167, 154, 152, 167, 182
};

template<>
int ContextInit<coeff_abs_level_greater2_flag>::initValue[18] =
{
        138, 153, 136, 167, 152, 152, 107, 167, 91, 122, 107, 167, 107, 167, 91, 107, 107, 167
};

template<>
int ContextInit<explicit_rdpcm_flag>::initValue[4] =
{
        139, 139, 139, 139
};

template<>
int ContextInit<explicit_rdpcm_dir_flag>::initValue[4] =
{
        139, 139, 139, 139
};

template<>
int ContextInit<log2_res_scale_abs_plus1>::initValue[24] =
{
        154, 154, 154, 154, 154, 154, 154, 154, 154,
        154, 154, 154, 154, 154, 154, 154, 154, 154,
        154, 154, 154, 154, 154, 154
};

template<>
int ContextInit<res_scale_sign_flag>::initValue[6] =
{
        154, 154, 154, 154, 154, 154
};
