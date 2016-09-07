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

// RBSP syntax functions.
// In general these are more-or-less copy/pasted from the HEVC standard text and should not be modified.

#ifndef INCLUDED_SyntaxRbsp_hpp
#define INCLUDED_SyntaxRbsp_hpp

#pragma once

#include "Syntax.h"
#include "SyntaxElements.h"
#include "HevcTypes.h"
#include "StateEncode.h"
#include "Global.h"
#include <cassert>


// Raw byte sequence payloads, trailing bits, and byte alignment syntax

template <class H>
void Syntax<video_parameter_set_rbsp>::go(const video_parameter_set_rbsp &, H &h)
{
    h(vps_video_parameter_set_id(), u(4));
    h(vps_reserved_three_2bits(), u(2));
    h(vps_max_layers_minus1(), u(6));
    h(vps_max_sub_layers_minus1(), u(3));
    h(vps_temporal_id_nesting_flag(), u(1));
    h(vps_reserved_0xffff_16bits(), u(16));

    h(profile_tier_level(h[vps_max_sub_layers_minus1()]));

    h(vps_sub_layer_ordering_info_present_flag(), u(1));

    for (int i = (h[vps_sub_layer_ordering_info_present_flag()] ?
            0 : h[vps_max_sub_layers_minus1()]);
            i <= h[vps_max_sub_layers_minus1()]; i++)
    {
        h(vps_max_dec_pic_buffering_minus1(i), ue(v));
        h(vps_max_num_reorder_pics(i), ue(v));
        h(vps_max_latency_increase_plus1(i), ue(v));
    }

    h(vps_max_layer_id(), u(6));
    h(vps_num_layer_sets_minus1(), ue(v));

    for (int i = 1; i <= h[vps_num_layer_sets_minus1()]; i++)
    {
        for (int j = 0; j <= h[vps_max_layer_id()]; j++)
        {
            h(layer_id_included_flag(), u(1));
        }
    }

    h(vps_timing_info_present_flag(), u(1));

    if (h[vps_timing_info_present_flag()])
    {
        h(vps_num_units_in_tick(), u(32));
        h(vps_time_scale(), u(32));
        h(vps_poc_proportional_to_timing_flag(), u(1));

        if (h[vps_poc_proportional_to_timing_flag()])
            h(vps_num_ticks_poc_diff_one_minus1(), ue(v));

        h(vps_num_hrd_parameters(), ue(v));

        for (int i = 0; i < h[vps_num_hrd_parameters()]; i++)
        {
            h(hrd_layer_set_idx(i), ue(v));
            if (i > 0)
            {
                h(cprms_present_flag(i), u(1));
            }
            h(hrd_parameters(h[cprms_present_flag(i)], h[vps_max_sub_layers_minus1()]));
        }
    }

    h(vps_extension_flag(), u(1));
    if (h[vps_extension_flag()])
    {
        while (h[more_rbsp_data()])
        {
            h(vps_extension_data_flag(), u(1));
        }
    }

    h(rbsp_trailing_bits());
}


template <class H>
void Syntax<seq_parameter_set_rbsp>::go(const seq_parameter_set_rbsp&, H &h)
{
    h(sps_video_parameter_set_id(), u(4));
    if (h[nuh_layer_id()] == 0)
        h(sps_max_sub_layers_minus1(), u(3));
    else
        h(sps_ext_or_max_sub_layers_minus1(), u(3));
    h[MultiLayerExtSpsFlag()] = (h[nuh_layer_id()] != 0 && h[sps_ext_or_max_sub_layers_minus1()] == 7);
    if (!h[MultiLayerExtSpsFlag()])
    {
        h(sps_temporal_id_nesting_flag(), u(1));
        h(profile_tier_level(h[sps_max_sub_layers_minus1()]));
    }
    h(sps_seq_parameter_set_id(), ue(v));
    if (h[MultiLayerExtSpsFlag()])
    {
        assert(!"not yet implemented");
    }
    else
    {
        h(chroma_format_idc(), ue(v));
        if (h[chroma_format_idc()] == 3)
            h(separate_colour_plane_flag(), u(1));
        h(pic_width_in_luma_samples(), ue(v));
        h(pic_height_in_luma_samples(), ue(v));
        h(conformance_window_flag(), u(1));
        if (h[conformance_window_flag()])
        {
            h(conf_win_left_offset(), ue(v));
            h(conf_win_right_offset(), ue(v));
            h(conf_win_top_offset(), ue(v));
            h(conf_win_bottom_offset(), ue(v));
        }
        h(bit_depth_luma_minus8(), ue(v));
        h(bit_depth_chroma_minus8(), ue(v));
    }
    h(log2_max_pic_order_cnt_lsb_minus4(), ue(v));
    if (!h[MultiLayerExtSpsFlag()])
    {
        h(sps_sub_layer_ordering_info_present_flag(), u(1));
        for (int i = (h[sps_sub_layer_ordering_info_present_flag()] ? 0 : h[sps_max_sub_layers_minus1()]); i <= h[sps_max_sub_layers_minus1()]; i++)
        {
            h(sps_max_dec_pic_buffering_minus1(i), ue(v));
            h(sps_max_num_reorder_pics(i), ue(v));
            h(sps_max_latency_increase_plus1(i), ue(v));
        }
    }
    h(log2_min_luma_coding_block_size_minus3(), ue(v));
    h(log2_diff_max_min_luma_coding_block_size(), ue(v));
    h(log2_min_luma_transform_block_size_minus2(), ue(v));
    h(log2_diff_max_min_luma_transform_block_size(), ue(v));
    h(max_transform_hierarchy_depth_inter(), ue(v));
    h(max_transform_hierarchy_depth_intra(), ue(v));
    h(scaling_list_enabled_flag(), u(1));
    if (h[scaling_list_enabled_flag()])
    {
        if (h[MultiLayerExtSpsFlag()])
        {
            assert(!"not yet implemented");
        }
        else
        {
            h(sps_scaling_list_data_present_flag(), u(1));
            if (h[sps_scaling_list_data_present_flag()])
                h(scaling_list_data());
        }
    }
    h(amp_enabled_flag(), u(1));
    h(sample_adaptive_offset_enabled_flag(), u(1));
    h(pcm_enabled_flag(), u(1));
    if (h[pcm_enabled_flag()])
    {
        h(pcm_sample_bit_depth_luma_minus1(), u(4));
        h(pcm_sample_bit_depth_chroma_minus1(), u(4));
        h(log2_min_pcm_luma_coding_block_size_minus3(), ue(v));
        h(log2_diff_max_min_pcm_luma_coding_block_size(), ue(v));
        h(pcm_loop_filter_disabled_flag(), u(1));
    }
    h(num_short_term_ref_pic_sets(), ue(v));
    for (int i = 0; i < h[num_short_term_ref_pic_sets()]; i++)
        h(short_term_ref_pic_set(i));

    h(long_term_ref_pics_present_flag(), u(1));
    if (h[long_term_ref_pics_present_flag()])
    {
        h(num_long_term_ref_pics_sps(), ue(v));
        for (int i = 0; i < h[num_long_term_ref_pics_sps()]; i++)
        {
            h(lt_ref_pic_poc_lsb_sps(i), uv());
            h(used_by_curr_pic_lt_sps_flag(i), u(1));
        }
    }
    h(sps_temporal_mvp_enabled_flag(), u(1));
    h(strong_intra_smoothing_enabled_flag(), u(1));

    h(vui_parameters_present_flag(), u(1));
    if (h[vui_parameters_present_flag()])
        h(vui_parameters());

    h(sps_extension_present_flag(), u(1));
    if (h[sps_extension_present_flag()])
    {
        h(sps_range_extension_flag(), u(1));
        h(sps_multilayer_extension_flag(), u(1));
        h(sps_extension_6bits(), u(6));
    }

    if (h[sps_range_extension_flag()])
        h(sps_range_extension());

    if (h[sps_multilayer_extension_flag()])
    {
        assert(!"not yet implemented");
        //sps_multilayer_extension()
    }

    if (h[sps_extension_6bits()])
        while (h[more_rbsp_data()])
            h(sps_extension_data_flag(), u(1));

    h(rbsp_trailing_bits());
}

NUMBER_OF_BITS_UV(lt_ref_pic_poc_lsb_sps, h[log2_max_pic_order_cnt_lsb_minus4()] + 4)



template <class H>
void Syntax<sps_range_extension>::go(const sps_range_extension&, H &h)
{
    h(transform_skip_rotation_enabled_flag(), u(1));
    h(transform_skip_context_enabled_flag(), u(1));
    h(implicit_rdpcm_enabled_flag(), u(1));
    h(explicit_rdpcm_enabled_flag(), u(1));
    h(extended_precision_processing_flag(), u(1));
    h(intra_smoothing_disabled_flag(), u(1));
    h(high_precision_offsets_enabled_flag(), u(1));
    h(persistent_rice_adaptation_enabled_flag(), u(1));
    h(cabac_bypass_alignment_enabled_flag(), u(1));
}


template <class H>
void Syntax<pic_parameter_set_rbsp>::go(const pic_parameter_set_rbsp&, H &h)
{
    h(pps_pic_parameter_set_id(), ue(v));
    h(pps_seq_parameter_set_id(), ue(v));
    h(dependent_slice_segments_enabled_flag(), u(1));
    h(output_flag_present_flag(), u(1));
    h(num_extra_slice_header_bits(), u(3));
    h(sign_data_hiding_enabled_flag(), u(1));
    h(cabac_init_present_flag(), u(1));
    h(num_ref_idx_l0_default_active_minus1(), ue(v));
    h(num_ref_idx_l1_default_active_minus1(), ue(v));
    h(init_qp_minus26(), se(v));
    h(constrained_intra_pred_flag(), u(1));
    h(transform_skip_enabled_flag(), u(1));
    h(cu_qp_delta_enabled_flag(), u(1));
    if (h[cu_qp_delta_enabled_flag()])
        h(diff_cu_qp_delta_depth(), ue(v));
    h(pps_cb_qp_offset(), se(v));
    h(pps_cr_qp_offset(), se(v));
    h(pps_slice_chroma_qp_offsets_present_flag(), u(1));
    h(weighted_pred_flag(), u(1));
    h(weighted_bipred_flag(), u(1));
    h(transquant_bypass_enabled_flag(), u(1));
    h(tiles_enabled_flag(), u(1));
    h(entropy_coding_sync_enabled_flag(), u(1));
    if (h[tiles_enabled_flag()])
    {
        h(num_tile_columns_minus1(), ue(v));
        h(num_tile_rows_minus1(), ue(v));
        h(uniform_spacing_flag(), u(1));
        if (!h[uniform_spacing_flag()])
        {
            for (int i = 0; i < h[num_tile_columns_minus1()]; i++)
                h(column_width_minus1(i), ue(v));
            for (int i = 0; i < h[num_tile_rows_minus1()]; i++)
                h(row_height_minus1(i), ue(v));
        }
        h(loop_filter_across_tiles_enabled_flag(), u(1));
    }
    h(pps_loop_filter_across_slices_enabled_flag(), u(1));
    h(deblocking_filter_control_present_flag(), u(1));

    if (h[deblocking_filter_control_present_flag()])
    {
        h(deblocking_filter_override_enabled_flag(), u(1));
        h(pps_deblocking_filter_disabled_flag(), u(1));
        if (!h[pps_deblocking_filter_disabled_flag()])
        {
            h(pps_beta_offset_div2(), se(v));
            h(pps_tc_offset_div2(), se(v));
        }
    }

    h(pps_scaling_list_data_present_flag(), u(1));
    if (h[pps_scaling_list_data_present_flag()])
        h(scaling_list_data());

    h(lists_modification_present_flag(), u(1));
    h(log2_parallel_merge_level_minus2(), ue(v));
    h(slice_segment_header_extension_present_flag(), u(1));

    h(pps_extension_present_flag(), u(1));
    if (h[pps_extension_present_flag()])
    {
        h(pps_range_extension_flag(), u(1));
        h(pps_multilayer_extension_flag(), u(1));
        h(pps_extension_6bits(), u(6));
    }

    if (h[pps_range_extension_flag()])
        h(pps_range_extension());

    if (h[pps_multilayer_extension_flag()])
    {
        assert(!"not yet implemented");
        //		h(pps_multilayer_extension());  /* specified in Annex F */
    }

    if (h[pps_extension_6bits()])
        while (h[more_rbsp_data()])
            h(pps_extension_data_flag(), u(1));

    h(rbsp_trailing_bits());
}


template <class H>
void Syntax<pps_range_extension>::go(const pps_range_extension&, H &h)
{
    if (h[transform_skip_enabled_flag()])
        h(log2_max_transform_skip_block_size_minus2(), ue(v));
    h(cross_component_prediction_enabled_flag(), u(1));
    h(chroma_qp_offset_list_enabled_flag(), u(1));
    if (h[chroma_qp_offset_list_enabled_flag()])
    {
        h(diff_cu_chroma_qp_offset_depth(), ue(v));
        h(chroma_qp_offset_list_len_minus1(), ue(v));
        for (int i = 0; i <= h[chroma_qp_offset_list_len_minus1()]; i++)
        {
            h(cb_qp_offset_list(i), se(v));
            h(cr_qp_offset_list(i), se(v));
        }
    }
    h(log2_sao_offset_scale_luma(), ue(v));
    h(log2_sao_offset_scale_chroma(), ue(v));
}


template <class H>
void Syntax<sei_rbsp>::go(const sei_rbsp&, H &h)
{
    do
    {
        h(sei_message());
    } while (h[more_rbsp_data()]);
    h(rbsp_trailing_bits());
}


template <class H>
void Syntax<access_unit_delimiter_rbsp>::go(const access_unit_delimiter_rbsp&, H &h)
{
    h(pic_type(), u(3));
    h(rbsp_trailing_bits());
}


template <class H>
void Syntax<end_of_seq_rbsp>::go(const end_of_seq_rbsp&, H &h)
{
}


template <class H>
void Syntax<end_of_bitstream_rbsp>::go(const end_of_bitstream_rbsp&, H &h)
{
}


template <class H>
void Syntax<filler_data_rbsp>::go(const filler_data_rbsp&, H &h)
{
    while (next_bits<int>(h, 8, true) == 0xFF)
        h(ff_byte() /* equal to 0xFF */, f(8));
    h(rbsp_trailing_bits());
}


template <class H>
void Syntax<slice_segment_layer_rbsp>::go(const slice_segment_layer_rbsp&, H &h)
{
    h(slice_segment_header());
    h(slice_segment_data());
    h(rbsp_slice_segment_trailing_bits());
}


template <class H>
void Syntax<rbsp_slice_segment_trailing_bits>::go(const rbsp_slice_segment_trailing_bits &, H &h)
{
    h(rbsp_trailing_bits());
    while (h[more_rbsp_trailing_data()])
        h(cabac_zero_word() /* equal to 0x0000 */, f(16));
}


template <class H>
void Syntax<rbsp_trailing_bits>::go(const rbsp_trailing_bits &, H &h)
{
    h(rbsp_stop_one_bit() /* equal to 1 */, f(1));
    while (!h[byte_aligned()])
        h(rbsp_alignment_zero_bit() /* equal to 0 */, f(1));
}


template <class H>
void  Syntax<byte_alignment>::go(const byte_alignment &, H &h)
{
    h(alignment_bit_equal_to_one() /* equal to 1 */, f(1));
    while (!h[byte_aligned()])
        h(alignment_bit_equal_to_zero() /* equal to 0 */, f(1));
}


template <class H>
void Syntax<profile_tier_level>::go(const profile_tier_level &ptl, H &h)
{
    h(general_profile_space(), u(2));
    h(general_tier_flag(), u(1));
    h(general_profile_idc(), u(5));
    for (int j = 0; j < 32; j++)
    {
        h(general_profile_compatibility_flag(j), u(1));
    }
    h(general_progressive_source_flag(), u(1));
    h(general_interlaced_source_flag(), u(1));
    h(general_non_packed_constraint_flag(), u(1));
    h(general_frame_only_constraint_flag(), u(1));
    if (h[general_profile_idc()] == 4 || h[general_profile_compatibility_flag(4)] ||
            h[general_profile_idc()] == 5 || h[general_profile_compatibility_flag(5)] ||
            h[general_profile_idc()] == 6 || h[general_profile_compatibility_flag(6)] ||
            h[general_profile_idc()] == 7 || h[general_profile_compatibility_flag(7)])
    {
        /* The number of bits in this syntax structure is not affected by this condition */
        h(general_max_12bit_constraint_flag(), u(1));
        h(general_max_10bit_constraint_flag(), u(1));
        h(general_max_8bit_constraint_flag(), u(1));
        h(general_max_422chroma_constraint_flag(), u(1));
        h(general_max_420chroma_constraint_flag(), u(1));
        h(general_max_monochrome_constraint_flag(), u(1));
        h(general_intra_constraint_flag(), u(1));
        h(general_one_picture_only_constraint_flag(), u(1));
        h(general_lower_bit_rate_constraint_flag(), u(1));
        h(general_reserved_zero_34bits(), u(34));
    }
    else
    {
        h(general_reserved_zero_43bits(), u(43));
    }

    if ((h[general_profile_idc()] >= 1 && h[general_profile_idc()] <= 5) ||
            h[general_profile_compatibility_flag(1)] || h[general_profile_compatibility_flag(2)] ||
            h[general_profile_compatibility_flag(3)] || h[general_profile_compatibility_flag(4)] ||
            h[general_profile_compatibility_flag(5)])
        /* The number of bits in this syntax structure is not affected by this condition */
    {
        h(general_inbld_flag(), u(1));
    }
    else
    {
        h(general_reserved_zero_bit(), u(1));
    }

    h(general_level_idc(), u(8));
    for (int i = 0; i < ptl.maxNumSubLayersMinus1; i++)
    {
        h(sub_layer_profile_present_flag(i), u(1));
        h(sub_layer_level_present_flag(i), u(1));
    }

    if (ptl.maxNumSubLayersMinus1 > 0)
        for (int i = ptl.maxNumSubLayersMinus1; i < 8; i++)
            h(reserved_zero_2bits(i), u(2));

    for (int i = 0; i < ptl.maxNumSubLayersMinus1; i++)
    {
        if (h[sub_layer_profile_present_flag(i)])
        {
            h(sub_layer_profile_space(i), u(2));
            h(sub_layer_tier_flag(i), u(1));
            h(sub_layer_profile_idc(i), u(5));

            for (int j = 0; j < 32; j++)
                h(sub_layer_profile_compatibility_flag(i, j), u(1));

            h(sub_layer_progressive_source_flag(i), u(1));
            h(sub_layer_interlaced_source_flag(i), u(1));
            h(sub_layer_non_packed_constraint_flag(i), u(1));
            h(sub_layer_frame_only_constraint_flag(i), u(1));

            if (h[sub_layer_profile_idc(i)] == 4 || h[sub_layer_profile_compatibility_flag(i, 4)] ||
                    h[sub_layer_profile_idc(i)] == 5 || h[sub_layer_profile_compatibility_flag(i, 5)] ||
                    h[sub_layer_profile_idc(i)] == 6 || h[sub_layer_profile_compatibility_flag(i, 6)] ||
                    h[sub_layer_profile_idc(i)] == 7 || h[sub_layer_profile_compatibility_flag(i, 7)])
                /* The number of bits in this syntax structure is not affected by this condition */
            {
                h(sub_layer_max_12bit_constraint_flag(i), u(1));
                h(sub_layer_max_10bit_constraint_flag(i), u(1));
                h(sub_layer_max_8bit_constraint_flag(i), u(1));
                h(sub_layer_max_422chroma_constraint_flag(i), u(1));
                h(sub_layer_max_420chroma_constraint_flag(i), u(1));
                h(sub_layer_max_monochrome_constraint_flag(i), u(1));
                h(sub_layer_intra_constraint_flag(i), u(1));
                h(sub_layer_one_picture_only_constraint_flag(i), u(1));
                h(sub_layer_lower_bit_rate_constraint_flag(i), u(1));
                h(sub_layer_reserved_zero_34bits(i), u(34));
            }
            else
            {
                h(sub_layer_reserved_zero_43bits(i), u(43));
            }

            if ((h[sub_layer_profile_idc(i)] >= 1 && h[sub_layer_profile_idc(i)] <= 5) ||
                    h[sub_layer_profile_compatibility_flag(1)] ||
                    h[sub_layer_profile_compatibility_flag(2)] ||
                    h[sub_layer_profile_compatibility_flag(3)] ||
                    h[sub_layer_profile_compatibility_flag(4)] ||
                    h[sub_layer_profile_compatibility_flag(5)])
                /* The number of bits in this syntax structure is not affected by this condition */
            {
                h(sub_layer_inbld_flag(i), u(1));
            }
            else
            {
                h(sub_layer_reserved_zero_bit(i), u(1));
            }
        }

        if (h[sub_layer_level_present_flag(i)])
        {
            h(sub_layer_level_idc(i), u(8));
        }
    }
}


template <class H>
void Syntax<scaling_list_data>::go(const scaling_list_data &, H &h)
{
    for (int sizeId = 0; sizeId < 4; sizeId++)
    {
        for (int matrixId = 0; matrixId < ((sizeId == 3) ? 2 : 6); matrixId++)
        {
            h(scaling_list_pred_mode_flag(sizeId, matrixId), u(1));
            if (!h[scaling_list_pred_mode_flag(sizeId, matrixId)])
                h(scaling_list_pred_matrix_id_delta(sizeId, matrixId), ue(v));
            else
            {
                int nextCoef = 8;
                const int coefNum = std::min(64, (1 << (4 + (sizeId << 1))));
                if (sizeId > 1)
                {
                    h(scaling_list_dc_coef_minus8(sizeId - 2, matrixId), se(v));
                    nextCoef = h[scaling_list_dc_coef_minus8(sizeId - 2, matrixId)] + 8;
                }
                for (int i = 0; i < coefNum; i++)
                {
                    h(scaling_list_delta_coef(), se(v));
                    nextCoef = (nextCoef + h[scaling_list_delta_coef()] + 256) % 256;
                    h[ScalingList(sizeId, matrixId, i)] = nextCoef;
                }
            }
        }
    }
}


template <class H>
void Syntax<slice_segment_header>::go(const slice_segment_header &, H &h)
{
    h(first_slice_segment_in_pic_flag(), u(1));
    if (h[nal_unit_type()] >= BLA_W_LP && h[nal_unit_type()] <= RSV_IRAP_VCL23)
        h(no_output_of_prior_pics_flag(), u(1));
    h(slice_pic_parameter_set_id(), ue(v));
    if (!h[first_slice_segment_in_pic_flag()])
    {
        if (h[dependent_slice_segments_enabled_flag()])
            h(dependent_slice_segment_flag(), u(1));
        h(slice_segment_address(), u(ceilLog2(h[PicSizeInCtbsY()])));
    }
    if (!h[dependent_slice_segment_flag()])
    {
        for (int i = 0; i < h[num_extra_slice_header_bits()]; i++)
        {
            h(slice_reserved_flag(i), u(1));
        }
        h(slice_type(), ue(v));
        if (h[output_flag_present_flag()])
            h(pic_output_flag(), u(1));
        if (h[separate_colour_plane_flag()] == 1)
            h(colour_plane_id(), u(2));
        if (h[nal_unit_type()] != IDR_W_RADL && h[nal_unit_type()] != IDR_N_LP)
        {
            h(slice_pic_order_cnt_lsb(), u(h[log2_max_pic_order_cnt_lsb_minus4()] + 4));
            h(short_term_ref_pic_set_sps_flag(), u(1));
            if (!h[short_term_ref_pic_set_sps_flag()])
                h(short_term_ref_pic_set(h[num_short_term_ref_pic_sets()]));
            else if (h[num_short_term_ref_pic_sets()] > 1)
                h(short_term_ref_pic_set_idx(), u(ceilLog2(h[num_short_term_ref_pic_sets()])));
            else
                h[short_term_ref_pic_set_idx()] = 0;
            if (h[long_term_ref_pics_present_flag()])
            {
                if (h[num_long_term_ref_pics_sps()] > 0)
                {
                    h(num_long_term_sps(), ue(v));
                }
                h(num_long_term_pics(), ue(v));
                for (int i = 0; i < h[num_long_term_sps()] + h[num_long_term_pics()]; i++)
                {
                    if (i < h[num_long_term_sps()])
                    {
                        if (h[num_long_term_ref_pics_sps()] > 1)
                            h(lt_idx_sps(i), u(ceilLog2(h[num_long_term_ref_pics_sps()])));
                    }
                    else
                    {
                        h(poc_lsb_lt(i), u(h[log2_max_pic_order_cnt_lsb_minus4()] + 4));
                        h(used_by_curr_pic_lt_flag(i), u(1));
                    }
                    h(delta_poc_msb_present_flag(i), u(1));
                    if (h[delta_poc_msb_present_flag(i)])
                        h(delta_poc_msb_cycle_lt(i), ue(v));

                    if (i < h[num_long_term_sps()])
                    {
                        h[PocLsbLt(i)] = h[lt_ref_pic_poc_lsb_sps(h[lt_idx_sps(i)])];
                        h[UsedByCurrPicLt(i)] = h[used_by_curr_pic_lt_sps_flag(h[lt_idx_sps(i)])];
                    }
                    else
                    {
                        h[PocLsbLt(i)] = h[poc_lsb_lt(i)];
                        h[UsedByCurrPicLt(i)] = h[used_by_curr_pic_lt_flag(i)];
                    }
                }
            }
            if (h[sps_temporal_mvp_enabled_flag()])
                h(slice_temporal_mvp_enabled_flag(), u(1));
        }
        if (h[sample_adaptive_offset_enabled_flag()])
        {
            h(slice_sao_luma_flag(), u(1));
            h(slice_sao_chroma_flag(), u(1));
        }
        if (h[slice_type()] == P || h[slice_type()] == B)
        {
            h(num_ref_idx_active_override_flag(), u(1));
            if (h[num_ref_idx_active_override_flag()])
            {
                h(num_ref_idx_l0_active_minus1(), ue(v));
                if (h[slice_type()] == B)
                    h(num_ref_idx_l1_active_minus1(), ue(v));
            }
            if (h[lists_modification_present_flag()] && h[NumPocTotalCurr()] > 1)
                h(ref_pic_lists_modification());
            if (h[slice_type()] == B)
                h(mvd_l1_zero_flag(), u(1));
            if (h[cabac_init_present_flag()])
                h(cabac_init_flag(), u(1));
            if (h[slice_temporal_mvp_enabled_flag()])
            {
                if (h[slice_type()] == B)
                    h(collocated_from_l0_flag(), u(1));
                if ((h[collocated_from_l0_flag()] && h[num_ref_idx_l0_active_minus1()] > 0) ||
                        (!h[collocated_from_l0_flag()] && h[num_ref_idx_l1_active_minus1()] > 0))
                    h(collocated_ref_idx(), ue(v));
            }
            if ((h[weighted_pred_flag()] && h[slice_type()] == P) ||
                    (h[weighted_bipred_flag()] && h[slice_type()] == B))
                h(pred_weight_table());
            h(five_minus_max_num_merge_cand(), ue(v));
        }
        h(slice_qp_delta(), se(v));
        if (h[pps_slice_chroma_qp_offsets_present_flag()])
        {
            h(slice_cb_qp_offset(), se(v));
            h(slice_cr_qp_offset(), se(v));
        }
        if (h[chroma_qp_offset_list_enabled_flag()])
            h(cu_chroma_qp_offset_enabled_flag(), u(1));
        if (h[deblocking_filter_override_enabled_flag()])
            h(deblocking_filter_override_flag(), u(1));
        h[slice_deblocking_filter_disabled_flag()] = h[pps_deblocking_filter_disabled_flag()]; // review: infer this?
        if (h[deblocking_filter_override_flag()])
        {
            h(slice_deblocking_filter_disabled_flag(), u(1));
            if (!h[slice_deblocking_filter_disabled_flag()])
            {
                h(slice_beta_offset_div2(), se(v));
                h(slice_tc_offset_div2(), se(v));
            }
        }
        if (h[pps_loop_filter_across_slices_enabled_flag()] &&
                (h[slice_sao_luma_flag()] || h[slice_sao_chroma_flag()] ||
                        !h[slice_deblocking_filter_disabled_flag()]))
            h(slice_loop_filter_across_slices_enabled_flag(), u(1));
    }
    if (h[tiles_enabled_flag()] || h[entropy_coding_sync_enabled_flag()])
    {
        h(num_entry_point_offsets(), ue(v));
        if (h[num_entry_point_offsets()] > 0)
        {
            h(offset_len_minus1(), ue(v));
            for (int i = 0; i < h[num_entry_point_offsets()]; i++)
                h(entry_point_offset_minus1(i), u(h[offset_len_minus1()] + 1));
        }
    }
    if (h[slice_segment_header_extension_present_flag()])
    {
        h(slice_segment_header_extension_length(), ue(v));
        for (int i = 0; i < h[slice_segment_header_extension_length()]; i++)
        {
            h(slice_segment_header_extension_data_byte(i), u(8));
        }
    }
    h(byte_alignment());
}


template <class H>
void Syntax<ref_pic_lists_modification>::go(const ref_pic_lists_modification&, H &h)
{
    h(ref_pic_list_modification_flag_l0(), u(1));
    if (h[ref_pic_list_modification_flag_l0()])
        for (int i = 0; i <= h[num_ref_idx_l0_active_minus1()]; i++)
            h(list_entry_l0(i), u(ceilLog2(h[NumPocTotalCurr()])));
    if (h[slice_type()] == B)
    {
        h(ref_pic_list_modification_flag_l1(), u(1));
        if (h[ref_pic_list_modification_flag_l1()])
            for (int i = 0; i <= h[num_ref_idx_l1_active_minus1()]; i++)
                h(list_entry_l1(i), u(ceilLog2(h[NumPocTotalCurr()])));
    }
}


template <class H>
void Syntax<pred_weight_table>::go(const pred_weight_table&, H &h)
{
    h(luma_log2_weight_denom(), ue(v));
    if (h[chroma_format_idc()] != 0)
        h(delta_chroma_log2_weight_denom(), se(v));
    for (int i = 0; i <= h[num_ref_idx_l0_active_minus1()]; i++)
        h(luma_weight_l0_flag(i), u(1));
    if (h[chroma_format_idc()] != 0)
        for (int i = 0; i <= h[num_ref_idx_l0_active_minus1()]; i++)
            h(chroma_weight_l0_flag(i), u(1));
    for (int i = 0; i <= h[num_ref_idx_l0_active_minus1()]; i++)
    {
        if (h[luma_weight_l0_flag(i)])
        {
            h(delta_luma_weight_l0(i), se(v));
            h(luma_offset_l0(i), se(v));
        }
        if (h[chroma_weight_l0_flag(i)])
            for (int j = 0; j < 2; j++)
            {
                h(delta_chroma_weight_l0(i, j), se(v));
                h(delta_chroma_offset_l0(i, j), se(v));
            }
    }
    if (h[slice_type()] == B)
    {
        for (int i = 0; i <= h[num_ref_idx_l1_active_minus1()]; i++)
            h(luma_weight_l1_flag(i), u(1));
        if (h[chroma_format_idc()] != 0)
            for (int i = 0; i <= h[num_ref_idx_l1_active_minus1()]; i++)
                h(chroma_weight_l1_flag(i), u(1));
        for (int i = 0; i <= h[num_ref_idx_l1_active_minus1()]; i++)
        {
            if (h[luma_weight_l1_flag(i)])
            {
                h(delta_luma_weight_l1(i), se(v));
                h(luma_offset_l1(i), se(v));
            }
            if (h[chroma_weight_l1_flag(i)])
                for (int j = 0; j < 2; j++)
                {
                    h(delta_chroma_weight_l1(i, j), se(v));
                    h(delta_chroma_offset_l1(i, j), se(v));
                }
        }
    }
}


template <class H>
void Syntax<short_term_ref_pic_set>::go(const short_term_ref_pic_set &fun, H &h)
{
    if (fun.stRpsIdx != 0)
        h(inter_ref_pic_set_prediction_flag(), u(1));
    if (h[inter_ref_pic_set_prediction_flag()])
    {
        if (fun.stRpsIdx == h[num_short_term_ref_pic_sets()])
            h(delta_idx_minus1(), ue(v));
        h(delta_rps_sign(), u(1));
        h(abs_delta_rps_minus1(), ue(v));
        const int RefRpsIdx = fun.stRpsIdx - (h[delta_idx_minus1()] + 1);
        for (int j = 0; j <= h[NumDeltaPocs(RefRpsIdx)]; j++)
        {
            h(used_by_curr_pic_flag(j), u(1));
            if (!h[used_by_curr_pic_flag(j)])
                h(use_delta_flag(j), u(1));
            else
                h[use_delta_flag(j)] = 1;
        }
    }
    else
    {
        h(num_negative_pics(), ue(v));
        h(num_positive_pics(), ue(v));
        for (int i = 0; i < h[num_negative_pics()]; i++)
        {
            h(delta_poc_s0_minus1(i), ue(v));
            h(used_by_curr_pic_s0_flag(i), u(1));
        }
        for (int i = 0; i < h[num_positive_pics()]; i++)
        {
            h(delta_poc_s1_minus1(i), ue(v));
            h(used_by_curr_pic_s1_flag(i), u(1));
        }
    }
}


// Slice segment data syntax

template <class H>
void Syntax<slice_segment_data>::iteration(H &h)
{
    h(coding_tree_unit());
    h(end_of_slice_segment_flag(), ae(v));
    h[CtbAddrInTs()]++;
    h[CtbAddrInRs()] = h[CtbAddrTsToRs(h[CtbAddrInTs()])];
    if (!h[end_of_slice_segment_flag()] &&
            ((h[tiles_enabled_flag()] && h[TileId(h[CtbAddrInTs()])] != h[TileId(h[CtbAddrInTs()] - 1)]) ||
                    (h[entropy_coding_sync_enabled_flag()] && h[CtbAddrInTs()] % h[PicWidthInCtbsY()] == 0)))
    {
        h(end_of_subset_one_bit() /* equal to 1 */, ae(v));
        h(byte_alignment());
    }
}


template <class H>
void Syntax<slice_segment_data>::go(const slice_segment_data &, H &h)
{
    do
    {
        iteration(h);
    }
    while (!h[end_of_slice_segment_flag()]);
}


template <class H>
void Syntax<vui_parameters>::go(const vui_parameters &, H &h)
{
    static const int EXTENDED_SAR = 255;

    h(aspect_ratio_info_present_flag(), u(1));
    if (h[aspect_ratio_info_present_flag()])
    {
        h(aspect_ratio_idc(), u(8));
        if (h[aspect_ratio_idc()] == EXTENDED_SAR)
        {
            h(sar_width(), u(16));
            h(sar_height(), u(16));
        }
    }
    h(overscan_info_present_flag(), u(1));
    if (h[overscan_info_present_flag()])
    {
        h(overscan_appropriate_flag(), u(1));
    }
    h(video_signal_type_present_flag(), u(1));
    if (h[video_signal_type_present_flag()])
    {
        h(video_format(), u(3));
        h(video_full_range_flag(), u(1));
        h(colour_description_present_flag(), u(1));
        if (h[colour_description_present_flag()])
        {
            h(colour_primaries(), u(8));
            h(transfer_characteristics(), u(8));
            h(matrix_coeffs(), u(8));
        }
    }
    h(chroma_loc_info_present_flag(), u(1));
    if (h[chroma_loc_info_present_flag()])
    {
        h(chroma_sample_loc_type_top_field(), ue(v));
        h(chroma_sample_loc_type_bottom_field(), ue(v));
    }
    h(neutral_chroma_indication_flag(), u(1));
    h(field_seq_flag(), u(1));
    h(frame_field_info_present_flag(), u(1));
    h(default_display_window_flag(), u(1));
    if (h[default_display_window_flag()])
    {
        h(def_disp_win_left_offset(), ue(v));
        h(def_disp_win_right_offset(), ue(v));
        h(def_disp_win_top_offset(), ue(v));
        h(def_disp_win_bottom_offset(), ue(v));
    }
    h(vui_timing_info_present_flag(), u(1));
    if (h[vui_timing_info_present_flag()])
    {
        h(vui_num_units_in_tick(), u(32));
        h(vui_time_scale(), u(32));
        h(vui_poc_proportional_to_timing_flag(), u(1));
        if (h[vui_poc_proportional_to_timing_flag()])
        {
            h(vui_num_ticks_poc_diff_one_minus1(), ue(v));
        }
        h(vui_hrd_parameters_present_flag(), u(1));
        if (h[vui_hrd_parameters_present_flag()])
            h(hrd_parameters(1, h[sps_max_sub_layers_minus1()]));
    }
    h(bitstream_restriction_flag(), u(1));
    if (h[bitstream_restriction_flag()])
    {
        h(tiles_fixed_structure_flag(), u(1));
        h(motion_vectors_over_pic_boundaries_flag(), u(1));
        h(restricted_ref_pic_lists_flag(), u(1));
        h(min_spatial_segmentation_idc(), ue(v));
        h(max_bytes_per_pic_denom(), ue(v));
        h(max_bits_per_min_cu_denom(), ue(v));
        h(log2_max_mv_length_horizontal(), ue(v));
        h(log2_max_mv_length_vertical(), ue(v));
    }
}


template <class H>
void Syntax<hrd_parameters>::go(const hrd_parameters &fun, H &h)
{
    if (fun.commonInfPresentFlag)
    {
        h(nal_hrd_parameters_present_flag(), u(1));
        h(vcl_hrd_parameters_present_flag(), u(1));
        if (h[nal_hrd_parameters_present_flag()] || h[vcl_hrd_parameters_present_flag()])
        {
            h(sub_pic_hrd_params_present_flag(), u(1));
            if (h[sub_pic_hrd_params_present_flag()])
            {
                h(tick_divisor_minus2(), u(8));
                h(du_cpb_removal_delay_increment_length_minus1(), u(5));
                h(sub_pic_cpb_params_in_pic_timing_sei_flag(), u(1));
                h(dpb_output_delay_du_length_minus1(), u(5));
            }
            h(bit_rate_scale(), u(4));
            h(cpb_size_scale(), u(4));
            if (h[sub_pic_hrd_params_present_flag()])
                h(cpb_size_du_scale(), u(4));
            h(initial_cpb_removal_delay_length_minus1(), u(5));
            h(au_cpb_removal_delay_length_minus1(), u(5));
            h(dpb_output_delay_length_minus1(), u(5));
        }
    }
    for (int i = 0; i <= fun.maxNumSubLayersMinus1; i++)
    {
        h(fixed_pic_rate_general_flag(i), u(1));
        if (!h[fixed_pic_rate_general_flag(i)])
            h(fixed_pic_rate_within_cvs_flag(i), u(1));
        if (h[fixed_pic_rate_within_cvs_flag(i)])
            h(elemental_duration_in_tc_minus1(i), ue(v));
        else
            h(low_delay_hrd_flag(i), u(1));
        if (!h[low_delay_hrd_flag(i)])
        {
            h(cpb_cnt_minus1(i), ue(v));
        }
        if (h[nal_hrd_parameters_present_flag()])
            h(sub_layer_hrd_parameters(i));
        if (h[vcl_hrd_parameters_present_flag()])
            h(sub_layer_hrd_parameters(i));
    }
}


template <class H>
void Syntax<sub_layer_hrd_parameters>::go(const sub_layer_hrd_parameters &fun, H &h)
{
    for (int i = 0; i <= h[CpbCnt(fun.subLayerId)]; i++)
    {
        h(bit_rate_value_minus1(i), ue(v));
        h(cpb_size_value_minus1(i), ue(v));
        if (h[sub_pic_hrd_params_present_flag()])
        {
            h(cpb_size_du_value_minus1(i), ue(v));
            h(bit_rate_du_value_minus1(i), ue(v));
        }
        h(cbr_flag(i), u(1));
    }
}

template <class H> void Syntax<sei_message>::go(const sei_message &fun, H &h)
{
    int payloadType = 0;
    while (next_bits<int>(h, 8, true) == 0xFF)
    {
        h(ff_byte() /* equal to 0xFF */, f(8));
        payloadType += 255;
    }
    h(last_payload_type_byte(), u(8));
    payloadType += h[last_payload_type_byte()];
    int payloadSize = 0;
    while (next_bits<int>(h, 8, true) == 0xFF)
    {
        h(ff_byte() /* equal to 0xFF */, f(8));
        payloadSize += 255;
    }
    h(last_payload_size_byte(), u(8));
    payloadSize += h[last_payload_size_byte()];
    h(sei_payload(payloadType, payloadSize));
}


struct alternative_transfer_characteristics;
struct reserved_sei_message;
struct buffering_period;
struct pic_timing;
struct pan_scan_rect;
struct filler_payload;
struct user_data_registered_itu_t_t35;
struct user_data_unregistered;
struct recovery_point;
struct scene_info;
struct picture_snapshot;
struct progressive_refinement_segment_start;
struct progressive_refinement_segment_end;
struct film_grain_characteristics;
struct post_filter_hint;
struct tone_mapping_info;
struct frame_packing_arrangement;
struct display_orientation;
struct structure_of_pictures_info;
struct active_parameter_sets;
struct decoding_unit_info;
struct temporal_sub_layer_zero_index;
struct decoded_picture_hash;
struct scalable_nesting;
struct region_refresh_info;
struct no_display;
struct time_code;
struct mastering_display_colour_volume;
struct segmented_rect_frame_packing_arrangement;
struct temporal_motion_constrained_tile_sets;
struct chroma_resampling_filter_hint;
struct knee_function_info;
struct colour_remapping_info;
struct deinterlaced_field_identification;
struct content_light_level;
struct layers_not_present;
struct inter_layer_constrained_tile_sets;
struct bsp_nesting;
struct bsp_initial_arrival_time;
struct sub_bitstream_property;
struct alpha_channel_info;
struct overlay_info;
struct temporal_mv_prediction_constraints;
struct frame_field_info;
struct three_dimensional_reference_displays_info;
struct depth_representation_info;
struct multiview_scene_info;
struct multiview_acquisition_info;
struct multiview_view_position;

struct reserved_payload_extension_data;
struct payload_bit_equal_to_one;
struct payload_bit_equal_to_zero;


template <class H> void Syntax<sei_payload>::go(sei_payload fun, H &h)
{
    StateEncode *stateEncode = h;
    if (h[nal_unit_type()] == 39 /*PREFIX_SEI_NUT*/)
    {
        if (fun.payloadType == 0)
            h(buffering_period(fun.payloadSize));
        else if (fun.payloadType == 1)
            h(pic_timing(fun.payloadSize));
        else if (fun.payloadType == 2)
            h(pan_scan_rect(fun.payloadSize));
        else if (fun.payloadType == 3)
            h(filler_payload(fun.payloadSize));
        else if (fun.payloadType == 4)
            h(user_data_registered_itu_t_t35(fun.payloadSize));
        else if (fun.payloadType == 5)
        {
            fun.payloadSize = 16 + stateEncode->userDataUnregMsgLen;
            h(user_data_unregistered(fun.payloadSize));
        }
        else if (fun.payloadType == 6)
            h(recovery_point(fun.payloadSize));
        else if (fun.payloadType == 9)
            h(scene_info(fun.payloadSize));
        else if (fun.payloadType == 15)
            h(picture_snapshot(fun.payloadSize));
        else if (fun.payloadType == 16)
            h(progressive_refinement_segment_start(fun.payloadSize));
        else if (fun.payloadType == 17)
            h(progressive_refinement_segment_end(fun.payloadSize));
        else if (fun.payloadType == 19)
            h(film_grain_characteristics(fun.payloadSize));
        else if (fun.payloadType == 22)
            h(post_filter_hint(fun.payloadSize));
        else if (fun.payloadType == 23)
            h(tone_mapping_info(fun.payloadSize));
        else if (fun.payloadType == 45)
            h(frame_packing_arrangement(fun.payloadSize));
        else if (fun.payloadType == 47)
            h(display_orientation(fun.payloadSize));
        else if (fun.payloadType == 128)
            h(structure_of_pictures_info(fun.payloadSize));
        else if (fun.payloadType == 129)
            h(active_parameter_sets(fun.payloadSize));
        else if (fun.payloadType == 130)
            h(decoding_unit_info(fun.payloadSize));
        else if (fun.payloadType == 131)
            h(temporal_sub_layer_zero_index(fun.payloadSize));
        else if (fun.payloadType == 133)
            h(scalable_nesting(fun.payloadSize));
        else if (fun.payloadType == 134)
            h(region_refresh_info(fun.payloadSize));
        else if (fun.payloadType == 135)
            h(no_display(fun.payloadSize));
        else if (fun.payloadType == 136)
            h(time_code(fun.payloadSize));
        else if (fun.payloadType == 137)
            h(mastering_display_colour_volume(fun.payloadSize));
        else if (fun.payloadType == 138)
            h(segmented_rect_frame_packing_arrangement(fun.payloadSize));
        else if (fun.payloadType == 139)
            h(temporal_motion_constrained_tile_sets(fun.payloadSize));
        else if (fun.payloadType == 140)
            h(chroma_resampling_filter_hint(fun.payloadSize));
        else if (fun.payloadType == 141)
            h(knee_function_info(fun.payloadSize));
        else if (fun.payloadType == 142)
            h(colour_remapping_info(fun.payloadSize));
        else if (fun.payloadType == 143)
            h(deinterlaced_field_identification(fun.payloadSize));
        else if (fun.payloadType == 144)
            h(content_light_level(fun.payloadSize));
        else if (fun.payloadType == 147)
            h(alternative_transfer_characteristics(fun.payloadSize));
        else if (fun.payloadType == 160)
            h(layers_not_present(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 161)
            h(inter_layer_constrained_tile_sets(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 162)
            h(bsp_nesting(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 163)
            h(bsp_initial_arrival_time(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 164)
            h(sub_bitstream_property(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 165)
            h(alpha_channel_info(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 166)
            h(overlay_info(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 167)
            h(temporal_mv_prediction_constraints(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 168)
            h(frame_field_info(fun.payloadSize)); /* specified in Annex F */
        else if (fun.payloadType == 176)
            h(three_dimensional_reference_displays_info(fun.payloadSize)); /* specified in Annex G */
        else if (fun.payloadType == 177)
            h(depth_representation_info(fun.payloadSize)); /* specified in Annex G */
        else if (fun.payloadType == 178)
            h(multiview_scene_info(fun.payloadSize)); /* specified in Annex G */
        else if (fun.payloadType == 179)
            h(multiview_acquisition_info(fun.payloadSize)); /* specified in Annex G */
        else if (fun.payloadType == 180)
            h(multiview_view_position(fun.payloadSize)); /* specified in Annex G */
        else
            h(reserved_sei_message(fun.payloadSize));
    }
    else /* nal_unit_type == SUFFIX_SEI_NUT */
    {
        if (fun.payloadType == 3)
            h(filler_payload(fun.payloadSize));
        else if (fun.payloadType == 4)
            h(user_data_registered_itu_t_t35(fun.payloadSize));
        else if (fun.payloadType == 5)
            h(user_data_unregistered(fun.payloadSize));
        else if (fun.payloadType == 17)
            h(progressive_refinement_segment_end(fun.payloadSize));
        else if (fun.payloadType == 22)
            h(post_filter_hint(fun.payloadSize));
        else if (fun.payloadType == 132)
            h(decoded_picture_hash(fun.payloadSize));
        else
            h(reserved_sei_message(fun.payloadSize));
    }

    if (h[more_data_in_payload()])
    {
        if (h[payload_extension_present()])
            h(reserved_payload_extension_data(), uv());
        h(payload_bit_equal_to_one(), f(1));
        while (!h[byte_aligned()])
            h(payload_bit_equal_to_zero(), f(1));
    }
}

#endif
