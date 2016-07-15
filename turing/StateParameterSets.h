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

#ifndef INCLUDED_StateParameterSets_h
#define INCLUDED_StateParameterSets_h

#pragma once

// Constructs to store and manipulate the state conveyed in VPS, SPS, PPS and slice_segment_header.

#include "StateValues.h"
#include "Global.h"
#include <boost/uuid/uuid.hpp>


struct Hrd :
    AccessOperators<Hrd>,
    ValueHolder<nal_hrd_parameters_present_flag>,
    ValueHolder<vcl_hrd_parameters_present_flag>,
    ValueHolder<sub_pic_hrd_params_present_flag>,
    ValueHolder<tick_divisor_minus2>,
    ValueHolder<du_cpb_removal_delay_increment_length_minus1>,
    ValueHolder<sub_pic_cpb_params_in_pic_timing_sei_flag>,
    ValueHolder<dpb_output_delay_du_length_minus1>,
    ValueHolder<bit_rate_scale>,
    ValueHolder<cpb_size_scale>,
    ValueHolder<cpb_size_du_scale>,
    ValueHolder<initial_cpb_removal_delay_length_minus1>,
    ValueHolder<au_cpb_removal_delay_length_minus1>,
    ValueHolder<dpb_output_delay_length_minus1>,
    ValueHolder<fixed_pic_rate_general_flag>,
    ValueHolder<fixed_pic_rate_within_cvs_flag>,
    ValueHolder<elemental_duration_in_tc_minus1>,
    ValueHolder<low_delay_hrd_flag>,
    ValueHolder<cpb_cnt_minus1>
    {
        // sub_layer_hrd_parameters
        struct SubLayer :
            ValueHolder<bit_rate_value_minus1>,
            ValueHolder<cpb_size_value_minus1>,
            ValueHolder<cpb_size_du_value_minus1>,
            ValueHolder<bit_rate_du_value_minus1>,
            ValueHolder<cbr_flag>
            {
            };

            std::vector<SubLayer> sublayers;
    };


struct HrdArray
{
    std::vector<Hrd> hrd;
};

template <class S> struct Active;
struct Sps;

template <class H>
Hrd *getHrd(H &h)
{
    auto activeSps = h[Active<Sps>()];

    if (activeSps)
    {
        if (activeSps->hrdArray.hrd.empty())
        {
            return 0;
        }
        else
        {
            return &activeSps->hrdArray.hrd.back();
        }
    }
    return 0;
}


struct ProfileTierLevel  :
    ValueHolder<general_profile_space>,
    ValueHolder<general_tier_flag>,
    ValueHolder<general_profile_idc>,
    ValueHolder<general_profile_compatibility_flag>,
    ValueHolder<general_progressive_source_flag>,
    ValueHolder<general_interlaced_source_flag>,
    ValueHolder<general_non_packed_constraint_flag>,
    ValueHolder<general_frame_only_constraint_flag>,
    ValueHolder<general_max_12bit_constraint_flag>,
    ValueHolder<general_max_10bit_constraint_flag>,
    ValueHolder<general_max_8bit_constraint_flag>,
    ValueHolder<general_max_422chroma_constraint_flag>,
    ValueHolder<general_max_420chroma_constraint_flag>,
    ValueHolder<general_max_monochrome_constraint_flag>,
    ValueHolder<general_intra_constraint_flag>,
    ValueHolder<general_one_picture_only_constraint_flag>,
    ValueHolder<general_lower_bit_rate_constraint_flag>,
    ValueHolder<general_reserved_zero_34bits>,
    ValueHolder<general_reserved_zero_43bits>,
    ValueHolder<general_inbld_flag>,
    ValueHolder<general_reserved_zero_bit>,
    ValueHolder<general_level_idc >,
    ValueHolder<sub_layer_profile_present_flag>,
    ValueHolder<sub_layer_level_present_flag>,
    ValueHolder<reserved_zero_2bits>,
    ValueHolder<sub_layer_profile_space>,
    ValueHolder<sub_layer_tier_flag>,
    ValueHolder<sub_layer_profile_idc>,
    ValueHolder<sub_layer_profile_compatibility_flag>,
    ValueHolder<sub_layer_progressive_source_flag>,
    ValueHolder<sub_layer_interlaced_source_flag>,
    ValueHolder<sub_layer_non_packed_constraint_flag>,
    ValueHolder<sub_layer_frame_only_constraint_flag>,
    ValueHolder<sub_layer_max_12bit_constraint_flag>,
    ValueHolder<sub_layer_max_10bit_constraint_flag>,
    ValueHolder<sub_layer_max_8bit_constraint_flag>,
    ValueHolder<sub_layer_max_422chroma_constraint_flag>,
    ValueHolder<sub_layer_max_420chroma_constraint_flag>,
    ValueHolder<sub_layer_max_monochrome_constraint_flag>,
    ValueHolder<sub_layer_intra_constraint_flag>,
    ValueHolder<sub_layer_one_picture_only_constraint_flag>,
    ValueHolder<sub_layer_lower_bit_rate_constraint_flag>,
    ValueHolder<sub_layer_reserved_zero_34bits>,
    ValueHolder<sub_layer_reserved_zero_43bits>,
    ValueHolder<sub_layer_inbld_flag>,
    ValueHolder<sub_layer_reserved_zero_bit>,
    ValueHolder<sub_layer_level_idc>
    {
    };


struct Vps :
    ValueHolder<vps_video_parameter_set_id>,
    ValueHolder<vps_reserved_three_2bits>,
    ValueHolder<vps_max_layers_minus1>,
    ValueHolder<vps_max_sub_layers_minus1>,
    ValueHolder<vps_temporal_id_nesting_flag>,
    ValueHolder<vps_reserved_0xffff_16bits>,
    ValueHolder<vps_sub_layer_ordering_info_present_flag>,
    ValueHolder<vps_max_dec_pic_buffering_minus1>,
    ValueHolder<vps_max_num_reorder_pics>,
    ValueHolder<vps_max_latency_increase_plus1>,
    ValueHolder<vps_max_layer_id>,
    ValueHolder<vps_num_layer_sets_minus1>,
    ValueHolder<layer_id_included_flag>,
    ValueHolder<vps_timing_info_present_flag>,
    ValueHolder<vps_num_units_in_tick>,
    ValueHolder<vps_time_scale>,
    ValueHolder<vps_poc_proportional_to_timing_flag>,
    ValueHolder<vps_num_ticks_poc_diff_one_minus1>,
    ValueHolder<hrd_layer_set_idx>,
    ValueHolder<vps_extension_flag>,
    ValueHolder<vps_extension_data_flag>
    {
        std::vector<int> cprms_present_flag;
        HrdArray hrdArray;
        ProfileTierLevel ptl;
    };


template <class ParameterSetData>
struct Table { };


template <class ParameterSetData>
struct ValueType<Table<ParameterSetData>>
{
    typedef std::map<int, std::shared_ptr<ParameterSetData>> Type;
};


template <class ParameterSetData>
struct Active { };


template <class ParameterSetData>
struct ValueType<Active<ParameterSetData>>
{
    typedef std::shared_ptr<ParameterSetData> Type;
};


struct ScalingListState :
    ValueHolder<scaling_list_pred_mode_flag>,
    ValueHolder<scaling_list_pred_matrix_id_delta>,
    ValueHolder<scaling_list_dc_coef_minus8>,
    ValueHolder<scaling_list_delta_coef>,
    ValueHolder<ScalingList>,
    AccessOperators<ScalingListState>
    {
    };


struct Strps :
    ValueHolder<inter_ref_pic_set_prediction_flag>,
    ValueHolder<delta_idx_minus1>,
    ValueHolder<delta_rps_sign>,
    ValueHolder<abs_delta_rps_minus1>,
    ValueHolder<used_by_curr_pic_flag>,
    ValueHolder<use_delta_flag>,
    ValueHolder<num_negative_pics>,
    ValueHolder<num_positive_pics>,
    ValueHolder<delta_poc_s0_minus1>,
    ValueHolder<used_by_curr_pic_s0_flag>,
    ValueHolder<delta_poc_s1_minus1>,
    ValueHolder<used_by_curr_pic_s1_flag>
    {
    };


struct VuiParameters :
    ValueHolder<aspect_ratio_info_present_flag>,
    ValueHolder<aspect_ratio_idc>,
    ValueHolder<sar_width>,
    ValueHolder<sar_height>,
    ValueHolder<overscan_info_present_flag>,
    ValueHolder<overscan_appropriate_flag>,
    ValueHolder<video_signal_type_present_flag>,
    ValueHolder<video_format>,
    ValueHolder<video_full_range_flag>,
    ValueHolder<colour_description_present_flag>,
    ValueHolder<colour_primaries>,
    ValueHolder<transfer_characteristics>,
    ValueHolder<matrix_coeffs>,
    ValueHolder<chroma_loc_info_present_flag>,
    ValueHolder<chroma_sample_loc_type_top_field>,
    ValueHolder<chroma_sample_loc_type_bottom_field>,
    ValueHolder<neutral_chroma_indication_flag>,
    ValueHolder<field_seq_flag>,
    ValueHolder<frame_field_info_present_flag>,
    ValueHolder<default_display_window_flag>,
    ValueHolder<def_disp_win_left_offset>,
    ValueHolder<def_disp_win_right_offset>,
    ValueHolder<def_disp_win_top_offset>,
    ValueHolder<def_disp_win_bottom_offset>,
    ValueHolder<vui_timing_info_present_flag>,
    ValueHolder<vui_num_units_in_tick>,
    ValueHolder<vui_time_scale>,
    ValueHolder<vui_poc_proportional_to_timing_flag>,
    ValueHolder<vui_num_ticks_poc_diff_one_minus1>,
    ValueHolder<vui_hrd_parameters_present_flag>,
    ValueHolder<bitstream_restriction_flag>,
    ValueHolder<tiles_fixed_structure_flag>,
    ValueHolder<motion_vectors_over_pic_boundaries_flag>,
    ValueHolder<restricted_ref_pic_lists_flag>,
    ValueHolder<min_spatial_segmentation_idc>,
    ValueHolder<max_bytes_per_pic_denom>,
    ValueHolder<max_bits_per_min_cu_denom>,
    ValueHolder<log2_max_mv_length_horizontal>,
    ValueHolder<log2_max_mv_length_vertical>
    {
        HrdArray hrdArray;
    };


struct Sps :
    ValueHolder<sps_video_parameter_set_id>,
    ValueHolder<sps_max_sub_layers_minus1>,
    ValueHolder<sps_ext_or_max_sub_layers_minus1>,
    ValueHolder<MultiLayerExtSpsFlag>,
    ValueHolder<sps_temporal_id_nesting_flag>,
    ValueHolder<sps_seq_parameter_set_id>,
    ValueHolder<chroma_format_idc>,
    ValueHolder<separate_colour_plane_flag>,
    ValueHolder<pic_width_in_luma_samples>,
    ValueHolder<pic_height_in_luma_samples>,
    ValueHolder<conformance_window_flag>,
    ValueHolder<conf_win_left_offset>,
    ValueHolder<conf_win_right_offset>,
    ValueHolder<conf_win_top_offset>,
    ValueHolder<conf_win_bottom_offset>,
    ValueHolder<bit_depth_luma_minus8>,
    ValueHolder<bit_depth_chroma_minus8>,
    ValueHolder<log2_max_pic_order_cnt_lsb_minus4>,
    ValueHolder<sps_sub_layer_ordering_info_present_flag>,
    ValueHolder<sps_max_dec_pic_buffering_minus1>,
    ValueHolder<sps_max_num_reorder_pics>,
    ValueHolder<sps_max_latency_increase_plus1>,
    ValueHolder<log2_min_luma_coding_block_size_minus3>,
    ValueHolder<log2_diff_max_min_luma_coding_block_size>,
    ValueHolder<log2_min_luma_transform_block_size_minus2>,
    ValueHolder<log2_diff_max_min_luma_transform_block_size>,
    ValueHolder<max_transform_hierarchy_depth_inter>,
    ValueHolder<max_transform_hierarchy_depth_intra>,
    ValueHolder<scaling_list_enabled_flag>,
    ValueHolder<sps_scaling_list_data_present_flag>,
    ValueHolder<amp_enabled_flag>,
    ValueHolder<sample_adaptive_offset_enabled_flag>,
    ValueHolder<pcm_enabled_flag>,
    ValueHolder<pcm_sample_bit_depth_luma_minus1>,
    ValueHolder<pcm_sample_bit_depth_chroma_minus1>,
    ValueHolder<log2_min_pcm_luma_coding_block_size_minus3>,
    ValueHolder<log2_diff_max_min_pcm_luma_coding_block_size>,
    ValueHolder<pcm_loop_filter_disabled_flag>,
    ValueHolder<num_short_term_ref_pic_sets>,
    ValueHolder<long_term_ref_pics_present_flag>,
    ValueHolder<num_long_term_ref_pics_sps>,
    ValueHolder<lt_ref_pic_poc_lsb_sps>,
    ValueHolder<used_by_curr_pic_lt_sps_flag>,
    ValueHolder<sps_temporal_mvp_enabled_flag>,
    ValueHolder<strong_intra_smoothing_enabled_flag>,
    ValueHolder<vui_parameters_present_flag>,
    ValueHolder<sps_extension_present_flag>,
    ValueHolder<sps_range_extension_flag>,
    ValueHolder<sps_multilayer_extension_flag>,
    ValueHolder<sps_extension_6bits>,
    ValueHolder<sps_extension_data_flag>,
    ValueHolder<transform_skip_rotation_enabled_flag>,
    ValueHolder<transform_skip_context_enabled_flag>,
    ValueHolder<implicit_rdpcm_enabled_flag>,
    ValueHolder<explicit_rdpcm_enabled_flag>,
    ValueHolder<extended_precision_processing_flag>,
    ValueHolder<intra_smoothing_disabled_flag>,
    ValueHolder<high_precision_offsets_enabled_flag>,
    ValueHolder<persistent_rice_adaptation_enabled_flag>,
    ValueHolder<cabac_bypass_alignment_enabled_flag>,
    ValueHolder<NumNegativePics>,
    ValueHolder<NumPositivePics>,
    ValueHolder<UsedByCurrPicS0>,
    ValueHolder<UsedByCurrPicS1>,
    ValueHolder<DeltaPocS0>,
    ValueHolder<DeltaPocS1>,
    ValueHolder<PocLsbLt>,
    ValueHolder<UsedByCurrPicLt>,
    VuiParameters
    {
        ScalingListState scalingListState;
        ProfileTierLevel ptl;
    };


struct Pps :
    ValueHolder<pps_pic_parameter_set_id>,
    ValueHolder<pps_seq_parameter_set_id>,
    ValueHolder<dependent_slice_segments_enabled_flag>,
    ValueHolder<output_flag_present_flag>,
    ValueHolder<num_extra_slice_header_bits>,
    ValueHolder<sign_data_hiding_enabled_flag>,
    ValueHolder<cabac_init_present_flag>,
    ValueHolder<num_ref_idx_l0_default_active_minus1>,
    ValueHolder<num_ref_idx_l1_default_active_minus1>,
    ValueHolder<init_qp_minus26>,
    ValueHolder<constrained_intra_pred_flag>,
    ValueHolder<transform_skip_enabled_flag>,
    ValueHolder<cu_qp_delta_enabled_flag>,
    ValueHolder<diff_cu_qp_delta_depth>,
    ValueHolder<pps_cb_qp_offset>,
    ValueHolder<pps_cr_qp_offset>,
    ValueHolder<pps_slice_chroma_qp_offsets_present_flag>,
    ValueHolder<weighted_pred_flag>,
    ValueHolder<weighted_bipred_flag>,
    ValueHolder<transquant_bypass_enabled_flag>,
    ValueHolder<tiles_enabled_flag>,
    ValueHolder<entropy_coding_sync_enabled_flag>,
    ValueHolder<uniform_spacing_flag>,
    ValueHolder<loop_filter_across_tiles_enabled_flag>,
    ValueHolder<pps_loop_filter_across_slices_enabled_flag>,
    ValueHolder<deblocking_filter_control_present_flag>,
    ValueHolder<deblocking_filter_override_enabled_flag>,
    ValueHolder<pps_deblocking_filter_disabled_flag>,
    ValueHolder<pps_beta_offset_div2>,
    ValueHolder<pps_tc_offset_div2>,
    ValueHolder<pps_scaling_list_data_present_flag>,
    ValueHolder<lists_modification_present_flag>,
    ValueHolder<log2_parallel_merge_level_minus2>,
    ValueHolder<slice_segment_header_extension_present_flag>,
    ValueHolder<log2_max_transform_skip_block_size_minus2>,
    ValueHolder<cross_component_prediction_enabled_flag>,
    ValueHolder<chroma_qp_offset_list_enabled_flag>,
    ValueHolder<diff_cu_chroma_qp_offset_depth>,
    ValueHolder<chroma_qp_offset_list_len_minus1>,
    ValueHolder<cb_qp_offset_list>,
    ValueHolder<cr_qp_offset_list>,
    ValueHolder<log2_sao_offset_scale_luma>,
    ValueHolder<log2_sao_offset_scale_chroma>,
    ValueHolder<pps_extension_present_flag>,
    ValueHolder<pps_range_extension_flag>,
    ValueHolder<pps_multilayer_extension_flag>,
    ValueHolder<pps_extension_6bits>,
    ValueHolder<pps_extension_data_flag>,
    AccessOperators<Pps>
    {
        std::vector<int> columnWidthsMinus1;
        std::vector<int> rowHeightsMinus1;
        ScalingListState scalingListState;
    };


template <> struct Access<num_tile_columns_minus1, Pps>
{
    typedef int Type;
    static Type get(num_tile_columns_minus1, Pps &s)
    {
        return Type(s.columnWidthsMinus1.size());
    }
    static void set(num_tile_columns_minus1, Type i, Pps &s)
    {
        s.columnWidthsMinus1.resize(i);
    }
};


template <> struct Access<vps_num_hrd_parameters, Vps>
{
    typedef int Type;
    static Type get(vps_num_hrd_parameters, Vps &s)
    {
        return Type(s.cprms_present_flag.size());
    }
    static void set(vps_num_hrd_parameters, Type i, Vps & s)
    {
        s.cprms_present_flag.resize(i);
        if (i > 0)
        {
            s.cprms_present_flag[0] = 1;
        }
    }
};

template <> struct Access<cprms_present_flag, Vps> :
ValueType<cprms_present_flag>
{
    static Type get(cprms_present_flag v, Vps &s)
    {
        return Type(s.cprms_present_flag[v.i]);
    }
    static void set(cprms_present_flag v, Type x, Vps &s)
    {
        s.cprms_present_flag[v.i] = x;
    }
};

template <> struct Access<num_tile_rows_minus1, Pps>
{
    typedef int Type;
    static Type get(num_tile_rows_minus1, Pps &s)
    {
        return Type(s.rowHeightsMinus1.size());
    }
    static void set(num_tile_rows_minus1, Type i, Pps &s)
    {
        s.rowHeightsMinus1.resize(i);
    }
};


template <> struct Access<column_width_minus1, Pps>
{
    typedef int Type;
    static Type get(column_width_minus1 v, Pps &s)
    {
        return s.columnWidthsMinus1[v.i];
    }
    static void set(column_width_minus1 v, Type i, Pps &s)
    {
        s.columnWidthsMinus1[v.i] = i;
    }
};

template <> struct Access<row_height_minus1, Pps>
{
    typedef int Type;
    static Type get(row_height_minus1 v, Pps &s)
    {
        return s.rowHeightsMinus1[v.i];
    }
    static void set(row_height_minus1 v, Type i, Pps &s)
    {
        s.rowHeightsMinus1[v.i] = i;
    }
};


struct PredWeightTable :
    ValueHolder<luma_log2_weight_denom>,
    ValueHolder<delta_chroma_log2_weight_denom>,
    ValueHolder<luma_weight_l0_flag>,
    ValueHolder<chroma_weight_l0_flag>,
    ValueHolder<delta_luma_weight_l0>,
    ValueHolder<luma_offset_l0>,
    ValueHolder<delta_chroma_weight_l0>,
    ValueHolder<delta_chroma_offset_l0>,
    ValueHolder<luma_weight_l1_flag>,
    ValueHolder<chroma_weight_l1_flag>,
    ValueHolder<delta_luma_weight_l1>,
    ValueHolder<luma_offset_l1>,
    ValueHolder<delta_chroma_weight_l1>,
    ValueHolder<delta_chroma_offset_l1>
    {
    };

// These items may have different values if a picture has multiple slice_segment_header() instances
struct SliceSegmentHeaderIndependent :
    ValueHolder<first_slice_segment_in_pic_flag>,
    ValueHolder<no_output_of_prior_pics_flag>,
    ValueHolder<slice_pic_parameter_set_id>,
    ValueHolder<dependent_slice_segment_flag>,
    ValueHolder<slice_segment_address>,
    ValueHolder<offset_len_minus1>,
    ValueHolder<slice_segment_header_extension_length>,
    ValueHolder<slice_segment_header_extension_data_byte>
    {
        std::vector<size_t> entry_point_offset_minus1;

        // Populate slice segment header's entry point list according to the actual size of encoded substreams.
        // Substreams is a container type, each element of which has a size() member function, e.g. vector<vector<uint8_t>>
        template <class Substreams>
        void populateEntryPoints(const Substreams &substreams)
        {
            auto &offset_len_minus1 = this->ValueHolder<::offset_len_minus1>::value;

            assert(!substreams.empty());

            // Populate slice segment header's entry point list and set offset_len_minus1 accordingly.
            for (size_t i = 0; i < substreams.size() - 1; ++i)
            {
                assert(!substreams[i].empty());

                this->entry_point_offset_minus1.push_back(substreams[i].size() - 1);

                while (this->entry_point_offset_minus1.back() >= (size_t(1) << (offset_len_minus1 + 1)))
                {
                    ++offset_len_minus1;
                }
            }
        }
    };


template <class H>
struct Access<num_entry_point_offsets, H, typename std::enable_if<std::is_base_of<SliceSegmentHeaderIndependent, H>::value>::type>
{
    typedef int Type;
    static Type get(num_entry_point_offsets, SliceSegmentHeaderIndependent &s)
    {
        return static_cast<Type>(s.entry_point_offset_minus1.size());
    }
    static void set(num_entry_point_offsets, Type i, SliceSegmentHeaderIndependent &s)
    {
        s.entry_point_offset_minus1.resize(i);
    }
};


template <class H>
struct Access<entry_point_offset_minus1, H, typename std::enable_if<std::is_base_of<SliceSegmentHeaderIndependent, H>::value>::type>
{
    typedef int Type;
    static Type get(entry_point_offset_minus1 e, SliceSegmentHeaderIndependent &s)
    {
        return static_cast<Type>(s.entry_point_offset_minus1[e.i]);
    }
    static void set(entry_point_offset_minus1 e, Type i, SliceSegmentHeaderIndependent &s)
    {
        s.entry_point_offset_minus1[e.i] = i;
    }
};


struct Rplm :
    ValueHolder<ref_pic_list_modification_flag_l0>,
    ValueHolder<list_entry_l0>,
    ValueHolder<ref_pic_list_modification_flag_l1>,
    ValueHolder<list_entry_l1>
    {
    };


// These items have the same value in each slice_segment_header() of any given picture
struct SliceSegmentHeaderDependent :
    ValueHolder<pic_output_flag>,
    ValueHolder<slice_pic_order_cnt_lsb>,
    ValueHolder<short_term_ref_pic_set_sps_flag>,
    ValueHolder<short_term_ref_pic_set_idx>,
    ValueHolder<num_long_term_sps>,
    ValueHolder<num_long_term_pics>,
    ValueHolder<slice_temporal_mvp_enabled_flag>,
    ValueHolder<lt_idx_sps>,
    ValueHolder<poc_lsb_lt>,
    ValueHolder<used_by_curr_pic_lt_flag>,
    ValueHolder<delta_poc_msb_present_flag>,
    ValueHolder<delta_poc_msb_cycle_lt>,
    ValueHolder<slice_reserved_flag>,
    ValueHolder<slice_type>,
    ValueHolder<colour_plane_id>,
    ValueHolder<slice_sao_luma_flag>,
    ValueHolder<slice_sao_chroma_flag>,
    ValueHolder<num_ref_idx_active_override_flag>,
    ValueHolder<num_ref_idx_l0_active_minus1>,
    ValueHolder<num_ref_idx_l1_active_minus1>,
    ValueHolder<mvd_l1_zero_flag>,
    ValueHolder<cabac_init_flag>,
    ValueHolder<collocated_from_l0_flag>,
    ValueHolder<collocated_ref_idx>,
    ValueHolder<five_minus_max_num_merge_cand>,
    ValueHolder<slice_qp_delta>,
    ValueHolder<slice_cb_qp_offset>,
    ValueHolder<slice_cr_qp_offset>,
    ValueHolder<cu_chroma_qp_offset_enabled_flag>,
    ValueHolder<deblocking_filter_override_flag>,
    ValueHolder<slice_deblocking_filter_disabled_flag>,
    ValueHolder<slice_beta_offset_div2>,
    ValueHolder<slice_tc_offset_div2>,
    ValueHolder<slice_loop_filter_across_slices_enabled_flag>,
    Rplm,
    PredWeightTable
    {
    };

struct SliceSegmentHeader :
    SliceSegmentHeaderIndependent,
    SliceSegmentHeaderDependent,
    ValueHolder<alignment_bit_equal_to_one>,
    ValueHolder<alignment_bit_equal_to_zero>,
    AccessOperators<SliceSegmentHeader>
    {
    };

struct FirstSliceSegmentConstants { };

template <>
struct ValueType<FirstSliceSegmentConstants>
{
    typedef SliceSegmentHeader Type;
};


struct bp_seq_parameter_set_id { };
struct irap_cpb_params_present_flag { };
struct cpb_delay_offset { };
struct dpb_delay_offset { };
struct concatenation_flag { };
struct au_cpb_removal_delay_delta_minus1 { };
DEFINE_VALUE_ARRAY_1(nal_initial_cpb_removal_delay, i, 32)
DEFINE_VALUE_ARRAY_1(nal_initial_cpb_removal_offset, i, 32)
DEFINE_VALUE_ARRAY_1(nal_initial_alt_cpb_removal_delay, i, 32)
DEFINE_VALUE_ARRAY_1(nal_initial_alt_cpb_removal_offset, i, 32)
DEFINE_VALUE_ARRAY_1(vcl_initial_cpb_removal_delay, i, 32)
DEFINE_VALUE_ARRAY_1(vcl_initial_cpb_removal_offset, i, 32)
DEFINE_VALUE_ARRAY_1(vcl_initial_alt_cpb_removal_delay, i, 32)
DEFINE_VALUE_ARRAY_1(vcl_initial_alt_cpb_removal_offset, i, 32)


struct BufferingPeriod :
    ValueHolder<bp_seq_parameter_set_id>,
    ValueHolder<irap_cpb_params_present_flag>,
    ValueHolder<cpb_delay_offset>,
    ValueHolder<dpb_delay_offset>,
    ValueHolder<concatenation_flag>,
    ValueHolder<au_cpb_removal_delay_delta_minus1>,
    ValueHolder<nal_initial_cpb_removal_delay>,
    ValueHolder<nal_initial_cpb_removal_offset>,
    ValueHolder<nal_initial_alt_cpb_removal_delay>,
    ValueHolder<nal_initial_alt_cpb_removal_offset>,
    ValueHolder<vcl_initial_cpb_removal_delay>,
    ValueHolder<vcl_initial_cpb_removal_offset>,
    ValueHolder<vcl_initial_alt_cpb_removal_delay>,
    ValueHolder<vcl_initial_alt_cpb_removal_offset>
    {
    };


struct StateParameterSets :
    ValueHolder<Table<Vps>>,
    ValueHolder<Table<Sps>>,
    ValueHolder<Table<Pps>>
    {
        // review: somewhere we need a dedicated state structure containing SEI message data and managing its persistence
        std::shared_ptr<struct MasteringDisplayColourVolume> masteringDisplayColourVolume;
        std::map<boost::uuids::uuid, std::shared_ptr<std::vector<uint8_t>>> unregisteredUserData;
    };


struct ActiveParameterSets :
    ValueHolder<Active<Vps>>,
    ValueHolder<Active<Sps>>,
    ValueHolder<Active<Pps>>,
    ValueHolder<Active<BufferingPeriod>>
    {
    };


template <class V, class S>
struct Access<V, S, typename  std::enable_if<std::is_base_of<ValueHolder<Active<Vps>>, S>::value && Accessible<V, Vps>::value>::type>
:
Access<V, Vps>
{
    static typename Access<V, Vps>::Type get(V v, ValueHolder<Active<Vps>> &s)
    {
        return Access<V, Vps>::get(v, *s.get(Active<Vps>()));
    }
};


template <class V, class S>
struct Access<V, S, typename  std::enable_if<std::is_base_of<ValueHolder<Active<Sps>>, S>::value && Accessible<V, Sps>::value>::type>
:
Access<V, Sps>
{
    static typename Access<V, Sps>::Type get(V v, ValueHolder<Active<Sps>> &s)
    {
        return Access<V, Sps>::get(v, *s.get(Active<Sps>()));
    }
};


template <class V, class S>
struct Access<V, S, typename  std::enable_if<std::is_base_of<ValueHolder<Active<Pps>>, S>::value && Accessible<V, Pps>::value>::type> :
Access<V, Pps>
{
    static typename Access<V, Pps>::Type get(V v, ValueHolder<Active<Pps>> &s)
    {
        return Access<V, Pps>::get(v, *s.get(Active<Pps>()));
    }
};


template <class V, class S>
struct Access<V, S, typename std::enable_if<std::is_base_of<ActiveParameterSets, S>::value && Accessible<V, BufferingPeriod>::value>::type> :
Access<V, BufferingPeriod>
{
    static typename Access<V, BufferingPeriod>::Type get(V v, ValueHolder<Active<BufferingPeriod>> &s)
    {
        return Access<V, BufferingPeriod>::get(v, *s.get(Active<BufferingPeriod>()));
    }
};


template <class H>
ScalingListState *getScalingListState(H &h)
{
    if (h[pps_scaling_list_data_present_flag()]) return &h[Active<Pps>()]->scalingListState;
    if (h[sps_scaling_list_data_present_flag()]) return &h[Active<Sps>()]->scalingListState;
    return 0;
}

#endif
