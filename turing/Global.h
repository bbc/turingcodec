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

// review: remove this header - organise this code better into smaller files

#ifndef INCLUDED_Global_h
#define INCLUDED_Global_h

#pragma once

#include "StructMacros.h"
#include "StateValues.h"
#include "Handlers.h"
#include "HevcTypes.h"
#include "HevcMath.h"
#include "MotionVector.h"
#include <cstdint>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <algorithm>
#include <array>
#include <cassert>


#if defined(WIN32)  && 0
// This is an assert that fires for Visual Studio Release configurations. Disable for best performance.
#define ASSERT(a) if (!(a)) { std::cout << #a << "\n"; __debugbreak(); }
#else
#define ASSERT(a) assert(a)
#endif


template <class V> struct Fixed;

template <typename F> struct Syntax;


struct Bitstream
{
    Bitstream(wchar_t const*name=0) : name(name) { }
    wchar_t const *name;
};


struct CodedVideoSequence { };

struct AccessUnit { };



DEFINE_STRUCT_ARITY_0(more_data_in_byte_stream);
DEFINE_STRUCT_ARITY_0(byte_aligned)
DEFINE_STRUCT_ARITY_0(more_data_in_payload)
DEFINE_STRUCT_ARITY_0(more_rbsp_data)
DEFINE_STRUCT_ARITY_0(more_rbsp_trailing_data)
DEFINE_STRUCT_ARITY_0(payload_extension_present)



// Annex B
DEFINE_STRUCT_ARITY_1(byte_stream_nal_unit, NumBytesInNalUnit);
struct leading_zero_8bits { };
template <> struct Fixed<leading_zero_8bits> { static const int value = 0x00; };
struct zero_byte { };
template <> struct Fixed<zero_byte> { static const int value = 0x00; };
struct start_code_prefix_one_3bytes { };
template <> struct Fixed<start_code_prefix_one_3bytes> { static const int value = 0x000001; };
struct trailing_zero_8bits { };
template <> struct Fixed<trailing_zero_8bits> { static const int value = 0x00; };


DEFINE_STRUCT_ARITY_1(nal_unit, NumBytesInNalUnit);
DEFINE_STRUCT_ARITY_1(rbsp_byte, k);
struct emulation_prevention_three_byte { };
template <> struct Fixed<emulation_prevention_three_byte> { static const int value = 0x03; };


DEFINE_STRUCT_ARITY_0(nal_unit_header)
struct forbidden_zero_bit { };
template <> struct Fixed<forbidden_zero_bit> { static const int value = 0; }; // CondCheck 7.4.2.2-A
struct nal_unit_type { };
struct nuh_layer_id { };
struct nuh_temporal_id_plus1 { };
DEFINE_DERIVED(TemporalId , h[nuh_temporal_id_plus1()] - 1);


DEFINE_STRUCT_ARITY_0(RbspReserved)


DEFINE_STRUCT_ARITY_0(RbspUnspecified)


DEFINE_STRUCT_ARITY_0(video_parameter_set_rbsp)
struct vps_video_parameter_set_id { };
struct vps_reserved_three_2bits { };
struct vps_max_layers_minus1 { };
DEFINE_DERIVED(MaxLayersMinus1, std::min(62, h[vps_max_layers_minus1()]))
struct vps_max_sub_layers_minus1 { };
struct vps_temporal_id_nesting_flag { };
struct vps_reserved_0xffff_16bits { };
struct vps_sub_layer_ordering_info_present_flag { };
DEFINE_STRUCT_ARITY_1(vps_max_dec_pic_buffering_minus1, i);
template <> struct IsValueArray<vps_max_dec_pic_buffering_minus1> : std::true_type { };
template <> struct ValueHolder<vps_max_dec_pic_buffering_minus1> : ValueArray1<vps_max_dec_pic_buffering_minus1, 7> { };
DEFINE_VALUE_ARRAY_1(vps_max_num_reorder_pics, i, 7);
template <> struct ValueType<struct vps_max_latency_increase_plus1> { typedef unsigned Type; };
DEFINE_VALUE_ARRAY_1(vps_max_latency_increase_plus1, i, 7);
struct vps_max_layer_id { };
struct vps_num_layer_sets_minus1 { };
struct layer_id_included_flag { };
struct vps_timing_info_present_flag { };
struct vps_num_units_in_tick { };
template <> struct ValueType<vps_num_units_in_tick> { typedef unsigned Type; };
struct vps_time_scale { };
template <> struct ValueType<vps_time_scale> { typedef unsigned Type; };
struct vps_poc_proportional_to_timing_flag { };
struct vps_num_ticks_poc_diff_one_minus1 { };
template <> struct ValueType<struct vps_num_ticks_poc_diff_one_minus1> { typedef unsigned Type; };
struct vps_num_hrd_parameters { };
DEFINE_VALUE_ARRAY_1(hrd_layer_set_idx, i, 1024);
DEFINE_VALUE_ARRAY_1(cprms_present_flag, i, 1024);
struct vps_extension_flag { };
struct vps_extension_data_flag { };



DEFINE_STRUCT_ARITY_0(seq_parameter_set_rbsp);
struct sps_video_parameter_set_id { };
struct sps_max_sub_layers_minus1 { };
struct sps_ext_or_max_sub_layers_minus1 { };
struct MultiLayerExtSpsFlag { };
struct sps_temporal_id_nesting_flag { };
struct sps_seq_parameter_set_id { };
struct chroma_format_idc { };
struct SubWidthC { };
DEFINE_DERIVED_LONG(SubWidthC)
{
    const int lookup[] = {1, 2, 2, 1, 1};
    return lookup[h[chroma_format_idc()]]; // review: should this be ChromaArrayType ?
}
};
struct SubHeightC { };
DEFINE_DERIVED_LONG(SubHeightC)
{
    const int lookup[] = {1, 2, 1, 1, 1};
    return lookup[h[chroma_format_idc()]];
}
};
struct separate_colour_plane_flag { };
DEFINE_DERIVED(ChromaArrayType, h[separate_colour_plane_flag()] ? 0 : h[chroma_format_idc()]);
struct pic_width_in_luma_samples { };
DEFINE_DERIVED(PicWidthInSamplesC, h[pic_width_in_luma_samples()] / h[SubWidthC()]);
struct pic_height_in_luma_samples { };
DEFINE_DERIVED(PicHeightInSamplesC, h[pic_height_in_luma_samples()] / h[SubHeightC()]);
DEFINE_DERIVED(PicSizeInSamplesY, h[pic_width_in_luma_samples()] * h[pic_height_in_luma_samples()]);
struct conformance_window_flag { };
struct conf_win_left_offset { };
struct conf_win_right_offset { };
struct conf_win_top_offset { };
struct conf_win_bottom_offset { };
struct bit_depth_luma_minus8 { };
DEFINE_DERIVED(BitDepthY, 8 + h[bit_depth_luma_minus8()])
DEFINE_DERIVED(QpBdOffsetY, 6 * h[bit_depth_luma_minus8()])
struct bit_depth_chroma_minus8 { };
DEFINE_DERIVED(BitDepthC, 8 + h[bit_depth_chroma_minus8()])
DEFINE_DERIVED(QpBdOffsetC, 6 * h[bit_depth_chroma_minus8()])
struct log2_max_pic_order_cnt_lsb_minus4 { };
DEFINE_DERIVED(MaxPicOrderCntLsb, 1 << (h[log2_max_pic_order_cnt_lsb_minus4()] + 4))
struct PicOrderCntVal { };
template <> struct ValueType < PicOrderCntVal > { typedef int32_t Type; };
struct sps_sub_layer_ordering_info_present_flag { };
DEFINE_VALUE_ARRAY_1(sps_max_dec_pic_buffering_minus1, i, 7);
DEFINE_VALUE_ARRAY_1(sps_max_num_reorder_pics, i, 7);
template <> struct ValueType<struct sps_max_latency_increase_plus1> { typedef unsigned Type; };
DEFINE_VALUE_ARRAY_1(sps_max_latency_increase_plus1, i, 7);
DEFINE_STRUCT_ARITY_1(SpsMaxLatencyPictures, i)
DEFINE_DERIVED_LONG(SpsMaxLatencyPictures)
{
    assert(h[sps_max_latency_increase_plus1(v.i) ] != 0);
    return h[sps_max_num_reorder_pics(v.i)] + h[sps_max_latency_increase_plus1(v.i)] - 1;
}
};
struct log2_min_luma_coding_block_size_minus3 { };
DEFINE_DERIVED(MinCbLog2SizeY, h[log2_min_luma_coding_block_size_minus3()] + 3)
DEFINE_DERIVED(MinCbSizeY, 1 << h[MinCbLog2SizeY()])
DEFINE_DERIVED(PicWidthInMinCbsY, h[pic_width_in_luma_samples()] /  h[MinCbSizeY()])
DEFINE_DERIVED(PicHeightInMinCbsY, h[pic_height_in_luma_samples()] /  h[MinCbSizeY()])
DEFINE_DERIVED(PicSizeInMinCbsY, h[PicWidthInMinCbsY()] * h[PicHeightInMinCbsY()])
DEFINE_DERIVED(RawMinCuBits, h[MinCbSizeY()] * h[MinCbSizeY()] * ( h[BitDepthY()] + h[BitDepthC()]/ 2 ))
struct log2_diff_max_min_luma_coding_block_size { };
DEFINE_DERIVED(CtbLog2SizeY , h[MinCbLog2SizeY()] + h[log2_diff_max_min_luma_coding_block_size()])
struct CtbSizeY;
DEFINE_DERIVED(CtbSizeY, 1  <<  h[CtbLog2SizeY()])
DEFINE_DERIVED(PicWidthInCtbsY, ceilDiv( h[pic_width_in_luma_samples()],  h[CtbSizeY()]))
DEFINE_DERIVED(PicHeightInCtbsY, ceilDiv( h[pic_height_in_luma_samples()],  h[CtbSizeY()]))
DEFINE_DERIVED(PicSizeInCtbsY, h[PicWidthInCtbsY()] * h[PicHeightInCtbsY()])
struct log2_min_luma_transform_block_size_minus2 { };
DEFINE_DERIVED(MinTbLog2SizeY , h[log2_min_luma_transform_block_size_minus2()] + 2)
struct log2_diff_max_min_luma_transform_block_size { };
DEFINE_DERIVED(MaxTbLog2SizeY, h[log2_min_luma_transform_block_size_minus2()] + 2 + h[ log2_diff_max_min_luma_transform_block_size()])
struct max_transform_hierarchy_depth_inter { };
struct max_transform_hierarchy_depth_intra { };
struct scaling_list_enabled_flag { };
struct sps_scaling_list_data_present_flag { };
struct amp_enabled_flag { };
struct sample_adaptive_offset_enabled_flag { };
struct pcm_enabled_flag { };
struct pcm_sample_bit_depth_luma_minus1 { };
DEFINE_DERIVED(PcmBitDepthY, h[pcm_sample_bit_depth_luma_minus1()] + 1)
struct pcm_sample_bit_depth_chroma_minus1 { };
DEFINE_DERIVED(PcmBitDepthC, h[pcm_sample_bit_depth_chroma_minus1()] + 1)
struct log2_min_pcm_luma_coding_block_size_minus3 { };
DEFINE_DERIVED(Log2MinIpcmCbSizeY , h[log2_min_pcm_luma_coding_block_size_minus3()] + 3)
struct log2_diff_max_min_pcm_luma_coding_block_size { };
DEFINE_DERIVED(Log2MaxIpcmCbSizeY, h[log2_diff_max_min_pcm_luma_coding_block_size ()] + h[Log2MinIpcmCbSizeY()])
struct pcm_loop_filter_disabled_flag { };
struct num_short_term_ref_pic_sets { };
struct long_term_ref_pics_present_flag { };
struct num_long_term_ref_pics_sps { };
DEFINE_VALUE_ARRAY_1(lt_ref_pic_poc_lsb_sps, i, 32);
DEFINE_VALUE_ARRAY_1(used_by_curr_pic_lt_sps_flag, i, 32);
struct sps_temporal_mvp_enabled_flag { };
struct strong_intra_smoothing_enabled_flag { };
struct vui_parameters_present_flag { };
struct sps_extension_present_flag { };
struct sps_range_extension_flag { };
struct sps_multilayer_extension_flag { };
struct sps_extension_6bits { };
struct sps_extension_data_flag { };
DEFINE_STRUCT_ARITY_2(MinTbAddrZs, x, y)

struct SliceAddrRs { };
struct CtbAddrInTs { };
struct CtbAddrInRs { };
DEFINE_STRUCT_ARITY_1(colBd, i)
DEFINE_STRUCT_ARITY_1(rowBd, i)
DEFINE_STRUCT_ARITY_1(TileId, ctbAddrTS)
DEFINE_STRUCT_ARITY_1(CtbAddrRsToTs, ctbAddrRS)
DEFINE_STRUCT_ARITY_1(CtbAddrTsToRs, ctbAddrTS)

DEFINE_DERIVED_LONG(MinTbAddrZs)
{
    const int tbX = ( v.x << h[MinTbLog2SizeY()] ) >> h[CtbLog2SizeY()];
    const int tbY = ( v.y << h[MinTbLog2SizeY()] ) >> h[CtbLog2SizeY()];
    const int ctbAddrRs = h[PicWidthInCtbsY()] * tbY + tbX;
    int val = h[CtbAddrRsToTs(ctbAddrRs) ] << ( ( h[CtbLog2SizeY()] - h[MinTbLog2SizeY()] ) * 2);
    int i, p;
    for( i = 0, p = 0; i < ( h[CtbLog2SizeY()] - h[MinTbLog2SizeY()] ); i++ )
    {
        int m = 1 << i;
        p += ( m & v.x ? m * m : 0 ) + ( m & v.y ? 2 * m * m : 0 );
    }
    val += p;
    return val;
}
};


template <typename Sample> struct StateReconstructedPicture;


// To-do - move availability code to separate file

// loop-free comparison of Z-order positions
// returns true if (xN, yN) precedes (xC, yC) in Z-order
static bool compareZ(int xC, int yC, int xN, int yN)
{
    int const xNot = ~xN;
    int const yNot = ~yN;

    int const yXor = yC ^ yNot;
    int const yAnd = yC & yNot;

    int const p = (yAnd | (xC & yXor));
    int const q = ~(yAnd | (xNot & yXor));

    return p > q;
}


// object designed to persist for lifetime of CTU
struct AvailabilityCtu
{
    template <class H>
    void init(H &h)
    {
        this->xCtb = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
        this->yCtb = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];
        this->width = h[pic_width_in_luma_samples()];
        this->height = h[pic_height_in_luma_samples()];
        this->left = initialiseNeighbour(h, -1, 0);
        this->aboveLeft = initialiseNeighbour(h, -1, -1);
        this->above = initialiseNeighbour(h, 0, -1);
        this->aboveRight = initialiseNeighbour(h, 1, -1);
    }

    template <char mode=0>
    bool available(int xCurr, int yCurr, int xN, int yN, int log2CtuSize) const
    {
        assert(xCurr >> log2CtuSize == this->xCtb);
        assert(yCurr >> log2CtuSize == this->yCtb);

        if (mode == 'L') assert(xN == xCurr - 1);
        if (mode == 'A') assert(yN == yCurr - 1);

        if (mode != 'L' && xN >= this->width) return false;
        if (mode != 'A' && yN >= this->height) return false;

        const int dx = (xN >> log2CtuSize) - this->xCtb;
        const int dy = (yN >> log2CtuSize) - this->yCtb;

        if (dy == 0)
        {
            if (dx == 0)
            {
                // (xN, yN) is within the current CTU
                if (mode != 'L' && xN >= this->width) return false;
                if (mode != 'A' && yN >= this->height) return false;
                return compareZ(xCurr, yCurr, xN, yN);
            }
            else if (mode != 'L' && dx == 1)
            {
                return false;
            }
            else
            {
                assert(dx == -1);
                return this->left;
            }
        }

        if (mode != 'A' && dy == 1)
        {
            return false;
        }
        else
        {
            assert(dy == -1);
            if (dx == 0)
            {
                return this->above;
            }
            else if (mode != 'L' && dx == 1)
            {
                return this->aboveRight;
            }
            else
            {
                assert(dx == -1);
                return this->aboveLeft;
            }
        }
    }
private:
    template <class H>
    bool initialiseNeighbour(H &h, int dx, int dy)
    {
        {
            const int xN = (h[CtbAddrInRs()] % h[PicWidthInCtbsY()] + dx) << h[CtbLog2SizeY()];
            const int yN = (h[CtbAddrInRs()] / h[PicWidthInCtbsY()] + dy) << h[CtbLog2SizeY()];

            if (xN < 0) return false;
            if (yN < 0) return false;
            if (xN >= h[pic_width_in_luma_samples()]) return false;
            if (yN >= h[pic_height_in_luma_samples()]) return false;
        }

        {
            const int CtbAddrInRsN = h[CtbAddrInRs()] + dx + h[PicWidthInCtbsY()] * dy;

            if (CtbAddrInRsN < h[SliceAddrRs()])
            {
                return false;
            }

            if (h[TileId(h[CtbAddrInTs()])] != h[TileId(h[CtbAddrRsToTs(CtbAddrInRsN)])])
            {
                return false;
            }
        }
        return true;
    }

    int xCtb;
    int yCtb;

    int width;
    int height;

    bool left;
    bool aboveLeft;
    bool above;
    bool aboveRight;
};


DEFINE_STRUCT_ARITY_0(sps_range_extension);
struct transform_skip_rotation_enabled_flag { };
struct transform_skip_context_enabled_flag { };
struct implicit_rdpcm_enabled_flag { };
struct explicit_rdpcm_enabled_flag { };
struct extended_precision_processing_flag { };
struct intra_smoothing_disabled_flag { };
struct high_precision_offsets_enabled_flag { };
struct persistent_rice_adaptation_enabled_flag { };
struct cabac_bypass_alignment_enabled_flag { };


DEFINE_STRUCT_ARITY_0(pic_parameter_set_rbsp);
struct pps_pic_parameter_set_id { };
struct pps_seq_parameter_set_id { };
struct dependent_slice_segments_enabled_flag { };
struct output_flag_present_flag { };
struct num_extra_slice_header_bits { };
struct sign_data_hiding_enabled_flag { };
struct cabac_init_present_flag { };
struct num_ref_idx_l0_default_active_minus1 { };
struct num_ref_idx_l1_default_active_minus1 { };
struct init_qp_minus26 { };
struct constrained_intra_pred_flag { };
struct transform_skip_enabled_flag { };
struct cu_qp_delta_enabled_flag { };
struct diff_cu_qp_delta_depth { };
DEFINE_DERIVED(Log2MinCuQpDeltaSize, h[CtbLog2SizeY()] - h[diff_cu_qp_delta_depth()]);
DEFINE_DERIVED(MinCuQpDeltaSize, 1 << h[Log2MinCuQpDeltaSize()]);
struct pps_cb_qp_offset { };
struct pps_cr_qp_offset { };
struct pps_slice_chroma_qp_offsets_present_flag { };
struct weighted_pred_flag { };
struct weighted_bipred_flag { };
struct slice_type;
DEFINE_DERIVED(weightedPredFlag, (h[slice_type()] == P) ? h[weighted_pred_flag()] : h[weighted_bipred_flag()]); // todo - ValueCache for this
struct transquant_bypass_enabled_flag { };
struct tiles_enabled_flag { };
struct entropy_coding_sync_enabled_flag { };
struct num_tile_columns_minus1 { };
struct num_tile_rows_minus1 { };
struct uniform_spacing_flag { };
DEFINE_STRUCT_ARITY_1(column_width_minus1, i);
DEFINE_STRUCT_ARITY_1(row_height_minus1, i);
struct loop_filter_across_tiles_enabled_flag { };
struct pps_loop_filter_across_slices_enabled_flag { };
struct deblocking_filter_control_present_flag { };
struct deblocking_filter_override_enabled_flag { };
struct pps_deblocking_filter_disabled_flag { };
struct pps_beta_offset_div2 { };
struct pps_tc_offset_div2 { };
struct pps_scaling_list_data_present_flag { };
struct lists_modification_present_flag { };
struct log2_parallel_merge_level_minus2 { };
DEFINE_DERIVED(Log2ParMrgLevel , h[log2_parallel_merge_level_minus2()] + 2);
struct slice_segment_header_extension_present_flag { };
struct pps_extension_present_flag { };
struct pps_range_extension_flag { };
struct pps_multilayer_extension_flag { };
struct pps_extension_6bits { };
struct pps_extension_data_flag { };


DEFINE_STRUCT_ARITY_0(pps_range_extension);
struct log2_max_transform_skip_block_size_minus2 { };
DEFINE_DERIVED(Log2MaxTransformSkipSize, h[log2_max_transform_skip_block_size_minus2()] + 2);
struct cross_component_prediction_enabled_flag { };
struct chroma_qp_offset_list_enabled_flag { };
struct diff_cu_chroma_qp_offset_depth { };
DEFINE_DERIVED(Log2MinCuChromaQpOffsetSize, h[CtbLog2SizeY()] - h[diff_cu_chroma_qp_offset_depth()])
struct chroma_qp_offset_list_len_minus1 { };
DEFINE_VALUE_ARRAY_1(cb_qp_offset_list, i, 32);
DEFINE_VALUE_ARRAY_1(cr_qp_offset_list, i, 32);
struct log2_sao_offset_scale_luma { };
struct log2_sao_offset_scale_chroma { };


DEFINE_STRUCT_ARITY_0(sei_rbsp);


DEFINE_STRUCT_ARITY_0(access_unit_delimiter_rbsp);
struct pic_type { };


DEFINE_STRUCT_ARITY_0(end_of_seq_rbsp);


DEFINE_STRUCT_ARITY_0(end_of_bitstream_rbsp);


DEFINE_STRUCT_ARITY_0(filler_data_rbsp);


DEFINE_STRUCT_ARITY_0(slice_segment_layer_rbsp);


DEFINE_STRUCT_ARITY_0(rbsp_slice_segment_trailing_bits);
struct cabac_zero_word { };
template <> struct Fixed<cabac_zero_word> { static const int value = 0x0000; };

DEFINE_STRUCT_ARITY_0(rbsp_trailing_bits);
struct rbsp_stop_one_bit { };
template <> struct Fixed<rbsp_stop_one_bit> { static const int value = 1; }; // CondCheck 7.4.3.11-A
struct rbsp_alignment_zero_bit { };
template <> struct Fixed<rbsp_alignment_zero_bit> { static const int value = 0; }; // CondCheck 7.4.3.11-B


DEFINE_STRUCT_ARITY_0(byte_alignment);
struct alignment_bit_equal_to_one { };
template <> struct Fixed<alignment_bit_equal_to_one> { static const int value = 1; }; // CondCheck 7.4.3.12-A
struct alignment_bit_equal_to_zero { };
template <> struct Fixed<alignment_bit_equal_to_zero> { static const int value = 0; }; // CondCheck 7.4.3.12-B


DEFINE_STRUCT_ARITY_1(profile_tier_level, maxNumSubLayersMinus1 );
struct general_profile_space { };
struct general_tier_flag { };
struct general_profile_idc { };
DEFINE_VALUE_ARRAY_1(general_profile_compatibility_flag, j, 32);
struct general_progressive_source_flag { };
struct general_interlaced_source_flag { };
struct general_non_packed_constraint_flag { };
struct general_frame_only_constraint_flag { };
struct general_max_12bit_constraint_flag { };
struct general_max_10bit_constraint_flag { };
struct general_max_8bit_constraint_flag { };
struct general_max_422chroma_constraint_flag { };
struct general_max_420chroma_constraint_flag { };
struct general_max_monochrome_constraint_flag { };
struct general_intra_constraint_flag { };
struct general_one_picture_only_constraint_flag { };
struct general_lower_bit_rate_constraint_flag { };
struct general_reserved_zero_34bits { };
template <> struct ValueType<general_reserved_zero_34bits> { typedef int64_t Type; };
struct general_reserved_zero_43bits { };
template <> struct ValueType<general_reserved_zero_43bits> { typedef int64_t Type; };
struct general_inbld_flag { };
struct general_reserved_zero_bit { };
struct general_level_idc  { };
DEFINE_VALUE_ARRAY_1(sub_layer_profile_present_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_level_present_flag, i, 8);
DEFINE_VALUE_ARRAY_1(reserved_zero_2bits, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_profile_space, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_tier_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_profile_idc, i, 8);
DEFINE_VALUE_ARRAY_2(sub_layer_profile_compatibility_flag, i, 8, j, 32);
DEFINE_VALUE_ARRAY_1(sub_layer_progressive_source_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_interlaced_source_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_non_packed_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_frame_only_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_max_12bit_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_max_10bit_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_max_8bit_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_max_422chroma_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_max_420chroma_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_max_monochrome_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_intra_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_one_picture_only_constraint_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_lower_bit_rate_constraint_flag, i, 8);
DEFINE_STRUCT_ARITY_1(sub_layer_reserved_zero_34bits, i);
template <> struct ValueType<sub_layer_reserved_zero_34bits> { typedef int64_t Type; };
DEFINE_STRUCT_ARITY_1(sub_layer_reserved_zero_43bits, i);
template <> struct ValueType<sub_layer_reserved_zero_43bits> { typedef int64_t Type; };
DEFINE_VALUE_ARRAY_1(sub_layer_inbld_flag, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_reserved_zero_bit, i, 8);
DEFINE_VALUE_ARRAY_1(sub_layer_level_idc, i, 8);


DEFINE_STRUCT_ARITY_0(scaling_list_data);
DEFINE_VALUE_ARRAY_2(scaling_list_pred_mode_flag, sizeId, 4, matrixId, 6);
DEFINE_VALUE_ARRAY_2(scaling_list_pred_matrix_id_delta, sizeId, 4, matrixId, 6);
DEFINE_VALUE_ARRAY_2(scaling_list_dc_coef_minus8, sizeId, 4, matrixId, 6);
//template <> struct ValueType<struct ScalingList> { typedef int16_t Type; };
DEFINE_VALUE_ARRAY_3(ScalingList, sizeId, 4, matrixId, 6, i, 64);
//DEFINE_STRUCT_ARITY_4(ScalingFactor, sizeId, matrixId, x, y);
//template <> struct ValueType<ScalingFactor> : ValueType<ScalingList> { };

struct scaling_list_delta_coef { };


// sei_message
DEFINE_STRUCT_ARITY_0(sei_message);
struct last_payload_type_byte { };
struct last_payload_size_byte { };
struct ff_byte { };
template <> struct Fixed<ff_byte> { static const int value = 0xff; }; // CondCheck D.3.5-A, 7.4.3.8-A


DEFINE_STRUCT_ARITY_2(sei_payload, payloadType, payloadSize)


DEFINE_STRUCT_ARITY_0(slice_segment_header);
struct first_slice_segment_in_pic_flag { };
struct no_output_of_prior_pics_flag { };
struct slice_pic_parameter_set_id { };
struct dependent_slice_segment_flag { };
struct slice_segment_address { };
struct slice_reserved_flag { slice_reserved_flag(int) { } };
struct slice_type { };
struct pic_output_flag { };
struct colour_plane_id { };
struct slice_pic_order_cnt_lsb { };
struct short_term_ref_pic_set_sps_flag { };
struct short_term_ref_pic_set_idx { };
struct num_long_term_sps { };
struct num_long_term_pics { };
struct CurrRpsIdx { };
DEFINE_DERIVED_LONG(CurrRpsIdx)
{
    if (h[short_term_ref_pic_set_sps_flag()] == 1)
    {
        return h[short_term_ref_pic_set_idx()];
    }
    return h[num_short_term_ref_pic_sets()];
}
};
struct NumPocTotalCurr { };
DEFINE_VALUE_ARRAY_1(UsedByCurrPicLt, i, 32);
DEFINE_VALUE_ARRAY_2(UsedByCurrPicS0, stRpsIdx , 65, i, 16);
DEFINE_VALUE_ARRAY_2(UsedByCurrPicS1, stRpsIdx , 65, i, 16);
DEFINE_VALUE_ARRAY_1(NumNegativePics, stRpsIdx, 65);
DEFINE_VALUE_ARRAY_1(NumPositivePics, stRpsIdx, 65);
DEFINE_DERIVED_LONG(NumPocTotalCurr)
{
    int NumPocTotalCurr = 0;
    for( int i = 0; i < h[NumNegativePics(h[CurrRpsIdx()])]; i++ )
        if( h[UsedByCurrPicS0(h[CurrRpsIdx()], i)] )
            NumPocTotalCurr++;
    for( int i = 0; i < h[NumPositivePics(h[CurrRpsIdx()])]; i++)
        if( h[UsedByCurrPicS1(h[CurrRpsIdx()], i)] )
            NumPocTotalCurr++;
    for( int i = 0; i < h[num_long_term_sps()] + h[num_long_term_pics()]; i++ )
        if( h[UsedByCurrPicLt(i)] )
            NumPocTotalCurr++;
    return NumPocTotalCurr;
}
};
DEFINE_STRUCT_ARITY_1(NumDeltaPocs, idx)
DEFINE_DERIVED_LONG(NumDeltaPocs)
{
    return h[NumNegativePics(v.idx)] + h[NumPositivePics(v.idx)];
}
};
DEFINE_VALUE_ARRAY_1(lt_idx_sps, i, 32);
DEFINE_VALUE_ARRAY_1(poc_lsb_lt, i, 32);
DEFINE_VALUE_ARRAY_1(used_by_curr_pic_lt_flag, i, 32);
DEFINE_VALUE_ARRAY_1(delta_poc_msb_present_flag, i, 32);
DEFINE_VALUE_ARRAY_1(delta_poc_msb_cycle_lt, i, 32);
struct slice_temporal_mvp_enabled_flag { };
struct slice_sao_luma_flag { };
struct slice_sao_chroma_flag { };
struct num_ref_idx_active_override_flag { };
struct num_ref_idx_l0_active_minus1 { };
struct num_ref_idx_l1_active_minus1 { };
struct mvd_l1_zero_flag { };
struct cabac_init_flag { };
struct initType { };
DEFINE_DERIVED_LONG(initType)
{

    if( h[slice_type()] == I )
    {
        return 0;
    }
    else if( h[slice_type()] == P )
    {
        return h[cabac_init_flag()] ? 2 : 1;
    }
    else
    {
        return h[cabac_init_flag()] ? 1 : 2;
    }
}
};
struct collocated_from_l0_flag { };
struct collocated_ref_idx { };
struct five_minus_max_num_merge_cand { };
DEFINE_DERIVED(MaxNumMergeCand, 5 - h[five_minus_max_num_merge_cand()]);
struct slice_qp_delta { };
DEFINE_DERIVED(SliceQpY, 26 + h[init_qp_minus26()] + h[slice_qp_delta()]);
struct QpY { };
struct slice_cb_qp_offset { };
struct slice_cr_qp_offset { };
struct cu_chroma_qp_offset_enabled_flag { };
struct deblocking_filter_override_flag { };
struct slice_deblocking_filter_disabled_flag { };
struct slice_beta_offset_div2 { };
struct slice_tc_offset_div2 { };
struct slice_loop_filter_across_slices_enabled_flag { };
struct num_entry_point_offsets { };
struct offset_len_minus1 { };
DEFINE_STRUCT_ARITY_1(entry_point_offset_minus1, i);
struct slice_segment_header_extension_length { };
struct slice_segment_header_extension_data_byte { slice_segment_header_extension_data_byte(int) { } };



DEFINE_STRUCT_ARITY_0(ref_pic_lists_modification);
struct ref_pic_list_modification_flag_l0 { };
DEFINE_VALUE_ARRAY_1(list_entry_l0, i, 16)
struct ref_pic_list_modification_flag_l1 { };
DEFINE_VALUE_ARRAY_1(list_entry_l1, i, 16)


DEFINE_STRUCT_ARITY_0(pred_weight_table)
struct luma_log2_weight_denom { };
struct delta_chroma_log2_weight_denom { };
DEFINE_DERIVED(ChromaLog2WeightDenom, h[luma_log2_weight_denom()] + h[delta_chroma_log2_weight_denom()])
DEFINE_VALUE_ARRAY_1(luma_weight_l0_flag, i, 32)
DEFINE_VALUE_ARRAY_1(chroma_weight_l0_flag, i, 32)
DEFINE_VALUE_ARRAY_1(delta_luma_weight_l0, i, 32)
DEFINE_VALUE_ARRAY_1(luma_offset_l0, i, 32)
DEFINE_VALUE_ARRAY_2(delta_chroma_weight_l0, i, 32, j, 2)
DEFINE_VALUE_ARRAY_2(delta_chroma_offset_l0, i, 32, j, 2)
DEFINE_VALUE_ARRAY_1(luma_weight_l1_flag, i, 32)
DEFINE_VALUE_ARRAY_1(chroma_weight_l1_flag, i, 32)
DEFINE_VALUE_ARRAY_1(delta_luma_weight_l1, i, 32)
DEFINE_VALUE_ARRAY_1(luma_offset_l1, i, 32)
DEFINE_VALUE_ARRAY_2(delta_chroma_weight_l1, i, 32, j, 2)
DEFINE_VALUE_ARRAY_2(delta_chroma_offset_l1, i, 32, j, 2)
DEFINE_STRUCT_ARITY_1(LumaWeightL0, i);
DEFINE_DERIVED_LONG(LumaWeightL0)
{
    return h[luma_weight_l0_flag(v.i)] ? (1 << h[luma_log2_weight_denom()]) + h[delta_luma_weight_l0(v.i)] : (1 << h[luma_log2_weight_denom()]);
}
};
DEFINE_STRUCT_ARITY_1(LumaWeightL1, i)
DEFINE_DERIVED_LONG(LumaWeightL1)
{
    return h[luma_weight_l1_flag(v.i)] ? (1 << h[luma_log2_weight_denom()]) + h[delta_luma_weight_l1(v.i)] : (1 << h[luma_log2_weight_denom()]);
}
};
DEFINE_STRUCT_ARITY_2(ChromaWeightL0, i, j)
DEFINE_DERIVED_LONG(ChromaWeightL0)
{
    return h[chroma_weight_l0_flag(v.i)] ? (1 << h[ChromaLog2WeightDenom()]) + h[delta_chroma_weight_l0(v.i, v.j)] : (1 << h[ChromaLog2WeightDenom()]);
}
};
DEFINE_STRUCT_ARITY_2(ChromaWeightL1, i, j)
DEFINE_DERIVED_LONG(ChromaWeightL1)
{
    return h[chroma_weight_l1_flag(v.i)] ? (1 << h[ChromaLog2WeightDenom()]) + h[delta_chroma_weight_l1(v.i, v.j)] : (1 << h[ChromaLog2WeightDenom()]);
}
};
DEFINE_STRUCT_ARITY_2(ChromaOffsetL0, i, j)
DEFINE_DERIVED_LONG(ChromaOffsetL0)
{
    if (h[chroma_weight_l0_flag(v.i)] == 0) return 0;
    return Clip3( -128, 127, ( h[delta_chroma_offset_l0(v.i, v.j)] -
            ( ( 128 * h[ChromaWeightL0(v.i, v.j)] )  >>  h[ChromaLog2WeightDenom()] ) + 128 ) );
}
};
DEFINE_STRUCT_ARITY_2(ChromaOffsetL1, i, j)
DEFINE_DERIVED_LONG(ChromaOffsetL1)
{
    if (h[chroma_weight_l1_flag(v.i)] == 0) return 0;
    return Clip3( -128, 127, ( h[delta_chroma_offset_l1(v.i, v.j)] -
            ( ( 128 * h[ChromaWeightL1(v.i, v.j)] )  >>  h[ChromaLog2WeightDenom()] ) + 128 ) );
}
};


DEFINE_STRUCT_ARITY_1(short_term_ref_pic_set, stRpsIdx )
static const int maxDecodedPictures = 16;
struct inter_ref_pic_set_prediction_flag { };
struct delta_idx_minus1 { };
struct delta_rps_sign { };
struct abs_delta_rps_minus1 { };
DEFINE_DERIVED(DeltaRPS, (1 - 2 * h[delta_rps_sign()]) * (h[abs_delta_rps_minus1()] + 1))
DEFINE_VALUE_ARRAY_1(used_by_curr_pic_flag, j, maxDecodedPictures)
DEFINE_VALUE_ARRAY_1(use_delta_flag, j, maxDecodedPictures)
struct num_negative_pics { };
struct num_positive_pics { };
DEFINE_VALUE_ARRAY_1(delta_poc_s0_minus1, j, maxDecodedPictures)
DEFINE_VALUE_ARRAY_1(used_by_curr_pic_s0_flag, j, maxDecodedPictures)
DEFINE_VALUE_ARRAY_1(delta_poc_s1_minus1, j, maxDecodedPictures)
DEFINE_VALUE_ARRAY_1(used_by_curr_pic_s1_flag, j, maxDecodedPictures)


DEFINE_STRUCT_ARITY_0(slice_segment_data)
struct end_of_slice_segment_flag { };
struct end_of_subset_one_bit { };
template <> struct Fixed<end_of_subset_one_bit> { static const int value = 1; }; // CondCheck 7.4.9.1-A
struct EndOfSubStream { };


DEFINE_STRUCT_ARITY_0(coding_tree_unit)
struct xCtb { };
struct yCtb { };


DEFINE_STRUCT_ARITY_2(sao, rx, ry)
struct sao_merge_left_flag { };
struct sao_merge_up_flag { };
struct sao_type_idx_luma { };
struct sao_type_idx_chroma { };
DEFINE_STRUCT_ARITY_3(SaoTypeIdx, cIdx, rx, ry)
template <> struct IsValueArray<SaoTypeIdx> : std::true_type { };
template <> struct ValueHolder<SaoTypeIdx> : ValueArray1<SaoTypeIdx, 3> { };
DEFINE_STRUCT_ARITY_3(SaoEoClass, cIdx, rx, ry)
template <> struct IsValueArray<SaoEoClass> : std::true_type { };
template <> struct ValueHolder<SaoEoClass> : ValueArray1<SaoEoClass, 3> { };
DEFINE_STRUCT_ARITY_4(sao_offset_abs, cIdx, rx, ry, i)
template <> struct IsValueArray<sao_offset_abs> : std::true_type { };
template <> struct ValueHolder<sao_offset_abs> : ValueArray2<sao_offset_abs, 3, 4>{};
template <> struct Index<sao_offset_abs, 0> { static const int value = 0; };
template <> struct Index<sao_offset_abs, 1> { static const int value = 3; };
DEFINE_STRUCT_ARITY_4(sao_offset_sign, cIdx, rx, ry, i)
template <> struct IsValueArray<sao_offset_sign> : std::true_type{};
template <> struct ValueHolder<sao_offset_sign> : ValueArray2<sao_offset_sign, 3, 4>{};
template <> struct Index<sao_offset_sign, 0> { static const int value = 0; };
template <> struct Index<sao_offset_sign, 1> { static const int value = 3; };
DEFINE_STRUCT_ARITY_3(sao_band_position, cIdx, rx, ry)
template <> struct IsValueArray<sao_band_position> : std::true_type{};
template <> struct ValueHolder<sao_band_position> : ValueArray1<sao_band_position, 3>{};
template <> struct Index<sao_band_position, 0> { static const int value = 0; };
struct sao_eo_class_luma { };
struct sao_eo_class_chroma { };


DEFINE_STRUCT_ARITY_4(coding_quadtree, x0, y0, log2CbSize, cqtDepth)


// A special syntax construct designed to represent spatial element V that are not walked by the syntax.
// This is used to invoke code for coding_quadtree()s that are not parsed because they fall outside of the picture.
// The Direction parameter allows more intelligent behaviour based on whether the deleted coding_quadtree was
// Right or Down of the picture boundary.
template <class V, class Direction>
struct Deleted :
    V
    {
        using V::V;
    };


#if defined(_MSC_VER) && _MSC_VER <= 1800
// This specialisation exists because the VS2013 compiler does not support inherited constructors (the "using V::V" above)
template <class Direction>
struct Deleted<coding_quadtree, Direction>
:
coding_quadtree
{
    Deleted(int x0, int y0, int log2CbSize, int cqtDepth)
:
    coding_quadtree(x0, y0, log2CbSize, cqtDepth)
    {
    };
};
#endif


DEFINE_STRUCT_ARITY_2(split_cu_flag, x0, y0)

DEFINE_STRUCT_ARITY_2(CtDepth, x, y)
template <> struct ValueType < CtDepth > { typedef int8_t Type; };

struct IsCuQpDeltaCoded { };
struct CuQpDeltaVal { };
struct IsCuChromaQpOffsetCoded { };
struct CuQpOffsetCb { };
struct CuQpOffsetCr { };



DEFINE_STRUCT_ARITY_3(coding_unit, x0, y0, log2CbSize)
struct cu_transquant_bypass_flag { };
DEFINE_STRUCT_ARITY_2(cu_skip_flag, x0, y0);
struct pred_mode_flag { };


// Direction cues for spatial memories and other purposes.
// Tag types used to specialise templates.
struct Corner {}; // above and left of the current block
struct Up {}; // above the current block
struct Left {}; // left of the current block
struct Down {}; // down of the current block
struct Right {}; // right of the current block
struct Current {}; // within the current block
struct CurrentCu {}; // within the current coding_unit
struct Anywhere {}; // wanted value could be anywhere


template <class V, class Direction>
struct Neighbouring
{
    int dx, dy;
    V v;
};

// hints to accessor code that the parameter is the currently applicable value (i.e. not that of a neighbour)
template <class V> Neighbouring<V, Current> current(V v)
{
    return Neighbouring<V, Current>{0, 0, v};
}

template <class V> Neighbouring<V, Left> left(int dx, int dy, V v)
{
    ASSERT(dx == -1);
    return Neighbouring < V, Left > {dx, dy, v};
}

template <class V> Neighbouring<V, Up> up(int dx, int dy, V v)
{
    ASSERT(dy == -1);
    return Neighbouring < V, Up > {dx, dy, v};
}


DEFINE_STRUCT_ARITY_2(CuPredMode, x0, y0);


struct part_mode { };

struct PartMode { };
template <> struct ValueType < PartMode > { typedef PartModeType Type; };
DEFINE_DERIVED_LONG(PartMode)
{
    const coding_quadtree *cqt = h;
    if (h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA)
    {
        switch (h[part_mode()])
        {
            default:
            case 0: 	return PART_2Nx2N;
            case 1:	return PART_NxN;
        }
    }
    else // CuPredMode == MODE_INTER
    {
        switch (h[part_mode()])
        {
            default:
            case 0: 	return PART_2Nx2N;
            case 1:	return PART_2NxN;
            case 2: 	return PART_Nx2N;
            case 3:	return PART_NxN;
            case 4: 	return PART_2NxnU;
            case 5:	return PART_2NxnD;
            case 6: 	return PART_nLx2N;
            case 7:	return PART_nRx2N;
        }
    }
}
};

struct IntraSplitFlag { };
DEFINE_DERIVED_LONG(IntraSplitFlag)
{
    const coding_quadtree *cqt = h;
    if (h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA)
    {
        return h[part_mode()];
    }
    else // CuPredMode == MODE_INTER
    {
        return 0;
    }
}
};
struct pcm_alignment_one_bit { };
template <> struct Fixed<pcm_alignment_one_bit> { static const int value = 1; };
struct pcm_alignment_zero_bit { };
template <> struct Fixed<pcm_alignment_zero_bit> { static const int value = 0; };
DEFINE_STRUCT_ARITY_2(pcm_flag, x0, y0);
DEFINE_STRUCT_ARITY_2(prev_intra_luma_pred_flag, x0, y0);
DEFINE_STRUCT_ARITY_2(mpm_idx, x0, y0);
DEFINE_STRUCT_ARITY_2(rem_intra_luma_pred_mode, x0, y0);
DEFINE_STRUCT_ARITY_2(intra_chroma_pred_mode, x0, y0);
struct rqt_root_cbf { };
struct MaxTrafoDepth { };

DEFINE_STRUCT_ARITY_2(IntraPredModeY, x0, y0);
template <> struct ValueType < IntraPredModeY > { typedef int8_t Type; };

DEFINE_STRUCT_ARITY_2(IntraPredModeC, xPb, yPb);
template <> struct ValueType < IntraPredModeC > { typedef int8_t Type; };
DEFINE_DERIVED_LONG(IntraPredModeC)
{
    const coding_quadtree *cqt = h;
    const coding_unit cu(cqt->x0, cqt->y0, cqt->log2CbSize);

    int intraPredModeY;
    if (h[ChromaArrayType()] == 3 && h[PartMode()] == PART_NxN)
    {
        intraPredModeY = h[IntraPredModeY(v.xPb, v.yPb)];
    }
    else
    {
        intraPredModeY = h[IntraPredModeY(cu.x0, cu.y0)];
    }

    assert(h[intra_chroma_pred_mode(v.xPb, v.yPb)] >= 0);
    assert(h[intra_chroma_pred_mode(v.xPb, v.yPb)] <= 4);

    int modeIdx;
    if (h[intra_chroma_pred_mode(v.xPb, v.yPb)] == 4)
    {
        modeIdx = intraPredModeY;
    }
    else
    {
        static const int lookup[4] = { 0, 26, 10, 1 };
        modeIdx = lookup[h[intra_chroma_pred_mode(v.xPb, v.yPb)]];
        if (modeIdx == intraPredModeY) modeIdx = 34;
    }

    if (h[ChromaArrayType()] == 2)
    {
        // note typo in HEVC standard text: intraPredModeC appears without capital 'I'

        static const uint8_t lookupChromaAngleFor422[] =
        { 0, 1, 2, 2, 2, 2, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 19, 20, 21, 22, 23, 23, 24, 24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 31 };

        modeIdx = lookupChromaAngleFor422[modeIdx];
    }

    return modeIdx;
}
};


DEFINE_STRUCT_ARITY_4(prediction_unit, x0, y0, nPbW, nPbH )
DEFINE_STRUCT_ARITY_2(merge_flag, x0, y0);
DEFINE_STRUCT_ARITY_2(merge_idx, x0, y0);
DEFINE_STRUCT_ARITY_2(inter_pred_idc, x0, y0);
DEFINE_STRUCT_ARITY_2(mvp_l0_flag, x0, y0);
DEFINE_STRUCT_ARITY_2(mvp_l1_flag, x0, y0);
DEFINE_STRUCT_ARITY_2(ref_idx_l0, x0, y0);
DEFINE_STRUCT_ARITY_2(ref_idx_l1, x0, y0);


DEFINE_STRUCT_ARITY_3(pcm_sample, x0, y0, log2CbSize )
DEFINE_STRUCT_ARITY_1(pcm_sample_luma, i )
DEFINE_STRUCT_ARITY_1(pcm_sample_chroma, i )


DEFINE_STRUCT_ARITY_7(transform_tree, x0, y0, xBase, yBase, log2TrafoSize, trafoDepth, blkIdx)
struct interSplitFlag { };
DEFINE_STRUCT_ARITY_3(split_transform_flag, x0, y0, trafoDepth )

template <class H>
inline int infer(split_transform_flag e, H &h)
{
    const transform_tree &tt = *static_cast<transform_tree *>(h);
    if (tt.log2TrafoSize > h[MaxTbLog2SizeY()])
    {
        return 1;
    }
    if (h[IntraSplitFlag()] == 1 && e.trafoDepth == 0)
    {
        return 1;
    }
    if (h[interSplitFlag()] == 1)
    {
        return 1;
    }
    return 0;
}

DEFINE_DERIVED_LONG(interSplitFlag)
    {
        const transform_tree &tt = *static_cast<transform_tree *>(h);
        return (h[current(CuPredMode(tt.x0, tt.y0))] == MODE_INTER
            && h[max_transform_hierarchy_depth_inter()] == 0
            && h[PartMode()] != PART_2Nx2N
            && tt.trafoDepth == 0);
    }
};

DEFINE_STRUCT_ARITY_3(cbf_cb, x0, y0, trafoDepth )
template <> struct IsValueArray<cbf_cb> : std::true_type { };
template <> struct ValueHolder<cbf_cb> : ValueArray1<cbf_cb, 6> { };
template <> struct Index<cbf_cb, 0> { static const int value = 2; };
DEFINE_STRUCT_ARITY_3(cbf_cr, x0, y0, trafoDepth )
template <> struct IsValueArray<cbf_cr> : std::true_type { };
template <> struct ValueHolder<cbf_cr> : ValueArray1<cbf_cr, 6> { };
template <> struct Index<cbf_cr, 0> { static const int value = 2; };
DEFINE_STRUCT_ARITY_3(cbf_luma, x0, y0, trafoDepth )
DEFINE_STRUCT_ARITY_3(mvd_coding, x0, y0, refList )
DEFINE_VALUE_ARRAY_1(abs_mvd_greater0_flag, compIdx, 2 )
DEFINE_VALUE_ARRAY_1(abs_mvd_greater1_flag, compIdx, 2 )
DEFINE_VALUE_ARRAY_1(abs_mvd_minus2, compIdx, 2 )
DEFINE_VALUE_ARRAY_1(mvd_sign_flag, compIdx, 2 )
DEFINE_STRUCT_ARITY_1(lMvd, compIdx)
DEFINE_DERIVED_LONG(lMvd)
{
    return h[abs_mvd_greater0_flag(v.compIdx)] *
            ( h[abs_mvd_minus2(v.compIdx)] + 2 ) * 	( 1 - 2 * h[mvd_sign_flag(v.compIdx)] );
}
};

DEFINE_STRUCT_ARITY_3(Mvd, refList, x0, y0)
template <> struct ValueType < Mvd > { typedef struct MotionVector Type; };
template <> struct IsValueArray<Mvd> : std::true_type { };
template <> struct ValueHolder<Mvd> : ValueArray1<Mvd, 2>{};
template <> struct Index<Mvd, 0> { static const int value = 0; };

struct refIdxL0 { };
struct refIdxL1 { };
struct predFlagL0 { };
struct predFlagL1 { };



struct mvpL0 { }; // unused?
template <> struct ValueType<mvpL0> { typedef struct MotionVector Type; };
struct mvpL1 { }; // unused?
template <> struct ValueType<mvpL1> { typedef struct MotionVector Type; };
struct mvL0 { };
template <> struct ValueType<mvL0> { typedef struct MotionVector Type; };
struct mvL1 { };
template <> struct ValueType<mvL1> { typedef struct MotionVector Type; };
struct mvL0A { };
template <> struct ValueType<mvL0A> { typedef struct MotionVector Type; };
struct mvL1A { };
template <> struct ValueType<mvL1A> { typedef struct MotionVector Type; };
struct mvL0B { };
template <> struct ValueType<mvL0B> { typedef struct MotionVector Type; };
struct mvL1B { };
template <> struct ValueType<mvL1B> { typedef struct MotionVector Type; };



DEFINE_STRUCT_ARITY_2(MvL0, x, y);
template <> struct ValueType<MvL0> { typedef struct MotionVector Type; };
DEFINE_STRUCT_ARITY_2(MvL1, x, y);
template <> struct ValueType<MvL1> { typedef struct MotionVector Type; };
DEFINE_STRUCT_ARITY_2(RefIdxL0, x, y);
DEFINE_STRUCT_ARITY_2(RefIdxL1, x, y);
DEFINE_STRUCT_ARITY_2(PredFlagL0, x, y);
DEFINE_STRUCT_ARITY_2(PredFlagL1, x, y);

DEFINE_STRUCT_ARITY_2(predFlagL0Col, x, y)
DEFINE_STRUCT_ARITY_2(predFlagL1Col, x, y)
DEFINE_STRUCT_ARITY_2(refIdxL0Col, x, y)
DEFINE_STRUCT_ARITY_2(refIdxL1Col, x, y)
DEFINE_STRUCT_ARITY_2(mvL0Col, x, y)
template <> struct ValueType<mvL0Col> { typedef struct MotionVector Type; };
DEFINE_STRUCT_ARITY_2(mvL1Col, x, y)
template <> struct ValueType<mvL1Col> { typedef struct MotionVector Type; };

DEFINE_STRUCT_ARITY_7(transform_unit, x0, y0, xBase, yBase, log2TrafoSize, trafoDepth, blkIdx );
struct cu_qp_delta_abs { };
struct cu_qp_delta_sign_flag { };
struct cu_chroma_qp_offset_flag { };
struct cu_chroma_qp_offset_idx { };


DEFINE_STRUCT_ARITY_4(residual_coding, x0, y0, log2TrafoSize, cIdx);
static std::ostream &operator<<(std::ostream &o, residual_coding rc)
{
    return o << "residual_coding(" << rc.x0 << ", " << rc.y0 << ", " << rc.log2TrafoSize << ", " << rc.cIdx << ")";
}

// Pseudo-element to allow specialisations to intecept deleted as well as coded residual_coding() and transform_tree()
template <class V, class F>
struct IfCbf
{
    V cbf;
    F f;
    operator F() const { return this->f; }
};

DEFINE_STRUCT_ARITY_3(transform_skip_flag, x0, y0, cIdx);
DEFINE_STRUCT_ARITY_3(explicit_rdpcm_flag, x0, y0, cIdx);
DEFINE_STRUCT_ARITY_3(explicit_rdpcm_dir_flag, x0, y0, cIdx);
struct last_sig_coeff_x_prefix { };
struct last_sig_coeff_y_prefix { };
struct last_sig_coeff_x_suffix { };
struct last_sig_coeff_y_suffix { };
DEFINE_VALUE_ARRAY_2(coded_sub_block_flag, xS, 8, yS, 8)
DEFINE_STRUCT_ARITY_2(sig_coeff_flag, xC, yC)
DEFINE_STRUCT_ARITY_1(coeff_abs_level_greater1_flag, n)
DEFINE_STRUCT_ARITY_1(coeff_abs_level_greater2_flag, n)
DEFINE_VALUE_ARRAY_1(coeff_sign_flag, n, 16)
DEFINE_STRUCT_ARITY_1(coeff_abs_level_remaining, n)
DEFINE_VALUE_ARRAY_1(StatCoeff, sbType, 4)
struct sbType { };
DEFINE_DERIVED_LONG(sbType)
{
    const residual_coding &rc = *static_cast<residual_coding *>(h);

    if (h[transform_skip_flag(rc.x0, rc.y0, rc.cIdx)] == 0 && h[cu_transquant_bypass_flag()] == 0)
    {
        return 2 * (rc.cIdx == 0 ? 1 : 0);
    }
    else
    {
        return 2 * (rc.cIdx == 0 ? 1 : 0) + 1;
    }
}
};
DEFINE_DERIVED(initRiceValue, h[persistent_rice_adaptation_enabled_flag()] ? h[StatCoeff(h[sbType()])] / 4 : 0)
struct transform_skip_flag_0 { };
struct transform_skip_flag_1 { };
struct transform_skip_flag_2 { };
struct scanIdx { };
DEFINE_DERIVED_LONG(scanIdx)
    {
        coding_quadtree const *cqt = h;
        residual_coding const *rc = h;
        if (h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA)
        {
            if (rc->log2TrafoSize == 2 || (rc->log2TrafoSize == 3 && rc->cIdx == 0) || (rc->log2TrafoSize == 3 && h[ChromaArrayType()] == 3))
            {
                int const predModeIntra = rc->cIdx ? h[IntraPredModeC(rc->x0, rc->y0)] : h[IntraPredModeY(rc->x0, rc->y0)];
                static uint8_t const lookup[35] = { 0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0 };
                return lookup[predModeIntra];
            }
        }
        return 0;
    }
};
struct LastSignificantCoeffX { };
struct LastSignificantCoeffY { };
DEFINE_STRUCT_ARITY_5(TransCoeffLevel, x0, y0, cIdx, xC, yC);
template <> struct IsValueArray<TransCoeffLevel> : std::true_type { };
template <> struct ValueHolder<TransCoeffLevel> : ValueArray2<TransCoeffLevel, 32, 32> { };
template <> struct Index<TransCoeffLevel, 0> { static const int value = 3; };
template <> struct Index<TransCoeffLevel, 1> { static const int value = 4; };

DEFINE_STRUCT_ARITY_3(cross_comp_pred, x0, y0, c);
DEFINE_VALUE_ARRAY_1(log2_res_scale_abs_plus1, c, 3)
DEFINE_VALUE_ARRAY_1(res_scale_sign_flag, c, 3)


DEFINE_STRUCT_ARITY_0(vui_parameters)
struct aspect_ratio_info_present_flag { };
struct aspect_ratio_idc { };
struct sar_width { };
struct sar_height { };
struct overscan_info_present_flag { };
struct overscan_appropriate_flag { };
struct video_signal_type_present_flag { };
struct video_format { };
struct video_full_range_flag { };
struct colour_description_present_flag { };
struct colour_primaries { };
struct transfer_characteristics { };
struct matrix_coeffs { };
struct chroma_loc_info_present_flag { };
struct chroma_sample_loc_type_top_field { };
struct chroma_sample_loc_type_bottom_field { };
struct neutral_chroma_indication_flag { };
struct field_seq_flag { };
struct frame_field_info_present_flag { };
struct default_display_window_flag { };
struct def_disp_win_left_offset { };
struct def_disp_win_right_offset { };
struct def_disp_win_top_offset { };
struct def_disp_win_bottom_offset { };
struct vui_timing_info_present_flag { };
struct vui_num_units_in_tick { };
template <> struct ValueType<vui_num_units_in_tick> { typedef unsigned Type; };
struct vui_time_scale { };
template <> struct ValueType<vui_time_scale> { typedef unsigned Type; };
struct vui_poc_proportional_to_timing_flag { };
struct vui_num_ticks_poc_diff_one_minus1 { };
template <> struct ValueType<vui_num_ticks_poc_diff_one_minus1> { typedef unsigned Type; };
struct vui_hrd_parameters_present_flag { };
struct bitstream_restriction_flag { };
struct tiles_fixed_structure_flag { };
struct motion_vectors_over_pic_boundaries_flag { };
struct restricted_ref_pic_lists_flag { };
struct min_spatial_segmentation_idc { };
struct max_bytes_per_pic_denom { };
struct max_bits_per_min_cu_denom { };
struct log2_max_mv_length_horizontal { };
struct log2_max_mv_length_vertical { };


DEFINE_STRUCT_ARITY_2(hrd_parameters, commonInfPresentFlag, maxNumSubLayersMinus1 );
struct nal_hrd_parameters_present_flag { };
struct vcl_hrd_parameters_present_flag { };
struct sub_pic_hrd_params_present_flag { };
struct tick_divisor_minus2 { };
struct du_cpb_removal_delay_increment_length_minus1 { };
struct sub_pic_cpb_params_in_pic_timing_sei_flag { };
struct dpb_output_delay_du_length_minus1 { };
struct bit_rate_scale { };
struct cpb_size_scale { };
struct cpb_size_du_scale { };
struct initial_cpb_removal_delay_length_minus1 { };
struct au_cpb_removal_delay_length_minus1 { };
struct dpb_output_delay_length_minus1 { };
DEFINE_VALUE_ARRAY_1(fixed_pic_rate_general_flag, i, 7)
DEFINE_VALUE_ARRAY_1(fixed_pic_rate_within_cvs_flag, i, 7)
DEFINE_VALUE_ARRAY_1(elemental_duration_in_tc_minus1, i, 7)
DEFINE_VALUE_ARRAY_1(low_delay_hrd_flag, i, 7)
DEFINE_VALUE_ARRAY_1(cpb_cnt_minus1, i, 7)
DEFINE_STRUCT_ARITY_1(CpbCnt, i)
DEFINE_DERIVED_LONG(CpbCnt)
{
    return h[cpb_cnt_minus1(v.i)];
}
};


DEFINE_STRUCT_ARITY_1(sub_layer_hrd_parameters, subLayerId)
template <> struct ValueType<struct bit_rate_value_minus1> { typedef unsigned Type; };
DEFINE_VALUE_ARRAY_1(bit_rate_value_minus1, i, 32)
template <> struct ValueType<struct cpb_size_value_minus1> { typedef unsigned Type; };
DEFINE_VALUE_ARRAY_1(cpb_size_value_minus1, i, 32)
template <> struct ValueType<struct cpb_size_du_value_minus1> { typedef unsigned Type; };
DEFINE_VALUE_ARRAY_1(cpb_size_du_value_minus1, i, 32)
template <> struct ValueType<struct bit_rate_du_value_minus1> { typedef unsigned Type; };
DEFINE_VALUE_ARRAY_1(bit_rate_du_value_minus1, i, 32)
DEFINE_VALUE_ARRAY_1(cbr_flag, i, 32)






DEFINE_VALUE_ARRAY_1(PocLsbLt, i, 32);

DEFINE_VALUE_ARRAY_2(DeltaPocS0, stRpsIdx , 65, i, 16);
DEFINE_VALUE_ARRAY_2(DeltaPocS1, stRpsIdx , 65, i, 16);





template <class V>
struct VectorOffset
{
    static const int value = 0;
};

template <class V>
struct ValueVector1
{
    std::vector<typename ValueType<V>::Type> value;
    typename ValueType<V>::Type &get(V v)
    {
        return this->value[v.argValue(0) + VectorOffset<V>::value];
    }
};

template <> struct ValueHolder<colBd> : ValueVector1<colBd> { };
template <> struct ValueHolder<rowBd> : ValueVector1<rowBd> { };
template <> struct ValueHolder<CtbAddrRsToTs> : ValueVector1<CtbAddrRsToTs> { };
template <> struct VectorOffset<CtbAddrRsToTs> { static const int value = 1; };
template <> struct ValueHolder<CtbAddrTsToRs> : ValueVector1<CtbAddrTsToRs> { };
template <> struct ValueHolder<TileId> : ValueVector1<TileId> { };
template <> struct VectorOffset<TileId> { static const int value = 1; };

template <class V>
struct ContainerOf { };

template <class V, class H>
struct Access<ContainerOf<V>, H, typename std::enable_if<std::is_base_of<ValueVector1<V>, H>::value>::type>
{
    typedef std::vector<typename ValueType<V>::Type> &Type;
    static Type get(ContainerOf<V>, ValueVector1<V> &h)
    {
        return h.value;
    }
};


// Thrown if parsing cannot continue
struct Abort { };

static int cartesianToZ(int x, int y)
{
    int z = 0;
    int shift = 0;
    while (x || y)
    {
        z |= (x & 1) << shift;
        x >>= 1;
        ++shift;
        z |= (y & 1) << shift;
        y >>= 1;
        ++shift;
    }
    return z;
}

struct Stop { };

template <class H>
int computeMaxDpbSize(H &h, int MaxLumaPs)
{
    const int maxDpbPicBuf = 6;
    if (h[PicSizeInSamplesY()] <= (MaxLumaPs >> 2))
        return std::min(4 * maxDpbPicBuf, 16);
    else if (h[PicSizeInSamplesY()] <= (MaxLumaPs >> 1))
        return std::min(2 * maxDpbPicBuf, 16);
    else if (h[PicSizeInSamplesY()] <= ((3 * MaxLumaPs) >> 2))
        return std::min((4 * maxDpbPicBuf) / 3, 16);
    else
        return maxDpbPicBuf;
}


DEFINE_STRUCT_ARITY_0(ContextsInitialize);
DEFINE_STRUCT_ARITY_1(ContextsSave, i);
DEFINE_STRUCT_ARITY_1(ContextsRestore, i);

static int QpC(int qPi)
{
    if (qPi < 30) return qPi;
    if (qPi > 42) return qPi - 6;
    static const int lookup[13] = { 29, 30, 31, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37 };
    return lookup[qPi - 30];
}




DEFINE_STRUCT_ARITY_1(RefPicList, x)

// This probably should just be member function of a state object.
DEFINE_STRUCT_ARITY_0(CabacRestart);

template <typename T, class H, class Enable=void>
struct NextBits;

template <typename T, class H>
T next_bits(H &h, int n, bool aligned=false)
{
    return NextBits<T, H>::go(h, n, aligned);
}

template <typename T, class H, class Enable=void>
struct ReadBytes;

template <typename T, class H>
T read_bytes(H &h, int n)
{
    return ReadBytes<T, H>::go(h, n);
}



// A pseudo syntax function which is a useful construct during cost measurement when encoding.
// Note: excludes rate of intra_chroma_pred_mode()
// There is one IntraPartitionPrediction per CU for PART_2Nx2N
// There are four IntraPartitionPredictions per CU for PART_NxN
// review: consider removing (transform_tree does this job)
DEFINE_STRUCT_ARITY_4(IntraPartitionPrediction, x0, y0, log2CbSize, blkIdx)


// A pseudo syntax function which is a useful construct during cost measurement when encoding.
// review: consider removing (transform_tree does this job)
DEFINE_STRUCT_ARITY_5(IntraPartition, x0, y0, log2CbSize, split, blkIdx)

// Namespace necessary here because Windows headers have a "Turing::Rectangle" too
// review: better fix - untangle dependencies so that no compilation unit has both Windows and Global.h
namespace Turing {
    struct Rectangle
    {
        int x0;
        int y0;
        int width;
        int height;
    };
}

static inline int xPositionOf(const Turing::Rectangle& rectangle) { return rectangle.x0; }
static inline int yPositionOf(const Turing::Rectangle& rectangle) { return rectangle.y0; }
static inline int widthOf(const Turing::Rectangle& rectangle) { return rectangle.width; }
static inline int heightOf(const Turing::Rectangle& rectangle) { return rectangle.height; }

static inline int xPositionOf(const coding_quadtree& cqt) { return cqt.x0; }
static inline int yPositionOf(const coding_quadtree& cqt) { return cqt.y0; }
static inline int widthOf(const coding_quadtree& cqt) { return 1 << cqt.log2CbSize; }
static inline int heightOf(const coding_quadtree& cqt) { return 1 << cqt.log2CbSize; }

static inline int xPositionOf(const coding_unit& cu) { return cu.x0; }
static inline int yPositionOf(const coding_unit& cu) { return cu.y0; }
static inline int widthOf(const coding_unit& cu) { return 1 << cu.log2CbSize; }
static inline int heightOf(const coding_unit& cu) { return 1 << cu.log2CbSize; }

static inline int xPositionOf(const prediction_unit& pu) { return pu.x0; }
static inline int yPositionOf(const prediction_unit& pu) { return pu.y0; }
static inline int widthOf(const prediction_unit& pu) { return pu.nPbW; }
static inline int heightOf(const prediction_unit& pu) { return pu.nPbH; }

static inline int xPositionOf(const transform_tree& tt) { return tt.x0; }
static inline int yPositionOf(const transform_tree& tt) { return tt.y0; }
static inline int widthOf(const transform_tree& tt) { return 1 << tt.log2TrafoSize; }
static inline int heightOf(const transform_tree& tt) { return 1 << tt.log2TrafoSize; }

static inline int xPositionOf(const IntraPartition& e) { return e.x0 + ((e.blkIdx & 1) << (e.log2CbSize - e.split)); }
static inline int yPositionOf(const IntraPartition& e) { return e.y0 + ((e.blkIdx & 2) << (e.log2CbSize - e.split - 1)); }
static inline int widthOf(const IntraPartition& e) { return 1 << (e.log2CbSize - e.split); }
static inline int heightOf(const IntraPartition& e) { return 1 << (e.log2CbSize - e.split); }

static inline int xPositionOf(const residual_coding& e) { return e.x0; }
static inline int yPositionOf(const residual_coding& e) { return e.y0; }
template <class E> static int zPositionOf(const E &e) { return cartesianToZ(xPositionOf(e), yPositionOf(e)); }

struct numGreater1Flag { };


struct ReconstructedSamples
{
    ReconstructedSamples(int x, int y, int cIdx) : x(x), y(y), cIdx(cIdx) {}
    int x; int y; int cIdx;
};


// review: these operator<< functions could be automatically generated

static std::ostream &operator<<(std::ostream &o, IntraPartition e)
{
    return o << "IntraPartition(" << e.x0 << ", " << e.y0 << ", " << e.log2CbSize << ", " << e.split << ", " << e.blkIdx << ")";
}

static std::ostream &operator<<(std::ostream &o, coding_unit cu)
{
    return o << "coding_unit(" << cu.x0 << ", " << cu.y0 << ", " << cu.log2CbSize << ")";
}

static std::ostream &operator<<(std::ostream &o, coding_quadtree cqt)
{
    return o << "coding_quadtree(" << cqt.x0 << ", " << cqt.y0 << ", " << cqt.log2CbSize << ", " << cqt.cqtDepth << ")";
}

static std::ostream &operator<<(std::ostream &o, prediction_unit pu)
{
    return o << "prediction_unit(" << pu.x0 << ", " << pu.y0 << ", " << pu.nPbW << ", " << pu.nPbH << ")";
}

static std::ostream &operator<<(std::ostream &o, transform_tree tt)
{
    return o << "transform_tree(" << tt.x0 << ", " << tt.y0 << ", " << tt.xBase << ", " << tt.yBase << ", " << tt.log2TrafoSize << ", " << tt.trafoDepth << ", " << tt.blkIdx << ")";
}

static std::ostream &operator<<(std::ostream &o, mvd_coding mvdc)
{
    return o << "mvd_coding(" << mvdc.x0 << ", " << mvdc.y0 << ", " << mvdc.refList << ")";
}

static std::ostream &operator<<(std::ostream &o, Turing::Rectangle rectangle)
{
    return o << "Rectangle(" << rectangle.x0 << ", " << rectangle.y0 << ", " << rectangle.width << ", " << rectangle.height << ")";
}


enum When { before = 0, after = 1, corner };

const int ctbWaitingOffsetX = 4;
const int ctbWaitingOffsetY = 3;

#endif
