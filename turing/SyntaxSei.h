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

// Review: remove parsing from here

#ifndef INCLUDED_SyntaxSei_h
#define INCLUDED_SyntaxSei_h

#pragma once

#include "Syntax.h"
#include "HevcTypes.h"
#include "Global.h"
#include "ScanOrder.h"
#include "SyntaxElements.h"
#include "Violation.h"
#include "StatePictures.h"
#include <cassert>


#define SEI_PAYLOAD_TYPES \
        BOTH(-1, reserved_sei_message) \
        PREFIX(0, buffering_period) \
        PREFIX(1, pic_timing) \
        PREFIX(2, pan_scan_rect) \
        BOTH(3, filler_payload) \
        BOTH(4, user_data_registered_itu_t_t35) \
        BOTH(5, user_data_unregistered) \
        PREFIX(6, recovery_point) \
        PREFIX(9, scene_info) \
        PREFIX(15, picture_snapshot) \
        PREFIX(16, progressive_refinement_segment_start) \
        BOTH(17, progressive_refinement_segment_end) \
        PREFIX(19, film_grain_characteristics) \
        BOTH(22, post_filter_hint) \
        PREFIX(23, tone_mapping_info) \
        PREFIX(45, frame_packing_arrangement) \
        PREFIX(47, display_orientation) \
        PREFIX(128, structure_of_pictures_info) \
        PREFIX(129, active_parameter_sets) \
        PREFIX(130, decoding_unit_info) \
        PREFIX(131, temporal_sub_layer_zero_index) \
        SUFFIX(132, decoded_picture_hash) \
        PREFIX(133, scalable_nesting) \
        PREFIX(134, region_refresh_info) \
        PREFIX(135, no_display) \
        PREFIX(136, time_code) \
        PREFIX(137, mastering_display_colour_volume) \
        PREFIX(138, segmented_rect_frame_packing_arrangement) \
        PREFIX(139, temporal_motion_constrained_tile_sets) \
        PREFIX(140, chroma_resampling_filter_hint) \
        PREFIX(141, knee_function_info) \
        PREFIX(142, colour_remapping_info) \
        PREFIX(143, deinterlaced_field_identification) \
        PREFIX(144, content_light_level) \
        PREFIX(147, alternative_transfer_characteristics) \
        PREFIX(160, layers_not_present) \
        PREFIX(161, inter_layer_constrained_tile_sets) \
        PREFIX(162, bsp_nesting) \
        PREFIX(163, bsp_initial_arrival_time) \
        PREFIX(164, sub_bitstream_property) \
        PREFIX(165, alpha_channel_info) \
        PREFIX(166, overlay_info) \
        PREFIX(167, temporal_mv_prediction_constraints) \
        PREFIX(168, frame_field_info) \
        PREFIX(176, three_dimensional_reference_displays_info) \
        PREFIX(177, depth_representation_info) \
        PREFIX(178, multiview_scene_info) \
        PREFIX(179, multiview_acquisition_info) \
        PREFIX(180, multiview_view_position) \


template <class Type> struct PayloadTypeOf;


#define PREFIX(n, Type) DEFINE_STRUCT_ARITY_1(Type, payloadSize) template <> struct PayloadTypeOf<Type> { static const int value=n; };
#define SUFFIX(n, Type)  DEFINE_STRUCT_ARITY_1(Type, payloadSize) template <> struct PayloadTypeOf<Type> { static const int value=n; };
#define BOTH(n, Type) DEFINE_STRUCT_ARITY_1(Type, payloadSize)  template <> struct PayloadTypeOf<Type> { static const int value=n; };

SEI_PAYLOAD_TYPES

#undef PREFIX
#undef SUFFIX
#undef BOTH


template <class Type>
int payloadTypeOf()
{
    const int payloadType = PayloadTypeOf<Type>::value;
    return payloadType;
}


static const char *payloadTypeName(int payloadType)
{
    switch (payloadType)
    {
#define PREFIX(n, Type) case n: return #Type;
#define SUFFIX(n, Type)  case n: return #Type;
#define BOTH(n, Type) case n: return #Type;

        SEI_PAYLOAD_TYPES

#undef PREFIX
#undef SUFFIX
#undef BOTH

        default: return "reserved_sei_message";
    }
}


struct reserved_payload_extension_data { };
struct payload_bit_equal_to_one { };
template <> struct Fixed<payload_bit_equal_to_one> { static const int value = 1; }; // CondCheck D.3.1-B
struct payload_bit_equal_to_zero { };
template <> struct Fixed<payload_bit_equal_to_zero> { static const int value = 0; }; // CondCheck D.3.1-C

template <class F> struct Read;

struct Stream;
struct ExceptionOverrun;

template <>
struct Read<Element<reserved_payload_extension_data, uv>>
{
    template <class H> static void go(Element<reserved_payload_extension_data, uv> e, H &h)
    {
        h(Violation("D.3.1", "SEI payload contains reserved_payload_extension_data")); // CondCheck D.3.1-A

        h[::Stream()].state.p = h[::Stream()].end - 1;
        h[::Stream()].state.mask = 0x01;
        while (!(h[::Stream()].state.mask & *h[::Stream()].state.p))
        {
            h[::Stream()].state.mask <<= 1;
            if (!(h[::Stream()].state.mask & 0xff))
            {
                h(Violation("", "payload_bit_equal_to_one not present"));
                throw Abort();
            }
        }
    }
};


struct bitstream_subset_flag { };
struct nesting_op_flag { };
struct default_op_flag { };
struct nesting_num_ops_minus1 { };
DEFINE_VALUE_ARRAY_1(nesting_max_temporal_id_plus1, i, 1024);
DEFINE_VALUE_ARRAY_1(nesting_op_idx, i, 1024);
struct all_layers_flag { };
struct nesting_no_op_max_temporal_id_plus1 { };
struct nesting_num_layers_minus1 { };
DEFINE_VALUE_ARRAY_1(nesting_layer_id, i, 64);
struct nesting_zero_bit { };
template <> struct Fixed<nesting_zero_bit> { static const int value = 0; }; // CondCheck


struct sei_ols_idx { };
struct sei_partitioning_scheme_idx { };
struct bsp_idx { };
struct bsp_nesting_zero_bit { };
struct num_seis_in_bsp_minus1 { };


struct ScalableNesting :
    ValueHolder<bitstream_subset_flag>,
    ValueHolder<nesting_op_flag>,
    ValueHolder<default_op_flag>,
    ValueHolder<nesting_num_ops_minus1>,
    ValueHolder<nesting_max_temporal_id_plus1>,
    ValueHolder<nesting_op_idx>,
    ValueHolder<all_layers_flag>,
    ValueHolder<nesting_no_op_max_temporal_id_plus1>,
    ValueHolder<nesting_num_layers_minus1>,
    ValueHolder<nesting_layer_id>,
    ValueHolder<nesting_zero_bit>
    {
        bool nested;
    };


struct BspNesting :
    ValueHolder<sei_ols_idx>,
    ValueHolder<sei_partitioning_scheme_idx>,
    ValueHolder<bsp_idx>,
    ValueHolder<bsp_nesting_zero_bit>,
    ValueHolder<num_seis_in_bsp_minus1>
    {
        bool nested;
    };


struct StateSei :
    ValueHolder<last_payload_type_byte>,
    ValueHolder<last_payload_size_byte>,
    sei_payload,
    ValueHolder<payload_bit_equal_to_one>,
    ValueHolder<payload_bit_equal_to_zero>,
    ValueHolder<reserved_payload_extension_data>,
    ScalableNesting,
    BspNesting
    {
    };


template <>
struct Read<sei_payload>
{
    template <class H> static void go(const sei_payload &f, H &h)
    {
        StateSei *stateSei = h;
        stateSei->ScalableNesting::nested = false;
        stateSei->BspNesting::nested = false;

        auto &stream = h[::Stream()];
        const int rbspBytes = static_cast<int>(stream.end - stream.state.p);
        if (f.payloadSize > rbspBytes)
        {
            h(Violation("7.4.6", "sei_payload()'s payloadSize{%1%} is greater than number of available RBSP bytes {%2%}") % f.payloadSize % rbspBytes);
            return;
        }

        auto streamCopy = h[::Stream()];

        h[::Stream()].begin = h[::Stream()].state.p;
        h[::Stream()].end = h[::Stream()].state.p + f.payloadSize;

        auto &h2 = h;
        *static_cast<sei_payload *>(h2) = f;

        const char *streamTypePrevious = static_cast<StatePicturesBase *>(h)->streamType;
        static_cast<StatePicturesBase *>(h)->streamType = StatePicturesBase::streamTypeSei();
        try
        {
            Syntax<sei_payload>::go(f, h2);
        }
        catch (ExceptionOverrun &)
        {
            h2(Violation("D.2","SEI parsing overrun - attempting to read data after payloadSize{%1%} bytes") % f.payloadSize);
        }
        catch (Abort &)
        {
        }

        static_cast<StatePicturesBase *>(h)->streamType = streamTypePrevious;

        streamCopy.state.p = h[::Stream()].end;
        streamCopy.state.mask = 0x80;

        h[::Stream()] = streamCopy;
    }
};


template <class V, class M> struct ReadU;


// Parsing of an element whose bit length is set according to another "_minus1" element's value.
template <class V, class L>
struct ReadUv
{
    template <class H> static void go(Element<V, u> &e, H &h)
    {
        e.m.n = h[L()] + 1;
        ReadU<V, u>::go(e, h);
    }
};


// Useful shorthand to impement parsing when the bit length of an element is set according to another "_minus1" element's value.
#define PARSE_UV(V, L) template <> struct Read<Element<V, u>> : ReadUv<V, L> { };


// Declare Read<F> for all SEI payload types
#define DECLARE_READ_SEI(F) template <> struct Read<F> { template <class H> static void go(F, H&); };
#define PREFIX(id, F) DECLARE_READ_SEI(F)
#define SUFFIX(id, F) DECLARE_READ_SEI(F)
#define BOTH(id, F) DECLARE_READ_SEI(F)

SEI_PAYLOAD_TYPES

#undef DECLARE_READ_SEI
#undef PREFIX
#undef SUFFIX
#undef BOTH


// review: the following needs to be implemented
#define UNIMPLEMENTED_SEI(F) \
template <class H> void Read<F>::go(F f, H &h) \
{ \
    seek(h, 8 * f.payloadSize); \
} \

UNIMPLEMENTED_SEI(bsp_nesting);
UNIMPLEMENTED_SEI(bsp_initial_arrival_time);
UNIMPLEMENTED_SEI(sub_bitstream_property);
UNIMPLEMENTED_SEI(three_dimensional_reference_displays_info);
UNIMPLEMENTED_SEI(depth_representation_info);
UNIMPLEMENTED_SEI(multiview_scene_info);
UNIMPLEMENTED_SEI(multiview_acquisition_info);
UNIMPLEMENTED_SEI(multiview_view_position);

#endif
