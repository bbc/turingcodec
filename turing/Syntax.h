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

// Definition of HEVC syntax, all in one place in a template form that can be used by encoders, decoders, parsers and analysers.
// These functions be very much like the syntax tables in the HEVC standard - with no need to modify unless syntax is changed.

#ifndef INCLUDED_Syntax_h
#define INCLUDED_Syntax_h

#pragma once

#include "HevcTypes.h"
#include "Global.h"
#include <cassert>


template <typename F=void> struct Syntax;

template <> struct Syntax<void> { };


// Raw byte sequence payloads, trailing bits, and byte alignment syntax

template <>
struct Syntax<video_parameter_set_rbsp>
{
    template <class H> static void go(const video_parameter_set_rbsp &fun, H &h);
};


template <>
struct Syntax<seq_parameter_set_rbsp>
{
    template <class H> static void go(const seq_parameter_set_rbsp &fun, H &h);
};


template <>
struct Syntax<sps_range_extension>
{
    template <class H> static void go(const sps_range_extension &fun, H &h);
};


template <>
struct Syntax<pic_parameter_set_rbsp>
{
    template <class H> static void go(const pic_parameter_set_rbsp &fun, H &h);
};


template <>
struct Syntax<pps_range_extension>
{
    template <class H> static void go(const pps_range_extension &fun, H &h);
};


template <>
struct Syntax<sei_rbsp>
{
    template <class H> static void go(const sei_rbsp &fun, H &h);
};


template <>
struct Syntax<access_unit_delimiter_rbsp>
{
    template <class H> static void go(const access_unit_delimiter_rbsp &fun, H &h);
};


template <>
struct Syntax<end_of_seq_rbsp>
{
    template <class H> static void go(const end_of_seq_rbsp &fun, H &h);
};


template <>
struct Syntax<end_of_bitstream_rbsp>
{
    template <class H> static void go(const end_of_bitstream_rbsp &fun, H &h);
};


template <>
struct Syntax<filler_data_rbsp>
{
    template <class H> static void go(const filler_data_rbsp &fun, H &h);
};


template <>
struct Syntax<slice_segment_layer_rbsp>
{
    template <class H> static void go(const slice_segment_layer_rbsp &, H&);
};


template <>
struct Syntax<rbsp_slice_segment_trailing_bits>
{
    template <class H> static void go(const rbsp_slice_segment_trailing_bits&, H&);
};


template <>
struct Syntax<rbsp_trailing_bits>
{
    template <class H> static void go(const rbsp_trailing_bits&, H&);
};


template <>
struct Syntax<byte_alignment>
{
    template <class H> static void go(const byte_alignment&, H&);
};


template <>
struct Syntax<profile_tier_level>
{
    template <class H> static void go(const profile_tier_level&, H&);
};


template <>
struct Syntax<scaling_list_data>
{
    template <class H> static void go(const scaling_list_data&, H&);
};


template <>
struct Syntax<slice_segment_header>
{
    template <class H> static void go(const slice_segment_header&, H&);
};


template <>
struct Syntax<ref_pic_lists_modification>
{
    template <class H> static void go(const ref_pic_lists_modification&, H&);
};


template <>
struct Syntax<pred_weight_table>
{
    template <class H> static void go(const pred_weight_table&, H&);
};


template <>
struct Syntax<short_term_ref_pic_set>
{
    template <class H> static void go(const short_term_ref_pic_set&, H&);
};


// Slice segment data syntax

template <>
struct Syntax<slice_segment_data>
{
    template <class H> static void go(const slice_segment_data&, H&);
    template <class H> static void iteration(H &h);
};


template <>
struct Syntax<coding_tree_unit>
{
    template <class H> static void go(const coding_tree_unit&, H&);
};


template <>
struct Syntax<sao>
{
    template <class H> static void go(const sao&, H&);
};


template <>
struct Syntax<coding_quadtree>
{
    template <class H> static void go(const coding_quadtree&, H&);
};


template <class Direction> struct Syntax<Deleted<coding_quadtree, Direction>> : Null<Deleted<coding_quadtree, Direction>>{};


template <>
struct Syntax < coding_unit >
{
    template <class H> static void go(const coding_unit &cu, H &h);
};


template <>
struct Syntax<prediction_unit>
{
    template <class H> static void go(const prediction_unit&, H&);
};


template <>
struct Syntax<pcm_sample>
{
    template <class H> static void go(const pcm_sample&, H&);
};


template <class V, class F>
struct Syntax<IfCbf<V, F>>
{
    template <class H> static void go(const IfCbf<V, F> &f, H &h)
    {
        if (h[f.cbf]) h(f.f);
    }
};


template <>
struct Syntax<transform_tree>
{
    template <class H> static void go(const transform_tree&, H&);
};


template <>
struct Syntax<mvd_coding>
{
    template <class H> static void go(const mvd_coding&, H&);
};


template <>
struct Syntax<transform_unit>
{
    template <class H> static void go(const transform_unit&, H&);
};


template <>
struct Syntax<residual_coding>
{
    template <class H> static void go(const residual_coding&, H&);
};


// Annex B - Byte stream format

template <>
struct Syntax<byte_stream_nal_unit>
{
    template <class H> static void go(const byte_stream_nal_unit&, H&);
};


template <>
struct Syntax<nal_unit>
{
    template <class H> static void go(const nal_unit&, H&);
};


template <>
struct Syntax<nal_unit_header>
{
    template <class H> static void go(const nal_unit_header&, H&);
};


// Annex E - Video usability information

template <>
struct Syntax<vui_parameters>
{
    template <class H> static void go(const vui_parameters&, H&);
};


template <>
struct Syntax<hrd_parameters>
{
    template <class H> static void go(const hrd_parameters&, H&);
};


template <>
struct Syntax<sub_layer_hrd_parameters>
{
    template <class H> static void go(const sub_layer_hrd_parameters&, H&);
};


template <> struct Syntax<struct RbspReserved> : Null<RbspReserved>{};

template <> struct Syntax<struct RbspUnspecified> : Null<RbspUnspecified>{};

template <> struct Syntax<struct Violation> : Null<Violation>{};

template <> struct Syntax<struct JustDecoded> : Null<JustDecoded>{};


template <>
struct Syntax<sei_message>
{
    template <class H> static void go(const sei_message &, H &);
};

template <>
struct Syntax<sei_payload>
{
    template <class H> static void go(sei_payload, H &);
};

#endif
