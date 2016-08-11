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

#pragma once

#include "../SyntaxSei.h"


struct colour_remap_id { };
struct colour_remap_cancel_flag { };
struct colour_remap_persistence_flag { };
struct colour_remap_video_signal_info_present_flag { };
struct colour_remap_full_range_flag { };
struct colour_remap_primaries { };
struct colour_remap_transfer_function { };
struct colour_remap_matrix_coefficients { };
struct colour_remap_input_bit_depth { };
struct colour_remap_bit_depth { };
struct colour_remap_matrix_present_flag { };
struct log2_matrix_denom { };
DEFINE_STRUCT_ARITY_1(pre_lut_num_val_minus1, c);
DEFINE_STRUCT_ARITY_1(post_lut_num_val_minus1, c);
DEFINE_STRUCT_ARITY_2(pre_lut_coded_value, c, i);
DEFINE_STRUCT_ARITY_2(pre_lut_target_value, c, i);
DEFINE_STRUCT_ARITY_2(colour_remap_coeffs, c, i);
DEFINE_STRUCT_ARITY_2(post_lut_coded_value, c, i);
DEFINE_STRUCT_ARITY_2(post_lut_target_value, c, i);


template <>
struct Syntax<colour_remapping_info>
{
    template <class H> static void go(colour_remapping_info const &fun, H &h)
    {
        h(colour_remap_id(), ue(v));
        h(colour_remap_cancel_flag(), u(1));
        if (!h[colour_remap_cancel_flag()])
        {
            h(colour_remap_persistence_flag(), u(1));
            h(colour_remap_video_signal_info_present_flag(), u(1));
            if (h[colour_remap_video_signal_info_present_flag()])
            {
                h(colour_remap_full_range_flag(), u(1));
                h(colour_remap_primaries(), u(8));
                h(colour_remap_transfer_function(), u(8));
                h(colour_remap_matrix_coefficients(), u(8));
            }
            h(colour_remap_input_bit_depth(), u(8));
            h(colour_remap_bit_depth(), u(8));
            for (int c = 0; c < 3; c++)
            {
                h(pre_lut_num_val_minus1(c), u(8));
                if (h[pre_lut_num_val_minus1(c)] > 0)
                {
                    for (int i = 0; i <= h[pre_lut_num_val_minus1(c)]; i++)
                    {
                        h(pre_lut_coded_value(c, i), uv());
                        h(pre_lut_target_value(c, i), uv());
                    }
                }
            }
            h(colour_remap_matrix_present_flag(), u(1));
            if (h[colour_remap_matrix_present_flag()])
            {
                h(log2_matrix_denom(), u(4));
                for (int c = 0; c < 3; c++)
                    for (int i = 0; i < 3; i++)
                        h(colour_remap_coeffs(c, i), se(v));
            }
            for (int c = 0; c < 3; c++)
            {
                h(post_lut_num_val_minus1(c), u(8));
                if (h[post_lut_num_val_minus1(c)] > 0)
                {
                    for (int i = 0; i <= h[post_lut_num_val_minus1(c)]; i++)
                    {
                        h(post_lut_coded_value(c, i), uv());
                        h(post_lut_target_value(c, i), uv());
                    }
                }
            }
        }
    }
};


NUMBER_OF_BITS_UV(pre_lut_coded_value, ((h[colour_remap_input_bit_depth()] + 7) >> 3) << 3)
NUMBER_OF_BITS_UV(pre_lut_target_value, ((h[colour_remap_bit_depth()] + 7) >> 3) << 3)
NUMBER_OF_BITS_UV(post_lut_coded_value, ((h[colour_remap_bit_depth()] + 7) >> 3) << 3)
NUMBER_OF_BITS_UV(post_lut_target_value, ((h[colour_remap_bit_depth()] + 7) >> 3) << 3)


struct ColourRemappingInfo :
    ValueHolder<colour_remap_id>,
    ValueHolder<colour_remap_cancel_flag>,
    ValueHolder<colour_remap_persistence_flag>,
    ValueHolder<colour_remap_video_signal_info_present_flag>,
    ValueHolder<colour_remap_full_range_flag>,
    ValueHolder<colour_remap_primaries>,
    ValueHolder<colour_remap_transfer_function>,
    ValueHolder<colour_remap_matrix_coefficients>,
    ValueHolder<colour_remap_input_bit_depth>,
    ValueHolder<colour_remap_bit_depth>,
    ValueHolder<colour_remap_matrix_present_flag>,
    ValueHolder<log2_matrix_denom>,
    ValueHolder<pre_lut_num_val_minus1>,
    ValueHolder<post_lut_num_val_minus1>,
    ValueHolder<pre_lut_coded_value>,
    ValueHolder<pre_lut_target_value>,
    ValueHolder<colour_remap_coeffs>,
    ValueHolder<post_lut_coded_value>,
    ValueHolder<post_lut_target_value>
    {
    };


template <class H> void Read<colour_remapping_info>::go(colour_remapping_info f, H &h)
{
    ColourRemappingInfo colourRemappingInfo;
    auto hh = h.extend(&colourRemappingInfo);
    Syntax<colour_remapping_info>::go(f, hh);
}


