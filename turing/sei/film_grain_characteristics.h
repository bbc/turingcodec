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







struct film_grain_characteristics_cancel_flag { };
struct film_grain_model_id { };
struct separate_colour_description_present_flag { };
struct film_grain_bit_depth_luma_minus8 { };
struct film_grain_bit_depth_chroma_minus8 { };
struct film_grain_full_range_flag { };
struct film_grain_colour_primaries { };
struct film_grain_transfer_characteristics { };
struct film_grain_matrix_coeffs { };
struct blending_mode_id { };
struct log2_scale_factor { };
DEFINE_VALUE_ARRAY_1(comp_model_present_flag, c, 3);
DEFINE_VALUE_ARRAY_1(num_intensity_intervals_minus1, c, 3);
DEFINE_VALUE_ARRAY_1(num_model_values_minus1, c, 3);
DEFINE_VALUE_ARRAY_2(intensity_interval_lower_bound, c, 3, i, 256);
DEFINE_VALUE_ARRAY_2(intensity_interval_upper_bound, c, 3, i, 256);
DEFINE_VALUE_ARRAY_3(comp_model_value, c, 3, i, 256, j, 8);
struct film_grain_characteristics_persistence_flag { };

struct FilmGrainCharacteristics :
    ValueHolder<film_grain_characteristics_cancel_flag>,
    ValueHolder<film_grain_model_id>,
    ValueHolder<separate_colour_description_present_flag>,
    ValueHolder<film_grain_bit_depth_luma_minus8>,
    ValueHolder<film_grain_bit_depth_chroma_minus8>,
    ValueHolder<film_grain_full_range_flag>,
    ValueHolder<film_grain_colour_primaries>,
    ValueHolder<film_grain_transfer_characteristics>,
    ValueHolder<film_grain_matrix_coeffs>,
    ValueHolder<blending_mode_id>,
    ValueHolder<log2_scale_factor>,
    ValueHolder<comp_model_present_flag>,
    ValueHolder<num_intensity_intervals_minus1>,
    ValueHolder<num_model_values_minus1>,
    ValueHolder<intensity_interval_lower_bound>,
    ValueHolder<intensity_interval_upper_bound>,
    ValueHolder<comp_model_value>,
    ValueHolder<film_grain_characteristics_persistence_flag>
    {
    };

template <>
struct Syntax<film_grain_characteristics>
{
    template <class H> static void go(film_grain_characteristics fun, H &h);
};


template <class H>
void Syntax<film_grain_characteristics>::go(film_grain_characteristics fun, H &h)
{
    h(film_grain_characteristics_cancel_flag(), u(1));
    if (!h[film_grain_characteristics_cancel_flag()])
    {
        h(film_grain_model_id(), u(2));
        h(separate_colour_description_present_flag(), u(1));
        if (h[separate_colour_description_present_flag()])
        {
            h(film_grain_bit_depth_luma_minus8(), u(3));
            h(film_grain_bit_depth_chroma_minus8(), u(3));
            h(film_grain_full_range_flag(), u(1));
            h(film_grain_colour_primaries(), u(8));
            h(film_grain_transfer_characteristics(), u(8));
            h(film_grain_matrix_coeffs(), u(8));
        }
        h(blending_mode_id(), u(2));
        h(log2_scale_factor(), u(4));

        for (int c = 0; c < 3; c++)
            h(comp_model_present_flag(c), u(1));

        for (int c = 0; c < 3; c++)
            if (h[comp_model_present_flag(c)])
            {
                h(num_intensity_intervals_minus1(c), u(8));
                h(num_model_values_minus1(c), u(3));
                for (int i = 0; i <= h[num_intensity_intervals_minus1(c)]; i++)
                {
                    h(intensity_interval_lower_bound(c, i), u(8));
                    h(intensity_interval_upper_bound(c, i), u(8));
                    for (int j = 0; j <= h[num_model_values_minus1(c)]; j++)
                    {
                        h(comp_model_value(c, i, j), se(v));
                    }
                }
            }

        h(film_grain_characteristics_persistence_flag(), u(1));
    }
}


template <class H> void Read<film_grain_characteristics>::go(film_grain_characteristics f, H &h)
{
    FilmGrainCharacteristics filmGrainCharacteristics;
    auto h3 = h.extend(&filmGrainCharacteristics);

    Syntax<film_grain_characteristics>::go(f, h3);
}
