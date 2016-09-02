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


struct tone_map_id { };
struct tone_map_cancel_flag { };
struct tone_map_persistence_flag { };
struct coded_data_bit_depth { };
struct target_bit_depth { };
struct tone_map_model_id { };
struct min_value { };
struct max_value { };
struct sigmoid_midpoint { };
struct sigmoid_width { };
DEFINE_VALUE_ARRAY_1(start_of_coded_interval, i, 65536);
struct num_pivots { };
DEFINE_VALUE_ARRAY_1(coded_pivot_value, i, 65536);
DEFINE_VALUE_ARRAY_1(target_pivot_value, i, 65536);
struct camera_iso_speed_idc { };
struct camera_iso_speed_value { };
struct exposure_index_idc { };
struct exposure_index_value { };
struct exposure_compensation_value_sign_flag { };
struct exposure_compensation_value_numerator { };
struct exposure_compensation_value_denom_idc { };
struct ref_screen_luminance_white { };
struct extended_range_white_level { };
struct nominal_black_level_code_value { };
struct nominal_white_level_code_value { };
struct extended_white_level_code_value { };

struct ToneMappingInfo :
    ValueHolder<tone_map_id>,
    ValueHolder<tone_map_cancel_flag>,
    ValueHolder<tone_map_persistence_flag>,
    ValueHolder<coded_data_bit_depth>,
    ValueHolder<target_bit_depth>,
    ValueHolder<tone_map_model_id>,
    ValueHolder<min_value>,
    ValueHolder<max_value>,
    ValueHolder<sigmoid_midpoint>,
    ValueHolder<sigmoid_width>,
    ValueHolder<start_of_coded_interval>,
    ValueHolder<num_pivots>,
    ValueHolder<coded_pivot_value>,
    ValueHolder<target_pivot_value>,
    ValueHolder<camera_iso_speed_idc>,
    ValueHolder<camera_iso_speed_value>,
    ValueHolder<exposure_index_idc>,
    ValueHolder<exposure_index_value>,
    ValueHolder<exposure_compensation_value_sign_flag>,
    ValueHolder<exposure_compensation_value_numerator>,
    ValueHolder<exposure_compensation_value_denom_idc>,
    ValueHolder<ref_screen_luminance_white>,
    ValueHolder<extended_range_white_level>,
    ValueHolder<nominal_black_level_code_value>,
    ValueHolder<nominal_white_level_code_value>,
    ValueHolder<extended_white_level_code_value>
    {
    };

template <>
struct Syntax<tone_mapping_info>
{
    template <class H> static void go(tone_mapping_info fun, H &h);
};


template <class H>
void Syntax<tone_mapping_info>::go(tone_mapping_info fun, H &h)
{
    h(tone_map_id(), ue(v));
    h(tone_map_cancel_flag(), u(1));
    if (!h[tone_map_cancel_flag()])
    {
        h(tone_map_persistence_flag(), u(1));
        h(coded_data_bit_depth(), u(8));
        h(target_bit_depth(), u(8));
        h(tone_map_model_id(), ue(v));
        if (h[tone_map_model_id()] == 0)
        {
            h(min_value(), u(32));
            h(max_value(), u(32));
        }
        else if (h[tone_map_model_id()] == 1)
        {
            h(sigmoid_midpoint(), u(32));
            h(sigmoid_width(), u(32));
        }
        else if (h[tone_map_model_id()] == 2)
        {
            for (int i = 0; i < (1 << h[target_bit_depth()]); i++)
                h(start_of_coded_interval(i), u(((h[coded_data_bit_depth()] + 7) >> 3) << 3));
        }
        else if (h[tone_map_model_id()] == 3)
        {
            h(num_pivots(), u(16));
            for (int i = 0; i < h[num_pivots()]; i++)
            {
                h(coded_pivot_value(i), u(((h[coded_data_bit_depth()] + 7) >> 3) << 3));
                h(target_pivot_value(i), u(((h[target_bit_depth()] + 7) >> 3) << 3));
            }
        }
        else if (h[tone_map_model_id()] == 4)
        {
#define EXTENDED_ISO 255
            h(camera_iso_speed_idc(), u(8));
            if (h[camera_iso_speed_idc()] == EXTENDED_ISO)
                h(camera_iso_speed_value(), u(32));
            h(exposure_index_idc(), u(8));
            if (h[exposure_index_idc()] == EXTENDED_ISO)
                h(exposure_index_value(), u(32));
            h(exposure_compensation_value_sign_flag(), u(1));
            h(exposure_compensation_value_numerator(), u(16));
            h(exposure_compensation_value_denom_idc(), u(16));
            h(ref_screen_luminance_white(), u(32));
            h(extended_range_white_level(), u(32));
            h(nominal_black_level_code_value(), u(16));
            h(nominal_white_level_code_value(), u(16));
            h(extended_white_level_code_value(), u(16));
        }
    }
}


template <class H> void Read<tone_mapping_info>::go(tone_mapping_info f, H &h)
{
    ToneMappingInfo toneMappingInfo;
    auto h3 = h.extend(&toneMappingInfo);

    Syntax<tone_mapping_info>::go(f, h3);
}
