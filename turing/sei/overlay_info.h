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


// st(v) : null-terminated string encoded as UTF - 8 characters as specified in ISO / IEC 10646. The parsing process is
// specified as follows : st(v) begins at a byte - aligned position in the bitstream and reads and returns a series of
// bytes from the bitstream, beginning at the current position and continuing up to but not including the next bytealigned
// byte that is equal to 0x00, and advances the bitstream pointer by(stringLength + 1) * 8 bit positions,
// where stringLength is equal to the number of bytes returned.
struct st
{
    st(int v) {}
};


struct overlay_info_cancel_flag { };
struct overlay_content_aux_id_minus128 { };
struct overlay_label_aux_id_minus128 { };
struct overlay_alpha_aux_id_minus128 { };
struct overlay_element_label_value_length_minus8 { }; // also referred to as "overlay_element_label_min_max_length_minus8" in the standard - typo?
struct num_overlays_minus1 { };
DEFINE_VALUE_ARRAY_1(overlay_idx, i, 16);
DEFINE_VALUE_ARRAY_1(language_overlay_present_flag, i, 16);
DEFINE_VALUE_ARRAY_1(overlay_content_layer_id, i, 16);
DEFINE_VALUE_ARRAY_1(overlay_label_present_flag, i, 16);
DEFINE_VALUE_ARRAY_1(overlay_label_layer_id, i, 16);
DEFINE_VALUE_ARRAY_1(overlay_alpha_present_flag, i, 16);
DEFINE_VALUE_ARRAY_1(overlay_alpha_layer_id, i, 16);
DEFINE_VALUE_ARRAY_1(num_overlay_elements_minus1, i, 16);
DEFINE_VALUE_ARRAY_2(overlay_element_label_min, i, 16, j, 256);
DEFINE_VALUE_ARRAY_2(overlay_element_label_max, i, 16, j, 256);
struct overlay_zero_bit { };
template <> struct Fixed<overlay_zero_bit> { static const int value = 0; };
template <> struct ValueType<struct overlay_language> { typedef std::string Type; };
DEFINE_VALUE_ARRAY_1(overlay_language, i, 16);
template <> struct ValueType<struct overlay_name> { typedef std::string Type; };
DEFINE_VALUE_ARRAY_1(overlay_name, i, 16);
template <> struct ValueType<struct overlay_element_name> { typedef std::string Type; };
DEFINE_VALUE_ARRAY_2(overlay_element_name, i, 16, j, 256);
struct overlay_info_persistence_flag { };


template <>
struct Syntax<overlay_info>
{
    template <class H> static void go(overlay_info const &fun, H &h)
    {
        h(overlay_info_cancel_flag(), u(1));
        if (!h[overlay_info_cancel_flag()])
        {
            h(overlay_content_aux_id_minus128(), ue(v));
            h(overlay_label_aux_id_minus128(), ue(v));
            h(overlay_alpha_aux_id_minus128(), ue(v));
            h(overlay_element_label_value_length_minus8(), ue(v));
            h(num_overlays_minus1(), ue(v));
            for (int i = 0; i <= h[num_overlays_minus1()]; i++)
            {
                h(overlay_idx(i), ue(v));
                h(language_overlay_present_flag(i), u(1));
                h(overlay_content_layer_id(i), u(6));
                h(overlay_label_present_flag(i), u(1));
                if (h[overlay_label_present_flag(i)])
                    h(overlay_label_layer_id(i), u(6));
                h(overlay_alpha_present_flag(i), u(1));
                if (h[overlay_alpha_present_flag(i)])
                    h(overlay_alpha_layer_id(i), u(6));
                if (h[overlay_label_present_flag(i)])
                {
                    h(num_overlay_elements_minus1(i), ue(v));
                    for (int j = 0; j <= h[num_overlay_elements_minus1(i)]; j++)
                    {
                        h(overlay_element_label_min(i, j), uv());
                        h(overlay_element_label_max(i, j), uv());
                    }
                }
            }
            while (!h[byte_aligned()])
                h(overlay_zero_bit() /* equal to 0 */, f(1));
            for (int i = 0; i <= h[num_overlays_minus1()]; i++)
            {
                if (h[language_overlay_present_flag(i)])
                    h(overlay_language(i), st(v));
                h(overlay_name(i), st(v));
                if (h[overlay_label_present_flag(i)])
                    for (int j = 0; j <= h[num_overlay_elements_minus1(i)]; j++)
                        h(overlay_element_name(i, j), st(v));
            }
            h(overlay_info_persistence_flag(), u(1));
        }
    }
};


NUMBER_OF_BITS_UV(overlay_element_label_min, h[overlay_element_label_value_length_minus8()] + 8)
NUMBER_OF_BITS_UV(overlay_element_label_max, h[overlay_element_label_value_length_minus8()] + 8)



struct OverlayInfo :
    ValueHolder<overlay_info_cancel_flag>,
    ValueHolder<overlay_content_aux_id_minus128>,
    ValueHolder<overlay_label_aux_id_minus128>,
    ValueHolder<overlay_alpha_aux_id_minus128>,
    ValueHolder<overlay_element_label_value_length_minus8>,
    ValueHolder<num_overlays_minus1>,
    ValueHolder<overlay_idx>,
    ValueHolder<language_overlay_present_flag>,
    ValueHolder<overlay_content_layer_id>,
    ValueHolder<overlay_label_present_flag>,
    ValueHolder<overlay_label_layer_id>,
    ValueHolder<overlay_alpha_present_flag>,
    ValueHolder<overlay_alpha_layer_id>,
    ValueHolder<num_overlay_elements_minus1>,
    ValueHolder<overlay_element_label_min>,
    ValueHolder<overlay_element_label_max>,
    ValueHolder<overlay_zero_bit>,
    ValueHolder<overlay_language>,
    ValueHolder<overlay_name>,
    ValueHolder<overlay_element_name>,
    ValueHolder<overlay_info_persistence_flag>
    {
    };


template <class H> void Read<overlay_info>::go(overlay_info f, H &h)
{
    OverlayInfo overlayInfo;
    auto hh = h.extend(&overlayInfo);
    Syntax<overlay_info>::go(f, hh);
}


template <>
struct Read<Element<overlay_element_label_min, u>>
{
    template <class H> static void go(Element<overlay_element_label_min, u> fun, H &h)
    {
        const int nBits = h[overlay_element_label_value_length_minus8()] + 8;
        h[fun.v] = read_bits<typename ValueType<overlay_element_label_min>::Type>(h, nBits);
    }
};

template <>
struct Read<Element<overlay_element_label_max, u>>
{
    template <class H> static void go(Element<overlay_element_label_max, u> fun, H &h)
    {
        const int nBits = h[overlay_element_label_value_length_minus8()] + 8;
        h[fun.v] = read_bits<typename ValueType<overlay_element_label_max>::Type>(h, nBits);
    }
};


template <class V>
struct Read<Element<V, st>>
{
    template <class H> static void go(Element<V, st> fun, H &h)
    {
        while (true)
        {
            char c = read_bytes<char>(h, 1);
            if (!c) 
                break;
            h[fun.v] += c;
        }
    }
};
