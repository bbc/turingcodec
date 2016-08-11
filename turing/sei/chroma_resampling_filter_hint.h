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


struct ver_chroma_filter_idc { };
struct hor_chroma_filter_idc { };
struct ver_filtering_field_processing_flag { };
struct target_format_idc { };
struct num_vertical_filters { };
DEFINE_STRUCT_ARITY_1(ver_tap_length_minus1, i);
DEFINE_STRUCT_ARITY_2(ver_filter_coeff, i, j);
struct num_horizontal_filters { };
DEFINE_STRUCT_ARITY_1(hor_tap_length_minus1, i);
DEFINE_STRUCT_ARITY_2(hor_filter_coeff, i, j);


template <>
struct Syntax<chroma_resampling_filter_hint>
{
    template <class H> static void go(chroma_resampling_filter_hint const &fun, H &h)
    {
        h(ver_chroma_filter_idc(), u(8));
        h(hor_chroma_filter_idc(), u(8));
        h(ver_filtering_field_processing_flag(), u(1));
        if (h[ver_chroma_filter_idc()] == 1 || h[hor_chroma_filter_idc()] == 1)
        {
            h(target_format_idc(), ue(v));
            if (h[ver_chroma_filter_idc()] == 1)
            {
                h(num_vertical_filters(), ue(v));
                for (int i = 0; i < h[num_vertical_filters()]; i++)
                {
                    h(ver_tap_length_minus1(i), ue(v));
                    for (int j = 0; j <= h[ver_tap_length_minus1(i)]; j++)
                        h(ver_filter_coeff(i, j), se(v));
                }
            }
            if (h[hor_chroma_filter_idc()] == 1)
            {
                h(num_horizontal_filters(), ue(v));
                for (int i = 0; i < h[num_horizontal_filters()]; i++)
                {
                    h(hor_tap_length_minus1(i), ue(v));
                    for (int j = 0; j <= h[hor_tap_length_minus1(i)]; j++)
                        h(hor_filter_coeff(i, j), se(v));
                }
            }
        }
    }
};


struct ChromaResamplingFilterHint :
    ValueHolder<ver_chroma_filter_idc>,
    ValueHolder<hor_chroma_filter_idc>,
    ValueHolder<ver_filtering_field_processing_flag>,
    ValueHolder<target_format_idc>,
    ValueHolder<num_vertical_filters>,
    ValueHolder<ver_tap_length_minus1>,
    ValueHolder<ver_filter_coeff>,
    ValueHolder<num_horizontal_filters>,
    ValueHolder<hor_tap_length_minus1>,
    ValueHolder<hor_filter_coeff>
    {
    };


template <class H> void Read<chroma_resampling_filter_hint>::go(chroma_resampling_filter_hint f, H &h)
{
    ChromaResamplingFilterHint chromaResamplingFilterHint;
    auto hh = h.extend(&chromaResamplingFilterHint);
    Syntax<chroma_resampling_filter_hint>::go(f, hh);
}


