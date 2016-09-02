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

struct pic_struct { };
struct source_scan_type { };
struct duplicate_flag { };
DEFINE_DERIVED(CpbDpbDelaysPresentFlag, h[nal_hrd_parameters_present_flag()] || h[vcl_hrd_parameters_present_flag()]);
struct au_cpb_removal_delay_minus1 { };
struct pic_dpb_output_delay { };
struct pic_dpb_output_du_delay { };
struct num_decoding_units_minus1 { };
struct du_common_cpb_removal_delay_flag { };
struct du_common_cpb_removal_delay_increment_minus1 { };
DEFINE_STRUCT_ARITY_1(num_nalus_in_du_minus1, i);
DEFINE_STRUCT_ARITY_1(du_cpb_removal_delay_increment_minus1, i);


template <>
struct Syntax<pic_timing>
{
    template <class H> static void go(pic_timing fun, H &h)
    {
        if (h[frame_field_info_present_flag()])
        {
            h(pic_struct(), u(4));
            h(source_scan_type(), u(2));
            h(duplicate_flag(), u(1));
        }
        if (h[CpbDpbDelaysPresentFlag()])
        {
            h(au_cpb_removal_delay_minus1(), uv());
            h(pic_dpb_output_delay(), uv());
            if (h[sub_pic_hrd_params_present_flag()])
                h(pic_dpb_output_du_delay(), uv());
            if (h[sub_pic_hrd_params_present_flag()] && h[sub_pic_cpb_params_in_pic_timing_sei_flag()])
            {
                h(num_decoding_units_minus1(), ue(v));
                h(du_common_cpb_removal_delay_flag(), u(1));
                if (h[du_common_cpb_removal_delay_flag()])
                    h(du_common_cpb_removal_delay_increment_minus1(), uv());
                for (int i = 0; i <= h[num_decoding_units_minus1()]; i++)
                {
                    h(num_nalus_in_du_minus1(i), ue(v));
                    if (!h[du_common_cpb_removal_delay_flag()] && i < h[num_decoding_units_minus1()])
                        h(du_cpb_removal_delay_increment_minus1(i), uv());
                }
            }
        }
    }
};

NUMBER_OF_BITS_MINUS1(au_cpb_removal_delay_minus1, au_cpb_removal_delay_length_minus1)
NUMBER_OF_BITS_MINUS1(pic_dpb_output_delay, dpb_output_delay_length_minus1)
NUMBER_OF_BITS_MINUS1(pic_dpb_output_du_delay, dpb_output_delay_du_length_minus1)
NUMBER_OF_BITS_MINUS1(du_common_cpb_removal_delay_increment_minus1, du_cpb_removal_delay_increment_length_minus1)
NUMBER_OF_BITS_MINUS1(du_cpb_removal_delay_increment_minus1, du_cpb_removal_delay_increment_length_minus1)


struct PicTiming :
    ValueHolder<pic_struct>,
    ValueHolder<source_scan_type>,
    ValueHolder<duplicate_flag>,
    ValueHolder<au_cpb_removal_delay_minus1>,
    ValueHolder<pic_dpb_output_delay>,
    ValueHolder<pic_dpb_output_du_delay>,
    ValueHolder<num_decoding_units_minus1>,
    ValueHolder<du_common_cpb_removal_delay_flag>,
    ValueHolder<du_common_cpb_removal_delay_increment_minus1>,
    ValueHolder<num_nalus_in_du_minus1>,
    ValueHolder<du_cpb_removal_delay_increment_minus1>
    {
    };

template <class H> void Read<pic_timing>::go(pic_timing f, H &h)
{
    Hrd *hrd = getHrd(h);

    if (hrd)
    {
        auto h2 = h.extend(hrd);
        PicTiming picTiming;
        auto h3 = h2.extend(&picTiming);

        Syntax<pic_timing>::go(f, h3);
    }
    else
    {
        // Seek to end of SEI payload so as not to trigger a further error
        seek(h, bitLen(h[::Stream()]));
    }
}

