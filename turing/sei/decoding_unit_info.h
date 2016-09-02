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

//

struct decoding_unit_idx { };
struct du_spt_cpb_removal_delay_increment { };
struct dpb_output_du_delay_present_flag { };
struct pic_spt_dpb_output_du_delay { };

struct DecodingUnitInfo :
    ValueHolder<decoding_unit_idx>,
    ValueHolder<du_spt_cpb_removal_delay_increment>,
    ValueHolder<dpb_output_du_delay_present_flag>,
    ValueHolder<pic_spt_dpb_output_du_delay>
    {
    };

template <>
struct Syntax<decoding_unit_info>
{
    template <class H> static void go(decoding_unit_info fun, H &h)
    {
        h(decoding_unit_idx(), ue(v));
        if (!h[sub_pic_cpb_params_in_pic_timing_sei_flag()])
            h(du_spt_cpb_removal_delay_increment(), uv());
        h(dpb_output_du_delay_present_flag(), u(1));
        if (h[dpb_output_du_delay_present_flag()])
            h(pic_spt_dpb_output_du_delay(), uv());
    }
};

NUMBER_OF_BITS_MINUS1(du_spt_cpb_removal_delay_increment, du_cpb_removal_delay_increment_length_minus1)
NUMBER_OF_BITS_MINUS1(pic_spt_dpb_output_du_delay, dpb_output_delay_du_length_minus1)


template <class H> void Read<decoding_unit_info>::go(decoding_unit_info f, H &h)
{
    Hrd *hrd = getHrd(h);

    if (hrd)
    {
        auto h2 = h.extend(hrd);

        DecodingUnitInfo decodingUnitInfo;
        auto h3 = h2.extend(&decodingUnitInfo);

        Syntax<decoding_unit_info>::go(f, h3);
    }
    else
    {
        // Seek to end of SEI payload so as not to trigger a further error
        seek(h, bitLen(h[::Stream()]));
    }
}


