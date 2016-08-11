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
#include "../Global.h"
#include "../StateParameterSets.h"


struct NalHrdBpPresentFlag {};
DEFINE_DERIVED_LONG(NalHrdBpPresentFlag)
    {
        Hrd *hrd = getHrd(h);
        if (!hrd)
            return 0;
        return (*hrd)[nal_hrd_parameters_present_flag()];
    }
};

struct VclHrdBpPresentFlag {};
DEFINE_DERIVED_LONG(VclHrdBpPresentFlag)
    {
        Hrd *hrd = getHrd(h);
        if (!hrd) 
            return 0;
        return (*hrd)[vcl_hrd_parameters_present_flag()];
    }
};


template <>
struct Syntax<buffering_period>
{
    template <class H> static void go(buffering_period fun, H &h)
    {
        h(bp_seq_parameter_set_id(), ue(v));
        if (!h[sub_pic_hrd_params_present_flag()])
            h(irap_cpb_params_present_flag(), u(1));
        if (h[irap_cpb_params_present_flag()])
        {
            h(cpb_delay_offset(), uv());
            h(dpb_delay_offset(), uv());
        }
        h(concatenation_flag(), u(1));
        h(au_cpb_removal_delay_delta_minus1(), uv());
        if (h[NalHrdBpPresentFlag()])
        {
            for (int i = 0; i <= h[CpbCnt()]; i++)
            {
                h(nal_initial_cpb_removal_delay(i), uv());
                h(nal_initial_cpb_removal_offset(i), uv());
                if (h[sub_pic_hrd_params_present_flag()] || h[irap_cpb_params_present_flag()])
                {
                    h(nal_initial_alt_cpb_removal_delay(i), uv());
                    h(nal_initial_alt_cpb_removal_offset(i), uv());
                }
            }
        }
        if (h[VclHrdBpPresentFlag()])
        {
            for (int i = 0; i <= h[CpbCnt()]; i++)
            {
                h(vcl_initial_cpb_removal_delay(i), uv());
                h(vcl_initial_cpb_removal_offset(i), uv());
                if (h[sub_pic_hrd_params_present_flag()] || h[irap_cpb_params_present_flag()])
                {
                    h(vcl_initial_alt_cpb_removal_delay(i), uv());
                    h(vcl_initial_alt_cpb_removal_offset(i), uv());
                }
            }
        }
    }
};


NUMBER_OF_BITS_MINUS1(cpb_delay_offset, au_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(dpb_delay_offset, dpb_output_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(au_cpb_removal_delay_delta_minus1, au_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(nal_initial_cpb_removal_delay, initial_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(nal_initial_cpb_removal_offset, initial_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(nal_initial_alt_cpb_removal_delay, initial_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(nal_initial_alt_cpb_removal_offset, initial_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(vcl_initial_cpb_removal_delay, initial_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(vcl_initial_cpb_removal_offset, initial_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(vcl_initial_alt_cpb_removal_delay, initial_cpb_removal_delay_length_minus1);
NUMBER_OF_BITS_MINUS1(vcl_initial_alt_cpb_removal_offset, initial_cpb_removal_delay_length_minus1);



template <class H> void Read<buffering_period>::go(buffering_period f, H &h)
{
    BufferingPeriod bufferingPeriod;
    auto h3 = h.extend(&bufferingPeriod);
    Hrd hrd;
    auto h2 = h3.extend(&hrd);

    Syntax<buffering_period>::go(f, h2);
}


template <> struct Read<Element<bp_seq_parameter_set_id, ue>>
{
    template <class H> static void go(Element<bp_seq_parameter_set_id, ue> f, H &h)
    {
        ReadUe<bp_seq_parameter_set_id>::go(f, h);

        auto found = h[Table<Sps>()].find(h[f.v]);

        if (found == h[Table<Sps>()].end())
        {
            h(Violation("D.3.2", "bp_seq_parameter_set_id does not indicate a valid SPS ID"));
            throw Abort();
        }
        else
        {
            h[Active<Sps>()] = found->second;
            auto &hrd = h[Active<Sps>()]->hrdArray.hrd;
            h.state = &h[Active<Sps>()]->hrdArray.hrd.front();
        }
    }
};
