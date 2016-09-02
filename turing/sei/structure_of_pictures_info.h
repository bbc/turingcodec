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
#include "../HevcTypes.h"


struct sop_seq_parameter_set_id { };
struct num_entries_in_sop_minus1 { };
DEFINE_VALUE_ARRAY_1(sop_vcl_nut, i, 1024);
DEFINE_VALUE_ARRAY_1(sop_temporal_id, i, 1024);
DEFINE_VALUE_ARRAY_1(sop_short_term_rps_idx, i, 1024);
DEFINE_VALUE_ARRAY_1(sop_poc_delta, i, 1024);


template <> struct Syntax<structure_of_pictures_info>
{
    template <class H> static void go(structure_of_pictures_info fun, H &h)
    {
        h(sop_seq_parameter_set_id(), ue(v));
        h(num_entries_in_sop_minus1(), ue(v));
        for (int i = 0; i <= h[num_entries_in_sop_minus1()]; i++)
        {
            h(sop_vcl_nut(i), u(6));
            h(sop_temporal_id(i), u(3));
            if (h[sop_vcl_nut(i)] != IDR_W_RADL  &&  h[sop_vcl_nut(i)] != IDR_N_LP)
                h(sop_short_term_rps_idx(i), ue(v));
            if (i > 0)
                h(sop_poc_delta(i), se(v));
        }
    }
};


struct StructureOfPicturesInfo :
    ValueHolder<sop_seq_parameter_set_id>,
    ValueHolder<num_entries_in_sop_minus1>,
    ValueHolder<sop_vcl_nut>,
    ValueHolder<sop_temporal_id>,
    ValueHolder<sop_short_term_rps_idx>,
    ValueHolder<sop_poc_delta>
    {
    };


template <class H> void Read<structure_of_pictures_info>::go(structure_of_pictures_info f, H &h)
{
    StructureOfPicturesInfo structureOfPicturesInfo;
    auto hh = h.extend(&structureOfPicturesInfo);
    Syntax<structure_of_pictures_info>::go(f, hh);
}
