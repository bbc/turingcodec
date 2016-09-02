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


struct il_all_tiles_exact_sample_value_match_flag { };
struct il_one_tile_per_tile_set_flag { };
struct il_num_sets_in_message_minus1 { };
struct skipped_tile_set_present_flag { };
struct all_tiles_ilc_idc { };
DEFINE_STRUCT_ARITY_1(ilcts_id, i);
DEFINE_STRUCT_ARITY_1(il_num_tile_rects_in_set_minus1, i);
DEFINE_STRUCT_ARITY_1(ilc_idc, i);
DEFINE_STRUCT_ARITY_1(il_exact_sample_value_match_flag, i);
DEFINE_STRUCT_ARITY_2(il_top_left_tile_index, i, j);
DEFINE_STRUCT_ARITY_2(il_bottom_right_tile_index, i, j);


template <>
struct Syntax<inter_layer_constrained_tile_sets>
{
    template <class H> static void go(inter_layer_constrained_tile_sets const &fun, H &h)
    {
        h(il_all_tiles_exact_sample_value_match_flag(), u(1));
        h(il_one_tile_per_tile_set_flag(), u(1));
        if (!h[il_one_tile_per_tile_set_flag()])
        {
            h(il_num_sets_in_message_minus1(), ue(v));
            if (h[il_num_sets_in_message_minus1()])
                h(skipped_tile_set_present_flag(), u(1));
            int const numSignificantSets = h[il_num_sets_in_message_minus1()] - h[skipped_tile_set_present_flag()] + 1;
            for (int i = 0; i < numSignificantSets; i++)
            {
                h(ilcts_id(i), ue(v));
                h(il_num_tile_rects_in_set_minus1(i), ue(v));
                for (int j = 0; j <= h[il_num_tile_rects_in_set_minus1(i)]; j++)
                {
                    h(il_top_left_tile_index(i, j), ue(v));
                    h(il_bottom_right_tile_index(i, j), ue(v));
                }
                h(ilc_idc(i), u(2));
                if (!h[il_all_tiles_exact_sample_value_match_flag()])
                    h(il_exact_sample_value_match_flag(i), u(1));
            }
        }
        else
            h(all_tiles_ilc_idc(), u(2));
    }
};


struct InterLayerConstrainedTileSets :
    ValueHolder<il_all_tiles_exact_sample_value_match_flag>,
    ValueHolder<il_one_tile_per_tile_set_flag>,
    ValueHolder<il_num_sets_in_message_minus1>,
    ValueHolder<skipped_tile_set_present_flag>,
    ValueHolder<all_tiles_ilc_idc>,
    ValueHolder<ilcts_id>,
    ValueHolder<il_num_tile_rects_in_set_minus1>,
    ValueHolder<ilc_idc>,
    ValueHolder<il_exact_sample_value_match_flag>,
    ValueHolder<il_top_left_tile_index>,
    ValueHolder<il_bottom_right_tile_index>
    {
    };


template <class H> void Read<inter_layer_constrained_tile_sets>::go(inter_layer_constrained_tile_sets f, H &h)
{
    InterLayerConstrainedTileSets interLayerConstrainedTileSets;
    auto hh = h.extend(&interLayerConstrainedTileSets);
    Syntax<inter_layer_constrained_tile_sets>::go(f, hh);
}
