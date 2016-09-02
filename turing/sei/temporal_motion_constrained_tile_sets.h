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


struct mc_all_tiles_exact_sample_value_match_flag { };
struct each_tile_one_tile_set_flag { };
struct limited_tile_set_display_flag { };
struct num_sets_in_message_minus1 { };
DEFINE_VALUE_ARRAY_1(mcts_id, i, 256)
DEFINE_VALUE_ARRAY_1(display_tile_set_flag, i, 256)
DEFINE_VALUE_ARRAY_1(num_tile_rects_in_set_minus1, i, 256)
DEFINE_VALUE_ARRAY_2(top_left_tile_index, i, 256, j, 256)
DEFINE_VALUE_ARRAY_2(bottom_right_tile_index, i, 256, j, 256)
DEFINE_VALUE_ARRAY_1(mc_exact_sample_value_match_flag, i, 256)
DEFINE_VALUE_ARRAY_1(mcts_tier_level_idc_present_flag, i, 256)
DEFINE_VALUE_ARRAY_1(mcts_tier_flag, i, 256)
DEFINE_VALUE_ARRAY_1(mcts_level_idc, i, 256)
struct mcts_max_tier_level_idc_present_flag { };
struct mcts_max_tier_flag { };
struct mcts_max_level_idc { };



template <> struct Syntax<temporal_motion_constrained_tile_sets>
{
    template <class H> static void go(temporal_motion_constrained_tile_sets fun, H &h);
};

template <class H>
void Syntax<temporal_motion_constrained_tile_sets>::go(temporal_motion_constrained_tile_sets fun, H &h)
{
    h(mc_all_tiles_exact_sample_value_match_flag(), u(1));
    h(each_tile_one_tile_set_flag(), u(1));
    if (!h[each_tile_one_tile_set_flag()])
    {
        h(limited_tile_set_display_flag(), u(1));
        h(num_sets_in_message_minus1(), ue(v));
        for (int i = 0; i <= h[num_sets_in_message_minus1()]; i++)
        {
            h(mcts_id(i), ue(v));
            if (h[limited_tile_set_display_flag()])
                h(display_tile_set_flag(i), u(1));
            h(num_tile_rects_in_set_minus1(i), ue(v));
            for (int j = 0; j <= h[num_tile_rects_in_set_minus1(i)]; j++)
            {
                h(top_left_tile_index(i, j), ue(v));
                h(bottom_right_tile_index(i, j), ue(v));
            }
            if (!h[mc_all_tiles_exact_sample_value_match_flag()])
                h(mc_exact_sample_value_match_flag(i), u(1));
            h(mcts_tier_level_idc_present_flag(i), u(1));
            if (h[mcts_tier_level_idc_present_flag(i)])
            {
                h(mcts_tier_flag(i), u(1));
                h(mcts_level_idc(i), u(8));
            }
        }
    }
    else
    {
        h(mcts_max_tier_level_idc_present_flag(), u(1));
        if (h[mcts_max_tier_level_idc_present_flag()])
        {
            h(mcts_max_tier_flag(), u(1));
            h(mcts_max_level_idc(), u(8));
        }
    }
}


struct TemporalMotionConstrainedTileSets :
    ValueHolder<mc_all_tiles_exact_sample_value_match_flag>,
    ValueHolder<each_tile_one_tile_set_flag>,
    ValueHolder<limited_tile_set_display_flag>,
    ValueHolder<num_sets_in_message_minus1>,
    ValueHolder<mcts_id>,
    ValueHolder<display_tile_set_flag>,
    ValueHolder<num_tile_rects_in_set_minus1>,
    ValueHolder<top_left_tile_index>,
    ValueHolder<bottom_right_tile_index>,
    ValueHolder<mc_exact_sample_value_match_flag>,
    ValueHolder<mcts_tier_level_idc_present_flag>,
    ValueHolder<mcts_tier_flag>,
    ValueHolder<mcts_level_idc>,
    ValueHolder<mcts_max_tier_level_idc_present_flag>,
    ValueHolder<mcts_max_tier_flag>,
    ValueHolder<mcts_max_level_idc>
    {
    };

template <class H> void Read<temporal_motion_constrained_tile_sets>::go(temporal_motion_constrained_tile_sets f, H &h)
{
    TemporalMotionConstrainedTileSets s;
    auto h2 = h.extend(&s);
    Syntax<temporal_motion_constrained_tile_sets>::go(f, h2);
}
