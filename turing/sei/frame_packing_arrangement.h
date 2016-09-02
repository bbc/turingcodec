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

struct frame_packing_arrangement_id { };
struct frame_packing_arrangement_cancel_flag { };
struct frame_packing_arrangement_type { };
struct quincunx_sampling_flag { };
struct content_interpretation_type { };
struct spatial_flipping_flag { };
struct frame0_flipped_flag { };
struct field_views_flag { };
struct current_frame_is_frame0_flag { };
struct frame0_self_contained_flag { };
struct frame1_self_contained_flag { };
struct frame0_grid_position_x { };
struct frame0_grid_position_y { };
struct frame1_grid_position_x { };
struct frame1_grid_position_y { };
struct frame_packing_arrangement_reserved_byte { };
struct frame_packing_arrangement_persistence_flag { };
struct upsampled_aspect_ratio_flag { };


struct FramePackingArrangement :
    ValueHolder<frame_packing_arrangement_id>,
    ValueHolder<frame_packing_arrangement_cancel_flag>,
    ValueHolder<frame_packing_arrangement_type>,
    ValueHolder<quincunx_sampling_flag>,
    ValueHolder<content_interpretation_type>,
    ValueHolder<spatial_flipping_flag>,
    ValueHolder<frame0_flipped_flag>,
    ValueHolder<field_views_flag>,
    ValueHolder<current_frame_is_frame0_flag>,
    ValueHolder<frame0_self_contained_flag>,
    ValueHolder<frame1_self_contained_flag>,
    ValueHolder<frame0_grid_position_x>,
    ValueHolder<frame0_grid_position_y>,
    ValueHolder<frame1_grid_position_x>,
    ValueHolder<frame1_grid_position_y>,
    ValueHolder<frame_packing_arrangement_reserved_byte>,
    ValueHolder<frame_packing_arrangement_persistence_flag>,
    ValueHolder<upsampled_aspect_ratio_flag>
    {
    };

template <>
struct Syntax<frame_packing_arrangement>
{
    template <class H> static void go(frame_packing_arrangement fun, H &h);
};


template <class H>
void Syntax<frame_packing_arrangement>::go(frame_packing_arrangement fun, H &h)
{
    h(frame_packing_arrangement_id(), ue(v));
    h(frame_packing_arrangement_cancel_flag(), u(1));
    if (!h[frame_packing_arrangement_cancel_flag()])
    {
        h(frame_packing_arrangement_type(), u(7));
        h(quincunx_sampling_flag(), u(1));
        h(content_interpretation_type(), u(6));
        h(spatial_flipping_flag(), u(1));
        h(frame0_flipped_flag(), u(1));
        h(field_views_flag(), u(1));
        h(current_frame_is_frame0_flag(), u(1));
        h(frame0_self_contained_flag(), u(1));
        h(frame1_self_contained_flag(), u(1));
        if (!h[quincunx_sampling_flag()] && h[frame_packing_arrangement_type()] != 5)
        {
            h(frame0_grid_position_x(), u(4));
            h(frame0_grid_position_y(), u(4));
            h(frame1_grid_position_x(), u(4));
            h(frame1_grid_position_y(), u(4));
        }
        h(frame_packing_arrangement_reserved_byte(), u(8));
        h(frame_packing_arrangement_persistence_flag(), u(1));
    }
    h(upsampled_aspect_ratio_flag(), u(1));
}


template <class H> void Read<frame_packing_arrangement>::go(frame_packing_arrangement f, H &h)
{
    FramePackingArrangement framePackingArrangement;
    auto h3 = h.extend(&framePackingArrangement);

    Syntax<frame_packing_arrangement>::go(f, h3);
}
