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


struct segmented_rect_frame_packing_arrangement_cancel_flag { };
struct segmented_rect_content_interpretation_type { };
struct segmented_rect_frame_packing_arrangement_persistence_flag { };


struct SegmentedRectFramePackingArrangement :
    ValueHolder<segmented_rect_frame_packing_arrangement_cancel_flag>,
    ValueHolder<segmented_rect_content_interpretation_type>,
    ValueHolder<segmented_rect_frame_packing_arrangement_persistence_flag>
    {
    };


template <> struct Syntax<segmented_rect_frame_packing_arrangement>
{
    template <class H> static void go(segmented_rect_frame_packing_arrangement fun, H &h);
};

template <class H>
void Syntax<segmented_rect_frame_packing_arrangement>::go(segmented_rect_frame_packing_arrangement fun, H &h)
{
    h(segmented_rect_frame_packing_arrangement_cancel_flag(), u(1));
    if (!h[segmented_rect_frame_packing_arrangement_cancel_flag()])
    {
        h(segmented_rect_content_interpretation_type(), u(2));
        h(segmented_rect_frame_packing_arrangement_persistence_flag(), u(1));
    }
}


template <class H> void Read<segmented_rect_frame_packing_arrangement>::go(segmented_rect_frame_packing_arrangement f, H &h)
{
    SegmentedRectFramePackingArrangement s;
    auto h2 = h.extend(&s);
    Syntax<segmented_rect_frame_packing_arrangement>::go(f, h2);
}
