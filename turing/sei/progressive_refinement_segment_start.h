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


struct progressive_refinement_id { };
struct pic_order_cnt_delta { };


template <>
struct Syntax<progressive_refinement_segment_start>
{
    template <class H> static void go(progressive_refinement_segment_start fun, H &h)
    {
        h(progressive_refinement_id(), ue(v));
        h(pic_order_cnt_delta(), ue(v));
    }
};


struct ProgressiveRefinementSegmentStart :
    ValueHolder<progressive_refinement_id>,
    ValueHolder<pic_order_cnt_delta>
    {
    };


template <class H> void Read<progressive_refinement_segment_start>::go(progressive_refinement_segment_start f, H &h)
{
    ProgressiveRefinementSegmentStart progressiveRefinementSegmentStart;
    auto h3 = h.extend(&progressiveRefinementSegmentStart);

    Syntax<progressive_refinement_segment_start>::go(f, h3);
}
