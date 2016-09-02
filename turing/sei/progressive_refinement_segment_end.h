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
#include "progressive_refinement_segment_start.h"


template <>
struct Syntax<progressive_refinement_segment_end>
{
    template <class H> static void go(progressive_refinement_segment_end fun, H &h);
};

template <class H>
void Syntax<progressive_refinement_segment_end>::go(progressive_refinement_segment_end fun, H &h)
{
    h(progressive_refinement_id(), ue(v));
}

struct ProgressiveRefinementSegmentEnd :
    ValueHolder<progressive_refinement_id>
    {
    };



template <class H> void Read<progressive_refinement_segment_end>::go(progressive_refinement_segment_end f, H &h)
{
    ProgressiveRefinementSegmentEnd progressiveRefinementSegmentEnd;
    auto h3 = h.extend(&progressiveRefinementSegmentEnd);

    Syntax<progressive_refinement_segment_end>::go(f, h3);
}

