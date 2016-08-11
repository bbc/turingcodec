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


struct display_orientation_cancel_flag { };
struct hor_flip { };
struct ver_flip { };
struct anticlockwise_rotation { };
struct display_orientation_persistence_flag { };

struct DisplayOrientation :
    ValueHolder<display_orientation_cancel_flag>,
    ValueHolder<hor_flip>,
    ValueHolder<ver_flip>,
    ValueHolder<anticlockwise_rotation>,
    ValueHolder<display_orientation_persistence_flag>
    {
    };


template <>
struct Syntax<display_orientation>
{
    template <class H> static void go(display_orientation fun, H &h);
};


template <class H>
void Syntax<display_orientation>::go(display_orientation fun, H &h)
{
    h(display_orientation_cancel_flag(), u(1));
    if (!h[display_orientation_cancel_flag()])
    {
        h(hor_flip(), u(1));
        h(ver_flip(), u(1));
        h(anticlockwise_rotation(), u(16));
        h(display_orientation_persistence_flag(), u(1));
    }
}


template <class H> void Read<display_orientation>::go(display_orientation f, H &h)
{
    DisplayOrientation displayOrientation;
    auto h3 = h.extend(&displayOrientation);

    Syntax<display_orientation>::go(f, h3);
}


