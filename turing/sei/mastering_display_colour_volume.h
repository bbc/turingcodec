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
#include "../StateParameterSets.h"


DEFINE_VALUE_ARRAY_1(display_primaries_x, c, 3)
DEFINE_VALUE_ARRAY_1(display_primaries_y, c, 3)
struct white_point_x { };
struct white_point_y { };
struct max_display_mastering_luminance { };
struct min_display_mastering_luminance { };

struct MasteringDisplayColourVolume :
    AccessOperators<MasteringDisplayColourVolume>,
    ValueHolder<display_primaries_x>,
    ValueHolder<display_primaries_y>,
    ValueHolder<white_point_x>,
    ValueHolder<white_point_y>,
    ValueHolder<max_display_mastering_luminance>,
    ValueHolder<min_display_mastering_luminance>
    {
    };


template <> struct Syntax<mastering_display_colour_volume>
{
    template <class H> static void go(mastering_display_colour_volume fun, H &h)
    {
        for (int c = 0; c < 3; c++)
        {
            h(display_primaries_x(c), u(16));
            h(display_primaries_y(c), u(16));
        }
        h(white_point_x(), u(16));
        h(white_point_y(), u(16));
        h(max_display_mastering_luminance(), u(32));
        h(min_display_mastering_luminance(), u(32));
    }
};


template <class H> void  Read<mastering_display_colour_volume>::go(mastering_display_colour_volume f, H &h)
{
    StateParameterSets *stateParameterSets = h;
    stateParameterSets->masteringDisplayColourVolume.reset(new MasteringDisplayColourVolume());
    auto h3 = h.extend(stateParameterSets->masteringDisplayColourVolume.get());
    Syntax<mastering_display_colour_volume>::go(f, h3);
}
