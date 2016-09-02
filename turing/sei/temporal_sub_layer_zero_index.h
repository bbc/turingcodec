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

//

struct temporal_sub_layer_zero_idx { };
struct irap_pic_id { };

struct TemporalSubLayerZeroIndex :
    ValueHolder<temporal_sub_layer_zero_idx>,
    ValueHolder<irap_pic_id>
    {
    };

template <>
struct Syntax<temporal_sub_layer_zero_index>
{
    template <class H> static void go(temporal_sub_layer_zero_index fun, H &h)
    {
        h(temporal_sub_layer_zero_idx(), u(8));
        h(irap_pic_id(), u(8));
    }
};


template <class H> void Read<temporal_sub_layer_zero_index>::go(temporal_sub_layer_zero_index f, H &h)
{
    TemporalSubLayerZeroIndex temporalSubLayerZeroIndex;
    auto h3 = h.extend(&temporalSubLayerZeroIndex);

    Syntax<temporal_sub_layer_zero_index>::go(f, h3);
}
