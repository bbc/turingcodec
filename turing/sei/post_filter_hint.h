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




struct filter_hint_size_y { };
struct filter_hint_size_x { };
struct filter_hint_type { };
DEFINE_VALUE_ARRAY_3(filter_hint_value, cIdx, 3, cy, 16, cx, 16);


template <>
struct Syntax<post_filter_hint>
{
    template <class H> static void go(post_filter_hint fun, H &h);
};


template <class H>
void Syntax<post_filter_hint>::go(post_filter_hint fun, H &h)
{
    h(filter_hint_size_y(), ue(v));
    h(filter_hint_size_x(), ue(v));
    h(filter_hint_type(), u(2));
    for (int cIdx = 0; cIdx < (h[chroma_format_idc()] == 0 ? 1 : 3); cIdx++)
        for (int cy = 0; cy < h[filter_hint_size_y()]; cy++)
            for (int cx = 0; cx < h[filter_hint_size_x()]; cx++)
                h(filter_hint_value(cIdx, cy, cx), se(v));
}


struct PostFilterHint :
    ValueHolder<filter_hint_size_y>,
    ValueHolder<filter_hint_size_x>,
    ValueHolder<filter_hint_type>,
    ValueHolder<filter_hint_value>
    {
    };


template <class H> void Read<post_filter_hint>::go(post_filter_hint f, H &h)
{
    PostFilterHint postFilterHint;
    auto h3 = h.extend(&postFilterHint);
    Syntax<post_filter_hint>::go(f, h3);
}
