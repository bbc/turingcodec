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

struct pan_scan_rect_id { };
struct pan_scan_rect_cancel_flag { };
struct pan_scan_cnt_minus1 { };
DEFINE_VALUE_ARRAY_1(pan_scan_rect_left_offset, i, 3);
DEFINE_VALUE_ARRAY_1(pan_scan_rect_right_offset, i, 3);
DEFINE_VALUE_ARRAY_1(pan_scan_rect_top_offset, i, 3);
DEFINE_VALUE_ARRAY_1(pan_scan_rect_bottom_offset, i, 3);
struct pan_scan_rect_persistence_flag { };

template <>
struct Syntax<pan_scan_rect>
{
    template <class H> static void go(pan_scan_rect fun, H &h)
    {
        h(pan_scan_rect_id(), ue(v));
        h(pan_scan_rect_cancel_flag(), u(1));
        if (!h[pan_scan_rect_cancel_flag()])
        {
            h(pan_scan_cnt_minus1(), ue(v));
            for (int i = 0; i <= h[pan_scan_cnt_minus1()]; i++)
            {
                h(pan_scan_rect_left_offset(i), se(v));
                h(pan_scan_rect_right_offset(i), se(v));
                h(pan_scan_rect_top_offset(i), se(v));
                h(pan_scan_rect_bottom_offset(i), se(v));
            }
            h(pan_scan_rect_persistence_flag(), u(1));
        }
    }
};

struct PanScanRect :
    ValueHolder<pan_scan_rect_id>,
    ValueHolder<pan_scan_rect_cancel_flag>,
    ValueHolder<pan_scan_cnt_minus1>,
    ValueHolder<pan_scan_rect_left_offset>,
    ValueHolder<pan_scan_rect_right_offset>,
    ValueHolder<pan_scan_rect_top_offset>,
    ValueHolder<pan_scan_rect_bottom_offset>,
    ValueHolder<pan_scan_rect_persistence_flag>
    {
    };

template <class H> void Read<pan_scan_rect>::go(pan_scan_rect f, H &h)
{
    PanScanRect panScanRect;
    auto h3 = h.extend(&panScanRect);

    Syntax<pan_scan_rect>::go(f, h3);
}


