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


struct ffinfo_pic_struct { };
struct ffinfo_source_scan_type { };
struct ffinfo_duplicate_flag { };


template <>
struct Syntax<frame_field_info>
{
    template <class H> static void go(frame_field_info const &fun, H &h)
    {
        h(ffinfo_pic_struct(), u(4));
        h(ffinfo_source_scan_type(), u(2));
        h(ffinfo_duplicate_flag(), u(1));
    }
};


struct FieldFrameInfo :
    ValueHolder<ffinfo_pic_struct>,
    ValueHolder<ffinfo_source_scan_type>,
    ValueHolder<ffinfo_duplicate_flag>
    {
    };


template <class H> void Read<frame_field_info>::go(frame_field_info f, H &h)
{
    FieldFrameInfo fieldFrameInfo;
    auto hh = h.extend(&fieldFrameInfo);
    Syntax<frame_field_info>::go(f, hh);
}
