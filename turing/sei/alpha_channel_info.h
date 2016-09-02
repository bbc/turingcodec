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


struct alpha_channel_cancel_flag { };
struct alpha_channel_use_idc { };
struct alpha_channel_bit_depth_minus8 { };
struct alpha_transparent_value { };
struct alpha_opaque_value { };
struct alpha_channel_incr_flag { };
struct alpha_channel_clip_flag { };
struct alpha_channel_clip_type_flag { };


template <>
struct Syntax<alpha_channel_info>
{
    template <class H> static void go(alpha_channel_info const &fun, H &h)
    {
        h(alpha_channel_cancel_flag(), u(1));
        if (!h[alpha_channel_cancel_flag()])
        {
            h(alpha_channel_use_idc(), u(3));
            h(alpha_channel_bit_depth_minus8(), u(3));
            h(alpha_transparent_value(), uv());
            h(alpha_opaque_value(), uv());
            h(alpha_channel_incr_flag(), u(1));
            h(alpha_channel_clip_flag(), u(1));
            if (h[alpha_channel_clip_flag()])
                h(alpha_channel_clip_type_flag(), u(1));
        }
    }
};


NUMBER_OF_BITS_UV(alpha_transparent_value, h[alpha_channel_bit_depth_minus8()] + 9)
NUMBER_OF_BITS_UV(alpha_opaque_value, h[alpha_channel_bit_depth_minus8()] + 9)


struct AlphaChannelInfo :
    ValueHolder<alpha_channel_cancel_flag>,
    ValueHolder<alpha_channel_use_idc>,
    ValueHolder<alpha_channel_bit_depth_minus8>,
    ValueHolder<alpha_transparent_value>,
    ValueHolder<alpha_opaque_value>,
    ValueHolder<alpha_channel_incr_flag>,
    ValueHolder<alpha_channel_clip_flag>,
    ValueHolder<alpha_channel_clip_type_flag>
    {
    };


template <class H> void Read<alpha_channel_info>::go(alpha_channel_info f, H &h)
{
    AlphaChannelInfo alphaChannelInfo;
    auto hh = h.extend(&alphaChannelInfo);
    Syntax<alpha_channel_info>::go(f, hh);
}
