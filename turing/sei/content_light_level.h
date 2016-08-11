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
struct max_content_light_level { };
struct max_pic_average_light_level { };


struct ContentLightLevel :
    ValueHolder<max_content_light_level>,
    ValueHolder<max_pic_average_light_level>
    {
    };


template <> struct Syntax<content_light_level>
{
    template <class H> static void go(content_light_level fun, H &h);
};


template <class H>
void Syntax<content_light_level>::go(content_light_level fun, H &h)
{
    h(max_content_light_level(), u(16));
    h(max_pic_average_light_level(), u(16));
}


template <class H> void Read<content_light_level>::go(content_light_level f, H &h)
{
    ContentLightLevel s;
    auto h3 = h.extend(&s);

    Syntax<content_light_level>::go(f, h3);
}


