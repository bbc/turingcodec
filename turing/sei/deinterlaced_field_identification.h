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


struct deinterlaced_picture_source_parity_flag { };


template <>
struct Syntax<deinterlaced_field_identification>
{
    template <class H> static void go(deinterlaced_field_identification const &fun, H &h)
    {
        h(deinterlaced_picture_source_parity_flag(), u(1));
    }
};


struct DeinterlacedFieldIdentification :
    ValueHolder<deinterlaced_picture_source_parity_flag>
    {
    };


template <class H> void Read<deinterlaced_field_identification>::go(deinterlaced_field_identification f, H &h)
{
    DeinterlacedFieldIdentification deinterlacedFieldIdentification;
    auto hh = h.extend(&deinterlacedFieldIdentification);
    Syntax<deinterlaced_field_identification>::go(f, hh);
}


