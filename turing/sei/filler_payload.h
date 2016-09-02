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


template <>
struct Syntax<filler_payload>
{
    template <class H> static void go(filler_payload fun, H &h)
    {
        for (int k = 0; k < fun.payloadSize; k++)
            h(ff_byte() /* equal to 0xFF */, f(8));
    }
};


template <class H> void Read<filler_payload>::go(filler_payload f, H &h)
{
    Syntax<filler_payload>::go(f, h);
}
