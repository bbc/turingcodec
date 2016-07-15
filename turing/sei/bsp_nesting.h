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
struct Syntax<bsp_nesting>
{
    template <class H> static void go(bsp_nesting const &fun, H &h)
    {
        h(sei_ols_idx(), ue(v));
        h(sei_partitioning_scheme_idx(), ue(v));
        h(bsp_idx(), ue(v));
        while (!h[byte_aligned()])
            h(bsp_nesting_zero_bit() /* equal to 0 */, u(1));
        h(num_seis_in_bsp_minus1(), ue(v));
        for (int i = 0; i <= h[num_seis_in_bsp_minus1()]; i++)
            h(sei_message());
    }
};



template <class H> void Read<bsp_nesting>::go(bsp_nesting const &f, H &h)
{
    BspNesting *bspNesting = h;

    if (bspNesting->nested)
    {
        h(Violation("F.14.3.2.7", "bsp_nesting() within bsp_nesting()")); // CondCheck F.14.3.2.7-?
    }
    else
    {
        bspNesting->nested = true;
        Syntax<bsp_nesting>::go(f, h);
        bspNesting->nested = false;
    }
}


#ifdef EXPLICIT_INSTANTIATION
    EXPLICIT_INSTANTIATION(bsp_nesting)
#endif
