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
struct Syntax<scalable_nesting>
{
    template <class H> static void go(scalable_nesting fun, H &h)
    {
        h(bitstream_subset_flag(), u(1));
        h(nesting_op_flag(), u(1));
        if (h[nesting_op_flag()])
        {
            h(default_op_flag(), u(1));
            h(nesting_num_ops_minus1(), ue(v));
            for (int i = h[default_op_flag()]; i <= h[nesting_num_ops_minus1()]; i++)
            {
                h(nesting_max_temporal_id_plus1(i), u(3));
                h(nesting_op_idx(i), ue(v));
            }
        }
        else
        {
            h(all_layers_flag(), u(1));
            if (!h[all_layers_flag()])
            {
                h(nesting_no_op_max_temporal_id_plus1(), u(3));
                h(nesting_num_layers_minus1(), ue(v));
                for (int i = 0; i <= h[nesting_num_layers_minus1()]; i++)
                    h(nesting_layer_id(i), u(6));
            }
        }

        while (!h[byte_aligned()])
            h(nesting_zero_bit() /* equal to 0 */, u(1));

        do
            h(sei_message());
        while (h[more_rbsp_data()]);
    }
};


template <class H> void Read<scalable_nesting>::go(scalable_nesting f, H &h)
{
    ScalableNesting *scalableNesting = h;

    if (scalableNesting->nested)
    {
        h(Violation("D.3.1", "scalable_nesting() within scalable_nesting()")); // CondCheck D.3.1-M
    }
    else
    {
        scalableNesting->nested = true;
        Syntax<scalable_nesting>::go(f, h);
        scalableNesting->nested = false;
    }
}
