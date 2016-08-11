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


struct knee_function_id { };
struct knee_function_cancel_flag { };
struct knee_function_persistence_flag { };
struct input_d_range { };
struct input_disp_luminance { };
struct output_d_range { };
struct output_disp_luminance { };
struct num_knee_points_minus1 { };
DEFINE_STRUCT_ARITY_1(input_knee_point, i);
DEFINE_STRUCT_ARITY_1(output_knee_point, i);



template <>
struct Syntax<knee_function_info>
{
    template <class H> static void go(knee_function_info const &fun, H &h)
    {
        h(knee_function_id(), ue(v));
        h(knee_function_cancel_flag(), u(1));
        if (!h[knee_function_cancel_flag()])
        {
            h(knee_function_persistence_flag(), u(1));
            h(input_d_range(), u(32));
            h(input_disp_luminance(), u(32));
            h(output_d_range(), u(32));
            h(output_disp_luminance(), u(32));
            h(num_knee_points_minus1(), ue(v));
            for (int i = 0; i <= h[num_knee_points_minus1()]; i++)
            {
                h(input_knee_point(i), u(10));
                h(output_knee_point(i), u(10));
            }
        }
    }
};


struct KneeFunctionInfo :
    ValueHolder<knee_function_id>,
    ValueHolder<knee_function_cancel_flag>,
    ValueHolder<knee_function_persistence_flag>,
    ValueHolder<input_d_range>,
    ValueHolder<input_disp_luminance>,
    ValueHolder<output_d_range>,
    ValueHolder<output_disp_luminance>,
    ValueHolder<num_knee_points_minus1>,
    ValueHolder<input_knee_point>,
    ValueHolder<output_knee_point>
    {
    };


template <class H> void Read<knee_function_info>::go(knee_function_info f, H &h)
{
    KneeFunctionInfo kneeFunctionInfo;
    auto hh = h.extend(&kneeFunctionInfo);
    Syntax<knee_function_info>::go(f, hh);
}
