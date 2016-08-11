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


template <> struct Syntax<time_code>
{
    template <class H> static void go(time_code fun, H &h);
};



struct num_clock_ts { };
DEFINE_VALUE_ARRAY_1(clock_timestamp_flag, i, 4)
DEFINE_VALUE_ARRAY_1(units_field_based_flag, i, 4)
DEFINE_VALUE_ARRAY_1(counting_type, i, 4)
DEFINE_VALUE_ARRAY_1(full_timestamp_flag, i, 4)
DEFINE_VALUE_ARRAY_1(discontinuity_flag, i, 4)
DEFINE_VALUE_ARRAY_1(cnt_dropped_flag, i, 4)
DEFINE_VALUE_ARRAY_1(n_frames, i, 4)
DEFINE_VALUE_ARRAY_1(seconds_value, i, 4)
DEFINE_VALUE_ARRAY_1(minutes_value, i, 4)
DEFINE_VALUE_ARRAY_1(hours_value, i, 4)
DEFINE_VALUE_ARRAY_1(seconds_flag, i, 4)
DEFINE_VALUE_ARRAY_1(minutes_flag, i, 4)
DEFINE_VALUE_ARRAY_1(hours_flag, i, 4)
DEFINE_VALUE_ARRAY_1(time_offset_length, i, 4)
DEFINE_VALUE_ARRAY_1(time_offset_value, i, 4)

struct TimeCode :
    ValueHolder<num_clock_ts>,
    ValueHolder<clock_timestamp_flag>,
    ValueHolder<units_field_based_flag>,
    ValueHolder<counting_type>,
    ValueHolder<full_timestamp_flag>,
    ValueHolder<discontinuity_flag>,
    ValueHolder<cnt_dropped_flag>,
    ValueHolder<n_frames>,
    ValueHolder<seconds_value>,
    ValueHolder<minutes_value>,
    ValueHolder<hours_value>,
    ValueHolder<seconds_flag>,
    ValueHolder<minutes_flag>,
    ValueHolder<hours_flag>,
    ValueHolder<time_offset_length>,
    ValueHolder<time_offset_value>
    {
    };


template <class H> void Read<time_code>::go(time_code  f, H &h)
{
    TimeCode s;
    auto h3 = h.extend(&s);
    Syntax<time_code>::go(f, h3);
}


template <class H>
void Syntax<time_code>::go(time_code fun, H &h)
{
    h(num_clock_ts(), u(2));
    for (int i = 0; i < h[num_clock_ts()]; i++)
    {
        h(clock_timestamp_flag(i), u(1));
        if (h[clock_timestamp_flag(i)])
        {
            h(units_field_based_flag(i), u(1));
            h(counting_type(i), u(5));
            h(full_timestamp_flag(i), u(1));
            h(discontinuity_flag(i), u(1));
            h(cnt_dropped_flag(i), u(1));
            h(n_frames(i), u(9));
            if (h[full_timestamp_flag(i)])
            {
                h(seconds_value(i) /* 0..59 */, u(6));
                h(minutes_value(i) /* 0..59 */, u(6));
                h(hours_value(i) /* 0..23 */, u(5));
            }
            else
            {
                h(seconds_flag(i), u(1));
                if (h[seconds_flag(i)])
                {
                    h(seconds_value(i) /* 0..59 */, u(6));
                    h(minutes_flag(i), u(1));
                    if (h[minutes_flag(i)])
                    {
                        h(minutes_value(i) /* 0..59 */, u(6));
                        h(hours_flag(i), u(1));
                        if (h[hours_flag(i)])
                            h(hours_value(i) /* 0..23 */, u(5));
                    }
                }
                h(time_offset_length(i), u(5));
                if (h[time_offset_length(i)] > 0)
                    h(time_offset_value(i), ::i(v));
            }
        }
    }
}
