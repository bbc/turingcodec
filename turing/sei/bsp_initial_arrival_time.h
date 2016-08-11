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


struct nal_initial_arrival_delay { };
struct vcl_initial_arrival_delay { };

DEFINED_DERIVED(BspSchedCnt, h[num_bsp_schedules_minus1(h, i, t)] + 1)

struct BspInitialArrivalTime :
    ValueHolder<nal_initial_arrival_delay>,
    ValueHolder<vcl_initial_arrival_delay>
    {
    };


template <>
struct Syntax<bsp_initial_arrival_time>
{
    template <class H> static void go(bsp_initial_arrival_time fun, H &h)
    {
        auto const psIdx = h[sei_partitioning_scheme_idx()];
        if (nalInitialArrivalDelayPresent)
            for (i = 0; i < BspSchedCnt[sei_ols_idx][psIdx][MaxTemporalId[0]]; i++)
                h(nal_initial_arrival_delay(i), uv());
        if (vclInitialArrivalDelayPresent)
            for (i = 0; i < BspSchedCnt[sei_ols_idx][psIdx][MaxTemporalId[0]]; i++)
                h(vcl_initial_arrival_delay(i), uv());
    }
};


template <class H> void Read<bsp_initial_arrival_time>::go(bsp_initial_arrival_time f, H &h)
{
    RecoveryPoint recoveryPoint;
    auto h3 = h.extend(&recoveryPoint);
    Syntax<recovery_point>::go(f, h3);
}


#ifdef EXPLICIT_INSTANTIATION
    EXPLICIT_INSTANTIATION(recovery_point)
#endif
