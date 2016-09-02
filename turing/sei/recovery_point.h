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


struct recovery_poc_cnt { };
struct exact_match_flag { };
struct broken_link_flag { };

struct RecoveryPoint :
    ValueHolder<recovery_poc_cnt>,
    ValueHolder<exact_match_flag>,
    ValueHolder<broken_link_flag>
    {
    };


template <>
struct Syntax<recovery_point>
{
    template <class H> static void go(recovery_point fun, H &h)
    {
        h(recovery_poc_cnt(), se(v));
        h(exact_match_flag(), u(1));
        h(broken_link_flag(), u(1));
    }
};


template <class H> void Read<recovery_point>::go(recovery_point f, H &h)
{
    RecoveryPoint recoveryPoint;
    auto h3 = h.extend(&recoveryPoint);
    Syntax<recovery_point>::go(f, h3);
}
