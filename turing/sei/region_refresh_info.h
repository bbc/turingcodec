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

struct refreshed_region_flag { };

struct RegionRefreshInfo :
    ValueHolder<refreshed_region_flag>
    {
    };

template <>
struct Syntax<region_refresh_info>
{
    template <class H> static void go(region_refresh_info fun, H &h)
    {
        h(refreshed_region_flag(), u(1));
    }
};


template <class H> void Read<region_refresh_info>::go(region_refresh_info f, H &h)
{
    RegionRefreshInfo regionRefreshInfo;
    auto h3 = h.extend(&regionRefreshInfo);

    Syntax<region_refresh_info>::go(f, h3);
}
