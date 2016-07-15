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



struct lnp_sei_active_vps_id { };
DEFINE_VALUE_ARRAY_1(layer_not_present_flag, i, 16)


struct LayersNotPreset :
    ValueHolder<lnp_sei_active_vps_id>,
    ValueHolder<layer_not_present_flag>,
    Active<Vps>
    {
    };


template <> struct Syntax<layers_not_present>
{
    template <class H> static void go(layers_not_present fun, H &h)
    {
        h(lnp_sei_active_vps_id(), u(4));
        for (int i = 0; i <= h[MaxLayersMinus1()]; i++)
        {
            h(layer_not_present_flag(i), u(1));
        }
    }
};


#ifdef EXPLICIT_INSTANTIATION
    EXPLICIT_INSTANTIATION(layers_not_present)
#endif
