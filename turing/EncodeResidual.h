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

#ifndef INCLUDED_EncodeResidual_h
#define INCLUDED_EncodeResidual_h

#include "Global.h"
#include "StateFunctionTables.h"

template <class> struct Write;
template <class> struct EstimateRate;


struct EncodeResidual
{
    template <class H>
    static inline typename std::enable_if<std::is_same<typename H::Tag, Write<void>>::value>::type encode(H &h)
    {
        encode2(h);
    }

    template <class H>
    static inline typename std::enable_if<!std::is_same<typename H::Tag, Write<void>>::value>::type encode(H &h)
    {
        residual_coding *rc = h;
        auto hh = h.template change<EstimateRate<void>>();
        encode2(hh);
    }

    template <class H>
    static inline void encode2(H &h)
    {
        StateFunctionTables *stateFunctionTables = h;
        bool const havePopcnt
            = (stateFunctionTables->instruction_set_support & HAVOC_LZCNT)
            && (stateFunctionTables->instruction_set_support & HAVOC_POPCNT);
            
        residual_coding *rc = h;
        if (rc->log2TrafoSize == 2)
        {
            if (havePopcnt)
                inner<true, true>(h);
            else
                inner<false, true>(h);
        }
        else
        {
            if (havePopcnt)
                inner<true, false>(h);
            else
                inner<false, false>(h);
        }
    }

    template <bool havePopcnt, bool is4x4, class H> 
    static void inner(H &hParent);
};

#endif