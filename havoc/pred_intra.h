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

#ifndef INCLUDED_havoc_prediction_intra_h
#define INCLUDED_havoc_prediction_intra_h

// HEVC intra prediction


#include "havoc.h"

namespace havoc { namespace intra { 
    
template <typename Sample>
using Function = void(Sample *dst, intptr_t dstStride, Sample const *neighbours, int predModeIntra);

template <typename Sample>
struct Table
{
    Function<Sample> *entries[3 * sizeof(Sample) - 2 /* bitDepth */][4 /* log2TrafoSize - 2 */][38 /* intraPredMode */];

    inline Function<Sample> *&lookup(int cIdx, int bitDepth, int log2TrafoSize, int predModeIntra)
    {
        auto const edge_flag = cIdx == 0 && log2TrafoSize < 5;
        if (edge_flag)
            if (predModeIntra == 1)
                predModeIntra = 35;
            else if (predModeIntra == 10)
                predModeIntra = 36;
            else if (predModeIntra == 26)
                predModeIntra = 37;

        auto const bd = sizeof(Sample) == 2 ? 10 - bitDepth : 0;
        return this->entries[bd][log2TrafoSize - 2][predModeIntra];
    }

    void populate(havoc_code code);
};

template <typename Sample>
havoc_test_function test;

}}

#endif
