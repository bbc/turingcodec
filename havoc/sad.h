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


#ifndef INCLUDED_sad_h
#define INCLUDED_sad_h

#include "havoc.h"


#define X_HEVC_PU_SIZES \
        X(64, 64) \
        X(64, 48) \
        X(64, 32) \
        X(64, 16) \
        X(48, 64) \
        X(32, 64) \
        X(32, 32) \
        X(32, 24) \
        X(32, 16) \
        X(32, 8) \
        X(24, 32) \
        X(16, 64) \
        X(16, 32) \
        X(16, 16) \
        X(16, 12) \
        X(16, 8) \
        X(16, 4) \
        X(12, 16) \
        X(8, 32) \
        X(8, 16) \
        X(8, 8) \
        X(8, 4) \
        X(4, 8) \


//typedef std::function<int(const uint8_t * /*src*/, intptr_t /*stride_src*/, const uint8_t * /*ref*/, intptr_t /*stride_ref*/, uint32_t /*rect*/)> havoc_sad;

/* Rectangular SAD (Sum of Absolute Differences) with single reference */
template <typename Sample>
using havoc_sad = int(const Sample *src, intptr_t stride_src, const Sample *ref, intptr_t stride_ref, uint32_t rect);




template <typename Sample>
struct havoc_table_sad
{
#define X(w, h) \
        havoc_sad<Sample> *sad ## w ## x ## h;

    X_HEVC_PU_SIZES
#undef X

    havoc_sad<Sample> *sadGeneric;
}
;

template <typename Sample>
static havoc_sad<Sample>** havoc_get_sad(havoc_table_sad<Sample> *table, int width, int height)
{
    switch (HAVOC_RECT(width, height))
    {
#define X(w, h) \
        case HAVOC_RECT(w, h): return &table->sad ## w ## x ## h;

        X_HEVC_PU_SIZES
#undef X
        default:;
    }
    return &table->sadGeneric;
}


template <typename Sample>
void havoc_populate_sad(havoc_table_sad<Sample> *table, havoc_code code);

havoc_test_function havoc_test_sad;


/* Rectangular SAD (Sum of Absolute Differences) with multiple references */
template <typename Sample>
using havoc_sad_multiref  = void(const Sample *src, intptr_t stride_src, const Sample *ref[], intptr_t stride_ref, int sad[], uint32_t rect);


template <typename Sample>
struct havoc_table_sad_multiref
{
    havoc_sad_multiref<Sample> *lookup[16][16];
    havoc_sad_multiref<Sample> *sadGeneric_4;
};

template <typename Sample>
havoc_sad_multiref<Sample>** havoc_get_sad_multiref(havoc_table_sad_multiref<Sample> *table, int ways, int width, int height)
{
    if (ways != 4) return 0;
    return &table->lookup[(width >> 2) - 1][(height >> 2) - 1];
}

template <typename Sample>
void havoc_populate_sad_multiref(havoc_table_sad_multiref<Sample> *table, havoc_code code);

havoc_test_function havoc_test_sad_multiref;

#endif
