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

/* Functions for decoding HEVC residual: inverse transform and add the result to predictor */


#ifndef INCLUDED_havoc_residual_decode_h
#define INCLUDED_havoc_residual_decode_h

#include "havoc.h"
#include <cassert>


using havoc_inverse_transform = void(int16_t dst[], int16_t const coeffs[], int bitDepth);


struct havoc_table_inverse_transform
{
    havoc_inverse_transform *sine;
    havoc_inverse_transform *cosine[4];
};


static inline havoc_inverse_transform** havoc_get_inverse_transform(havoc_table_inverse_transform *table, int trType, int log2TrafoSize)
{
    if (trType)
    {
        assert(log2TrafoSize == 2);
        return &table->sine;
    }
    else
    {
        return &table->cosine[log2TrafoSize - 2];
    }
}


void havoc_populate_inverse_transform(havoc_table_inverse_transform *table, havoc_code code, int encoder);


template <typename Sample>
using havoc_inverse_transform_add = void(Sample *dst, intptr_t stride_dst, Sample const* pred, intptr_t stride_pred, int16_t const coeffs[], int bitDepth);


template <typename Sample>
struct havoc_table_inverse_transform_add
{
    havoc_inverse_transform_add<Sample> *sine;
    havoc_inverse_transform_add<Sample> *cosine[4];
};


template <typename Sample>
static inline havoc_inverse_transform_add<Sample>** havoc_get_inverse_transform_add(havoc_table_inverse_transform_add<Sample> *table, int trType, int log2TrafoSize)
{
    if (trType)
    {
        assert(log2TrafoSize == 2);
        return &table->sine;
    }
    else
    {
        return &table->cosine[log2TrafoSize - 2];
    }
}


template <typename Sample>
void havoc_populate_inverse_transform_add(havoc_table_inverse_transform_add<Sample> *table, havoc_code code, int encoder);


template <typename Sample>
void havoc_test_inverse_transform_add(int *error_count, havoc_instruction_set mask);


static int havoc_clip(int x, int bit_depth)
{
    if (x < 0) return 0;
    int const max = (int)(1 << bit_depth) - 1;
    if (x > max) return max;
    return x;
}


template <typename Sample>
void havoc_add_residual(int n, Sample* dst, intptr_t stride_dst, Sample const* pred, intptr_t stride_pred, int16_t *residual, int bitDepth)
{
    for (int y = 0; y < n; ++y)
    {
        for (int x = 0; x < n; ++x)
        {
            dst[x + y * stride_dst] = havoc_clip(pred[x + y * stride_pred] + residual[x + y * n], bitDepth);
        }
    }
}


// Review: this is an encode function in a file called "residual_decode.h"
typedef void havoc_transform(int16_t *coeffs, const int16_t *src, intptr_t src_stride);


template <int bitDepth>
struct havoc_table_transform
{
    havoc_transform *dst;
    havoc_transform *dct[4];
};


template <int bitDepth>
static havoc_transform** havoc_get_transform(havoc_table_transform<bitDepth> *table, int trType, int log2TrafoSize)
{
    if (trType)
    {
        assert(log2TrafoSize == 2);
        return &table->dst;
    }
    else
    {
        return &table->dct[log2TrafoSize - 2];
    }
}

template <int bitDepth>
void havoc_populate_transform(havoc_table_transform<bitDepth> *table, havoc_code code);


void havoc_test_transform(int *error_count, havoc_instruction_set mask);


#endif
