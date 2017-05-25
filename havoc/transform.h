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


#ifndef INCLUDED_residual_decode_h
#define INCLUDED_residual_decode_h

#include "havoc.h"
#include <cassert>


namespace havoc {

using inverse_transform = void(int16_t dst[], int16_t const coeffs[], int bitDepth);


struct table_inverse_transform
{
    inverse_transform *sine;
    inverse_transform *cosine[4];
};


static inline inverse_transform** get_inverse_transform(table_inverse_transform *table, int trType, int log2TrafoSize)
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


void populate_inverse_transform(table_inverse_transform *table, havoc_code code, int encoder);


template <typename Sample>
using inverse_transform_add = void(Sample *dst, intptr_t stride_dst, Sample const* pred, intptr_t stride_pred, int16_t const coeffs[], int bitDepth);


template <typename Sample>
struct table_inverse_transform_add
{
    inverse_transform_add<Sample> *sine;
    inverse_transform_add<Sample> *cosine[4];
};


template <typename Sample>
static inline inverse_transform_add<Sample>** get_inverse_transform_add(table_inverse_transform_add<Sample> *table, int trType, int log2TrafoSize)
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
void populate_inverse_transform_add(table_inverse_transform_add<Sample> *table, havoc_code code, int encoder);


template <typename Sample>
void test_inverse_transform_add(int *error_count, havoc_instruction_set mask);


static int clip(int x, int bit_depth)
{
    if (x < 0) return 0;
    int const max = (int)(1 << bit_depth) - 1;
    if (x > max) return max;
    return x;
}


template <typename Sample>
void add_residual(int n, Sample* dst, intptr_t stride_dst, Sample const* pred, intptr_t stride_pred, int16_t *residual, int bitDepth)
{
    for (int y = 0; y < n; ++y)
    {
        for (int x = 0; x < n; ++x)
        {
            dst[x + y * stride_dst] = clip(pred[x + y * stride_pred] + residual[x + y * n], bitDepth);
        }
    }
}


typedef void Transform(int16_t *coeffs, const int16_t *src, intptr_t src_stride);


template <int bitDepth>
struct table_transform
{
    Transform *dst;
    Transform *dct[4];
};


template <int bitDepth>
static Transform** get_transform(table_transform<bitDepth> *table, int trType, int log2TrafoSize)
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
void populate_transform(table_transform<bitDepth> *table, havoc_code code);


void test_transform(int *error_count, havoc_instruction_set mask);

}

#endif
