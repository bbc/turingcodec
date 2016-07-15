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

#include "sao.h"
#include <cstdint>
#include <cstddef>
#include <cassert>


static int Clip3(int min, int max, int value)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}


template <typename Sample>
void sao_filter_band(Sample *__restrict dst, intptr_t dst_stride, const Sample *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict offset_table, int bitDepth)
{
    int y;
    for (y = 0; y < h; ++y)
    {
        int x;
        for (x = 0; x < w; ++x)
        {
            const int bandShift = bitDepth - 5;
            const int index = (int)src[x + y * src_stride] >> bandShift;

            int value = Clip3(0, (1 << bitDepth) - 1, src[x + y * src_stride] + offset_table[index]);

            dst[x + y * dst_stride] = value;
        }
    }
}

static int Sign(int x)
{
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

template <typename Sample>
void sao_filter_edge(Sample *__restrict dst, intptr_t dst_stride, const Sample *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict SaoOffsetVal, int eoClass, int bitDepth)
{
    const int hLookup[4] = { -1, 0, -1, 1 };
    const int vLookup[4] = { 0, -1, -1, -1 };

    const int hPos[2] = { hLookup[eoClass], -hLookup[eoClass] };
    const int vPos[2] = { vLookup[eoClass], -vLookup[eoClass] };

    const Sample *neighbours[2] =
    {
            &src[hPos[0] + vPos[0] * src_stride],
            &src[hPos[1] + vPos[1] * src_stride]
    };

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            int edgeIdx = 2
                    + Sign((int)src[x + y * src_stride] - (int)neighbours[0][x + y * src_stride])
                    + Sign((int)src[x + y * src_stride] - (int)neighbours[1][x + y * src_stride]);

            if (edgeIdx == 0 || edgeIdx == 1 || edgeIdx == 2)
            {
                edgeIdx = (edgeIdx == 2) ? 0 : (edgeIdx + 1);
            }

            int value = Clip3(0, (1 << bitDepth) - 1, src[x + y * src_stride] + SaoOffsetVal[edgeIdx]);

            dst[x + y * dst_stride] = value;
        }
    }
}

template void sao_filter_band<uint16_t>(uint16_t *__restrict dst, intptr_t dst_stride, const uint16_t *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict offset_table, int bitDepth);
template void sao_filter_band<uint8_t>(uint8_t *__restrict dst, intptr_t dst_stride, const uint8_t *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict offset_table, int bitDepth);
template void sao_filter_edge<uint16_t>(uint16_t *__restrict dst, intptr_t dst_stride, const uint16_t *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict SaoOffsetVal, int eoClass, int bitDepth);
template void sao_filter_edge<uint8_t>(uint8_t *__restrict dst, intptr_t dst_stride, const uint8_t *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict SaoOffsetVal, int eoClass, int bitDepth);
