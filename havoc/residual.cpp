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

#include "residual.h"
#include "havoc_test.h"
#include "Jit.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <array>
#include <map>
#include <vector>

template <typename Sample>
static void havoc_residual_c_ref(int16_t *pRes, intptr_t strideRes, const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB, int w, int h)
{
    //subtract prediction from input
    //review: SIMD optimisations, perhaps integrate into forward transform
    //for (int y = 0; y < nTbS; ++y)
    //    for (int x = 0; x < nTbS; ++x)
    //        resSamplesRaw(x, y) = sourceSamples(x, y) - predSamples(x, y);
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            pRes[x + y * strideRes] = pA[x + y * strideA] - pB[x + y * strideB];
        }
    }
}


#define ORDER(a, b, c, d) ((a << 6) | (b << 4) | (c << 2) | d)

template <typename Sample>
struct ResidualSSE :
    Jit::Function
    {
        ResidualSSE(Jit::Buffer *buffer, int width, int height) :
        Jit::Function(buffer, Jit::CountArguments<havoc_residual<Sample>>::value),
        width(width),
        height(height)
        {
            this->build();
        }

        int width, height;

        void assemble() override
        {
            auto &r0 = arg64(0); // Sample   *pRes
            auto &r1 = arg64(1); // intptr_t strideRes
            auto &r2 = arg64(2); // Sample   *pA
            auto &r3 = arg64(3); // intptr_t strideA
            auto &r4 = arg64(4); // Sample   *pB
            auto &r5 = arg64(5); // intptr_t strideB
            auto &r6 = arg64(6); // int w
            auto &r7 = arg64(7); // int h

            int const widthBytes = width * sizeof(Sample);

            auto &m0 = regXmm(0);
            auto &m1 = regXmm(1);
            auto &m2 = regXmm(2);
            auto &m3 = regXmm(3);
            auto &m5 = regXmm(5);

            pxor(m0, m0);
            pxor(m1, m1);
            pxor(m2, m2);
            pxor(m3, m3);
            pxor(m5, m5);

            if (sizeof(Sample) == 1)
            {
#define INSERT_LINE_4(last)                         \
                {                                   \
                     movd(m0, ptr[r2]);             \
                     movd(m1, ptr[r4]);             \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     psubw(m0, m1);                 \
                     movq(ptr[r0], m0);             \
                     if (!last)                     \
                     {                              \
                         lea(r0, ptr[r0 + r1 * 2]); \
                         lea(r2, ptr[r2 + r3]);     \
                         lea(r4, ptr[r4 + r5]);     \
                     }                              \
                }
#define INSERT_LINE_8(last)                         \
                {                                   \
                     movq(m0, ptr[r2]);             \
                     movq(m1, ptr[r4]);             \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     psubw(m0, m1);                 \
                     movdqa(ptr[r0], m0);           \
                     if (!last)                     \
                     {                              \
                         lea(r0, ptr[r0 + r1 * 2]); \
                         lea(r2, ptr[r2 + r3]);     \
                         lea(r4, ptr[r4 + r5]);     \
                     }                              \
                }
#define INSERT_LINE_16(last)                        \
                {                                   \
                     movdqa(m0, ptr[r2]);           \
                     movdqa(m1, ptr[r4]);           \
                     movdqa(m2, m0);                \
                     movdqa(m3, m1);                \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     punpckhbw(m2, m5);             \
                     punpckhbw(m3, m5);             \
                     psubw(m0, m1);                 \
                     psubw(m2, m3);                 \
                     movdqa(ptr[r0],        m0);    \
                     movdqa(ptr[r0 + 0x10], m2);    \
                     if (!last)                     \
                     {                              \
                         lea(r0, ptr[r0 + r1 * 2]); \
                         lea(r2, ptr[r2 + r3]);     \
                         lea(r4, ptr[r4 + r5]);     \
                     }                              \
                }
#define INSERT_LINE_32(last)                        \
                {                                   \
                     movdqa(m0, ptr[r2]);           \
                     movdqa(m1, ptr[r4]);           \
                     movdqa(m2, m0);                \
                     movdqa(m3, m1);                \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     punpckhbw(m2, m5);             \
                     punpckhbw(m3, m5);             \
                     psubw(m0, m1);                 \
                     psubw(m2, m3);                 \
                     movdqa(ptr[r0],        m0);    \
                     movdqa(ptr[r0 + 0x10], m2);    \
                     movdqa(m0, ptr[r2 + 0x10]);    \
                     movdqa(m1, ptr[r4 + 0x10]);    \
                     movdqa(m2, m0);                \
                     movdqa(m3, m1);                \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     punpckhbw(m2, m5);             \
                     punpckhbw(m3, m5);             \
                     psubw(m0, m1);                 \
                     psubw(m2, m3);                 \
                     movdqa(ptr[r0 + 0x20], m0);    \
                     movdqa(ptr[r0 + 0x30], m2);    \
                     if (!last)                     \
                     {                              \
                         lea(r0, ptr[r0 + r1 * 2]); \
                         lea(r2, ptr[r2 + r3]);     \
                         lea(r4, ptr[r4 + r5]);     \
                     }                              \
                }
#define INSERT_LINE_64(last)                        \
                {                                   \
                     movdqa(m0, ptr[r2]);           \
                     movdqa(m1, ptr[r4]);           \
                     movdqa(m2, m0);                \
                     movdqa(m3, m1);                \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     punpckhbw(m2, m5);             \
                     punpckhbw(m3, m5);             \
                     psubw(m0, m1);                 \
                     psubw(m2, m3);                 \
                     movdqa(ptr[r0],        m0);    \
                     movdqa(ptr[r0 + 0x10], m2);    \
                     movdqa(m0, ptr[r2 + 0x10]);    \
                     movdqa(m1, ptr[r4 + 0x10]);    \
                     movdqa(m2, m0);                \
                     movdqa(m3, m1);                \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     punpckhbw(m2, m5);             \
                     punpckhbw(m3, m5);             \
                     psubw(m0, m1);                 \
                     psubw(m2, m3);                 \
                     movdqa(ptr[r0 + 0x20], m0);    \
                     movdqa(ptr[r0 + 0x30], m2);    \
                     movdqa(m0, ptr[r2 + 020]);     \
                     movdqa(m1, ptr[r4 + 020]);     \
                     movdqa(m2, m0);                \
                     movdqa(m3, m1);                \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     punpckhbw(m2, m5);             \
                     punpckhbw(m3, m5);             \
                     psubw(m0, m1);                 \
                     psubw(m2, m3);                 \
                     movdqa(ptr[r0 + 0x40], m0);    \
                     movdqa(ptr[r0 + 0x50], m2);    \
                     movdqa(m0, ptr[r2 + 0x30]);    \
                     movdqa(m1, ptr[r4 + 0x30]);    \
                     movdqa(m2, m0);                \
                     movdqa(m3, m1);                \
                     punpcklbw(m0, m5);             \
                     punpcklbw(m1, m5);             \
                     punpckhbw(m2, m5);             \
                     punpckhbw(m3, m5);             \
                     psubw(m0, m1);                 \
                     psubw(m2, m3);                 \
                     movdqa(ptr[r0 + 0x60], m0);    \
                     movdqa(ptr[r0 + 0x70], m2);    \
                     if (!last)                     \
                     {                              \
                         lea(r0, ptr[r0 + r1 * 2]); \
                         lea(r2, ptr[r2 + r3]);     \
                         lea(r4, ptr[r4 + r5]);     \
                     }                              \
                }
                if (widthBytes == 4)
                {
                    //db({ 0xcc });
                    INSERT_LINE_4(0); INSERT_LINE_4(0); INSERT_LINE_4(0); INSERT_LINE_4(1);
                }
                if (widthBytes == 8)
                {
                    INSERT_LINE_8(0); INSERT_LINE_8(0);
                    INSERT_LINE_8(0); INSERT_LINE_8(0);
                    INSERT_LINE_8(0); INSERT_LINE_8(0);
                    INSERT_LINE_8(0); INSERT_LINE_8(1);
                }
                if (widthBytes == 16)
                {
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(1);
                }
                if (widthBytes == 32)
                {
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(1);
                }
                if (widthBytes == 64)
                {
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                    INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(1);
               }
#undef INSERT_LINE_4
#undef INSERT_LINE_8 
#undef INSERT_LINE_16
#undef INSERT_LINE_32
#undef INSERT_LINE_64
            } 
            else
            {
#define INSERT_LINE_4(last)                    \
                {                              \
                    movq(m0, ptr[r2]);         \
                    movq(m1, ptr[r4]);         \
                    psubw(m0, m1);             \
                    movq(ptr[r0], m0);         \
                    if (!last)                 \
                    {                          \
                        lea(r0, ptr[r0 + r1 * sizeof(Sample)]); \
                        lea(r2, ptr[r2 + r3 * sizeof(Sample)]); \
                        lea(r4, ptr[r4 + r5 * sizeof(Sample)]); \
                    }                          \
                }
#define INSERT_LINE_8(last)                    \
                {                              \
                    movdqa(m0, ptr[r2]);       \
                    movdqa(m1, ptr[r4]);       \
                    psubw(m0, m1);             \
                    movdqa(ptr[r0], m0);       \
                    if (!last)                 \
                    {                          \
                        lea(r0, ptr[r0 + r1 * sizeof(Sample)]); \
                        lea(r2, ptr[r2 + r3 * sizeof(Sample)]); \
                        lea(r4, ptr[r4 + r5 * sizeof(Sample)]); \
                    }                          \
                }
#define INSERT_LINE_16(last)                    \
                {                               \
                    movdqa(m0, ptr[r2]);        \
                    movdqa(m1, ptr[r4]);        \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0], m0);        \
                    movdqa(m0, ptr[r2 + 0x10]); \
                    movdqa(m1, ptr[r4 + 0x10]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x10], m0); \
                    if (!last)                  \
                    {                           \
                        lea(r0, ptr[r0 + r1 * sizeof(Sample)]);  \
                        lea(r2, ptr[r2 + r3 * sizeof(Sample)]);  \
                        lea(r4, ptr[r4 + r5 * sizeof(Sample)]);  \
                    }                           \
                }
#define INSERT_LINE_32(last)                    \
                {                               \
                    movdqa(m0, ptr[r2]);        \
                    movdqa(m1, ptr[r4]);        \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0], m0);        \
                    movdqa(m0, ptr[r2 + 0x10]); \
                    movdqa(m1, ptr[r4 + 0x10]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x10], m0); \
                    movdqa(m0, ptr[r2 + 0x20]); \
                    movdqa(m1, ptr[r4 + 0x20]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x20], m0); \
                    movdqa(m0, ptr[r2 + 0x30]); \
                    movdqa(m1, ptr[r4 + 0x30]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x30], m0); \
                    if (!last)                  \
                    {                           \
                        lea(r0, ptr[r0 + r1 * sizeof(Sample)]);  \
                        lea(r2, ptr[r2 + r3 * sizeof(Sample)]);  \
                        lea(r4, ptr[r4 + r5 * sizeof(Sample)]);  \
                    }                           \
                }
#define INSERT_LINE_64(last)                    \
                {                               \
                    movdqa(m0, ptr[r2]);        \
                    movdqa(m1, ptr[r4]);        \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0], m0);        \
                    movdqa(m0, ptr[r2 + 0x10]); \
                    movdqa(m1, ptr[r4 + 0x10]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x10], m0); \
                    movdqa(m0, ptr[r2 + 0x20]); \
                    movdqa(m1, ptr[r4 + 0x20]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x20], m0); \
                    movdqa(m0, ptr[r2 + 0x30]); \
                    movdqa(m1, ptr[r4 + 0x30]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x30], m0); \
                    movdqa(m0, ptr[r2 + 0x40]); \
                    movdqa(m1, ptr[r4 + 0x40]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x40], m0); \
                    movdqa(m0, ptr[r2 + 0x50]); \
                    movdqa(m1, ptr[r4 + 0x50]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x50], m0); \
                    movdqa(m0, ptr[r2 + 0x60]); \
                    movdqa(m1, ptr[r4 + 0x60]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x60], m0); \
                    movdqa(m0, ptr[r2 + 0x70]); \
                    movdqa(m1, ptr[r4 + 0x70]); \
                    psubw(m0, m1);              \
                    movdqa(ptr[r0 + 0x70], m0); \
                     if (!last)                 \
                     {                          \
                         lea(r0, ptr[r0 + r1 * sizeof(Sample)]); \
                         lea(r2, ptr[r2 + r3 * sizeof(Sample)]); \
                         lea(r4, ptr[r4 + r5 * sizeof(Sample)]); \
                     }                          \
                }
                if (widthBytes == 8)
                {
                    INSERT_LINE_4(0); INSERT_LINE_4(0); INSERT_LINE_4(0); INSERT_LINE_4(1);
                }
                if (widthBytes == 16)
                {
                    INSERT_LINE_8(0); INSERT_LINE_8(0);
                    INSERT_LINE_8(0); INSERT_LINE_8(0);
                    INSERT_LINE_8(0); INSERT_LINE_8(0);
                    INSERT_LINE_8(0); INSERT_LINE_8(1);
                }
                if (widthBytes == 32)
                {
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(0);
                    INSERT_LINE_16(0); INSERT_LINE_16(1);
                }
                if (widthBytes == 64)
                {
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(0);
                    INSERT_LINE_32(0); INSERT_LINE_32(1);
                }
                if (widthBytes == 128)
                {
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0);
                   INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(0); INSERT_LINE_64(1);
               }
#undef INSERT_LINE_4
#undef INSERT_LINE_8 
#undef INSERT_LINE_16
#undef INSERT_LINE_32
#undef INSERT_LINE_64
            }
        }
    };

template <typename Sample>
struct ResidualAVX :
    Jit::Function
    {
        ResidualAVX(Jit::Buffer *buffer, int width, int height) :
        Jit::Function(buffer, Jit::CountArguments<havoc_residual<uint8_t>>::value),
        width(width),
        height(height)
        {
            this->build();
        }

        int width, height;

        void assemble() override
        {
/*
            auto &r0 = arg64(0);// Sample   *pA
            auto &r1 = arg64(1);// intptr_t strideA
            auto &r2 = arg64(2);// Sample   *pB
            auto &r3 = arg64(3);// intptr_t strideB
            auto &r4 = arg64(4);// int w
            auto &r5 = arg64(5);// int h
*/

            auto &r0 = arg64(0); // Sample   *pRes
            auto &r1 = arg64(1); // intptr_t strideRes
            auto &r2 = arg64(2); // Sample   *pA
            auto &r3 = arg64(3); // intptr_t strideA
            auto &r4 = arg64(4); // Sample   *pB
            auto &r5 = arg64(5); // intptr_t strideB
            auto &r6 = arg64(6); // int w
            auto &r7 = arg64(7); // int h

            if (sizeof(Sample) == 2)
            {
                shl(r3, 1);
                shl(r5, 1);
            }

            int const widthBytes = width * sizeof(Sample);

            if (widthBytes % 16 == 0)
            {
                auto &m0 = regXmm(0);
                auto &m1 = regXmm(1);
                auto &m2 = regXmm(2);
                auto &m3 = regXmm(3);
                auto &m4 = regXmm(4);
                auto &m5 = regXmm(5);

                vpxor(m0, m0);
                vpxor(m5, m5);

                L("loop");
                {
                    for (int x = 0; x < widthBytes; x += 16)
                    {
                        vmovdqa(m1, ptr[r2 + x]);
                        vmovdqa(m2, ptr[r4 + x]);
                        if (sizeof(Sample) == 1)
                        {
                            vpunpckhbw(m3, m1, m5);
                            vpunpckhbw(m4, m2, m5);
                            vpunpcklbw(m1, m5);
                            vpunpcklbw(m2, m5);
                            vpsubw(m1, m2);
                            vmovdqa(ptr[r0 + (2*x)], m1);
                            vpsubw(m3, m4);
                            vmovdqa(ptr[r0 + (2*x) + 0x10], m3);
                        }
                        else
                        {
                            vpsubw(m1, m2);
                            vmovdqa(ptr[r0 + x], m1);
                        }
                    }
                    lea(r0, ptr[r0 + r1 * 2]);
                    lea(r2, ptr[r2 + r3]);
                    lea(r4, ptr[r4 + r5]);
                }
                dec(Xbyak::Reg32(r7.getIdx()));
                jg("loop");
            }
        }
    };


template <typename Sample>
void havoc_populate_residual(havoc_table_residual<Sample> *table, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    *havoc_get_residual(table, 2) = 0;
    *havoc_get_residual(table, 3) = 0;
    *havoc_get_residual(table, 4) = 0;
    *havoc_get_residual(table, 5) = 0;
    *havoc_get_residual(table, 6) = 0;

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
    {
        *havoc_get_residual(table, 2) = havoc_residual_c_ref<Sample>;
        *havoc_get_residual(table, 3) = havoc_residual_c_ref<Sample>;
        *havoc_get_residual(table, 4) = havoc_residual_c_ref<Sample>;
        *havoc_get_residual(table, 5) = havoc_residual_c_ref<Sample>;
        *havoc_get_residual(table, 6) = havoc_residual_c_ref<Sample>;
    }

#if 1
    if (buffer.isa & HAVOC_SSE2)
    {
        {
            ResidualSSE<Sample> residual(&buffer, 4, 4);
            *havoc_get_residual(table, 2) = residual;
        }
        {
            ResidualSSE<Sample> residual(&buffer, 8, 8);
            *havoc_get_residual(table, 3) = residual;
        }
        {
            ResidualSSE<Sample> residual(&buffer, 16, 16);
            *havoc_get_residual(table, 4) = residual;
        }
        {
            ResidualSSE<Sample> residual(&buffer, 32, 32);
            *havoc_get_residual(table, 5) = residual;
        }
        {
            ResidualSSE<Sample> residual(&buffer, 64, 64);
            *havoc_get_residual(table, 6) = residual;
        }
    }

    if (buffer.isa & HAVOC_AVX2)
    {
        if (sizeof(Sample) == 2)
        {
            ResidualAVX<Sample> residual(&buffer, 8, 8);
            *havoc_get_residual(table, 3) = residual;
        }
        {
            ResidualAVX<Sample> residual(&buffer, 16, 16);
            *havoc_get_residual(table, 4) = residual;
        }
        {
            ResidualAVX<Sample> residual(&buffer, 32, 32);
            *havoc_get_residual(table, 5) = residual;
        }
        {
            ResidualAVX<Sample> residual(&buffer, 64, 64);
            *havoc_get_residual(table, 6) = residual;
        }
    }
#endif
}

struct BoundResidualBase
{
    int    log2TrafoSize;
    int    bits;
};

template <typename Sample>
struct BoundResidual :
    BoundResidualBase
    {
        Sample  *srcA, *srcB;
        int16_t *res;
        havoc_residual<Sample> *f;

        int init(void *p, havoc_code code)
        {
            auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

            auto s = this;

            havoc_table_residual<Sample> table;

            havoc_populate_residual(&table, code);

            s->f = *havoc_get_residual(&table, s->log2TrafoSize);

            if (buffer.isa == HAVOC_C_REF)
            {
                const int nCbS = 1 << s->log2TrafoSize;
                printf("\t%d bits %dx%d : ", sizeof(Sample) == 2 ? 10 : 8, nCbS, nCbS);
            }

            return !!s->f;
        }
    };

int init_residual(void *p, havoc_code code)
{
    BoundResidualBase *b = (BoundResidualBase *)p;
    if (b->bits == 8)
        return static_cast<BoundResidual<uint8_t> *>(b)->init(p, code);
    else
        return static_cast<BoundResidual<uint16_t> *>(b)->init(p, code);
};

template <typename Sample>
void invokeResidual(BoundResidualBase *b, int n)
{
    BoundResidual<Sample> *s = static_cast<BoundResidual<Sample> *>(b);
    const int nCbS = 1 << s->log2TrafoSize;
    while (n--)
    {
        s->f(s->res, 2*nCbS, s->srcA, 2*nCbS, s->srcB, 2*nCbS, nCbS, nCbS);
    }
}

void invoke_residual(void *p, int n)
{
    BoundResidualBase *b = (BoundResidualBase *)p;
    if (b->bits == 8)
        invokeResidual<uint8_t>(b, n);
    else
        invokeResidual<uint16_t>(b, n);
};

template <typename Sample>
int mismatch_residual(void *boundRef, void *boundTest)
{
    BoundResidualBase     *ref     = (BoundResidualBase *)boundRef;
    BoundResidualBase     *test    = (BoundResidualBase *)boundTest;
    BoundResidual<Sample> *sref    = static_cast<BoundResidual<Sample> *>(boundRef);
    BoundResidual<Sample> *stest   = static_cast<BoundResidual<Sample> *>(boundTest);
    auto const            nTbSref  = 1 << sref->log2TrafoSize;
    auto const            nTbStest = 1 << stest->log2TrafoSize;

    if (nTbStest == nTbSref)
    {
        if (memcmp(sref->res, stest->res, 32 * sizeof(Sample) * nTbSref))
            return 1;
        else
            return 0;
    }
    else
        return 0;
}

template <typename Sample>
void testResidual(int *error_count, havoc_instruction_set mask)
{
    int16_t res [64 * 64 * 2];
    Sample  srcA[64 * 64 * 2];
    Sample  srcB[64 * 64 * 2];

    for (int i = 0; i < 64 * 64 * 2; ++i)
    {
        srcA[i] = rand() & 0x3ff;
        srcB[i] = rand() & 0x3ff;
    }

    BoundResidual<Sample> b[2];

    b[0].bits = sizeof(Sample) == 1 ? 8 : 10;
    b[0].res  = res;
    b[0].srcA = srcA;
    b[0].srcB = srcB;

    for (b[0].log2TrafoSize = 2; b[0].log2TrafoSize <= 6; ++b[0].log2TrafoSize)
    {
        b[1] = b[0];
        *error_count += havoc_test(&b[0], &b[1], init_residual, invoke_residual, mismatch_residual<Sample>, mask, 10);
    }
}

void havoc_test_residual(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_residual - Calculate residuals\n");
    testResidual<uint8_t>(error_count, mask);
    testResidual<uint16_t>(error_count, mask);
}
