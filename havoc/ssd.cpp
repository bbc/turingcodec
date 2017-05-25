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

#include "ssd.h"
#include "havoc_test.h"
#include "Jit.h"
#include <stdlib.h>
#include <assert.h>


template <typename Sample>
static uint32_t havoc_ssd_c_ref(const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB, int w, int h)
{
    uint32_t ssd = 0;
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            const int diff = pA[x + y * strideA] - pB[x + y * strideB];
            ssd += diff * diff;
        }
    }
    if (sizeof(Sample) == 2)
        ssd >>= 4;
    return ssd;
}


#if USE_HM_DERIVED

template <typename Sample>
static uint32_t havoc_ssd_c_opt_4x4(const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB, int w, int h)
{
    uint32_t ssd = 0;
    for (int y = 0; y < h; ++y)
    {
        int diff;
        diff = pA[0] - pB[0];  ssd += diff * diff;
        diff = pA[1] - pB[1];  ssd += diff * diff;
        diff = pA[2] - pB[2];  ssd += diff * diff;
        diff = pA[3] - pB[3];  ssd += diff * diff;
        pA += strideA;
        pB += strideB;
    }

    if (sizeof(Sample) == 2)
        ssd >>= 4;

    return ssd;
}

#endif


#define ORDER(a, b, c, d) ((a << 6) | (b << 4) | (c << 2) | d)


template <typename Sample>
struct Ssd :
    Jit::Function
    {
        Ssd(Jit::Buffer *buffer, int width, int height) :
        Jit::Function(buffer, Jit::CountArguments<havoc_ssd<uint8_t>>::value),
        width(width),
        height(height)
        {
            this->build();
        }

        int width, height;

        void assemble() override
                {
            auto &r0 = arg64(0);
            auto &r1 = arg64(1);
            auto &r2 = arg64(2);
            auto &r3 = arg64(3);
            auto &r4 = arg64(4);
            auto &r5 = arg64(5);

            if (sizeof(Sample) == 2)
            {
                shl(r1, 1);
                shl(r3, 1);
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
                        vmovdqa(m1, ptr[r0 + x]);
                        vmovdqa(m2, ptr[r2 + x]);
                        if (sizeof(Sample) == 1)
                        {
                            vpunpckhbw(m3, m1, m5);
                            vpunpckhbw(m4, m2, m5);
                            vpunpcklbw(m1, m5);
                            vpunpcklbw(m2, m5);
                        }
                        vpsubw(m1, m2);
                        if (sizeof(Sample) == 1) vpsubw(m3, m4);
                        vpmaddwd(m1, m1);
                        if (sizeof(Sample) == 1) vpmaddwd(m3, m3);
                        vpaddd(m0, m1);
                        if (sizeof(Sample) == 1) vpaddd(m0, m3);
                    }
                    lea(r0, ptr[r0 + r1]);
                    lea(r2, ptr[r2 + r3]);
                }
                dec(Xbyak::Reg32(r5.getIdx()));
                jg("loop");

                vpshufd(m1, m0, ORDER(3, 2, 3, 2));
                vpaddd(m0, m1); // Review: phaddd might be better here
                vpshufd(m1, m0, ORDER(3, 2, 0, 1));
                vpaddd(m0, m1);
                vmovd(eax, m0);
                if (sizeof(Sample) == 2)
                    shr(eax, 4);
            }
                }
    };


template <typename Sample>
void havoc_populate_ssd(havoc_table_ssd<Sample> *table, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    *havoc_get_ssd(table, 2) = 0;
    *havoc_get_ssd(table, 3) = 0;
    *havoc_get_ssd(table, 4) = 0;
    *havoc_get_ssd(table, 5) = 0;
    *havoc_get_ssd(table, 6) = 0;

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
    {
        *havoc_get_ssd(table, 2) = havoc_ssd_c_ref<Sample>;
        *havoc_get_ssd(table, 3) = havoc_ssd_c_ref<Sample>;
        *havoc_get_ssd(table, 4) = havoc_ssd_c_ref<Sample>;
        *havoc_get_ssd(table, 5) = havoc_ssd_c_ref<Sample>;
        *havoc_get_ssd(table, 6) = havoc_ssd_c_ref<Sample>;
    }

#if USE_HM_DERIVED
    if (buffer.isa & HAVOC_C_OPT)
    {
        *havoc_get_ssd(table, 2) = havoc_ssd_c_opt_4x4<Sample>;
    }
#endif

    if (buffer.isa & HAVOC_AVX)
    {
        if (sizeof(Sample) == 2)
        {
            Ssd<Sample> ssd(&buffer, 8, 8);
            *havoc_get_ssd(table, 3) = ssd;
        }
        {
            Ssd<Sample> ssd(&buffer, 16, 16);
            *havoc_get_ssd(table, 4) = ssd;
        }
        {
            Ssd<Sample> ssd(&buffer, 32, 32);
            *havoc_get_ssd(table, 5) = ssd;
        }
        {
            Ssd<Sample> ssd(&buffer, 64, 64);
            *havoc_get_ssd(table, 6) = ssd;
        }
    }
}



struct BoundSsdBase
{
    int log2TrafoSize;
    uint32_t ssd;
    int bits;
};


template <typename Sample>
struct BoundSsd :
    BoundSsdBase
    {
        Sample *srcA, *srcB;
        havoc_ssd<Sample> *f;

        int init(void *p, havoc_code code)
        {
            auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

            auto s = this;

            havoc_table_ssd<Sample> table;

            havoc_populate_ssd(&table, code);

            s->f = *havoc_get_ssd(&table, s->log2TrafoSize);

            if (buffer.isa == HAVOC_C_REF)
            {
                const int nCbS = 1 << s->log2TrafoSize;
                printf("\t%d bits %dx%d : ", sizeof(Sample) == 2 ? 10 : 8, nCbS, nCbS);
            }

            return !!s->f;
        }
    };

int init_ssd(void *p, havoc_code code)
{
    BoundSsdBase *b = (BoundSsdBase *)p;
    if (b->bits == 8)
        return static_cast<BoundSsd<uint8_t> *>(b)->init(p, code);
    else
        return static_cast<BoundSsd<uint16_t> *>(b)->init(p, code);
};



template <typename Sample>
void invokeSsd(BoundSsdBase *b, int n)
{
    BoundSsd<Sample> *s = static_cast<BoundSsd<Sample> *>(b);
    const int nCbS = 1 << s->log2TrafoSize;
    while (n--)
    {
        s->ssd = s->f(s->srcA, 2*nCbS, s->srcB, 2*nCbS, nCbS, nCbS);
    }
}


void invoke_ssd(void *p, int n)
{
    BoundSsdBase *b = (BoundSsdBase *)p;
    if (b->bits == 8)
        invokeSsd<uint8_t>(b, n);
    else
        invokeSsd<uint16_t>(b, n);
};



int mismatch_ssd(void *boundRef, void *boundTest)
{
    BoundSsdBase *ref = (BoundSsdBase *)boundRef;
    BoundSsdBase *test = (BoundSsdBase *)boundTest;

    return ref->ssd != test->ssd;
}


template <typename Sample>
void testSsd(int *error_count, havoc_instruction_set mask)
{
    Sample  srcA[64 * 64 * 2];
    Sample  srcB[64 * 64 * 2];

    for (int i = 0; i < 64 * 64 * 2; ++i)
    {
        srcA[i] = rand() & 0x3ff;
        srcB[i] = rand() & 0x3ff;
    }

    BoundSsd<Sample> b[2];

    b[0].bits = sizeof(Sample) == 1 ? 8 : 10;
    b[0].srcA = srcA;
    b[0].srcB = srcB;

    for (b[0].log2TrafoSize = 2; b[0].log2TrafoSize <= 6; ++b[0].log2TrafoSize)
    {
        b[1] = b[0];
        *error_count += havoc_test(&b[0], &b[1], init_ssd, invoke_ssd, mismatch_ssd, mask, 10);
    }
}


void havoc_test_ssd(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_ssd - Sum of Square Differences\n");
    testSsd<uint8_t>(error_count, mask);
    testSsd<uint16_t>(error_count, mask);
}
