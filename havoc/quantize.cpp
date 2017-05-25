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

#include "quantize.h"
#include "havoc_test.h"
#include "Jit.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>


static int Clip3(int min, int max, int value)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}


static void havoc_quantize_inverse_c_ref(int16_t *dst, const int16_t *src, int scale, int shift, int n)
{
    while (n--)
    {
        *dst++ = (int16_t)Clip3(
                -32768,
                32767,
                ((*src++ * scale) + (1 << (shift - 1))) >> shift);
    }
}

struct InverseQuantise :
    Jit::Function
    {
        InverseQuantise(Jit::Buffer *buffer, bool round) :
            Jit::Function(buffer, Jit::CountArguments<havoc_quantize_inverse>::value),
            round(round)
            {
                this->build();
            }

        bool round;
        Xbyak::Label ones_w;

        void data()
        {
            if (this->round)
            {
                align();
                L(ones_w);
                dw({ 1 }, 8);
            }
        }

        void assemble()
        {
            auto &r0 = arg64(0);
            auto &r1 = arg64(1);
            auto &r2 = arg64(2);
            auto &r3 = arg64(3);
            auto &r4 = arg64(4);

            auto &m0 = regXmm(0);
            auto &m1 = regXmm(1);
            auto &m2 = regXmm(2);
            auto &m3 = regXmm(3);
            auto &m4 = regXmm(4);
            auto &m5 = regXmm(5);

            Xbyak::Reg32 r2d(r2.getIdx());
            Xbyak::Reg32 r3d(r3.getIdx());
            Xbyak::Reg8 r3b(r3.getIdx());
            Xbyak::Reg32 r4d(r4.getIdx());

            if (this->round)
            {
                movdqa(m0, ptr[rip + ones_w]);

                //  m1 = shift
                movd(m1, r3d);
                add(r3b, 15);
                bts(r2d, r3d);
                // r2d = (0x10000 << (shift - 1)) + scale
            }
            else
            {

                pxor(m0, m0);

                // scale >>= shift
                auto &r5 = reg64(5);
                mov(r5, rcx);
                mov(cl, r3d);
                shr(r2d, cl);
                mov(rcx, r5);
            }

            movd(m2, r2d);
            pshufd(m2, m2, 0);

            shr(r4d, 4);
            // r4 = n / 16

            L("loop");
            {
                for (int offset = 0; offset < 32; offset += 16)
                {
                    movdqa(m4, ptr[r1 + offset]);
                    // m4 = src[7], src[6], src[5], src[4], src[3], src[2], src[1], src[0]

                    movdqa(m5, m4);

                    punpcklwd(m5, m0);
                    // m5 = 1, src[3], 1, src[2], 1, src[1], 1, src[0]

                    punpckhwd(m4, m0);
                    // m4 = 1, src[7], 1, src[6], 1, src[5], 1, src[4]

                    pmaddwd(m4, m2);
                    // m4 = (1 << (shift - 1)) + src[7] * scale, (1 << (shift - 1)) + src[6] * scale, (1 << (shift - 1)) + src[5] * scale, (1 << (shift - 1)) + src[4] * scale

                    pmaddwd(m5, m2);
                    // m5 = (1 << (shift - 1)) + src[3] * scale, (1 << (shift - 1)) + src[2] * scale, (1 << (shift - 1)) + src[1] * scale, (1 << (shift - 1)) + src[0] * scale

                    if (this->round)
                    {
                        psrad(m4, m1);
                        // m4 = ((1 << (shift - 1)) + src[7] * scale) >> shift, ((1 << (shift - 1)) + src[6] * scale) >> shift, ((1 << (shift - 1)) + src[5] * scale) >> shift, ((1 << (shift - 1)) + src[4] * scale) >> shift

                        psrad(m5, m1);
                        // m5 = ((1 << (shift - 1)) + src[3] * scale) >> shift, ((1 << (shift - 1)) + src[2] * scale) >> shift, ((1 << (shift - 1)) + src[1] * scale) >> shift, ((1 << (shift - 1)) + src[0] * scale) >> shift
                    }

                    packssdw(m5, m4);
                    // m5 = ((1 << (shift - 1)) + src[7] * scale) >> shift, ((1 << (shift - 1)) + src[6] * scale) >> shift, ((1 << (shift - 1)) + src[5] * scale) >> shift, ((1 << (shift - 1)) + src[4] * scale) >> shift, ((1 << (shift - 1)) + src[3] * scale) >> shift, ((1 << (shift - 1)) + src[2] * scale) >> shift, ((1 << (shift - 1)) + src[1] * scale) >> shift, ((1 << (shift - 1)) + src[0] * scale) >> shift

                    movdqa(ptr[r0 + offset], m5);
                }

                add(r1, 32);
                add(r0, 32);
            }
            dec(r4d);
            jg("loop");
        }
    };


static havoc_quantize_inverse * get_quantize_inverse(havoc_code code)
{
}


void havoc_populate_quantize_inverse(havoc_table_quantize_inverse *table, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    table->p[0] = 0;
    table->p[1] = 0;

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
    {
        table->p[0] = havoc_quantize_inverse_c_ref;
        table->p[1] = havoc_quantize_inverse_c_ref;
    }

    if (buffer.isa & HAVOC_SSE41)
    {
        InverseQuantise inverseQuantise0(&buffer, false);
        table->p[0] = inverseQuantise0;
        InverseQuantise inverseQuantise1(&buffer, true);
        table->p[1] = inverseQuantise1;
    }
}


typedef struct
{
    int16_t *src;
    HAVOC_ALIGN(32, int16_t, dst[32*32]);
    havoc_quantize_inverse *f;
    int scale;
    int shift;
    int log2TrafoSize;
}
havoc_bound_quantize_inverse;


int init_quantize_inverse(void *p, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    havoc_bound_quantize_inverse *s = (havoc_bound_quantize_inverse *)p;

    havoc_table_quantize_inverse table;

    havoc_populate_quantize_inverse(&table, code);

    s->f = *havoc_get_quantize_inverse(&table, s->scale, s->shift);

    if (s->f && buffer.isa == HAVOC_C_REF)
    {
        const int nCbS = 1 << s->log2TrafoSize;
        printf("\t%s%dx%d : ", s->scale < 1000 ? "round " : "", nCbS, nCbS);
    }

    return !!s->f;
}


void invoke_quantize_inverse(void *p, int count)
{
    havoc_bound_quantize_inverse *s = (havoc_bound_quantize_inverse *)p;
    s->dst[0] = rand();
    while (count--)
    {
        const int n = 1 << (2 * s->log2TrafoSize);
        s->f(s->dst, s->src, s->scale, s->shift, n);
    }
}


int mismatch_quantize_inverse(void *boundRef, void *boundTest)
{
    havoc_bound_quantize_inverse *ref = (havoc_bound_quantize_inverse *)boundRef;
    havoc_bound_quantize_inverse *test = (havoc_bound_quantize_inverse *)boundTest;

    const int nCbS = 1 << ref->log2TrafoSize;

    return memcmp(ref->dst, test->dst, nCbS * nCbS * sizeof(int16_t));
}


void havoc_test_quantize_inverse(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_quantize_inverse - Inverse Quantization (\"scaling\")\n");

    HAVOC_ALIGN(32, int16_t, src[32 * 32]);

    for (int x = 0; x < 32 * 32; ++x) src[x] = (rand() & 0xff) - 0x100;

    havoc_bound_quantize_inverse b[2];
    b[0].src = src;

    int const log2TrafoSize = 3;
    int const bitDepth = 8;

    for (int i = 0; i < 2; ++i)
    {
        b[0].scale = i ? 52224 : 51;

        for (b[0].log2TrafoSize = 2; b[0].log2TrafoSize <= 5; ++b[0].log2TrafoSize)
        {
            b[0].shift = b[0].log2TrafoSize - 1 + bitDepth - 8;
            b[1] = b[0];
            *error_count += havoc_test(&b[0], &b[1], init_quantize_inverse, invoke_quantize_inverse, mismatch_quantize_inverse, mask, 10);
        }
    }
}


static int havoc_quantize_c_ref(int16_t *dst, const int16_t *src, int scale, int shift, int offset, int n)
{
    assert(scale < 0x8000);
    assert(offset < 0x8000);
    assert(shift >= 16);
    assert(shift <= 27);

    offset <<= (shift - 16);
    assert(offset < 0x4000000);

    int cbf = 0;
    while (n--)
    {
        int x = *src++;
        int sign = x < 0 ? -1 : 1;

        x = abs(x);
        x = ((x * scale) + offset) >> shift;
        x *= sign;
        x = Clip3(-32768, 32767, x);

        cbf |= x;

        *dst++ = x;
    }
    return cbf;
}


struct Quantise :
    Jit::Function
{
    Quantise(Jit::Buffer *buffer) :
        Jit::Function(buffer, Jit::CountArguments<havoc_quantize>::value)
        {
            this->build();
        }

    void assemble()
    {
        auto &r0 = arg64(0);
        auto &r1 = arg64(1);
        auto &r2 = arg64(2);
        auto &r3 = arg64(3);
        auto &r4 = arg64(4);
        auto &r5 = arg64(5);

        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        auto &m4 = regXmm(4);
        auto &m5 = regXmm(5);
        auto &m6 = regXmm(6);
        auto &m7 = regXmm(7);

        Xbyak::Reg32 r2d(r2.getIdx());
        Xbyak::Reg32 r3d(r3.getIdx());
        Xbyak::Reg32 r4d(r4.getIdx());
        Xbyak::Reg32 r5d(r5.getIdx());

        movd(m1, r3d);
        // m1 = shift

        bts(r2d, r3d);
        // r2d = (1 << shift) + scale

        movd(m2, r2d);
        pshufd(m2, m2, 0);
        // m2 = 1 << (shift - 16), scale, 1 << (shift - 16), scale, 1 << (shift - 16), scale, 1 << (shift - 16), scale

        movd(m3, r4d);
        pshuflw(m3, m3, 0);
        pshufd(m3, m3, 0);
        // m3 = offset, offset, offset, offset, offset, offset, offset, offset

        pxor(m0, m0);
        // m0 = 0

        shr(r5d, 4);
        // r5 = n / 16

        L("loop");
        {
            for (int offset = 0; offset < 32; offset += 16)
            {
                movdqa(m4, ptr[r1 + offset]);
                // m4 = src[7], src[6], src[5], src[4], src[3], src[2], src[1], src[0]

                pabsw(m5, m4);
                // m5 = abs(src[7]), abs(src[6]), abs(src[5]), abs(src[4]), abs(src[3]), abs(src[2]), abs(src[1]), abs(src[0])

                movdqa(m6, m5);
                punpcklwd(m6, m3);
                // m6 = offset, abs(src[3]), offset, abs(src[2]), offset, abs(src[1]), offset, abs(src[0])

                punpckhwd(m5, m3);
                // m5 = offset, abs(src[7]), offset, abs(src[6]), offset, abs(src[5]), offset, abs(src[4])

                pmaddwd(m6, m2);
                // m6 = (offset << (shift - 16)) + abs(src[3])*scale, (offset << (shift - 16)) + abs(src[2])*scale, (offset << (shift - 16)) + abs(src[1])*scale, (offset << (shift - 16)) + abs(src[0])*scale

                psrad(m6, m1);
                // m6 = (offset << (shift - 16)) + abs(src[3])*scale >> shift, (offset << (shift - 16)) + abs(src[2])*scale >> shift, (offset << (shift - 16)) + abs(src[1])*scale >> shift, (offset << (shift - 16)) + abs(src[0])*scale >> shift

                pmaddwd(m5, m2);
                // m5 = (offset << (shift - 16)) + abs(src[7])*scale, (offset << (shift - 16)) + abs(src[6])*scale, (offset << (shift - 16)) + abs(src[5])*scale, (offset << (shift - 16)) + abs(src[4])*scale,

                psrad(m5, m1);
                // m5 = (offset << (shift - 16)) + abs(src[7])*scale >> shift, (offset << (shift - 16)) + abs(src[6])*scale >> shift, (offset << (shift - 16)) + abs(src[5])*scale >> shift, (offset << (shift - 16)) + abs(src[4])*scale >> shift

                punpcklwd(m7, m4);
                // m7 = (src[3] << 16) + 0x ? ? ? ? , (src[2] << 16) + 0x ? ? ? ? , (src[1] << 16) + 0x ? ? ? ? , (src[0] << 16) + 0x ? ? ? ?

                psignd(m6, m7);
                // m6 = ((offset << (shift - 16)) + abs(src[3])*scale >> shift)*sign(src[3]), ((offset << (shift - 16)) + abs(src[2])*scale >> shift)*sign(src[2]), ((offset << (shift - 16)) + abs(src[1])*scale >> shift)*sign(src[1]), ((offset << (shift - 16)) + abs(src[0])*scale >> shift)*sign(src[0])
                // m6 = dst[3], dst[2], dst[1], dst[0]

                punpckhwd(m4, m4);
                // m4 = (src[7] << 16) + 0x ? ? ? ? , (src[6] << 16) + 0x ? ? ? ? , (src[5] << 16) + 0x ? ? ? ? , (src[4] << 16) + 0x ? ? ? ?

                psignd(m5, m4);
                // m5 = ((offset << (shift - 16)) + abs(src[7])*scale >> shift)*sign(src[7]), ((offset << (shift - 16)) + abs(src[6])*scale >> shift)*sign(src[6]), ((offset << (shift - 16)) + abs(src[5])*scale >> shift)*sign(src[5]), ((offset << (shift - 16)) + abs(src[4])*scale >> shift)*sign(src[4])
                // m5 = dst[7], dst[6], dst[5], dst[4]

                packssdw(m6, m5);
                // m6 = dst[7], dst[6], dst[5], dst[4], dst[3], dst[2], dst[1], dst[0]

                por(m0, m6);
                // m0 is non - zero if we have seen any non - zero quantized coefficients

                movdqa(ptr[r0 + offset], m6);
            }

            add(r1, 32);
            add(r0, 32);
        }
        dec(r5d);
        jg("loop");

        // return zero only if m0 is zero - no non - zero quantized coefficients seen(cbf = 0)
        packsswb(m0, m0);
        packsswb(m0, m0);
        movd(eax, m0);
    }
};



static havoc_quantize * get_quantize(havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    havoc_quantize *f = 0;

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT)) 
        f = havoc_quantize_c_ref;

    Quantise quantise(&buffer);

    if (buffer.isa & HAVOC_SSE41) 
        f = quantise;

    return f;
}


void havoc_populate_quantize(havoc_table_quantize *table, havoc_code code)
{
    table->p = get_quantize(code);
}


typedef struct
{
    int16_t *src;
    HAVOC_ALIGN(32, int16_t, dst[32 * 32]);
    havoc_quantize *f;
    int scale;
    int shift;
    int offset;
    int log2TrafoSize;
    int cbf;
}
havoc_bound_quantize;


int init_quantize(void *p, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
    auto mask = buffer.isa;

    havoc_bound_quantize *s = (havoc_bound_quantize *)p;
    havoc_table_quantize table;
    havoc_populate_quantize(&table, code);
    s->f = *havoc_get_quantize(&table);

    if (buffer.isa == HAVOC_C_REF)
    {
        const int nCbS = 1 << s->log2TrafoSize;
        printf("\t%dx%d : ", nCbS, nCbS);
    }
    return !!s->f;
}


void invoke_quantize(void *p, int iterations)
{
    havoc_bound_quantize *s = (havoc_bound_quantize *)p;
    s->dst[0] = rand();
    while (iterations--)
    {
        const int n = 1 << (2 * s->log2TrafoSize);
        s->cbf = s->f(s->dst, s->src, s->scale, s->shift, s->offset, n);
    }
}


int mismatch_quantize(void *boundRef, void *boundTest)
{
    havoc_bound_quantize *ref = (havoc_bound_quantize *)boundRef;
    havoc_bound_quantize *test = (havoc_bound_quantize *)boundTest;

    const int n = 1 << (2 * ref->log2TrafoSize);

    if (!!ref->cbf != !!test->cbf) return 1;

    return memcmp(ref->dst, test->dst, n * sizeof(int16_t));
}


void havoc_test_quantize(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_quantize - Quantization\n");

    HAVOC_ALIGN(32, int16_t, src[32 * 32]);

    for (int x = 0; x < 32 * 32; ++x)
    {
        src[x] = rand() - rand();
    }

    havoc_bound_quantize b[2];

    b[0].src = src;
    b[0].scale = 51;
    b[0].shift = 20;
    b[0].offset = 14;

    for (b[0].log2TrafoSize = 2; b[0].log2TrafoSize <= 5; ++b[0].log2TrafoSize)
    {
        b[1] = b[0];
        *error_count += havoc_test(&b[0], &b[1], init_quantize, invoke_quantize, mismatch_quantize, mask, 10);
    }
}





static void havoc_quantize_reconstruct_c_ref(uint8_t *rec, intptr_t stride_rec, const uint8_t *predSamples, intptr_t stride_pred, const int16_t *resSamples, int n)
{
    for (int y = 0; y < n; ++y)
    {
        for (int x = 0; x < n; ++x)
        {
            rec[x + y * stride_rec] = (uint8_t)Clip3(0, 255, predSamples[x + y * stride_pred] + resSamples[x + y * n]);
        }
    }

}


#define ORDER(a, b, c, d) ((a << 6) | (b << 4) | (c << 2) | d)

struct QuantiseReconstruct :
    Jit::Function
    {
        QuantiseReconstruct(Jit::Buffer *buffer, int nCbS) :
            Jit::Function(buffer, 6),
            nCbS(nCbS)
        {
            this->build();
        }

        int const nCbS;

        void assemble()
        {
            auto &r0 = arg64(0);
            auto &r1 = arg64(1);
            auto &r2 = arg64(2);
            auto &r3 = arg64(3);
            auto &r4 = arg64(4);
            auto &r5 = arg64(5);

            auto &m0 = regXmm(0);
            auto &m1 = regXmm(1);
            auto &m2 = regXmm(2);
            auto &m3 = regXmm(3);

            Xbyak::Reg32 r5d(r5.getIdx());

            pxor(m0, m0);

            if (nCbS == 4)
            {
                mov(r5d, 2);

                L("loop");
                {
                    movd(m1, ptr[r2]);
                    movd(m2, ptr[r2 + r3]);

                    lea(r2, ptr[r2 + r3 * 2]);

                    punpckldq(m1, m2);
                    punpcklbw(m1, m0);

                    movdqu(m3, ptr[r4]);
                    paddw(m1, m3);
                    lea(r4, ptr[r4 + 16]);
                    packuswb(m1, m1);

                    movd(ptr[r0], m1);
                    pshufd(m1, m1, ORDER(0, 0, 0, 1));
                    movd(ptr[r0 + r1], m1);
                    lea(r0, ptr[r0 + r1 * 2]);
                }
                dec(r5d);
                jg("loop");
            }
            else
            {
                mov(r5d, nCbS);

                L("loop");
                {
                    movq(m1, ptr[r2]);
                    lea(r2, ptr[r2 + r3]);
                    if (nCbS >= 16)
                    {
                        movdqa(m2, m1);
                        punpckhbw(m2, m0);
                    }
                    punpcklbw(m1, m0);

                    paddw(m1, ptr[r4]);
                    lea(r4, ptr[r4 + 2*nCbS]);

                    packuswb(m1, m1);

                    movq(ptr[r0], m1);
                    lea(r0, ptr[r0 + r1]);
                }
                dec(r5d);
                jg("loop");
            }
        }
    };


havoc_quantize_reconstruct * get_quantize_reconstruct(int log2TrafoSize, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
    auto mask = buffer.isa;

    havoc_quantize_reconstruct *f = 0;

    if (mask & (HAVOC_C_REF | HAVOC_C_OPT)) f = havoc_quantize_reconstruct_c_ref;

    if (mask & HAVOC_SSE41)
    {
        const int nCbS = 1 << log2TrafoSize;
        if (nCbS == 4)
        {
            QuantiseReconstruct qr(&buffer, nCbS);
            f = qr;
        }
        if (nCbS == 8)
        {
            QuantiseReconstruct qr(&buffer, nCbS);
            f = qr;
        }
        if (nCbS == 16)
        {
            QuantiseReconstruct qr(&buffer, nCbS);
            f = qr;
        }
        if (nCbS == 32)
        {
            QuantiseReconstruct qr(&buffer, nCbS);
            f = qr;
        }
    }

    return f;
}


void havoc_populate_quantize_reconstruct(havoc_table_quantize_reconstruct *table, havoc_code code)
{
    for (int log2TrafoSize = 2; log2TrafoSize < 6; ++log2TrafoSize)
    {
        *havoc_get_quantize_reconstruct(table, log2TrafoSize) = get_quantize_reconstruct(log2TrafoSize, code);
    }
}


typedef struct
{
    HAVOC_ALIGN(32, uint8_t, rec[32 * 32]);
    intptr_t stride_rec;
    const uint8_t *pred;
    intptr_t stride_pred;
    const int16_t *res;
    int log2TrafoSize;
    havoc_quantize_reconstruct *f;
}
bound_quantize_reconstruct;


int init_quantize_reconstruct(void *p, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    bound_quantize_reconstruct *s = (bound_quantize_reconstruct *)p;

    havoc_table_quantize_reconstruct table;

    havoc_populate_quantize_reconstruct(&table, code);

    s->f = *havoc_get_quantize_reconstruct(&table, s->log2TrafoSize);

    if (buffer.isa == HAVOC_C_REF)
    {
        const int nCbS = 1 << s->log2TrafoSize;
        printf("\t%dx%d : ", nCbS, nCbS);
    }

    return !!s->f;
}


void invoke_quantize_reconstruct(void *p, int n)
{
    bound_quantize_reconstruct *s = (bound_quantize_reconstruct *)p;
    while (n--)
    {
        const int nCbS = 1 << s->log2TrafoSize;
        s->f(s->rec, s->stride_rec, s->pred, s->stride_pred, s->res, nCbS);
    }
}

int mismatch_quantize_reconstruct(void *boundRef, void *boundTest)
{
    bound_quantize_reconstruct *ref = (bound_quantize_reconstruct *)boundRef;
    bound_quantize_reconstruct *test = (bound_quantize_reconstruct *)boundTest;

    const int nCbS = 1 << ref->log2TrafoSize;

    int mismatch = 0;
    for (int y = 0; y < nCbS; ++y)
    {
        mismatch |= memcmp(
                &ref->rec[y * ref->stride_rec],
                &test->rec[y * test->stride_rec],
                nCbS);
    }

    return mismatch;
}



void havoc_test_quantize_reconstruct(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_quantize_reconstruct - Reconstruction\n");

    HAVOC_ALIGN(32, uint8_t, pred[32 * 32]);
    HAVOC_ALIGN(32, int16_t, res[32 * 32]);

    for (int x = 0; x < 32 * 32; ++x)
    {
        pred[x] = rand() & 0xff;
        res[x] = (rand() & 0x1ff) - 0x100;
    }

    bound_quantize_reconstruct b[2];

    b[0].stride_rec = 32;
    b[0].pred = pred;
    b[0].stride_pred = 32;
    b[0].res = res;

    for (b[0].log2TrafoSize = 2; b[0].log2TrafoSize <= 5; ++b[0].log2TrafoSize)
    {
        b[1] = b[0];
        *error_count += havoc_test(&b[0], &b[1], init_quantize_reconstruct, invoke_quantize_reconstruct, mismatch_quantize_reconstruct, mask, 10);
    }
}
