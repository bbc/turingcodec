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


#include "sad.h"
#include "havoc_test.h"
#include "Jit.h"
#include <array>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


template <typename Sample>
struct SadSse2 :
    Jit::Function
{
    SadSse2(Jit::Buffer *buffer, int width, int height) :
        Jit::Function(buffer, Jit::CountArguments<havoc_sad<uint8_t>>::value),
        width(width),
        height(height)
    {
        this->build();
    }

    int width, height;

    template <class A, class B>
    void sad(A const &a, B const &b)
    {
        if (sizeof(Sample) == 2)
        {
            psubw(a, b);
            pabsw(a, a);
        }
        else
            psadbw(a, b);
    }

    void assemble() override
    {
        int widthBytes = width * sizeof(Sample);

        if ((widthBytes % 16) && widthBytes != 8) return;

        auto &src = arg64(0);
        auto &stride_src = arg64(1);
        auto &ref = arg64(2);
        auto &stride_ref = arg64(3);

        if (sizeof(Sample) == 2)
        {
            shl(stride_src, 1);
            shl(stride_ref, 1);
        }

        auto &n = reg64(4);

        const Xbyak::Reg64 *stride_ref_x3 = widthBytes == 16 ? &reg64(5) : 0;
        const Xbyak::Reg64 *stride_src_x3 = widthBytes == 16 ? &reg64(6) : 0;

        if (widthBytes == 8)
        {
            mov(n, height / 2);
        }
        else if (widthBytes == 16)
        {
            mov(n, height / 4);
            lea(*stride_ref_x3, ptr[stride_ref + stride_ref * 2]);
            lea(*stride_src_x3, ptr[stride_src + stride_src * 2]);
        }
        else if (widthBytes == 32)
        {
            mov(n, height / 2);
        }
        else
        {
            mov(n, height);
        }

        auto &xmm0 = regXmm(0);
        auto &xmm1 = regXmm(1);
        auto &xmm2 = regXmm(2);
        auto &xmm3 = regXmm(3);
        auto &xmm4 = regXmm(4);
        auto &xmm5 = regXmm(5);

        pxor(xmm0, xmm0);
        pxor(xmm5, xmm5);

        L("loop");
        {
            if (widthBytes == 8)
            {
                movq(xmm1, ptr[ref]);
                movhps(xmm1, ptr[ref + stride_ref]);
                movq(xmm2, ptr[src]);
                movhps(xmm2, ptr[src + stride_src]);
                sad(xmm1, xmm2);

                lea(ref, ptr[ref + stride_ref * 2]);
                if (sizeof(Sample) == 2)
                {
                    phaddw(xmm1, xmm5);
                    phaddw(xmm1, xmm5);
                    phaddw(xmm1, xmm5);
                }
                paddd(xmm0, xmm1);
                lea(src, ptr[src + stride_src * 2]);
            }
            else if (widthBytes == 16)
            {
                movdqu(xmm1, ptr[ref]);
                movdqu(xmm2, ptr[ref + stride_ref]);
                movdqu(xmm3, ptr[ref + stride_ref * 2]);
                movdqu(xmm4, ptr[ref + *stride_ref_x3]);

                sad(xmm1, ptr[src]);
                sad(xmm2, ptr[src + stride_src]);
                sad(xmm3, ptr[src + stride_src * 2]);
                sad(xmm4, ptr[src + *stride_src_x3]);

                lea(src, ptr[src + stride_src * 4]);
                paddd(xmm1, xmm2);
                lea(ref, ptr[ref + stride_ref * 4]);
                paddd(xmm3, xmm4);
                if (sizeof(Sample) == 2)
                {
                    phaddw(xmm1, xmm5);
                    phaddw(xmm1, xmm5);
                    phaddw(xmm1, xmm5);
                    phaddw(xmm3, xmm5);
                    phaddw(xmm3, xmm5);
                    phaddw(xmm3, xmm5);
                }
                paddd(xmm0, xmm1);
                paddd(xmm0, xmm3);
            }
            else if (widthBytes == 32)
            {
                movdqu(xmm1, ptr[ref]);
                movdqu(xmm2, ptr[ref + 16]);
                movdqu(xmm3, ptr[ref + stride_ref]);
                movdqu(xmm4, ptr[ref + stride_ref + 16]);

                sad(xmm1, ptr[src]);
                sad(xmm2, ptr[src + 16]);
                sad(xmm3, ptr[src + stride_src]);
                sad(xmm4, ptr[src + stride_src + 16]);

                lea(src, ptr[src + stride_src * 2]);
                paddd(xmm1, xmm2);
                lea(ref, ptr[ref + stride_ref * 2]);
                paddd(xmm3, xmm4);
                if (sizeof(Sample) == 2)
                {
                    phaddw(xmm1, xmm5);
                    phaddw(xmm1, xmm5);
                    phaddw(xmm1, xmm5);
                    phaddw(xmm3, xmm5);
                    phaddw(xmm3, xmm5);
                    phaddw(xmm3, xmm5);
                }
                paddd(xmm0, xmm1);
                paddd(xmm0, xmm3);
            }
            else
            {
                for (int dx = 0; dx < widthBytes; dx += 64)
                {
                    movdqu(xmm1, ptr[ref + dx]);
                    movdqu(xmm2, ptr[ref + dx + 16]);
                    if (widthBytes - dx > 32) movdqu(xmm3, ptr[ref + dx + 32]);
                    if (widthBytes - dx > 48) movdqu(xmm4, ptr[ref + dx + 48]);
                    sad(xmm1, ptr[src + dx + 0]);
                    sad(xmm2, ptr[src + dx + 16]);
                    if (widthBytes - dx > 32) sad(xmm3, ptr[src + dx + 32]);
                    if (widthBytes - dx > 48) sad(xmm4, ptr[src + dx + 48]);
                    paddd(xmm1, xmm2);

                    if (widthBytes > 48) paddd(xmm3, xmm4);
                    if (sizeof(Sample) == 2)
                    {
                        phaddw(xmm1, xmm5);
                        phaddw(xmm1, xmm5);
                        phaddw(xmm1, xmm5);
                        if (widthBytes - dx > 32) phaddw(xmm3, xmm5);
                        if (widthBytes - dx > 32) phaddw(xmm3, xmm5);
                        if (widthBytes - dx > 32) phaddw(xmm3, xmm5);
                    }
                    paddd(xmm0, xmm1);
                    if (widthBytes - dx > 32) paddd(xmm0, xmm3);
                }
                add(ref, stride_ref);
                add(src, stride_src);
            }
        }
        dec(n);
        jg("loop");

        movhlps(xmm1, xmm0);
        paddd(xmm0, xmm1);
        movd(eax, xmm0);

        if (sizeof(Sample) == 2)
            shr(eax, 2); // 10-bit - divide SAD by 4
    }
};

template <typename Sample>
static int havoc_sad_c_opt_4(const Sample *src, intptr_t stride_src, const Sample *ref, intptr_t stride_ref, uint32_t rect)
{
    const int width = rect >> 8;
    const int height = rect & 0xff;
    int sad = 0;
    for (int y = 0; y < height; ++y)
    {
        sad += abs((int)src[0] - (int)ref[0]);
        sad += abs((int)src[1] - (int)ref[1]);
        sad += abs((int)src[2] - (int)ref[2]);
        sad += abs((int)src[3] - (int)ref[3]);
        src += stride_src;
        ref += stride_ref;
    }
    if (sizeof(Sample) == 2)
        sad >>= 2;
    return sad;
}

template <typename Sample>
static int havoc_sad_c_opt_8(const Sample *src, intptr_t stride_src, const Sample *ref, intptr_t stride_ref, uint32_t rect)
{
    const int width = rect >> 8;
    const int height = rect & 0xff;
    int sad = 0;
    for (int y = 0; y < height; ++y)
    {
        sad += abs((int)src[0] - (int)ref[0]);
        sad += abs((int)src[1] - (int)ref[1]);
        sad += abs((int)src[2] - (int)ref[2]);
        sad += abs((int)src[3] - (int)ref[3]);
        sad += abs((int)src[4] - (int)ref[4]);
        sad += abs((int)src[5] - (int)ref[5]);
        sad += abs((int)src[6] - (int)ref[6]);
        sad += abs((int)src[7] - (int)ref[7]);
        src += stride_src;
        ref += stride_ref;
    }
    if (sizeof(Sample) == 2)
        sad >>= 2;
    return sad;
}

template <typename Sample>
static int havoc_sad_c_opt_16(const Sample *src, intptr_t stride_src, const Sample *ref, intptr_t stride_ref, uint32_t rect)
{
    const int width = rect >> 8;
    const int height = rect & 0xff;
    int sad = 0;
    for (int y = 0; y < height; ++y)
    {
        sad += abs((int)src[0] - (int)ref[0]);
        sad += abs((int)src[1] - (int)ref[1]);
        sad += abs((int)src[2] - (int)ref[2]);
        sad += abs((int)src[3] - (int)ref[3]);
        sad += abs((int)src[4] - (int)ref[4]);
        sad += abs((int)src[5] - (int)ref[5]);
        sad += abs((int)src[6] - (int)ref[6]);
        sad += abs((int)src[7] - (int)ref[7]);
        sad += abs((int)src[8] - (int)ref[8]);
        sad += abs((int)src[9] - (int)ref[9]);
        sad += abs((int)src[10] - (int)ref[10]);
        sad += abs((int)src[11] - (int)ref[11]);
        sad += abs((int)src[12] - (int)ref[12]);
        sad += abs((int)src[13] - (int)ref[13]);
        sad += abs((int)src[14] - (int)ref[14]);
        sad += abs((int)src[15] - (int)ref[15]);
        src += stride_src;
        ref += stride_ref;
    }
    if (sizeof(Sample) == 2)
        sad >>= 2;
    return sad;
}

template <typename Sample>
static int havoc_sad_c_opt_32(const Sample *src, intptr_t stride_src, const Sample *ref, intptr_t stride_ref, uint32_t rect)
{
    const int width = rect >> 8;
    const int height = rect & 0xff;
    int sad = 0;
    for (int y = 0; y < height; ++y)
    {
        sad += abs((int)src[0] - (int)ref[0]);
        sad += abs((int)src[1] - (int)ref[1]);
        sad += abs((int)src[2] - (int)ref[2]);
        sad += abs((int)src[3] - (int)ref[3]);
        sad += abs((int)src[4] - (int)ref[4]);
        sad += abs((int)src[5] - (int)ref[5]);
        sad += abs((int)src[6] - (int)ref[6]);
        sad += abs((int)src[7] - (int)ref[7]);
        sad += abs((int)src[8] - (int)ref[8]);
        sad += abs((int)src[9] - (int)ref[9]);
        sad += abs((int)src[10] - (int)ref[10]);
        sad += abs((int)src[11] - (int)ref[11]);
        sad += abs((int)src[12] - (int)ref[12]);
        sad += abs((int)src[13] - (int)ref[13]);
        sad += abs((int)src[14] - (int)ref[14]);
        sad += abs((int)src[15] - (int)ref[15]);
        sad += abs((int)src[16] - (int)ref[16]);
        sad += abs((int)src[17] - (int)ref[17]);
        sad += abs((int)src[18] - (int)ref[18]);
        sad += abs((int)src[19] - (int)ref[19]);
        sad += abs((int)src[20] - (int)ref[20]);
        sad += abs((int)src[21] - (int)ref[21]);
        sad += abs((int)src[22] - (int)ref[22]);
        sad += abs((int)src[23] - (int)ref[23]);
        sad += abs((int)src[24] - (int)ref[24]);
        sad += abs((int)src[25] - (int)ref[25]);
        sad += abs((int)src[26] - (int)ref[26]);
        sad += abs((int)src[27] - (int)ref[27]);
        sad += abs((int)src[28] - (int)ref[28]);
        sad += abs((int)src[29] - (int)ref[29]);
        sad += abs((int)src[30] - (int)ref[30]);
        sad += abs((int)src[31] - (int)ref[31]);
        src += stride_src;
        ref += stride_ref;
    }
    if (sizeof(Sample) == 2)
        sad >>= 2;
    return sad;
}

template <typename Sample>
static int havoc_sad_c_opt_64(const Sample *src, intptr_t stride_src, const Sample *ref, intptr_t stride_ref, uint32_t rect)
{
    const int width = rect >> 8;
    const int height = rect & 0xff;
    int sad = 0;
    for (int y = 0; y < height; ++y)
    {
        sad += abs((int)src[0] - (int)ref[0]);
        sad += abs((int)src[1] - (int)ref[1]);
        sad += abs((int)src[2] - (int)ref[2]);
        sad += abs((int)src[3] - (int)ref[3]);
        sad += abs((int)src[4] - (int)ref[4]);
        sad += abs((int)src[5] - (int)ref[5]);
        sad += abs((int)src[6] - (int)ref[6]);
        sad += abs((int)src[7] - (int)ref[7]);
        sad += abs((int)src[8] - (int)ref[8]);
        sad += abs((int)src[9] - (int)ref[9]);
        sad += abs((int)src[10] - (int)ref[10]);
        sad += abs((int)src[11] - (int)ref[11]);
        sad += abs((int)src[12] - (int)ref[12]);
        sad += abs((int)src[13] - (int)ref[13]);
        sad += abs((int)src[14] - (int)ref[14]);
        sad += abs((int)src[15] - (int)ref[15]);
        sad += abs((int)src[16] - (int)ref[16]);
        sad += abs((int)src[17] - (int)ref[17]);
        sad += abs((int)src[18] - (int)ref[18]);
        sad += abs((int)src[19] - (int)ref[19]);
        sad += abs((int)src[20] - (int)ref[20]);
        sad += abs((int)src[21] - (int)ref[21]);
        sad += abs((int)src[22] - (int)ref[22]);
        sad += abs((int)src[23] - (int)ref[23]);
        sad += abs((int)src[24] - (int)ref[24]);
        sad += abs((int)src[25] - (int)ref[25]);
        sad += abs((int)src[26] - (int)ref[26]);
        sad += abs((int)src[27] - (int)ref[27]);
        sad += abs((int)src[28] - (int)ref[28]);
        sad += abs((int)src[29] - (int)ref[29]);
        sad += abs((int)src[30] - (int)ref[30]);
        sad += abs((int)src[31] - (int)ref[31]);
        sad += abs((int)src[32] - (int)ref[32]);
        sad += abs((int)src[33] - (int)ref[33]);
        sad += abs((int)src[34] - (int)ref[34]);
        sad += abs((int)src[35] - (int)ref[35]);
        sad += abs((int)src[36] - (int)ref[36]);
        sad += abs((int)src[37] - (int)ref[37]);
        sad += abs((int)src[38] - (int)ref[38]);
        sad += abs((int)src[39] - (int)ref[39]);
        sad += abs((int)src[40] - (int)ref[40]);
        sad += abs((int)src[41] - (int)ref[41]);
        sad += abs((int)src[42] - (int)ref[42]);
        sad += abs((int)src[43] - (int)ref[43]);
        sad += abs((int)src[44] - (int)ref[44]);
        sad += abs((int)src[45] - (int)ref[45]);
        sad += abs((int)src[46] - (int)ref[46]);
        sad += abs((int)src[47] - (int)ref[47]);
        sad += abs((int)src[48] - (int)ref[48]);
        sad += abs((int)src[49] - (int)ref[49]);
        sad += abs((int)src[50] - (int)ref[50]);
        sad += abs((int)src[51] - (int)ref[51]);
        sad += abs((int)src[52] - (int)ref[52]);
        sad += abs((int)src[53] - (int)ref[53]);
        sad += abs((int)src[54] - (int)ref[54]);
        sad += abs((int)src[55] - (int)ref[55]);
        sad += abs((int)src[56] - (int)ref[56]);
        sad += abs((int)src[57] - (int)ref[57]);
        sad += abs((int)src[58] - (int)ref[58]);
        sad += abs((int)src[59] - (int)ref[59]);
        sad += abs((int)src[60] - (int)ref[60]);
        sad += abs((int)src[61] - (int)ref[61]);
        sad += abs((int)src[62] - (int)ref[62]);
        sad += abs((int)src[63] - (int)ref[63]);
        src += stride_src;
        ref += stride_ref;
    }
    if (sizeof(Sample) == 2)
        sad >>= 2;
    return sad;
}


template <typename Sample>
static int havoc_sad_c_ref(const Sample *src, intptr_t stride_src, const Sample *ref, intptr_t stride_ref, uint32_t rect)
{
    const int width = rect >> 8;
    const int height = rect & 0xff;

    int sad = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            sad += abs((int)src[x + y * stride_src] - (int)ref[x + y * stride_ref]);
        }
    }
    if (sizeof(Sample) == 2)
        sad >>= 2;
    return sad;
}


template <typename Sample>
havoc_sad<Sample>* get_sad(int width, int height, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    if ((sizeof(Sample) == 1 && buffer.isa & HAVOC_SSE2) ||
        (sizeof(Sample) == 2 && buffer.isa & HAVOC_SSSE3))
    {
#define X(w, h) \
        { \
            SadSse2<Sample> sadSse2(&buffer, w, h); \
            if (w==width && h==height) \
            { \
                havoc_sad<Sample> *f = sadSse2; \
                if (f) return f; \
            } \
        }

        X_HEVC_PU_SIZES;
#undef X
    }

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
    {
        return (havoc_sad<Sample>*)&havoc_sad_c_ref<Sample>;
    }

#if USE_HM_DERIVED
    if (buffer.isa & HAVOC_C_REF)
    {
        if (width == 4) return (havoc_sad<Sample>*)&havoc_sad_c_opt_4;
        if (width == 8) return (havoc_sad<Sample>*)&havoc_sad_c_opt_8;
        if (width == 16) return (havoc_sad<Sample>*)&havoc_sad_c_opt_16;
        if (width == 32) return (havoc_sad<Sample>*)&havoc_sad_c_opt_32;
        if (width == 64) return (havoc_sad<Sample>*)&havoc_sad_c_opt_64;
    }
#endif

    return 0;
}


template <typename Sample>
void havoc_populate_sad(havoc_table_sad<Sample> *table, havoc_code code)
{
    for (int height = 4; height <= 64; height += 4)
    {
        for (int width = 4; width <= 64; width += 4)
        {
            *havoc_get_sad(table, width, height) = get_sad<Sample>(width, height, code);
        }
    }
}


template void havoc_populate_sad<uint8_t>(havoc_table_sad<uint8_t> *table, havoc_code code);


template void havoc_populate_sad<uint16_t>(havoc_table_sad<uint16_t> *table, havoc_code code);


template <typename Sample>
static void havoc_sad_multiref_4_c_ref(const Sample *src, intptr_t stride_src, const Sample *ref[], intptr_t stride_ref, int sad[], uint32_t rect)
{
    const int width = rect >> 8;
    const int height = rect & 0xff;

    sad[0] = 0;
    sad[1] = 0;
    sad[2] = 0;
    sad[3] = 0;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            for (int way = 0; way < 4; ++way)
            {
                sad[way] += abs((int)src[x + y * stride_src] - (int)ref[way][x + y * stride_ref]);
            }
        }
    }

    if (sizeof(Sample) == 2)
    {
        sad[0] >>= 2;
        sad[1] >>= 2;
        sad[2] >>= 2;
        sad[3] >>= 2;
    }
}


#if USE_WEBM_DERIVED

template <typename Sample>
struct Sad4Avx2 :
    Jit::Function
{
    Sad4Avx2(Jit::Buffer *buffer, int width, int height) :
        Jit::Function(buffer, Jit::CountArguments<havoc_sad_multiref<uint8_t>>::value),
        width(width),
        height(height)
    {
        this->build();
    }

    int width, height;

    Xbyak::Label mask_24_32;
    Xbyak::Label mask_12_16;

    void data() override
    {
        align();

        int widthBytes = width * sizeof(Sample);

        if (widthBytes == 24)
        {
            L(mask_24_32);
            db({ 0xff }, 24);
            db({ 0 }, 8);
        }

        if (widthBytes == 12)
        {
            L(mask_12_16);
            db({ 0xff }, 12);
            db({ 0 }, 4);
            db({ 0xff }, 12);
            db({ 0 }, 4);
        }
    }

    template <class A, class B>
    void sad(A const &a, B const &b)
    {
        if (sizeof(Sample) == 2)
        {
            vpsubw(a, b);
            vpabsw(a, a);
            vphaddw(a, ymm9);
            vpunpcklwd(a, ymm9);
            vphaddd(a, ymm9);
            vpunpckldq(a, ymm9);
        }
        else
            vpsadbw(a, b);
    }

    void assemble() override
    {
        auto &src = arg64(0);
        auto &src_stride = arg64(1);
        auto &ref = arg64(2);
        auto &ref_stride = arg64(3);
        auto &sads = arg64(4);
        auto &rect = arg64(5);

        if (sizeof(Sample) == 2)
        {
            shl(src_stride, 1);
            shl(ref_stride, 1);
        }

        auto &xmm0 = regXmm(0);
        auto &xmm1 = regXmm(1);
        auto &xmm2 = regXmm(2);
        auto &xmm3 = regXmm(3);
        auto &xmm4 = regXmm(4);
        auto &xmm5 = regXmm(5);
        auto &xmm6 = regXmm(6);
        auto &xmm7 = regXmm(7);
        auto &xmm8 = regXmm(8);
        auto *xmm9 = sizeof(Sample) == 2 ? &regXmm(9) : (Xbyak::Xmm const *)0;

        if (sizeof(Sample) == 2)
            vpxor(ymm9, ymm9);

        int widthBytes = width * sizeof(Sample);

        if (widthBytes == 8)
        {
            auto &ref0 = reg64(6);
            auto &ref1 = reg64(7);
            auto &ref2 = reg64(8);
            auto &ref3 = ref;

            mov(ref0, ptr[ref]);
            mov(ref1, ptr[ref + 0x8]);
            mov(ref2, ptr[ref + 0x10]);
            mov(ref3, ptr[ref + 0x18]);

            and (rect, 0xff);
            auto &n = rect;
            shr(n, 1);

            vmovq(xmm4, ptr[src]);
            vmovq(xmm0, ptr[ref0]);
            vmovq(xmm1, ptr[ref1]);
            vmovq(xmm2, ptr[ref2]);
            vmovq(xmm3, ptr[ref3]);

            sub(ref1, ref0);
            sub(ref2, ref0);
            sub(ref3, ref0);
            lea(ref0, ptr[ref0 + ref_stride]);

            vmovhps(xmm4, ptr[src + src_stride]);
            vmovhps(xmm0, ptr[ref0]);
            vmovhps(xmm1, ptr[ref1 + ref0]);
            vmovhps(xmm2, ptr[ref2 + ref0]);
            vmovhps(xmm3, ptr[ref3 + ref0]);

            lea(ref0, ptr[ref0 + ref_stride]);
            lea(src, ptr[src + src_stride * 2]);

            sad(xmm0, xmm4);
            sad(xmm1, xmm4);
            sad(xmm2, xmm4);
            sad(xmm3, xmm4);

            dec(n);

            L("loop");
            {
                vmovq(xmm4, ptr[src]);
                vmovq(xmm5, ptr[ref0]);
                vmovq(xmm6, ptr[ref1 + ref0]);
                vmovq(xmm7, ptr[ref2 + ref0]);
                vmovq(xmm8, ptr[ref3 + ref0]);

                lea(ref0, ptr[ref0 + ref_stride]);

                vmovhps(xmm4, ptr[src + src_stride]);
                vmovhps(xmm5, ptr[ref0]);
                vmovhps(xmm6, ptr[ref1 + ref0]);
                vmovhps(xmm7, ptr[ref2 + ref0]);
                vmovhps(xmm8, ptr[ref3 + ref0]);

                lea(ref0, ptr[ref0 + ref_stride]);
                lea(src, ptr[src + src_stride * 2]);

                sad(xmm5, xmm4);
                sad(xmm6, xmm4);
                sad(xmm7, xmm4);
                sad(xmm8, xmm4);

                vpaddd(xmm0, xmm5);
                vpaddd(xmm1, xmm6);
                vpaddd(xmm2, xmm7);
                vpaddd(xmm3, xmm8);
            }
            dec(n);
            jg("loop");

            vpslldq(xmm1, 4);
            vpslldq(xmm3, 4);
            vpor(xmm0, xmm1);
            vpor(xmm2, xmm3);
            vmovdqa(xmm1, xmm0);
            vmovdqa(xmm3, xmm2);
            vpunpcklqdq(xmm0, xmm2);
            vpunpckhqdq(xmm1, xmm3);
            vpaddd(xmm0, xmm1);
            if (sizeof(Sample) == 2)
                vpsrld(xmm0, 2);
            vmovdqu(ptr[sads], xmm0);
        }
        else if (widthBytes == 4)
        {
            auto &ref0 = reg64(6);
            auto &ref1 = reg64(7);
            auto &ref2 = reg64(8);
            auto &ref3 = ref;

            mov(ref0, ptr[ref]);
            mov(ref1, ptr[ref + 0x8]);
            mov(ref2, ptr[ref + 0x10]);
            mov(ref3, ptr[ref + 0x18]);

            vpxor(xmm0, xmm0);
            vpxor(xmm1, xmm1);
            vpxor(xmm2, xmm2);
            vpxor(xmm3, xmm3);

            sub(ref1, ref0);
            sub(ref2, ref0);
            sub(ref3, ref0);

            and (rect, 0xFF);
            auto &n = rect;
            shr(n, 1);

            L("loop");
            {
                vmovd(xmm4, ptr[src]);
                vmovd(xmm5, ptr[ref0]);
                vmovd(xmm6, ptr[ref1 + ref0]);
                vmovd(xmm7, ptr[ref2 + ref0]);
                vmovd(xmm8, ptr[ref3 + ref0]);
                lea(ref0, ptr[ref0 + ref_stride]);
                vpunpckldq(xmm4, ptr[src + src_stride]);
                vpunpckldq(xmm5, ptr[ref0]);
                vpunpckldq(xmm6, ptr[ref1 + ref0]);
                vpunpckldq(xmm7, ptr[ref2 + ref0]);
                vpunpckldq(xmm8, ptr[ref3 + ref0]);
                lea(ref0, ptr[ref0 + ref_stride]);
                lea(src, ptr[src + src_stride * 2]);
                sad(xmm5, xmm4);
                sad(xmm6, xmm4);
                sad(xmm7, xmm4);
                sad(xmm8, xmm4);
                vpunpckldq(xmm5, xmm6);
                vpunpckldq(xmm7, xmm8);
                vpaddd(xmm0, xmm5);
                vpaddd(xmm2, xmm7);
            }
            dec(n);
            jg("loop");

            vpslldq(xmm1, xmm1, 4);
            vpslldq(xmm3, xmm3, 4);
            vpor(xmm0, xmm1);
            vpor(xmm2, xmm3);
            vmovdqa(xmm1, xmm0);
            vmovdqa(xmm3, xmm2);
            vpunpcklqdq(xmm0, xmm2);
            vpunpckhqdq(xmm1, xmm3);
            vpaddd(xmm0, xmm1);
            if (sizeof(Sample) == 2)
                vpsrld(xmm0, 2);
            vmovdqu(ptr[sads], xmm0);
        }
        else/* if (widthBytes == 32 || widthBytes == 48 || widthBytes == 64)*/
        {
            auto &r0 = arg64(0);
            auto &r1 = arg64(1);
            auto &r2 = arg64(2);
            auto &r3 = arg64(3);
            auto &r4 = arg64(4);
            auto &r5 = arg64(5);

            auto &r6 = reg64(6);
            auto &r7 = reg64(7);
            auto &r8 = reg64(8);

            mov(r6, r1); //stride_src
            mov(r7, r3); // stride_ref

            and (rect, 0xFF);
            auto &n = rect;

            if (widthBytes <= 16)
            {
                shr(n, 1);
            }

            mov(r8, r4); // sad[]

            mov(r1, ptr[ref + 0 * 8]);
            mov(r3, ptr[ref + 2 * 8]);
            mov(r4, ptr[ref + 3 * 8]);
            mov(r2, ptr[ref + 1 * 8]);

            vpxor(ymm5, ymm5);
            vpxor(ymm6, ymm6);
            vpxor(ymm7, ymm7);
            vpxor(ymm8, ymm8);

            L("loop");
            {
                if (widthBytes <= 16)
                {
                    vmovdqu(xmm0, ptr[r0]);
                    vmovdqu(xmm1, ptr[r1]);
                    vmovdqu(xmm2, ptr[r2]);
                    vmovdqu(xmm3, ptr[r3]);
                    vmovdqu(xmm4, ptr[r4]);

                    vinserti128(ymm0, ymm0, ptr[r0 + r6], 1);
                    vinserti128(ymm1, ymm1, ptr[r1 + r7], 1);
                    vinserti128(ymm2, ymm2, ptr[r2 + r7], 1);
                    vinserti128(ymm3, ymm3, ptr[r3 + r7], 1);
                    vinserti128(ymm4, ymm4, ptr[r4 + r7], 1);

                    if (widthBytes == 12)
                    {
                        vpand(ymm0, ptr[rip + mask_12_16]);
                        vpand(ymm1, ptr[rip + mask_12_16]);
                        vpand(ymm2, ptr[rip + mask_12_16]);
                        vpand(ymm3, ptr[rip + mask_12_16]);
                        vpand(ymm4, ptr[rip + mask_12_16]);
                    }
                }
                else
                {
                    vmovdqu(ymm0, ptr[r0]);
                    vmovdqu(ymm1, ptr[r1]);
                    vmovdqu(ymm2, ptr[r2]);
                    vmovdqu(ymm3, ptr[r3]);
                    vmovdqu(ymm4, ptr[r4]);
                }

                sad(ymm1, ymm0);
                sad(ymm2, ymm0);
                sad(ymm3, ymm0);
                sad(ymm4, ymm0);

                if (widthBytes == 24)
                {
                    vpand(ymm1, ptr[rip + mask_24_32]);
                    vpand(ymm2, ptr[rip + mask_24_32]);
                    vpand(ymm3, ptr[rip + mask_24_32]);
                    vpand(ymm4, ptr[rip + mask_24_32]);
                }

                vpaddd(ymm5, ymm1);
                vpaddd(ymm6, ymm2);
                vpaddd(ymm7, ymm3);
                vpaddd(ymm8, ymm4);

                for (int i = 1; i < (widthBytes + 31) / 32; ++i)
                {
                    vmovdqu(ymm0, ptr[r0 + 32 * i]);
                    vmovdqu(ymm1, ptr[r1 + 32 * i]);
                    vmovdqu(ymm2, ptr[r2 + 32 * i]);
                    vmovdqu(ymm3, ptr[r3 + 32 * i]);
                    vmovdqu(ymm4, ptr[r4 + 32 * i]);

                    if (widthBytes == 48)
                    {
                        sad(xmm1, xmm0);
                        sad(xmm2, xmm0);
                        sad(xmm3, xmm0);
                        sad(xmm4, xmm0);
                    }
                    else
                    {
                        sad(ymm1, ymm0);
                        sad(ymm2, ymm0);
                        sad(ymm3, ymm0);
                        sad(ymm4, ymm0);
                    }

                    vpaddd(ymm5, ymm1);
                    vpaddd(ymm6, ymm2);
                    vpaddd(ymm7, ymm3);
                    vpaddd(ymm8, ymm4);
                }

                if (widthBytes <= 16)
                {
                    lea(r0, ptr[r0 + r6 * 2]);
                    lea(r1, ptr[r1 + r7 * 2]);
                    lea(r2, ptr[r2 + r7 * 2]);
                    lea(r3, ptr[r3 + r7 * 2]);
                    lea(r4, ptr[r4 + r7 * 2]);
                }
                else
                {
                    lea(r0, ptr[r0 + r6]);
                    lea(r1, ptr[r1 + r7]);
                    lea(r2, ptr[r2 + r7]);
                    lea(r3, ptr[r3 + r7]);
                    lea(r4, ptr[r4 + r7]);
                }
            }
            dec(n);
            jg("loop");

            vextracti128(xmm0, ymm5, 1);
            vextracti128(xmm1, ymm6, 1);
            vextracti128(xmm2, ymm7, 1);
            vextracti128(xmm3, ymm8, 1);

            vpaddd(xmm5, xmm0);
            vpaddd(xmm6, xmm1);
            vpaddd(xmm7, xmm2);
            vpaddd(xmm8, xmm3);

            //; two different ways of achieving the same thing - both seem to take similar number of cycles
            if (0)
            {
                //	pshufd xm0, xm5, ORDER(0, 0, 0, 2)
                //	pshufd xm1, xm6, ORDER(0, 0, 0, 2)
                //	pshufd xm2, xm7, ORDER(0, 0, 0, 2)
                //	pshufd xm3, xm8, ORDER(0, 0, 0, 2)
                //	paddd xm5, xm0
                //	paddd xm6, xm1
                //	paddd xm7, xm2
                //	paddd xm8, xm3
                //	movd[r8 + 0 * 4], xm5
                //	movd[r8 + 1 * 4], xm6
                //	movd[r8 + 2 * 4], xm7
                //	movd[r8 + 3 * 4], xm8
            }
            else
            {
                vpslldq(xmm6, 4);
                vpslldq(xmm8, 4);
                vpor(xmm5, xmm6);
                vpor(xmm7, xmm8);
                vmovdqa(xmm6, xmm5);
                vmovdqa(xmm8, xmm7);
                vpunpcklqdq(xmm5, xmm7);
                vpunpckhqdq(xmm6, xmm8);
                vpaddd(xmm5, xmm6);
                if (sizeof(Sample) == 2)
                    vpsrld(xmm5, 2);
                vmovdqu(ptr[r8], xmm5);
            }
        }
    }
};

#endif


template <typename Sample>
havoc_sad_multiref<Sample>* get_sad_multiref(int ways, int width, int height, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
    havoc_sad_multiref<Sample>* f = 0;

    if (ways != 4) return 0;

#if USE_WEBM_DERIVED
    if (buffer.isa & HAVOC_AVX2)
    {
#define X(w, h) \
        { \
            Sad4Avx2<Sample> sad4Avx2(&buffer, w, h); \
            if (w==width && h==height) \
            { \
                havoc_sad_multiref<Sample> *f = sad4Avx2; \
                if (f) return f; \
            } \
        }

        X_HEVC_PU_SIZES;
#undef X
    }
#endif

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
    {
        if (!f) f = &havoc_sad_multiref_4_c_ref<Sample>;
    }

    return f;
}

template <typename Sample>
void havoc_populate_sad_multiref(havoc_table_sad_multiref<Sample> *table, havoc_code code)
{
    for (int height = 4; height <= 64; height += 4)
    {
        for (int width = 4; width <= 64; width += 4)
        {
            *havoc_get_sad_multiref(table, 4, width, height) = get_sad_multiref<Sample>(4, width, height, code);
        }
    }

}

template void havoc_populate_sad_multiref<uint8_t>(havoc_table_sad_multiref<uint8_t> *table, havoc_code code);
template void havoc_populate_sad_multiref<uint16_t>(havoc_table_sad_multiref<uint16_t> *table, havoc_code code);


struct BoundSadBase
{
    int width;
    int height;
    int sad;
    int bitDepth;
    virtual int init(havoc_code code) = 0;
};


template <typename Sample>
struct BoundSad :
    BoundSadBase
{
    HAVOC_ALIGN(32, Sample, src[128 * 128]);
    HAVOC_ALIGN(32, Sample, ref[128 * 128]);
    havoc_sad<Sample> *f;
    int init(havoc_code code) override
    {
        this->bitDepth = sizeof(Sample) == 2 ? 10 : 8;

        auto s = this;

        havoc_table_sad<Sample> table;
        havoc_populate_sad<Sample>(&table, code);

        s->f = *havoc_get_sad<Sample>(&table, s->width, s->height);

        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
        if (buffer.isa == HAVOC_C_REF) printf("\t%d bits %dx%d:", s->bitDepth, s->width, s->height);

        return !!s->f;
    }
};

int init_sad(void *p, havoc_code code)
{
    BoundSadBase *s = (BoundSadBase *)p;
    return s->init(code);
}


typedef BoundSad<uint8_t> bound_sad;


template <typename Sample>
void invokeSad(void *p, int n)
{
    BoundSad<Sample> *s = (BoundSad<Sample> *)p;
    const Sample *unaligned_ref = &s->ref[1 + 1 * 128];
    while (n--)
    {
        s->sad = s->f(s->src, 64, unaligned_ref, 64, HAVOC_RECT(s->width, s->height));
    }
}


void invoke_sad(void *p, int n)
{
    BoundSadBase *s = (BoundSadBase *)p;

    if (s->bitDepth == 8)
        invokeSad<uint8_t>(s, n);
    else
        invokeSad<uint16_t>(s, n);
}


int mismatch_sad(void *boundRef, void *boundTest)
{
    BoundSadBase *ref = (BoundSadBase *)boundRef;
    BoundSadBase *test = (BoundSadBase *)boundTest;

    return  ref->sad != test->sad;
}


static const int partitions[][2] = {
    { 64, 64 },{ 64, 48 },{ 64, 32 },{ 64, 16 },
    { 48, 64 },
    { 32, 64 },{ 32, 32 },{ 32, 24 },{ 32, 16 },{ 32, 8 },
    { 24, 32 },
    { 16, 64 },{ 16, 32 },{ 16, 16 },{ 16, 12 },{ 16, 8 },{ 16, 4 },
    { 12, 16 },
    { 8, 32 },{ 8, 16 },{ 8, 8 },{ 8, 4 },
    { 4, 8 },
    { 0, 0 } };


template <typename Sample>
void testSad(int *error_count, havoc_instruction_set mask)
{
    BoundSad<Sample> b[2];

    for (int x = 0; x < 128 * 128; x++) b[0].src[x] = rand() & 0x3ff;
    for (int x = 0; x < 128 * 128; x++) b[0].ref[x] = rand() & 0x3ff;

    for (int i = 0; partitions[i][0]; ++i)
    {
        b[0].width = partitions[i][0];
        b[0].height = partitions[i][1];
        b[1] = b[0];
        *error_count += havoc_test(&b[0], &b[1], init_sad, invoke_sad, mismatch_sad, mask, 10);
    }
}


void havoc_test_sad(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_sad - Sum of Absolute Differences\n");

    //	testSad<uint8_t>(error_count, mask);
    testSad<uint16_t>(error_count, mask);
}


struct BoundsSadMultirefBase
{
    int bits;
    int ways;
    int width;
    int height;
    int sad[4];
    virtual int init(havoc_code code) = 0;
};


template <typename Sample>
struct BoundSadMultiref :
    BoundsSadMultirefBase
{
    havoc_sad_multiref<Sample> *f;
    HAVOC_ALIGN(32, Sample, src[128 * 128]);
    HAVOC_ALIGN(32, Sample, ref[128 * 128]);
    const Sample *ref_array[4];

    int init(havoc_code code) override
    {
        this->bits = sizeof(Sample) == 2 ? 10 : 8;

        auto s = this;

        havoc_table_sad_multiref<Sample> table;
        havoc_populate_sad_multiref(&table, code);
        s->f = *havoc_get_sad_multiref(&table, s->ways, s->width, s->height);

        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
        if (s->f && buffer.isa == HAVOC_C_REF)
        {
            printf("\t%d bits %d-way %dx%d : ", s->bits, s->ways, s->width, s->height);
        }

        return !!s->f;
    }
};


int init_sad_multiref(void *p, havoc_code code)
{
    BoundsSadMultirefBase *s = (BoundsSadMultirefBase *)p;
    return s->init(code);
}


template <typename Sample>
void invokeSadMultiref(void *p, int n)
{
    BoundSadMultiref<Sample> *s = (BoundSadMultiref<Sample> *)p;

    while (n--)
    {
        s->f(s->src, 64, s->ref_array, 64, s->sad, HAVOC_RECT(s->width, s->height));
    }
}


void invoke_sad_multiref(void *p, int n)
{
    BoundsSadMultirefBase *s = (BoundsSadMultirefBase *)p;

    s->sad[0] = 0;
    s->sad[1] = 0;
    s->sad[2] = 0;
    s->sad[3] = 0;

    if (s->bits == 8)
        invokeSadMultiref<uint8_t>(s, n);
    else
        invokeSadMultiref<uint16_t>(s, n);
}


int mismatch_sad_multiref(void *boundRef, void *boundTest)
{
    BoundsSadMultirefBase *ref = (BoundsSadMultirefBase *)boundRef;
    BoundsSadMultirefBase *test = (BoundsSadMultirefBase *)boundTest;

    assert(ref->ways);

    for (int i = 0; i < ref->ways; ++i)
    {
        if (ref->sad[i] != test->sad[i]) return 1;
    }

    return 0;
}


template <typename Sample>
void testSadMultiref(int *error_count, havoc_instruction_set mask)
{
    BoundSadMultiref<Sample> b[2];

    b[0].ways = 4;

    for (int x = 0; x < 128 * 128; x++) b[0].src[x] = rand() & 0x3ff;
    for (int x = 0; x < 128 * 128; x++) b[0].ref[x] = rand() & 0x3ff;

    b[0].ref_array[0] = &b[0].ref[1 + 2 * 128];
    b[0].ref_array[1] = &b[0].ref[2 + 1 * 128];
    b[0].ref_array[2] = &b[0].ref[3 + 2 * 128];
    b[0].ref_array[3] = &b[0].ref[2 + 3 * 128];

    for (int i = 0; partitions[i][0]; ++i)
    {
        b[0].width = partitions[i][0];
        b[0].height = partitions[i][1];
        b[1] = b[0];
        *error_count += havoc_test(&b[0], &b[1], init_sad_multiref, invoke_sad_multiref, mismatch_sad_multiref, mask, 10);
    }
}


void havoc_test_sad_multiref(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_sad_multiref - Sum Of Absolute Differences with multiple references (%d candidate references)\n", 4);
    testSadMultiref<uint16_t>(error_count, mask);
    testSadMultiref<uint8_t>(error_count, mask);
}
