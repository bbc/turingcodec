/*
Copyright(C) 2016 British Broadcasting Corporation, Parabola Research
and Queen Mary University of London.

This file is part of the Turing codec.

The Turing codec is free software; you can redistribute it and / or modify
it under the terms of version 2 of the GNU General Public License as
published by the Free Software Foundation.

The Turing codec is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

Commercial support and intellectual property rights for
the Turing codec are also available under a proprietary license.
For more information, contact us at info @ turingcodec.org.
 */


#include "pred_inter.h"
#include "havoc_test.h"
#include "Jit.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <type_traits>


static int Clip3(int min, int max, int value)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}


int havoc_pred_coefficient(int n, int fractionalPosition, int k)
{
    if (n == 8)
    {
        static const int kernel[4][8] =
        {
            { 0, 0, 0, 64, 0, 0, 0, 0 },
            { -1, 4, -10, 58, 17, -5, 1, 0 },
            { -1, 4, -11, 40, 40, -11, 4, -1 },
            { 0, 1, -5, 17, 58, -10, 4, -1 }
        };

        return kernel[fractionalPosition][k];
    }
    else
    {
        static const int kernel[8][4] =
        {
            { 0, 64, 0, 0 },
            { -2, 58, 10, -2 },
            { -4, 54, 16, -2 },
            { -6, 46, 28, -4 },
            { -4, 36, 36, -4 },
            { -4, 28, 46, -6 },
            { -2, 16, 54, -4 },
            { -2, 10, 58, -2 },
        };

        return kernel[fractionalPosition][k];
    }
};


/*
Generic function to handle all luma, chroma, 8 and 16-bit interpolation.
Consider this a reference implementation: it will not run fast.
 */
template <typename Dst, typename Src>
static void havoc_pred_uni_generic(
    Dst *dst, intptr_t stride_dst,
    Src const *src, intptr_t stride_src,
    int w, int h,
    intptr_t stride_tap,
    int n, int fractionalPosition, int shift, int add, int bitDepth)
{
    while (h--)
    {
        for (int x = 0; x < w; ++x)
        {
            int a = add << shift >> 1;

            for (int k = 0; k < n; ++k)
            {
                const intptr_t src_offset = x + (k - n / 2 + 1) * stride_tap;

                a += havoc_pred_coefficient(n, fractionalPosition, k) * src[src_offset];
            }

            a >>= shift;

            if (bitDepth)
            {
                a = Clip3(0, (1 << bitDepth) - 1, a);
            }

            dst[x] = a;
        }

        dst += stride_dst;
        src += stride_src;
    }
};


template <typename Sample>
void havoc_pred_uni_copy_block(Sample *dst, intptr_t stride_dst, const Sample *ref, intptr_t stride_ref, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    assert(!xFrac);
    assert(!yFrac);
    while (h--)
    {
        memcpy(dst, ref, w * sizeof(Sample));
        dst += stride_dst;
        ref += stride_ref;
    }
}


template <typename Sample>
void havoc_pred_uni_8tap_h(Sample *dst, intptr_t stride_dst, const Sample *ref, intptr_t stride_ref, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    assert(xFrac);
    assert(!yFrac);
    havoc_pred_uni_generic(dst, stride_dst, ref, stride_ref, w, h, 1, 8, xFrac, 6, 1, bitDepth);
}


template <typename Sample>
void havoc_pred_uni_8tap_v(Sample *dst, intptr_t stride_dst, const Sample *ref, intptr_t stride_ref, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    assert(!xFrac);
    assert(yFrac);
    havoc_pred_uni_generic(dst, stride_dst, ref, stride_ref, w, h, stride_ref, 8, yFrac, 6, 1, bitDepth);
}


template <typename Sample>
void havoc_pred_uni_8tap_hv(Sample *dst, intptr_t stride_dst, const Sample *ref, intptr_t stride_ref, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    int intermediate[(64 + 7) * 64];

    int shift1 = bitDepth - 8;
    if (shift1 > 4) shift1 = 4;

    int shift2 = 6;

    int shift3 = 14 - bitDepth;
    if (shift3 < 2) shift3 = 2;

    /* Horizontal filter */
    havoc_pred_uni_generic(intermediate, 64, ref - 3 * stride_ref, stride_ref, w, h + 7, 1, 8, xFrac, shift1, 0, 0);

    /* Vertical filter */
    havoc_pred_uni_generic(dst, stride_dst, intermediate + 3 * 64, 64, w, h, 64, 8, yFrac, shift2 + shift3, 1, bitDepth);
}


template <typename Sample>
void havoc_pred_uni_4tap_h(Sample *dst, intptr_t stride_dst, const Sample *ref, intptr_t stride_ref, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    assert(xFrac);
    assert(!yFrac);
    havoc_pred_uni_generic(dst, stride_dst, ref, stride_ref, w, h, 1, 4, xFrac, 6, 1, bitDepth);
}


template <typename Sample>
void havoc_pred_uni_4tap_v(Sample *dst, intptr_t stride_dst, const Sample *ref, intptr_t stride_ref, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    assert(!xFrac);
    assert(yFrac);
    havoc_pred_uni_generic(dst, stride_dst, ref, stride_ref, w, h, stride_ref, 4, yFrac, 6, 1, bitDepth);
}


template <typename Sample>
void havoc_pred_uni_4tap_hv(Sample *dst, intptr_t stride_dst, const Sample *ref, intptr_t stride_ref, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    int intermediate[(64 + 3) * 64];

    int shift1 = bitDepth - 8;
    if (shift1 > 4) shift1 = 4;

    int shift2 = 6;

    int shift3 = 14 - bitDepth;
    if (shift3 < 2) shift3 = 2;

    /* Horizontal filter */
    havoc_pred_uni_generic(intermediate, 64, ref - stride_ref, stride_ref, w, h + 3, 1, 4, xFrac, shift1, 0, 0);

    /* Vertical filter */
    havoc_pred_uni_generic(dst, stride_dst, intermediate + 64, 64, w, h, 64, 4, yFrac, shift2 + shift3, 1, bitDepth);
}


template <typename Sample>
struct PredUniCopy
    :
    Jit::Function
{
    int width;

    PredUniCopy(Jit::Buffer *buffer, int width)
        :
        Jit::Function(buffer, Jit::CountArguments<typename std::remove_pointer<HavocPredUni<Sample>>::type>::value),
        width(width)
    {
        this->build();
    }

    void assemble() override
    {
        auto &r0 = arg64(0);
        auto &r1 = arg64(1);
        auto &r2 = arg64(2);
        auto &r3 = arg64(3);
        auto height = Xbyak::Reg32(arg64(5).getIdx());

        L("loop");
        {
            int const n = sizeof(Sample) * width / 16;

            Xbyak::Xmm const *m[8];
            for (int i = 0; i < n; ++i)
            {
                m[i] = &regXmm(i);
                movdqu(*m[i], ptr[r2 + 16 * i]);
            }

            for (int i = 0; i < n; ++i)
            {
                movdqu(ptr[r0 + 16 * i], *m[i]);
            }

            lea(r2, ptr[r2 + r3 * sizeof(Sample)]);
            lea(r0, ptr[r0 + r1 * sizeof(Sample)]);
        }
        dec(height);
        jg("loop");
    }
};

typedef void havoc_pred_uni_8to16(int16_t *dst, intptr_t stride_dst, const uint8_t *ref, intptr_t stride_ref, int nPbW, int nPbH, int xFrac, int yFrac);
typedef void havoc_pred_uni_16to8(uint8_t *dst, intptr_t stride_dst, const int16_t *ref, intptr_t stride_ref, int nPbW, int nPbH, int xFrac, int yFrac);
typedef void havoc_pred_bi_v_16to16(uint8_t *dst, intptr_t stride_dst, const int16_t *refAtop, const int16_t *refBtop, intptr_t stride_ref, int nPbW, int nPbH, int yFracA, int yFracB);
typedef void havoc_pred_bi_copy(uint8_t *dst, intptr_t stride_dst, const uint8_t *ref0, const uint8_t *ref1, intptr_t stride_ref, int nPbW, int nPbH, int xFrac0, int yFrac0, int xFrac1, int yFrac1, int bitDepth);


template <typename Sample>
struct PredInter
    :
    Jit::Function
{
    PredInter(Jit::Buffer *buffer, int params, int taps, int width, int xFrac, int yFrac, int inputBitDepth, int inputSize)
        :
        Jit::Function(buffer, params),
        taps(taps),
        width(width),
        xFrac(xFrac),
        yFrac(yFrac),
        inputBitDepth(inputBitDepth),
        inputSize(inputSize),
        refs(1)
    {
        this->build();
    }

    PredInter(Jit::Buffer *buffer, int taps, int width, int inputBitDepth)
        :
        Jit::Function(buffer, Jit::CountArguments<havoc_pred_bi_v_16to16>::value),
        taps(taps),
        width(width),
        xFrac(xFrac),
        yFrac(yFrac),
        inputBitDepth(inputBitDepth),
        refs(2)
    {
        assert(0);
        this->build();
    }

    int refs;
    int taps;
    int width;
    int xFrac;
    int yFrac;
    int inputBitDepth;
    int inputSize;

    Xbyak::Label coefficients;
    Xbyak::Label coefficients16;
    Xbyak::Label constant_times_4_dd_0x800;
    Xbyak::Label constant_times_4_dd_0x200;
    Xbyak::Label constant_times_4_dd_0x20;
    Xbyak::Label constant_times_8_dw_0x20;
    Xbyak::Label constant_times_8_dw_0x40;

    void data() override
    {
        align();

        L(coefficients);
        for (int frac = 0; frac < (12 - taps); ++frac)
        {
            for (int k = 0; k < taps; k += 2)
            {
                int coeff0 = havoc_pred_coefficient(taps, frac, k);
                int coeff1 = havoc_pred_coefficient(taps, frac, k + 1);

                db({ coeff0, coeff1 }, 8);
            }
        }

        L(coefficients16);
        for (int frac = 0; frac < (12 - taps); ++frac)
        {
            for (int k = 0; k < taps; k += 2)
            {
                int coeff0 = havoc_pred_coefficient(taps, frac, k);
                int coeff1 = havoc_pred_coefficient(taps, frac, k + 1);

                dw({ coeff0, coeff1 }, 4);
            }
        }

        L(constant_times_4_dd_0x800);
        dd({ 0x800 }, 4);

        L(constant_times_4_dd_0x200);
        dd({ 0x200 }, 4);

        L(constant_times_8_dw_0x20);
        dw({ 0x20 }, 8);

        L(constant_times_4_dd_0x20);
        dd({ 0x20 }, 4);

        L(constant_times_8_dw_0x40);
        dw({ 0x40 }, 8);
    }

    void assemble() override
    {
        if (this->refs == 2)
            this->assembleBi();
        else
            this->assembleUni();
    }

    void assembleBi()
    {
        this->PRED_BI_V_8NxH(taps, 16, width);
    }

    int outputBitDepth = 8;

    void assembleUni()
    {
        auto &r0 = arg64(0);
        auto &r1 = arg64(1);
        auto &r2 = arg64(2);
        auto &r3 = arg64(3);
        auto &r4 = arg64(4);
        auto &r5 = arg64(5);
        auto &r6 = arg64(6);
        auto &r7 = arg64(7);

        if (xFrac && !yFrac)
        {
            this->PRED_UNI_H_16NxH(taps, outputBitDepth, width, false);
            return;
        }

        if (xFrac && yFrac)
        {
            this->stackSize = 64 * (64 + taps - 1) * 2 + 16;

            push(r0);
            lea(r0, ptr[rsp + 8]);
            push(r1);
            push(r5);
            push(r0);

            mov(r1, 64);

            for (int i = 0; i < sizeof(Sample); ++i)
            {
                if (taps == 8)
                {
                    sub(r2, r3);
                    sub(r2, r3);
                }
                sub(r2, r3);
            }

            add(r5, taps - 1);

            // first H then continue to V
            this->PRED_UNI_H_16NxH(taps, 16, width, true);

            pop(r0);
            lea(r2, ptr[r0 + 64 * 2 * (taps / 2 - 1)]);
            mov(r3, 64);

            pop(r5);
            pop(r1);
            pop(r0);

            inputBitDepth = 16;
        }

        Xbyak::Reg32 r5d(r5.getIdx());
        Xbyak::Reg32 r6d(r6.getIdx());

        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        //auto &m4 = regXmm(4);
        auto &m5 = regXmm(4); // !

        int const inputBits = (inputBitDepth == 16 || sizeof(Sample) == 2) ? 16 : 8;
        if (inputBits == 16) shl(r3, 1);

        if (sizeof(Sample) == 2) shl(r1, 1);

        // adjust input pointer (subtract (taps/2-1) * stride)
        for (int i = 0; i < taps / 2 - 1; ++i) sub(r2, r3);

        if (inputBits == 16)
        {
            lea(r4, ptr[rip + coefficients16]);
        }
        else
        {
            lea(r4, ptr[rip + coefficients]);
        }

        mov(r6, r7);
        shl(r6d, 4 + taps / 4);
        lea(r4, ptr[r4 + r6]);

        L("loopV");
        {
            for (int i = 0; i < width / 8; ++i)
            {
                // each iteration operates on a row of 8-samples
                if (inputBitDepth == 16)
                {
                    if (sizeof(Sample) == 2)
                        movdqa(m3, ptr[rip + constant_times_4_dd_0x200]);
                    else
                        movdqa(m3, ptr[rip + constant_times_4_dd_0x800]);
                    movdqa(m5, m3);
                }
                else if (sizeof(Sample) == 2)
                {
                    movdqa(m3, ptr[rip + constant_times_4_dd_0x20]);
                    movdqa(m5, m3);
                }
                else
                {
                    movdqa(m3, ptr[rip + constant_times_8_dw_0x20]);
                }

                for (int j = 0; j < taps / 2; ++j)
                {
                    // each iteration of this loop performs two filter taps

                    if (inputBits == 16)
                    {
                        movdqu(m0, ptr[r2]);
                        movdqu(m1, ptr[r2 + r3]);
                        movdqa(m2, m0);
                        punpckhwd(m2, m1);
                        punpcklwd(m0, m1);
                        pmaddwd(m0, ptr[r4]);
                        pmaddwd(m2, ptr[r4]);

                        paddd(m3, m0);
                        paddd(m5, m2);
                    }
                    else
                    {
                        movq(m0, ptr[r2]);
                        movq(m1, ptr[r2 + r3]);
                        punpcklbw(m0, m1);
                        pmaddubsw(m0, ptr[r4]);

                        paddw(m3, m0);
                    }

                    lea(r2, ptr[r2 + r3 * 2]);

                    add(r4, 16);
                }

                neg(r3);
                lea(r2, ptr[r2 + r3 * taps]);
                neg(r3);

                sub(r4, taps * 8);

                if (sizeof(Sample) == 2)
                {
                    if (inputBitDepth == 16)
                    {
                        psrad(m3, 4);
                        psrad(m5, 4);
                        packusdw(m3, m5); // should not saturate - just a shuffle
                        psrlw(m3, 6);
                    }
                    else
                    {
                        packusdw(m3, m5);
                        psrlw(m3, 6);
                    }
                    movdqu(ptr[r0], m3);
                }
                else
                {
                    if (inputBitDepth == 16)
                    {
                        psrad(m3, 12);
                        psrad(m5, 12);
                        packssdw(m3, m5);
                        packuswb(m3, m3);
                    }
                    else
                    {
                        psraw(m3, 6);
                        packuswb(m3, m3);
                    }
                    movq(ptr[r0], m3);
                }

                // review: these unnecessary in final loop iteration
                lea(r2, ptr[r2 + inputBits]);
                lea(r0, ptr[r0 + 8 * sizeof(Sample)]);
            }

            sub(r0, 8 * sizeof(Sample) * width / 8);
            sub(r2, inputBits * width / 8);

            add(r0, r1);
            add(r2, r3);
        }
        dec(r5d);
        jg("loopV");
    }

    void packedMultiplyAdd(Xbyak::Xmm const &dst, Xbyak::Xmm const &src)
    {
        if (sizeof(Sample) == 2) pmaddwd(dst, src); else pmaddubsw(dst, src);
    }

    void packedAdd(Xbyak::Xmm const &dst, Xbyak::Xmm const &src)
    {
        if (sizeof(Sample) == 2) paddd(dst, src); else paddw(dst, src);
    }

    // taps is number of filter taps (4 or 8)
    // outputTypeBits is size of output type (8 for uint8_t rounded, 16 for int16_t right shifted 6)
    // dx horizontal offset as integer number of bytes
    void filterHorizontal16Bytes(int taps, int, int dx, bool filteringHV)
    {
        // note: comments apply when Sample is uint8_t

        auto &r0 = reg64(0);
        auto &r2 = reg64(2);

        auto &m0 = xmm0;
        auto &m1 = xmm1;
        auto &m2 = xmm2;
        auto &m3 = xmm3;
        auto &m4 = xmm4;
        auto &m5 = xmm5;
        auto &m6 = xmm6;
        auto &m7 = xmm7;

        movdqu(m2, ptr[r2 + sizeof(Sample) * (1 - (taps / 2)) + dx]);
        packedMultiplyAdd(m2, m4);
        // m2 = eca86420 (even positions, first two taps)

        movdqu(m1, ptr[r2 + sizeof(Sample) * (2 - (taps / 2)) + dx]);
        packedMultiplyAdd(m1, m4);
        // m1 = fdb97531 (odd positions, first two taps)

        movdqu(m0, ptr[r2 + sizeof(Sample) * (3 - (taps / 2)) + dx]);
        packedMultiplyAdd(m0, m5);
        // m0 = eca86420 (even positions, two taps)

        packedAdd(m2, m0);
        // m2 = eca86420 (even positions)

        movdqu(m0, ptr[r2 + sizeof(Sample) * (4 - (taps / 2)) + dx]);
        packedMultiplyAdd(m0, m5);
        // m0 = fdb97531 (odd positions, two taps)

        packedAdd(m1, m0);
        // m1 = fdb97531 (odd positions)

        if (taps == 8)
        {
            // need four more taps...

            movdqu(m0, ptr[r2 + sizeof(Sample) * (5 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m6);
            // m0 = eca86420(even positions, two taps)

            packedAdd(m2, m0);
            // m2 = eca86420(even positions)

            movdqu(m0, ptr[r2 + sizeof(Sample) * (6 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m6);
            // m0 = fdb97531(odd positions, two taps)

            packedAdd(m1, m0);
            // m1 = fdb97531(odd positions)

            movdqu(m0, ptr[r2 + sizeof(Sample) * (7 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m7);
            // m0 = eca86420(even positions, two taps)

            packedAdd(m2, m0);
            // m2 = eca86420(even positions)

            movdqu(m0, ptr[r2 + sizeof(Sample) * (8 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m7);
            // m0 = fdb97531(odd positions, two taps)

            packedAdd(m1, m0);
            // m1 = fdb97531(odd positions)
        }

        movq(m0, m2);
        if (sizeof(Sample) == 2) punpckldq(m0, m1); else punpcklwd(m0, m1);
        // m0 = 76543210

        if (sizeof(Sample) == 2) punpckhdq(m2, m1); else punpckhwd(m2, m1);
        // m2 = fedcba98

        if (sizeof(Sample) == 2)
        {
            if (filteringHV)
            {
                psrad(m0, 2);
                psrad(m2, 2);
                packssdw(m0, m2);
                movdqu(ptr[r0 + dx], m0);
            }
            else
            {
                paddd(m0, ptr[rip + constant_times_4_dd_0x20]);
                paddd(m2, ptr[rip + constant_times_4_dd_0x20]);
                packusdw(m0, m2);
                psrlw(m0, 6);
                movdqu(ptr[r0 + dx], m0);
            }
        }
        else
        {
            if (filteringHV)
            {
                movdqu(ptr[r0 + 2 * dx], m0);
                movdqu(ptr[r0 + 2 * dx + 16], m2);
            }
            else
            {
                paddw(m0, ptr[rip + constant_times_8_dw_0x20]);
                paddw(m2, ptr[rip + constant_times_8_dw_0x20]);
                psraw(m0, 6);
                psraw(m2, 6);
                packuswb(m0, m2);
                movdqu(ptr[r0 + dx], m0);
            }
        }
    }


    // taps is number of filter taps (4 or 8)
    // width is block width (number of samples, multiple of 16)
    void PRED_UNI_H_16NxH(int taps, int, int width, bool filteringHV)
    {
        auto &r0 = reg64(0);
        auto &r1 = reg64(1);
        auto &r2 = reg64(2);
        auto &r3 = reg64(3);
        auto &r4 = reg64(4);
        auto &r5 = reg64(5);
        Xbyak::Reg32 r5d(r5.getIdx());
        auto &r6 = reg64(6);
        Xbyak::Reg32 r6d(r6.getIdx());

        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        auto &m4 = regXmm(4);
        auto &m5 = regXmm(5);
        auto &m6 = regXmm(6);
        auto &m7 = regXmm(7);

        if (sizeof(Sample) == 2)
            lea(r4, ptr[rip + coefficients16]);
        else
            lea(r4, ptr[rip + coefficients]);

        if (taps == 8)
        {
            shl(r6d, 6);  // frac *= 4 * 16
        }
        else
        {
            shl(r6d, 5);  // frac *= 2 * 16
        }

        movdqa(m4, ptr[r4 + r6]);
        movdqa(m5, ptr[r4 + r6 + 1 * 16]);
        if (taps == 8)
        {
            movdqa(m6, ptr[r4 + r6 + 2 * 16]);
            movdqa(m7, ptr[r4 + r6 + 3 * 16]);
        }

        int shiftDstStride = 0;
        if (sizeof(Sample) == 2 || filteringHV) ++shiftDstStride;

        if (shiftDstStride) shl(r1, shiftDstStride);

        if (sizeof(Sample) == 2)
        {
            // src is uint16_t * so need to double the stride
            shl(r3, 1);
        }

        L("loop");
        {
            for (int i = 0; i < width * sizeof(Sample) / 16; ++i)
            {
                int const dx = 16 * i;
                filterHorizontal16Bytes(taps, 0, dx, filteringHV);
            }

            add(r0, r1);
            add(r2, r3);
        }
        dec(r5d);
        jg("loop");
    }

    // void havoc_pred_bi_v_%1tap_16to16_%3xh_sse4(uint8_t *dst, intptr_t stride_dst, const int16_t *refAtop, const int16_t *refBtop, intptr_t stride_ref, int nPbW, int nPbH, int yFracA, int yFracB)//
    void PRED_BI_V_8NxH(int taps, int inputTypeSize, int width)
    {
        // taps is number of filter taps (4 or 8)//
        // inputTypeSize is size of input type (8 for uint8_t, 16 for int16_t right shifted 6)
        // width is block width (number of samples, multiple of 8)
        assert(inputTypeSize == 16);

        // // void havoc_pred_bi_v_%1tap_16to16_%3xh_sse4(uint8_t *dst, intptr_t stride_dst, const int16_t *refAtop, const int16_t *refBtop, intptr_t stride_ref, int nPbW, int nPbH, int yFracA, int yFracB)//
        // INIT_XMM sse4
        // cglobal pred_bi_v_%1tap_16to16_%3xh, 9, 9, 8

        auto &r0 = arg64(0); // dst
        auto &r1 = arg64(1); // stride_dst
        auto &r2 = arg64(2); // refA
        auto &r3 = arg64(3); // refB
        auto &r4 = arg64(4); // stride_ref
        auto &r5 = arg64(5); // width
        auto &r6 = arg64(6);
        auto r6d = Xbyak::Reg32(r6.getIdx()); // height
        auto &r7 = arg64(7);
        auto r7d = Xbyak::Reg32(r7.getIdx()); // yFracA
        auto &r8 = arg64(8);
        auto r8d = Xbyak::Reg32(r8.getIdx()); // yFracB

        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        auto &m4 = regXmm(4);
        auto &m5 = regXmm(5);
        auto &m6 = regXmm(6);
        auto &m7 = regXmm(7);

        shl(r4, 1);

#if 1 // ARCH_X86_64
        shl(r7d, 4 + taps / 4);
        shl(r8d, 4 + taps / 4);
        lea(r5, ptr[rip + coefficients16]);
        lea(r7, ptr[r5 + r7]);
        lea(r8, ptr[r5 + r8]);
#define coeffA r7
#define coeffB r8
#define stride_dst r1
#else
        mov r5, r7m
            shl r5, 4 + taps / 4
            lea r1, [pred_inter_ % 1tap_coefficient_pairs_4_dw + r5]
            mov r5, r8m
            shl r5, 4 + taps / 4
            lea r5, [pred_inter_ % 1tap_coefficient_pairs_4_dw + r5]
            % define coeffA r1
            %define coeffB r5
            %define stride_dst r1m
#endif
            L("loopBiV");
        {
            for (int i = 0; i < width / 8; ++i)
            {
                // each iteration of this loop operates on a row of 8-samples

                pxor(m3, m3);
                movdqa(m5, m3);
                movdqa(m6, m3);
                movdqa(m7, m3);

                for (int j = 0; j < taps / 2; ++j)
                {
                    // each iteration of this loop performs two filter taps

                    // reference picture A
                    movdqu(m0, ptr[r2]);
                    movdqu(m1, ptr[r2 + r4]);
                    movdqa(m2, m0);
                    punpckhwd(m2, m1);
                    punpcklwd(m0, m1);
                    pmaddwd(m0, ptr[coeffA]);
                    pmaddwd(m2, ptr[coeffA]);
                    paddd(m3, m0);
                    paddd(m5, m2);
                    lea(r2, ptr[r2 + r4 * 2]);

                    // reference picture B
                    movdqu(m0, ptr[r3]);
                    movdqu(m1, ptr[r3 + r4]);
                    movdqa(m2, m0);
                    punpckhwd(m2, m1);
                    punpcklwd(m0, m1);
                    pmaddwd(m0, ptr[coeffB]);
                    pmaddwd(m2, ptr[coeffB]);
                    paddd(m6, m0);
                    paddd(m7, m2);
                    lea(r3, ptr[r3 + r4 * 2]);

                    add(coeffA, 16);
                    add(coeffB, 16);
                }

                neg(r4);
                lea(r2, ptr[r2 + r4 * taps]);
                lea(r3, ptr[r3 + r4 * taps]);
                neg(r4);

                sub(coeffA, 8 * taps);
                sub(coeffB, 8 * taps);

                psrad(m3, 6);
                psrad(m5, 6);
                psrad(m6, 6);
                psrad(m7, 6);

                packssdw(m3, m5);
                packssdw(m6, m7);

                paddsw(m3, m6);
                paddsw(m3, ptr[rip + constant_times_8_dw_0x40]);
                psraw(m3, 7);

                packuswb(m3, m3);

                movdqu(ptr[r0], m3);

                lea(r2, ptr[r2 + 16]);
                lea(r3, ptr[r3 + 16]);
                lea(r0, ptr[r0 + 8]);
            }
            sub(r2, 2 * width);
            sub(r3, 2 * width);
            sub(r0, width);

            add(r2, r4);
            add(r3, r4);
            add(r0, stride_dst);
        }
        dec(r6d);
        jg("loopBiV");
#undef coeffA
#undef coeffB
#undef stride_dst
    }
};


// returns true if width is a PU width used in HEVC
bool legalWidth(int width)
{
    while ((width & 1) == 0) width >>= 1;
    return width == 1 || width == 3;
}


template <typename Sample>
void havocPopulatePredUni(HavocTablePredUni<Sample> *table, havoc_code code)
{
    int const maxBitDepth = 6 + 2 * sizeof(Sample);

    for (int taps = 4; taps <= 8; taps += 4)
        for (int w = taps; w <= 64; w += taps)
        {
            for (int xFrac = 0; xFrac < 2; ++xFrac)
                for (int yFrac = 0; yFrac < 2; ++yFrac)
                    for (int bitDepth = 8; bitDepth <= maxBitDepth; ++bitDepth)
                        *havocGetPredUni<Sample>(table, taps, w, 0, xFrac, yFrac, bitDepth) = 0;
        }

    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
    auto const n = Jit::CountArguments<typename std::remove_pointer<HavocPredUni<Sample>>::type>::value;

    if (buffer.isa & HAVOC_C_REF)
        for (int taps = 4; taps <= 8; taps += 4)
            for (int w = 0; w <= 64; w += taps)
                for (int xFrac = 0; xFrac < 2; ++xFrac)
                    for (int yFrac = 0; yFrac < 2; ++yFrac)
                        for (int bitDepth = 8; bitDepth <= maxBitDepth; ++bitDepth)
                            *havocGetPredUni<Sample>(table, taps, w, 0, xFrac, yFrac, bitDepth)
                            = taps == 8 ? havoc_pred_uni_8tap_hv<Sample> : havoc_pred_uni_4tap_hv<Sample>;

    if (buffer.isa & HAVOC_C_OPT)
        for (int bitDepth = 8; bitDepth <= maxBitDepth; ++bitDepth)
        {
            for (int w = 1; w <= 64; ++w)
            {
                if (!legalWidth(w)) continue;

                if (w >= 4 && w != 6)
                {
                    *havocGetPredUni(table, 8, w, 0, 0, 0, bitDepth) = havoc_pred_uni_copy_block<Sample>;
                    *havocGetPredUni(table, 8, w, 0, 1, 0, bitDepth) = havoc_pred_uni_8tap_h<Sample>;
                    *havocGetPredUni(table, 8, w, 0, 0, 1, bitDepth) = havoc_pred_uni_8tap_v<Sample>;
                    *havocGetPredUni(table, 8, w, 0, 1, 1, bitDepth) = havoc_pred_uni_8tap_hv<Sample>;
                }

                if (w >= 2 && w != 3)
                {
                    *havocGetPredUni(table, 4, w, 0, 0, 0, bitDepth) = havoc_pred_uni_copy_block<Sample>;
                    *havocGetPredUni(table, 4, w, 0, 1, 0, bitDepth) = havoc_pred_uni_4tap_h<Sample>;
                    *havocGetPredUni(table, 4, w, 0, 0, 1, bitDepth) = havoc_pred_uni_4tap_v<Sample>;
                    *havocGetPredUni(table, 4, w, 0, 1, 1, bitDepth) = havoc_pred_uni_4tap_hv<Sample>;
                }
            }
        }

    if (buffer.isa & HAVOC_SSE2)
    {
        int const step = 16 / sizeof(Sample);
        for (int bitDepth = 8; bitDepth <= std::min(10, 8 * (int)sizeof(Sample)); ++bitDepth)
            for (int w = step; w <= 64; w += step)
            {
                if (!legalWidth(w)) continue;

                PredUniCopy<Sample> a(&buffer, w);

                *havocGetPredUni(table, 8, w - 8, 0, 0, 0, bitDepth) = a;
                *havocGetPredUni(table, 8, w, 0, 0, 0, bitDepth) = a;
                *havocGetPredUni(table, 4, w - 12, 0, 0, 0, bitDepth) = a;
                *havocGetPredUni(table, 4, w - 8, 0, 0, 0, bitDepth) = a;
                *havocGetPredUni(table, 4, w - 4, 0, 0, 0, bitDepth) = a;
                *havocGetPredUni(table, 4, w, 0, 0, 0, bitDepth) = a;
            }
    }

    if (buffer.isa & HAVOC_SSE41)
        for (int taps = 4; taps <= 8; taps += 4)
            for (int w = 16; w <= 64; w += 16)
            {
                if (!legalWidth(w)) continue;

                int bitDepth = sizeof(Sample) == 2 ? 10 : 8;

                PredInter<Sample> aH(&buffer, n, taps, w, 1, 0, 8, 8);
                PredInter<Sample> aV(&buffer, n, taps, w, 0, 1, 8, 8);
                PredInter<Sample> aHV(&buffer, n, taps, w, 1, 1, 8, 8);

                *havocGetPredUni(table, taps, w, 0, 1, 0, bitDepth) = aH;
                *havocGetPredUni(table, taps, w - 8, 0, 1, 0, bitDepth) = aH;
                *havocGetPredUni(table, taps, w, 0, 0, 1, bitDepth) = aV;
                *havocGetPredUni(table, taps, w, 0, 1, 1, bitDepth) = aHV;
                *havocGetPredUni(table, taps, w - 8, 0, 1, 1, bitDepth) = aHV;

                *havocGetPredUni(table, taps, w - 4, 0, 1, 0, bitDepth) = aH;
                *havocGetPredUni(table, taps, w - 4 - 8, 0, 1, 0, bitDepth) = aH;
                *havocGetPredUni(table, taps, w - 4, 0, 0, 1, bitDepth) = aV;
                *havocGetPredUni(table, taps, w - 4, 0, 1, 1, bitDepth) = aHV;
                *havocGetPredUni(table, taps, w - 4 - 8, 0, 1, 1, bitDepth) = aHV;

                if (legalWidth(w - 8))
                {
                    PredInter<Sample> aV2(&buffer, n, taps, w - 8, 0, 1, 8, 8);
                    *havocGetPredUni(table, taps, w - 8, 0, 0, 1, bitDepth) = aV2;
                    if (taps == 4)
                        *havocGetPredUni(table, taps, w - 4 - 8, 0, 0, 1, bitDepth) = aV2;
                }
            }

}



#define STRIDE_DST 192

typedef struct
{
    HAVOC_ALIGN(32, uint8_t, dst8[64 * STRIDE_DST]);
    HAVOC_ALIGN(32, uint16_t, dst16[64 * STRIDE_DST]);
    havoc_pred_uni_8to8 *f8;
    havoc_pred_uni_16to16 *f16;
    intptr_t stride_dst;
    const uint8_t *ref8;
    const uint16_t *ref16;
    intptr_t stride_ref;
    int w;
    int h;
    int xFrac;
    int yFrac;
    int taps;
    int bitDepth;
}
bound_pred_uni;


static int get_pred_uni(void *p, havoc_code code)
{
    bound_pred_uni *s = (bound_pred_uni *)p;
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    s->f8 = 0;
    s->f16 = 0;

    if (s->bitDepth > 8)
    {
        HavocTablePredUni<uint16_t> table;
        havocPopulatePredUni<uint16_t>(&table, code);
        s->f16 = *havocGetPredUni<uint16_t>(&table, s->taps, s->w, s->h, s->xFrac, s->yFrac, s->bitDepth);
        memset(s->dst16, 0, 2 * 64 * s->stride_dst);
    }
    else
    {
        HavocTablePredUni<uint8_t> table;
        havocPopulatePredUni<uint8_t>(&table, code);
        s->f8 = *havocGetPredUni<uint8_t>(&table, s->taps, s->w, s->h, s->xFrac, s->yFrac, s->bitDepth);
        memset(s->dst8, 0, 64 * s->stride_dst);
    }

    if ((s->f8 || s->f16) && buffer.isa == HAVOC_C_REF)
    {
        printf("\t%d-bit %d-tap %dx%d %s%s : ", s->bitDepth, s->taps, s->w, s->h, s->xFrac ? "H" : "", s->yFrac ? "V" : "");
    }

    return s->f8 || s->f16;
}


void invoke_pred_uni(void *p, int n)
{
    bound_pred_uni *s = (bound_pred_uni *)p;
    if (s->bitDepth > 8) while (n--)
    {
        s->f16(s->dst16, s->stride_dst, s->ref16, s->stride_ref, s->w, s->h, s->xFrac, s->yFrac, s->bitDepth);
    }
    else while (n--)
    {
        s->f8(s->dst8, s->stride_dst, s->ref8, s->stride_ref, s->w, s->h, s->xFrac, s->yFrac, s->bitDepth);
    }
}


int mismatch_pred_uni(void *boundRef, void *boundTest)
{
    bound_pred_uni *ref = (bound_pred_uni *)boundRef;
    bound_pred_uni *test = (bound_pred_uni *)boundTest;

    for (int y = 0; y < ref->h; ++y)
    {
        if (ref->bitDepth > 8)
        {
            if (memcmp(&ref->dst16[y*ref->stride_dst], &test->dst16[y*test->stride_dst], 2 * ref->w))
                return 1;
        }
        else
        {
            if (memcmp(&ref->dst8[y*ref->stride_dst], &test->dst8[y*test->stride_dst], ref->w))
                return 1;
        }
    }

    return 0;
}


static void test_partitions(int *error_count, bound_pred_uni *b, havoc_instruction_set mask)
{
    const int partitions[24][2] =
    {
        { 8, 4 }, { 8, 8 }, { 4, 8 },
        { 16, 4 }, { 16, 8 }, { 16, 12 }, { 16, 16 }, { 12, 16 }, { 8, 16 }, { 4, 16 },
        { 32, 8 }, { 32, 16 }, { 32, 24 }, { 32, 32 }, { 24, 32 }, { 16, 32 }, { 8, 32 },
        { 64, 16 }, { 64, 32 }, { 64, 48 }, { 64, 64 }, { 48, 64 }, { 32, 64 }, { 16, 64 },
    };

    for (int k = 0; k < 24; ++k)
    {
        const int nPbW = partitions[k][0];
        const int nPbH = partitions[k][1];

        b[0].w = nPbW * b[0].taps / 8;
        b[0].h = nPbH * b[0].taps / 8;

        b[1] = b[0];

        *error_count += havoc_test(&b[0], &b[1], get_pred_uni, invoke_pred_uni, mismatch_pred_uni, mask, 10);
    }

    if (b[0].taps == 4)
    {
        // additional tests for large chroma (used with 4:2:2 and 4:4:4)

        for (int k = 17; k < 24; ++k)
        {
            const int nPbW = partitions[k][0];
            const int nPbH = partitions[k][1];

            b[0].w = nPbW;
            b[0].h = nPbH;

            b[1] = b[0];

            *error_count += havoc_test(&b[0], &b[1], get_pred_uni, invoke_pred_uni, mismatch_pred_uni, mask, 10);
        }
    }
}


void havoc_test_pred_uni(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_pred_uni - Unireference Inter Prediction (single-reference motion compensation)\n");

    bound_pred_uni b[2];

#define STRIDE_REF 192
    HAVOC_ALIGN(32, uint8_t, ref8[80 * STRIDE_REF]);
    HAVOC_ALIGN(32, uint16_t, ref16[80 * STRIDE_REF]);
    b[0].stride_dst = STRIDE_DST;
    b[0].stride_ref = STRIDE_REF;
    b[0].ref8 = ref8 + 8 * b[0].stride_ref;
    b[0].ref16 = ref16 + 8 * b[0].stride_ref;
#undef STRIDE_DST
#undef STRIDE_REF

    for (b[0].bitDepth = 10; b[0].bitDepth >= 8; --b[0].bitDepth)
    {
        for (int x = 0; x < 80 * b[0].stride_ref; x++)
        {
            ref8[x] = rand() & 0xff;
            ref16[x] = rand() % (1 << b[0].bitDepth);
        }

        for (b[0].taps = 8; b[0].taps >= 4; b[0].taps -= 4)
        {
            for (b[0].yFrac = 0; b[0].yFrac < 2; ++b[0].yFrac)
            {
                for (b[0].xFrac = 0; b[0].xFrac < 2; ++b[0].xFrac)
                {
                    test_partitions(error_count, b, mask);
                }
            }
        }
    }
}


void havoc_pred_bi_mean_32and32to8_c_ref(uint8_t *dst, intptr_t stride_dst, const int *ref0, const int *ref1, intptr_t stride_ref, int w, int h, int bits, int shift)
{
    assert(bits == 8);
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            int v = ((int)ref0[x + y * stride_ref] + (int)ref1[x + y * stride_ref] + (1 << (shift - 1))) >> shift;
            v = Clip3(0, 255, v);
            dst[x + y * stride_dst] = v;
        }
    }
}


template <typename Sample>
void havoc_pred_bi_mean_c_ref(Sample *dst, intptr_t stride_dst, const int *ref0, const int *ref1, intptr_t stride_ref, int w, int h, int bits, int shift)
{
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            int v = ((int)ref0[x + y * stride_ref] + (int)ref1[x + y * stride_ref] + (1 << (shift - 1))) >> shift;
            v = Clip3(0, (1 << bits) - 1, v);
            dst[x + y * stride_dst] = v;
        }
    }
}


template <typename Sample, int taps>
void havocPredBi_c_ref(Sample *dst, intptr_t stride_dst, const Sample *ref0, const Sample *ref1, intptr_t stride_ref, int nPbW, int nPbH, int xFrac0, int yFrac0, int xFrac1, int yFrac1, int bitDepth)
{
    int intermediate[4][(64 + taps - 1) * 64];

    const int w = nPbW;
    const int h = nPbH;

    int shift1 = bitDepth - 8;
    if (shift1 > 4) shift1 = 4;

    int shift2 = 6;

    int shift3 = 14 - bitDepth;
    if (shift3 < 2) shift3 = 2;

    /* Horizontal filter */
    havoc_pred_uni_generic(intermediate[2], 64, ref0 - (taps / 2 - 1) * stride_ref, stride_ref, w, h + taps - 1, 1, taps, xFrac0, shift1, 0, 0);

    /* Vertical filter */
    havoc_pred_uni_generic(intermediate[0], 64, intermediate[2] + (taps / 2 - 1) * 64, 64, w, h, 64, taps, yFrac0, shift2, 0, 0);

    /* Horizontal filter */
    havoc_pred_uni_generic(intermediate[3], 64, ref1 - (taps / 2 - 1) * stride_ref, stride_ref, w, h + taps - 1, 1, taps, xFrac1, shift1, 0, 0);

    /* Vertical filter */
    havoc_pred_uni_generic(intermediate[1], 64, intermediate[3] + (taps / 2 - 1) * 64, 64, w, h, 64, taps, yFrac1, shift2, 0, 0);

    /* Combine two references for bi pred */
    havoc_pred_bi_mean_c_ref(dst, stride_dst, intermediate[0], intermediate[1], 64, nPbW, nPbH, bitDepth, shift3 + 1);
}





template <typename Sample>
struct PredBi
    :
    Jit::Function
{
    PredBi(Jit::Buffer *buffer, int width, int taps = 0)
        :
        Jit::Function(buffer, Jit::CountArguments<HavocPredBi<uint8_t>>::value),
        width(width),
        taps(taps)
    {
        static int n = 0;
        this->build();
    }

    int width, taps;

    void assembleInterp()
    {
        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        auto &m4 = regXmm(4);
        auto &m5 = regXmm(5);
        auto &m6 = regXmm(6);
        auto &m7 = regXmm(7);

        auto &r0 = arg64(0); // dst
        auto &r1 = arg64(1);
        auto &r2 = arg64(2); //ref0
        auto &r3 = arg64(3); //ref1
        auto &r4 = arg64(4);
        auto &r5 = arg64(5);// w
        auto &r6 = arg64(6);// h
        auto &r7 = arg64(7);// xFrac0
        auto &r8 = arg64(8);// yFrac0
        auto &r9 = arg64(9);// xFrac1
        auto &r10 = arg64(10);// yFrac1

        this->stackSize = (64 + taps - 1) * 64 * 2 * 2 + 12 * 8;

        // save all arguments in allocated stack memory
        for (int i = 0; i < 11; ++i)
        {
            mov(ptr[rsp + 8 * i], reg64(i));
        }

        // setup arguments to call H filter (ref0)
        lea(r0, ptr[rsp + 12 * 8]);
        mov(r1, 64);
        mov(r3, r4);
        neg(r4);
        if (taps == 8)
        {
            lea(r2, ptr[r2 + r4 * (2 * sizeof(Sample))]);
        }
        lea(r2, ptr[r2 + r4 * sizeof(Sample)]);
        mov(r4, r5); // w
        mov(r5, r6); // h
        add(r5, taps - 1);
        mov(r6, r7);// xFrac
        xor (r7, r7); // yFrac

        // naked call (registers preserved)
        call(labelPRED_UNI_H_16NxH);

        // restore arguments from stack memory
        for (int i = 0; i < 11; ++i)
        {
            mov(reg64(i), ptr[rsp + 8 * i]);
        }

        // setup arguments to call H filter (ref1)
        lea(r0, ptr[rsp + 12 * 8 + (64 + taps - 1) * 64 * 2]);
        mov(r1, 64);
        mov(r2, r3);
        mov(r3, r4);
        neg(r4);
        if (taps == 8)
        {
            lea(r2, ptr[r2 + r4 * (2 * sizeof(Sample))]);
        }
        lea(r2, ptr[r2 + r4 * sizeof(Sample)]);
        mov(r4, r5); // w
        mov(r5, r6); // h
        add(r5, taps - 1);
        mov(r6, r9);// xFrac
        xor (r7, r7); // yFrac

        // naked call (registers preserved)
        call(labelPRED_UNI_H_16NxH);

        // restore arguments from stack memory
        for (int i = 0; i < 11; ++i)
        {
            mov(reg64(i), ptr[rsp + 8 * i]);
        }


        lea(r2, ptr[rsp + 12 * 8]); // refAtop
        lea(r3, ptr[rsp + 12 * 8 + (64 + taps - 1) * 64 * 2]); // refBtop
        mov(r4, 64); // stride_ref
        mov(r7, r8); //yFrac0
        mov(r8, r10); //yFrac1

        // inline "function"
        // void havoc_pred_bi_v_%1tap_16to16_%3xh_sse4(uint8_t *dst, intptr_t stride_dst, const int16_t *refAtop, const int16_t *refBtop, intptr_t stride_ref, int nPbW, int nPbH, int yFracA, int yFracB)//
        this->PRED_BI_V_8NxH(this->taps, 16, this->width);



        //p16(intermediate[2], 64, ref0 - (taps / 2 - 1) * stride_ref, stride_ref, w, h + taps - 1, xFrac0, 0); \


    }

    void assemble() override
    {
        if (taps)
        {
            assembleInterp();
            return;
        }

        auto &r0 = arg64(0);
        auto &r1 = arg64(1);
        auto &r2 = arg64(2);
        auto &r3 = arg64(3);
        auto &r4 = arg64(4);
        auto &r5 = arg64(5);
        auto r6d = Xbyak::Reg32(arg64(6).getIdx());

        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);

        if (sizeof(Sample) == 2)
        {
            shl(r1, 1);
            shl(r4, 1);
        }

        L("loop");
        {
            for (int dx = 0; dx < width * sizeof(Sample); dx += 16)
            {
                movdqu(m0, ptr[r2 + dx]);
                movdqu(m1, ptr[r3 + dx]);
                if (sizeof(Sample) == 2)
                    pavgw(m0, m1);
                else
                    pavgb(m0, m1);
                movdqu(ptr[r0 + dx], m0);
            }
            lea(r2, ptr[r2 + r4]);
            lea(r3, ptr[r3 + r4]);
            lea(r0, ptr[r0 + r1]);
        }
        dec(r6d);
        jg("loop");
    }

    Xbyak::Label coefficients;
    Xbyak::Label coefficients16;
    Xbyak::Label constant_times_4_dd_0x800;
    Xbyak::Label constant_times_4_dd_0x200;
    Xbyak::Label constant_times_4_dd_0x20;
    Xbyak::Label constant_times_8_dw_0x20;
    Xbyak::Label constant_times_8_dw_0x40;
    Xbyak::Label constant_times_8_dw_0x10;
    Xbyak::Label labelPRED_UNI_H_16NxH;
    Xbyak::Label labelPRED_BI_V_8NxH;

    void data() override
    {
        align();

        L(coefficients);
        for (int frac = 0; frac < (12 - taps); ++frac)
        {
            for (int k = 0; k < taps; k += 2)
            {
                int coeff0 = havoc_pred_coefficient(taps, frac, k);
                int coeff1 = havoc_pred_coefficient(taps, frac, k + 1);

                db({ coeff0, coeff1 }, 8);

            }
        }

        L(coefficients16);
        for (int frac = 0; frac < (12 - taps); ++frac)
        {
            for (int k = 0; k < taps; k += 2)
            {
                int coeff0 = havoc_pred_coefficient(taps, frac, k);
                int coeff1 = havoc_pred_coefficient(taps, frac, k + 1);

                dw({ coeff0, coeff1 }, 4);
            }
        }

        L(constant_times_4_dd_0x800);
        dd({ 0x800 }, 4);

        L(constant_times_4_dd_0x20);
        dd({ 0x20 }, 4);

        L(constant_times_8_dw_0x20);
        dw({ 0x20 }, 8);

        L(constant_times_8_dw_0x40);
        dw({ 0x40 }, 8);

        L(constant_times_8_dw_0x10);
        dw({ 0x10 }, 8);

        if (this->taps)
        {
            L(labelPRED_UNI_H_16NxH);
            this->PRED_UNI_H_16NxH(this->taps, 16, this->width, true);
            ret();
        }
    }

    void packedMultiplyAdd(Xbyak::Xmm const &dst, Xbyak::Xmm const &src)
    {
        if (sizeof(Sample) == 2) pmaddwd(dst, src); else pmaddubsw(dst, src);
    }

    void packedAdd(Xbyak::Xmm const &dst, Xbyak::Xmm const &src)
    {
        if (sizeof(Sample) == 2) paddd(dst, src); else paddw(dst, src);
    }

    // taps is number of filter taps (4 or 8)
    // outputTypeBits is size of output type (8 for uint8_t rounded, 16 for int16_t right shifted 6)
    // dx horizontal offset as integer number of bytes
    void filterHorizontal16Bytes(int taps, int, int dx, bool filteringHV)
    {
        // note: comments apply when Sample is uint8_t

        auto &r0 = reg64(0);
        auto &r2 = reg64(2);

        auto &m0 = xmm0;
        auto &m1 = xmm1;
        auto &m2 = xmm2;
        auto &m3 = xmm3;
        auto &m4 = xmm4;
        auto &m5 = xmm5;
        auto &m6 = xmm6;
        auto &m7 = xmm7;

        movdqu(m2, ptr[r2 + sizeof(Sample) * (1 - (taps / 2)) + dx]);
        packedMultiplyAdd(m2, m4);
        // m2 = eca86420 (even positions, first two taps)

        movdqu(m1, ptr[r2 + sizeof(Sample) * (2 - (taps / 2)) + dx]);
        packedMultiplyAdd(m1, m4);
        // m1 = fdb97531 (odd positions, first two taps)

        movdqu(m0, ptr[r2 + sizeof(Sample) * (3 - (taps / 2)) + dx]);
        packedMultiplyAdd(m0, m5);
        // m0 = eca86420 (even positions, two taps)

        packedAdd(m2, m0);
        // m2 = eca86420 (even positions)

        movdqu(m0, ptr[r2 + sizeof(Sample) * (4 - (taps / 2)) + dx]);
        packedMultiplyAdd(m0, m5);
        // m0 = fdb97531 (odd positions, two taps)

        packedAdd(m1, m0);
        // m1 = fdb97531 (odd positions)

        if (taps == 8)
        {
            // need four more taps...

            movdqu(m0, ptr[r2 + sizeof(Sample) * (5 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m6);
            // m0 = eca86420(even positions, two taps)

            packedAdd(m2, m0);
            // m2 = eca86420(even positions)

            movdqu(m0, ptr[r2 + sizeof(Sample) * (6 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m6);
            // m0 = fdb97531(odd positions, two taps)

            packedAdd(m1, m0);
            // m1 = fdb97531(odd positions)

            movdqu(m0, ptr[r2 + sizeof(Sample) * (7 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m7);
            // m0 = eca86420(even positions, two taps)

            packedAdd(m2, m0);
            // m2 = eca86420(even positions)

            movdqu(m0, ptr[r2 + sizeof(Sample) * (8 - (taps / 2)) + dx]);
            packedMultiplyAdd(m0, m7);
            // m0 = fdb97531(odd positions, two taps)

            packedAdd(m1, m0);
            // m1 = fdb97531(odd positions)
        }

        movq(m0, m2);
        if (sizeof(Sample) == 2) punpckldq(m0, m1); else punpcklwd(m0, m1);
        // m0 = 76543210

        if (sizeof(Sample) == 2) punpckhdq(m2, m1); else punpckhwd(m2, m1);
        // m2 = fedcba98

        if (sizeof(Sample) == 2)
        {
            if (filteringHV)
            {
                psrad(m0, 2);
                psrad(m2, 2);
                packssdw(m0, m2);
                movdqu(ptr[r0 + dx], m0);
            }
            else
            {
                paddd(m0, ptr[rip + constant_times_4_dd_0x20]);
                paddd(m2, ptr[rip + constant_times_4_dd_0x20]);
                packusdw(m0, m2);
                psrlw(m0, 6);
                movdqu(ptr[r0 + dx], m0);
            }
        }
        else
        {
            if (filteringHV)
            {
                movdqu(ptr[r0 + 2 * dx], m0);
                movdqu(ptr[r0 + 2 * dx + 16], m2);
            }
            else
            {
                paddw(m0, ptr[rip + constant_times_8_dw_0x20]);
                paddw(m2, ptr[rip + constant_times_8_dw_0x20]);
                psraw(m0, 6);
                psraw(m2, 6);
                packuswb(m0, m2);
                movdqu(ptr[r0 + dx], m0);
            }
        }
    }


    // taps is number of filter taps (4 or 8)
    // width is block width (number of samples, multiple of 16)
    void PRED_UNI_H_16NxH(int taps, int, int width, bool filteringHV)
    {
        auto &r0 = reg64(0);
        auto &r1 = reg64(1);
        auto &r2 = reg64(2);
        auto &r3 = reg64(3);
        auto &r4 = reg64(4);
        auto &r5 = reg64(5);
        Xbyak::Reg32 r5d(r5.getIdx());
        auto &r6 = reg64(6);
        Xbyak::Reg32 r6d(r6.getIdx());

        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        auto &m4 = regXmm(4);
        auto &m5 = regXmm(5);
        auto &m6 = regXmm(6);
        auto &m7 = regXmm(7);

        if (sizeof(Sample) == 2)
            lea(r4, ptr[rip + coefficients16]);
        else
            lea(r4, ptr[rip + coefficients]);

        if (taps == 8)
        {
            shl(r6d, 6);  // frac *= 4 * 16
        }
        else
        {
            shl(r6d, 5);  // frac *= 2 * 16
        }

        movdqa(m4, ptr[r4 + r6]);
        movdqa(m5, ptr[r4 + r6 + 1 * 16]);
        if (taps == 8)
        {
            movdqa(m6, ptr[r4 + r6 + 2 * 16]);
            movdqa(m7, ptr[r4 + r6 + 3 * 16]);
        }

        int shiftDstStride = 0;
        if (sizeof(Sample) == 2 || filteringHV) ++shiftDstStride;

        if (shiftDstStride) shl(r1, shiftDstStride);

        if (sizeof(Sample) == 2)
        {
            // src is uint16_t * so need to double the stride
            shl(r3, 1);
        }

        L("loop");
        {
            for (int i = 0; i < width * sizeof(Sample) / 16; ++i)
            {
                int const dx = 16 * i;
                filterHorizontal16Bytes(taps, 0, dx, filteringHV);
            }

            add(r0, r1);
            add(r2, r3);
        }
        dec(r5d);
        jg("loop");
    }


    // void havoc_pred_bi_v_%1tap_16to16_%3xh_sse4(uint8_t *dst, intptr_t stride_dst, const int16_t *refAtop, const int16_t *refBtop, intptr_t stride_ref, int nPbW, int nPbH, int yFracA, int yFracB)//
    void PRED_BI_V_8NxH(int taps, int inputTypeSize, int width)
    {
        // taps is number of filter taps (4 or 8)//
        // inputTypeSize is size of input type (8 for uint8_t, 16 for int16_t right shifted 6)
        // width is block width (number of samples, multiple of 8)
        assert(inputTypeSize == 16);

        // // void havoc_pred_bi_v_%1tap_16to16_%3xh_sse4(uint8_t *dst, intptr_t stride_dst, const int16_t *refAtop, const int16_t *refBtop, intptr_t stride_ref, int nPbW, int nPbH, int yFracA, int yFracB)//
        // INIT_XMM sse4
        // cglobal pred_bi_v_%1tap_16to16_%3xh, 9, 9, 8

        auto &r0 = arg64(0); // dst
        auto &r1 = arg64(1); // stride_dst
        auto &r2 = arg64(2); // refA
        auto &r3 = arg64(3); // refB
        auto &r4 = arg64(4); // stride_ref
        auto &r5 = arg64(5); // width
        auto &r6 = arg64(6);
        auto r6d = Xbyak::Reg32(r6.getIdx()); // height
        auto &r7 = arg64(7);
        auto r7d = Xbyak::Reg32(r7.getIdx()); // yFracA
        auto &r8 = arg64(8);
        auto r8d = Xbyak::Reg32(r8.getIdx()); // yFracB

        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        auto &m4 = regXmm(4);
        auto &m5 = regXmm(5);
        auto &m6 = regXmm(6);
        auto &m7 = regXmm(7);

        shl(r4, 1);

        //		db({ 0xcc });


        shl(r7d, 4 + taps / 4);
        shl(r8d, 4 + taps / 4);
        lea(r5, ptr[rip + coefficients16]);
        lea(r7, ptr[r5 + r7]);
        lea(r8, ptr[r5 + r8]);
#define coeffA r7
#define coeffB r8
#define stride_dst r1
        L("loopBiV");
        {
            int offset = 0;
            for (int i = 0; i < (width + 7) / 8; ++i)
            {
                // each iteration of this loop operates on a row of 8-samples

                pxor(m3, m3);
                movdqa(m5, m3);
                movdqa(m6, m3);
                movdqa(m7, m3);

                for (int j = 0; j < taps / 2; ++j)
                {
                    // each iteration of this loop performs two filter taps

                    // reference picture A
                    movdqu(m0, ptr[r2]);
                    movdqu(m1, ptr[r2 + r4]);
                    movdqa(m2, m0);
                    punpckhwd(m2, m1);
                    punpcklwd(m0, m1);
                    pmaddwd(m0, ptr[coeffA]);
                    pmaddwd(m2, ptr[coeffA]);
                    paddd(m3, m0);
                    paddd(m5, m2);
                    lea(r2, ptr[r2 + r4 * 2]);

                    // reference picture B
                    movdqu(m0, ptr[r3]);
                    movdqu(m1, ptr[r3 + r4]);
                    movdqa(m2, m0);
                    punpckhwd(m2, m1);
                    punpcklwd(m0, m1);
                    pmaddwd(m0, ptr[coeffB]);
                    pmaddwd(m2, ptr[coeffB]);
                    paddd(m6, m0);
                    paddd(m7, m2);
                    lea(r3, ptr[r3 + r4 * 2]);

                    add(coeffA, 16);
                    add(coeffB, 16);
                }

                neg(r4);
                lea(r2, ptr[r2 + r4 * taps]);
                lea(r3, ptr[r3 + r4 * taps]);
                neg(r4);

                sub(coeffA, 8 * taps);
                sub(coeffB, 8 * taps);

                if (sizeof(Sample) == 1)
                {
                    psrad(m3, 6);
                    psrad(m5, 6);
                    psrad(m6, 6);
                    psrad(m7, 6);

                    packssdw(m3, m5);
                    packssdw(m6, m7);

                    paddsw(m3, m6);

                    paddsw(m3, ptr[rip + constant_times_8_dw_0x40]);
                    psraw(m3, 7);
                    packuswb(m3, m3);
                }
                else
                {
                    psrad(m3, 6);
                    psrad(m5, 6);
                    psrad(m6, 6);
                    psrad(m7, 6);

                    paddd(m3, m6);
                    paddd(m5, m7);

                    pslld(m3, 1);
                    pslld(m5, 1);

                    paddd(m3, ptr[rip + constant_times_4_dd_0x20]);
                    paddd(m5, ptr[rip + constant_times_4_dd_0x20]);
                    packusdw(m3, m5);

                    psrlw(m3, 6);
                }

                movdqu(ptr[r0], m3);

                lea(r2, ptr[r2 + 16]);
                lea(r3, ptr[r3 + 16]);
                lea(r0, ptr[r0 + 8 * sizeof(Sample)]);

                offset += 8;
            }
            sub(r2, 2 * offset);
            sub(r3, 2 * offset);
            sub(r0, offset * sizeof(Sample));

            add(r2, r4);
            add(r3, r4);
            lea(r0, ptr[r0 + stride_dst * sizeof(Sample)]);
        }
        dec(r6d);
        jg("loopBiV");
#undef coeffA
#undef coeffB
#undef stride_dst
    }

};


template <typename Sample>
void havocPopulatePredBi(HavocTablePredBi<Sample> *table, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    for (int bitDepth = 8; bitDepth <= 10; ++bitDepth)
        for (int taps = 4; taps <= 8; taps += 4)
            for (int w = 0; w <= 64; w += 2 * taps)
                for (int frac = 0; frac < 2; ++frac)
                    *havocGetPredBi<Sample>(table, taps, w, 0, frac, frac, frac, frac, bitDepth) = 0;

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
        for (int bitDepth = 8; bitDepth <= 10; ++bitDepth)
            for (int taps = 4; taps <= 8; taps += 4)
                for (int w = 0; w <= 64; w += 2 * taps)
                    for (int frac = 0; frac < 2; ++frac)
                        if (taps == 4)
                            *havocGetPredBi<Sample>(table, taps, w, 0, frac, frac, frac, frac, bitDepth) = havocPredBi_c_ref<Sample, 4>;
                        else
                            *havocGetPredBi<Sample>(table, taps, w, 0, frac, frac, frac, frac, bitDepth) = havocPredBi_c_ref<Sample, 8>;

    if (buffer.isa & HAVOC_SSE2)
        for (int width = 64; width >= 16; width -= 16)
        {
            int bitDepth = sizeof(Sample) == 2 ? 10 : 8;

            PredBi<Sample> a(&buffer, width);
            for (int taps = 4; taps <= 8; taps += 4)
                for (int w = 0; w <= 64 && w <= width; w += 2 * taps)
                    *havocGetPredBi<Sample>(table, taps, w, 0, 0, 0, 0, 0, bitDepth) = a;

            for (int taps = 4; taps <= 8; taps += 4)
            {
                PredBi<Sample> a(&buffer, width, taps);
                for (int w = 0; w <= 64 && w <= width; w += 2 * taps)
                    *havocGetPredBi<Sample>(table, taps, w, 0, 1, 1, 1, 1, bitDepth) = a;
            }
        }
}



#define STRIDE_DST 192

typedef struct
{
    HavocPredBi<uint8_t> *f8;
    HavocPredBi<uint16_t> *f16;
    HAVOC_ALIGN(32, uint8_t, dst8[64 * STRIDE_DST]);
    HAVOC_ALIGN(32, uint16_t, dst16[64 * STRIDE_DST]);
    intptr_t stride_dst;
    const uint8_t *ref8[2];
    const uint16_t *ref16[2];
    intptr_t stride_ref;
    int w;
    int h;
    int xFracA;
    int yFracA;
    int xFracB;
    int yFracB;
    int taps;
    int bitDepth;
}
bound_pred_bi;


int init_pred_bi(void *p, havoc_code code)
{
    bound_pred_bi *s = (bound_pred_bi *)p;

    s->f8 = 0;
    s->f16 = 0;

    if (s->bitDepth > 8)
    {
        HavocTablePredBi<uint16_t> table;
        havocPopulatePredBi<uint16_t>(&table, code);
        s->f16 = *havocGetPredBi<uint16_t>(&table, s->taps, s->w, s->h, s->xFracA, s->yFracA, s->xFracB, s->yFracB, s->bitDepth);
        memset(s->dst16, 0, 2 * 64 * s->stride_dst);
    }
    else
    {
        HavocTablePredBi<uint8_t> table;
        havocPopulatePredBi<uint8_t>(&table, code);
        s->f8 = *havocGetPredBi<uint8_t>(&table, s->taps, s->w, s->h, s->xFracA, s->yFracA, s->xFracB, s->yFracB, s->bitDepth);
        memset(s->dst8, 0, 64 * s->stride_dst);
    }

    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    if ((s->f8 || s->f16) && buffer.isa == HAVOC_C_REF)
    {
        printf("\t%d-bit  %d-tap %dx%d %s%s %s%s : ", s->bitDepth, s->taps, s->w, s->h, s->xFracA ? "H" : "", s->yFracA ? "V" : "", s->xFracB ? "H" : "", s->yFracB ? "V" : "");
    }

    return s->f8 || s->f16;
}


void invoke_pred_bi(void *p, int n)
{
    bound_pred_bi *s = (bound_pred_bi *)p;
    if (s->bitDepth > 8) while (n--)
    {
        s->f16(s->dst16, s->stride_dst, s->ref16[0], s->ref16[1], s->stride_ref, s->w, s->h, s->xFracA, s->yFracA, s->xFracB, s->yFracB, s->bitDepth);
    }
    else while (n--)
    {
        s->f8(s->dst8, s->stride_dst, s->ref8[0], s->ref8[1], s->stride_ref, s->w, s->h, s->xFracA, s->yFracA, s->xFracB, s->yFracB, s->bitDepth);
    }
}


int mismatch_pred_bi(void *boundRef, void *boundTest)
{
    bound_pred_bi *ref = (bound_pred_bi *)boundRef;
    bound_pred_bi *test = (bound_pred_bi *)boundTest;

    for (int y = 0; y < ref->h; ++y)
    {
        if (ref->bitDepth > 8)
        {
            if (memcmp(&ref->dst16[y*ref->stride_dst], &test->dst16[y*test->stride_dst], 2 * ref->w))
                return 1;
        }
        else
        {
            if (memcmp(&ref->dst8[y*ref->stride_dst], &test->dst8[y*test->stride_dst], ref->w))
                return 1;
        }
    }

    return 0;
}


static void test_partitions_bi(int *error_count, bound_pred_bi *b, havoc_instruction_set mask)
{
    const int partitions[22][2] =
    {
        { 8, 8 },
        { 16, 4 }, { 16, 8 }, { 16, 12 }, { 16, 16 }, { 12, 16 }, { 8, 16 }, { 4, 16 },
        { 32, 8 }, { 32, 16 }, { 32, 24 }, { 32, 32 }, { 24, 32 }, { 16, 32 }, { 8, 32 },
        { 64, 16 }, { 64, 32 }, { 64, 48 }, { 64, 64 }, { 48, 64 }, { 32, 64 }, { 16, 64 },
    };

    for (int k = 0; k < 22; ++k)
    {
        const int nPbW = partitions[k][0];
        const int nPbH = partitions[k][1];

        b[0].w = nPbW;
        b[0].h = nPbH;

        b[1] = b[0];

        *error_count += havoc_test(&b[0], &b[1], init_pred_bi, invoke_pred_bi, mismatch_pred_bi, mask, 10);
    }

    if (b[0].taps == 4)
        for (int k = 15; k < 22; ++k)
        {
            const int nPbW = partitions[k][0];
            const int nPbH = partitions[k][1];

            b[0].w = nPbW * b[0].taps / 8;
            b[0].h = nPbH * b[0].taps / 8;

            b[1] = b[0];

            *error_count += havoc_test(&b[0], &b[1], init_pred_bi, invoke_pred_bi, mismatch_pred_bi, mask, 10);
        }
}


void havoc_test_pred_bi(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_pred_bi - Bireference Inter Prediction (two-reference motion compensation)\n");

    bound_pred_bi b[2];

#define STRIDE_REF 192
    HAVOC_ALIGN(32, uint8_t, ref8[2][80 * STRIDE_REF]);
    HAVOC_ALIGN(32, uint16_t, ref16[2][80 * STRIDE_REF]);
    b[0].stride_dst = STRIDE_DST;
    b[0].stride_ref = STRIDE_REF;
#undef STRIDE_DST
#undef STRIDE_REF

    for (int x = 0; x < 80 * b[0].stride_ref; x++)
    {
        ref8[0][x] = rand() & 0xff;
        ref8[1][x] = rand() & 0xff;
        ref16[0][x] = rand() % (1 << 10);
        ref16[1][x] = rand() % (1 << 10);
    }

    b[0].ref8[0] = ref8[0] + 8 * b[0].stride_ref;
    b[0].ref8[1] = ref8[1] + 8 * b[0].stride_ref;
    b[0].ref16[0] = ref16[0] + 8 * b[0].stride_ref;
    b[0].ref16[1] = ref16[1] + 8 * b[0].stride_ref;

    for (b[0].bitDepth = 10; b[0].bitDepth >= 8; --b[0].bitDepth)
        for (b[0].taps = 8; b[0].taps >= 4; b[0].taps -= 4)
        {
            b[0].xFracA = b[0].yFracA = b[0].xFracB = b[0].yFracB = 0;
            test_partitions_bi(error_count, b, mask);

            b[0].xFracA = 1;
            b[0].yFracA = 2;
            b[0].xFracB = 3;
            b[0].yFracB = 0;
            test_partitions_bi(error_count, b, mask);
        }
}

namespace havoc {

    template <typename Sample>
    void subtractBi_c_ref(Sample *dst, intptr_t stride_dst, const Sample *pred, intptr_t stride_pred, const Sample *src, intptr_t stride_src, int nPbW, int nPbH, int bitDepth)
    {
        assert(bitDepth <= 8 * sizeof(Sample));

        int const max = (1 << bitDepth) - 1;

        for (int y = 0; y < nPbH; ++y)
        {
            for (int x = 0; x < nPbW; ++x)
            {
                int a = 2 * src[x + y * stride_src] - pred[x + y * stride_pred];
                if (a < 0) a = 0;
                if (a > max) a = max;
                dst[x + y * stride_dst] = a;
            }
        }
    }

    template <typename Sample>
    void populateSubtractBi(TableSubtractBi<Sample> *table, havoc_code code, int bitDepth)
    {
        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

        table->get() = 0;

        if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
            table->get() = subtractBi_c_ref<Sample>;
    }

    template void populateSubtractBi<uint8_t>(TableSubtractBi<uint8_t> *table, havoc_code code, int bitDepth);
    template void populateSubtractBi<uint16_t>(TableSubtractBi<uint16_t> *table, havoc_code code, int bitDepth);


    struct BoundSubtractBiBase
    {
        intptr_t stride_dst;
        intptr_t stride_pred;
        intptr_t stride_src;
        int bitDepth;
        int nPbW;
        int nPbH;
    };

    template <typename Sample>
    struct BoundSubtractBi
        :
        BoundSubtractBiBase
    {
        Sample *dst;
        Sample *pred;
        Sample *src;

        SubtractBi<Sample> *f;

        int init(havoc_code code)
        {
            TableSubtractBi<Sample> table;
            populateSubtractBi(&table, code);
            this->f = *table.get();

            memset(this->dst, 0, 64 * this->stride_dst);

            auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

            if (this->f && buffer.isa == HAVOC_C_REF)
            {
                printf("\t%d-bit %dx%d: ", this->bitDepth, this->nPbW, this->nPbH);
            }

            return !!this->f;
        }

        void invoke(int n)
        {
            while (n--)
                this->f(this->dst, this->stride_dst, this->pred, this->stride_pred, this->src, this->stride_src, this->nPbW, this->nPbH, this->bitDepth);
        }

        int mismatch(BoundSubtractBi *ref)
        {
            for (int y = 0; y < ref->nPbH; ++y)
            {
                if (memcmp(&ref->dst[y*ref->stride_dst], &this->dst[y*this->stride_dst], ref->nPbW * sizeof(Sample)))
                    return 1;
            }
            return 0;
        }
    };

    int initSubtractBi(void *p, havoc_code code)
    {
        auto *base = reinterpret_cast<BoundSubtractBiBase *>(p);
        if (base->bitDepth == 8)
        {
            auto *s = static_cast<BoundSubtractBi<uint8_t> *>(base);
            return s->init(code);
        }
        else
        {
            auto *s = static_cast<BoundSubtractBi<uint16_t> *>(base);
            return s->init(code);
        }
    }


    void invokeSubtractBi(void *p, int n)
    {
        auto *base = reinterpret_cast<BoundSubtractBiBase *>(p);
        if (base->bitDepth == 8)
        {
            auto *s = static_cast<BoundSubtractBi<uint8_t> *>(base);
            s->invoke(n);
        }
        else
        {
            auto *s = static_cast<BoundSubtractBi<uint16_t> *>(base);
            s->invoke(n);
        }
    }


    int mismatchSubtractBi(void *boundRef, void *boundTest)
    {
        auto *base = reinterpret_cast<BoundSubtractBiBase *>(boundRef);
        if (base->bitDepth == 8)
        {
            auto *ref = static_cast<BoundSubtractBi<uint8_t> *>(boundRef);
            auto *test = static_cast<BoundSubtractBi<uint8_t> *>(boundTest);
            return test->mismatch(ref);
        }
        else
        {
            auto *ref = static_cast<BoundSubtractBi<uint16_t> *>(boundRef);
            auto *test = static_cast<BoundSubtractBi<uint16_t> *>(boundTest);
            return test->mismatch(ref);
        }
    }



    template <typename Sample>
    void testSubtractBi(int *error_count, havoc_instruction_set mask)
    {
        int constexpr bitDepth = 6 + 2 * sizeof(Sample);

        printf("\nhavoc::SubtractBi %d-bit - Bireference subtraction (compute ideal second prediction)\n", bitDepth);

        BoundSubtractBi<Sample> b[2];

        size_t constexpr stride = 192;
        HAVOC_ALIGN(32, Sample, src[64 * stride]);
        HAVOC_ALIGN(32, Sample, pred[64 * stride]);
        HAVOC_ALIGN(32, Sample, dst[64 * stride]);

        for (int x = 0; x < 64 * stride; x++)
        {
            src[x] = rand() % (1 << bitDepth);
            pred[x] = rand() % (1 << bitDepth);
        }

        b[0].bitDepth = bitDepth;
        b[0].src = src;
        b[0].pred = pred;
        b[0].dst = dst;
        b[0].stride_src = stride;
        b[0].stride_pred = stride;
        b[0].stride_dst = stride;

        for (int log2 = 3; log2 <= 6; ++log2)
        {
            b[0].nPbW = 1 << log2;
            b[0].nPbH = 1 << log2;

            b[1] = b[0];

            havoc_test(&b[0], &b[1], initSubtractBi, invokeSubtractBi, mismatchSubtractBi, mask, 10);
        }
    }

    template void testSubtractBi<uint8_t>(int *error_count, havoc_instruction_set mask);
    template void testSubtractBi<uint16_t>(int *error_count, havoc_instruction_set mask);


}
