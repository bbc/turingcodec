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

#include "hadamard.h"
#include "havoc_test.h"
#include "Jit.h"
#include <stdlib.h>
#include <assert.h>




static void hadamard_iteration(int m, int n, int *dst, int *src, intptr_t stride)
{
    for (int i = 0; i < m; i += 2 * n)
        for (int j = 0; j < n; ++j)
        {
            int a = src[(i + j)*stride];
            int b = src[(i + n + j)*stride];
            dst[i + j] = a + b;
            dst[i + n + j] = a - b;
        }
}

static void hadamard_transform(int m, int n, int *dst, int *src, intptr_t stride)
{
    if (n == 1)
    {
        hadamard_iteration(m, n, dst, src, stride);
    }
    else
    {
        assert(m <= 8);
        int temp[8];
        hadamard_iteration(m, n, temp, src, stride);
        hadamard_transform(m, n / 2, dst, temp, 1);
    }
}


template <int n, typename Sample>
static int compute_satd_c_ref(const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB)
{
    assert(n <= 8);

    int intermediate[8][8];

    // subtraction and horizontal transform
    for (int y = 0; y < n; ++y)
    {
        int diff[8];
        for (int x = 0; x < n; ++x)
        {
            diff[x] = pA[x] - pB[x];
        }

        hadamard_transform(n, n / 2, intermediate[y], diff, 1);

        pA += strideA;
        pB += strideB;
    }

    // vertical transform and sum of absolutes
    const int roundingOffset = n / 4;
    int sad = roundingOffset;
    for (int x = 0; x < n; ++x)
    {
        int transformed[8];

        hadamard_transform(n, n / 2, transformed, &intermediate[0][x], 8);

        for (int y = 0; y < n; ++y)
        {
            sad += abs(transformed[y]);
        }
    }
    sad /= n / 2;
    if (sizeof(Sample) == 2)
        sad >>= 2;
    return sad;
}


#if USE_HM_DERIVED

template <typename Sample>
static int compute_satd_c_opt_2x2(const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB)
{
    int satd = 0, diff[4], m[4];

    diff[0] = pA[0] - pB[0];
    diff[1] = pA[1] - pB[1];
    diff[2] = pA[strideA] - pB[0 + strideB];
    diff[3] = pA[strideA + 1] - pB[1 + strideB];
    m[0] = diff[0] + diff[2];
    m[1] = diff[1] + diff[3];
    m[2] = diff[0] - diff[2];
    m[3] = diff[1] - diff[3];

    satd += abs(m[0] + m[1]);
    satd += abs(m[0] - m[1]);
    satd += abs(m[2] + m[3]);
    satd += abs(m[2] - m[3]);

    if (sizeof(Sample) == 2)
        satd >>= 2;

    return satd;
}


template <typename Sample>
static int compute_satd_c_opt_4x4(const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB)
{
    int k, satd = 0, diff[16], m[16], d[16];

    for (k = 0; k < 16; k += 4)
    {
        diff[k + 0] = pA[0] - pB[0];
        diff[k + 1] = pA[1] - pB[1];
        diff[k + 2] = pA[2] - pB[2];
        diff[k + 3] = pA[3] - pB[3];

        pB += strideB;
        pA += strideA;
    }

    /*===== hadamard transform =====*/
    m[0] = diff[0] + diff[12];
    m[1] = diff[1] + diff[13];
    m[2] = diff[2] + diff[14];
    m[3] = diff[3] + diff[15];
    m[4] = diff[4] + diff[8];
    m[5] = diff[5] + diff[9];
    m[6] = diff[6] + diff[10];
    m[7] = diff[7] + diff[11];
    m[8] = diff[4] - diff[8];
    m[9] = diff[5] - diff[9];
    m[10] = diff[6] - diff[10];
    m[11] = diff[7] - diff[11];
    m[12] = diff[0] - diff[12];
    m[13] = diff[1] - diff[13];
    m[14] = diff[2] - diff[14];
    m[15] = diff[3] - diff[15];

    d[0] = m[0] + m[4];
    d[1] = m[1] + m[5];
    d[2] = m[2] + m[6];
    d[3] = m[3] + m[7];
    d[4] = m[8] + m[12];
    d[5] = m[9] + m[13];
    d[6] = m[10] + m[14];
    d[7] = m[11] + m[15];
    d[8] = m[0] - m[4];
    d[9] = m[1] - m[5];
    d[10] = m[2] - m[6];
    d[11] = m[3] - m[7];
    d[12] = m[12] - m[8];
    d[13] = m[13] - m[9];
    d[14] = m[14] - m[10];
    d[15] = m[15] - m[11];

    m[0] = d[0] + d[3];
    m[1] = d[1] + d[2];
    m[2] = d[1] - d[2];
    m[3] = d[0] - d[3];
    m[4] = d[4] + d[7];
    m[5] = d[5] + d[6];
    m[6] = d[5] - d[6];
    m[7] = d[4] - d[7];
    m[8] = d[8] + d[11];
    m[9] = d[9] + d[10];
    m[10] = d[9] - d[10];
    m[11] = d[8] - d[11];
    m[12] = d[12] + d[15];
    m[13] = d[13] + d[14];
    m[14] = d[13] - d[14];
    m[15] = d[12] - d[15];

    d[0] = m[0] + m[1];
    d[1] = m[0] - m[1];
    d[2] = m[2] + m[3];
    d[3] = m[3] - m[2];
    d[4] = m[4] + m[5];
    d[5] = m[4] - m[5];
    d[6] = m[6] + m[7];
    d[7] = m[7] - m[6];
    d[8] = m[8] + m[9];
    d[9] = m[8] - m[9];
    d[10] = m[10] + m[11];
    d[11] = m[11] - m[10];
    d[12] = m[12] + m[13];
    d[13] = m[12] - m[13];
    d[14] = m[14] + m[15];
    d[15] = m[15] - m[14];

    for (k = 0; k < 16; ++k)
    {
        satd += abs(d[k]);
    }

    //satd /= 4 / 2;
    satd = ((satd + 1) >> 1);

    if (sizeof(Sample) == 2)
        satd >>= 2;
    return satd;
}


template <typename Sample>
static int compute_satd_c_opt_8x8(const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB)
{
    int k, i, j, jj, sad = 0;
    int diff[64], m1[8][8], m2[8][8], m3[8][8];

    for (k = 0; k < 64; k += 8)
    {
        diff[k + 0] = pA[0] - pB[0];
        diff[k + 1] = pA[1] - pB[1];
        diff[k + 2] = pA[2] - pB[2];
        diff[k + 3] = pA[3] - pB[3];
        diff[k + 4] = pA[4] - pB[4];
        diff[k + 5] = pA[5] - pB[5];
        diff[k + 6] = pA[6] - pB[6];
        diff[k + 7] = pA[7] - pB[7];

        pA += strideA;
        pB += strideB;
    }

    //horizontal
    for (j = 0; j < 8; j++)
    {
        jj = j << 3;
        m2[j][0] = diff[jj] + diff[jj + 4];
        m2[j][1] = diff[jj + 1] + diff[jj + 5];
        m2[j][2] = diff[jj + 2] + diff[jj + 6];
        m2[j][3] = diff[jj + 3] + diff[jj + 7];
        m2[j][4] = diff[jj] - diff[jj + 4];
        m2[j][5] = diff[jj + 1] - diff[jj + 5];
        m2[j][6] = diff[jj + 2] - diff[jj + 6];
        m2[j][7] = diff[jj + 3] - diff[jj + 7];

        m1[j][0] = m2[j][0] + m2[j][2];
        m1[j][1] = m2[j][1] + m2[j][3];
        m1[j][2] = m2[j][0] - m2[j][2];
        m1[j][3] = m2[j][1] - m2[j][3];
        m1[j][4] = m2[j][4] + m2[j][6];
        m1[j][5] = m2[j][5] + m2[j][7];
        m1[j][6] = m2[j][4] - m2[j][6];
        m1[j][7] = m2[j][5] - m2[j][7];

        m2[j][0] = m1[j][0] + m1[j][1];
        m2[j][1] = m1[j][0] - m1[j][1];
        m2[j][2] = m1[j][2] + m1[j][3];
        m2[j][3] = m1[j][2] - m1[j][3];
        m2[j][4] = m1[j][4] + m1[j][5];
        m2[j][5] = m1[j][4] - m1[j][5];
        m2[j][6] = m1[j][6] + m1[j][7];
        m2[j][7] = m1[j][6] - m1[j][7];
    }

    //vertical
    for (i = 0; i < 8; i++)
    {
        m3[0][i] = m2[0][i] + m2[4][i];
        m3[1][i] = m2[1][i] + m2[5][i];
        m3[2][i] = m2[2][i] + m2[6][i];
        m3[3][i] = m2[3][i] + m2[7][i];
        m3[4][i] = m2[0][i] - m2[4][i];
        m3[5][i] = m2[1][i] - m2[5][i];
        m3[6][i] = m2[2][i] - m2[6][i];
        m3[7][i] = m2[3][i] - m2[7][i];

        m1[0][i] = m3[0][i] + m3[2][i];
        m1[1][i] = m3[1][i] + m3[3][i];
        m1[2][i] = m3[0][i] - m3[2][i];
        m1[3][i] = m3[1][i] - m3[3][i];
        m1[4][i] = m3[4][i] + m3[6][i];
        m1[5][i] = m3[5][i] + m3[7][i];
        m1[6][i] = m3[4][i] - m3[6][i];
        m1[7][i] = m3[5][i] - m3[7][i];

        m2[0][i] = m1[0][i] + m1[1][i];
        m2[1][i] = m1[0][i] - m1[1][i];
        m2[2][i] = m1[2][i] + m1[3][i];
        m2[3][i] = m1[2][i] - m1[3][i];
        m2[4][i] = m1[4][i] + m1[5][i];
        m2[5][i] = m1[4][i] - m1[5][i];
        m2[6][i] = m1[6][i] + m1[7][i];
        m2[7][i] = m1[6][i] - m1[7][i];
    }

    for (i = 0; i < 8; i++)
    {
        for (j = 0; j < 8; j++)
        {
            sad += abs(m2[i][j]);
        }
    }
    sad = ((sad + 2) >> 2);
    if (sizeof(Sample) == 2)
        sad >>= 2;

    return sad;
}

#endif


#define ORDER(a, b, c, d) ((a << 6) | (b << 4) | (c << 2) | d)




#define BUTTERFLY_HORIZONTAL_4(a1, a2) \
        pshufd(a2, a1, ORDER(1, 0, 3, 2)); \
        pxor(a1, m6); \
        psubw(a1, m6); \
        paddw(a1, a2);

#define BUTTERFLY_HORIZONTAL_2(a1, a2) \
        pshufd(a2, a1, ORDER(2, 3, 0, 1)); \
        pxor(a1, m6); \
        psubw(a1, m6); \
        paddw(a1, a2);

#define BUTTERFLY_HORIZONTAL_1(a1, a2) \
        movdqa(a2, a1); \
        phaddw(a2, a1); \
        phsubw(a1, a1); \
        punpcklwd(a1, a2);



template <typename Sample>
struct Satd4
    :
    Jit::Function
{
    Satd4(Jit::Buffer *buffer)
        :
        Jit::Function(buffer, Jit::CountArguments<havoc_hadamard_satd<Sample>>::value)
    {
        this->build();
    }

    Xbyak::Label constant_ffffffff00000000ffffffff00000000;

    void data()
    {
        align();

        L(constant_ffffffff00000000ffffffff00000000);
        dd({ -1, 0 }, 2);
    }

    void loadDiff(Xbyak::Xmm const &dst, Xbyak::Address srcA, Xbyak::Address srcB, Xbyak::Xmm const &temp, Xbyak::Xmm const &zero)
    {
        movdqu(dst, srcA);
        movdqu(temp, srcB);
        if (sizeof(Sample) == 1)
        {
            punpcklbw(dst, zero);
            punpcklbw(temp, zero);
        }
        psubw(dst, temp);
    }

    void assemble()
    {
        auto &m0 = regXmm(0);
        auto &m1 = regXmm(1);
        auto &m2 = regXmm(2);
        auto &m3 = regXmm(3);
        auto &m4 = regXmm(4);
        auto &m5 = regXmm(5);
        auto &m6 = regXmm(6);
        auto &m7 = regXmm(7);

        auto &r0 = arg64(0);
        auto &r1 = arg64(1);
        auto &r2 = arg64(2);
        auto &r3 = arg64(3);

        pxor(m7, m7);

        loadDiff(m0, ptr[r0], ptr[r2], m5, m7);
        loadDiff(m1, ptr[r0 + r1 * sizeof(Sample)], ptr[r2 + r3 * sizeof(Sample)], m5, m7);
        lea(r0, ptr[r0 + r1 * (2 * sizeof(Sample))]);
        lea(r2, ptr[r2 + r3 * (2 * sizeof(Sample))]);
        loadDiff(m2, ptr[r0], ptr[r2], m5, m7);
        loadDiff(m3, ptr[r0 + r1 * sizeof(Sample)], ptr[r2 + r3 * sizeof(Sample)], m5, m7);

        punpcklqdq(m0, m2);
        punpcklqdq(m1, m3);

        // diff in m0, m1, m2, m3

        movaps(m6, ptr[rip + constant_ffffffff00000000ffffffff00000000]);
        BUTTERFLY_HORIZONTAL_2(m0, m4);
        BUTTERFLY_HORIZONTAL_2(m1, m4);

        BUTTERFLY_HORIZONTAL_1(m0, m4);
        BUTTERFLY_HORIZONTAL_1(m1, m4);

        // rows in m0, m1, m2, m3

        // vertical butterfly 2
        // mova m6, [constant_ffffffffffffffff0000000000000000]
        pshufd(m6, m6, ORDER(3, 1, 2, 0));

        BUTTERFLY_HORIZONTAL_4(m0, m4);
        BUTTERFLY_HORIZONTAL_4(m1, m4);

        // rows in m0, m1

        // vertical butterfly 1
        movdqa(m2, m0);
        psubw(m2, m1);
        paddw(m0, m1);

        // rows in m0, m2

        pabsw(m0, m0);
        pabsw(m2, m2);

        paddw(m0, m2);
        phaddw(m0, m0);
        phaddw(m0, m0);
        phaddw(m0, m0);

        movd(eax, m0);
        and (eax, 0xffff);
        add(eax, 1);
        shr(eax, 1 + (sizeof(Sample) == 2 ? 2 : 0));
    }
};


#define	PUNPCKHDQQQ(a1, a2, a3) \
        vperm2i128(m ## a1, m ## a2, m ## a3, 0x31);


#define	PUNPCKLDQQQ(a1, a2, a3) \
        vinserti128(m ## a1, m ## a2, xmm ## a3, 0x1);


template <typename Sample>
struct Satd8
    :
    Jit::Function
{
    Satd8(Jit::Buffer *buffer)
        :
        Jit::Function(buffer, Jit::CountArguments<havoc_hadamard_satd<Sample>>::value)
    {
        this->build();
    }

    Xbyak::Label constant_010101010101010101ff01ff01ff01ff010101010101010101ff01ff01ff01ff;
    Xbyak::Label constant_times_16_dw_1;

    void data()
    {
        align(32);
        L(constant_010101010101010101ff01ff01ff01ff010101010101010101ff01ff01ff01ff);
        db({ 1, 1 }, 4);
        db({ 1, -1 }, 4);
        db({ 1, 1 }, 4);
        db({ 1, -1 }, 4);

        L(constant_times_16_dw_1);
        dw({ 1 }, 16);
    }

    void hadamardHorizontal8(int a1, int a2)
    {
        auto &r0 = arg64(0);
        auto &r1 = arg64(1);
        auto &r2 = arg64(2);
        auto &r3 = arg64(3);

        if (sizeof(Sample) == 1)
        {
            vmovq(Xbyak::Xmm(a1), ptr[r0]);  /*srcA(7..0, y) */
            vmovq(Xbyak::Xmm(a2), ptr[r2]);  /* srcB(7..0, y) */
            vmovq(xmm10, ptr[r0 + r1]);  /* srcA(7..0, y+2) */
            vmovq(xmm11, ptr[r2 + r3]);  /* srcB(7..0, y+2) */

            lea(r0, ptr[r0 + r1*(2 * sizeof(Sample))]);
            lea(r2, ptr[r2 + r3*(2 * sizeof(Sample))]);

            vinserti128(Xbyak::Ymm(a1), Xbyak::Ymm(a1), ptr[r0], 1);
            vinserti128(Xbyak::Ymm(a2), Xbyak::Ymm(a2), ptr[r2], 1);
            vinserti128(ymm10, ymm10, ptr[r0 + r1 * sizeof(Sample)], 1);
            vinserti128(ymm11, ymm11, ptr[r2 + r3 * sizeof(Sample)], 1);

            vpunpcklqdq(Xbyak::Ymm(a1), Xbyak::Ymm(a1));  /* b  srcA(7..0, y+1),  srcA(7..0, y+1),  srcA(7..0, y+0),  srcA(7..0, y+0) */
            vpunpcklqdq(Xbyak::Ymm(a2), Xbyak::Ymm(a2));  /* b  srcB(7..0, y+1),  srcB(7..0, y+1),  srcB(7..0, y+0),  srcB(7..0, y+0) */
            vpunpcklqdq(ymm10, ymm10);  /* b  srcA(7..0, y+3),  srcA(7..0, y+3),  srcA(7..0, y+2),  srcA(7..0, y+2) */
            vpunpcklqdq(ymm11, ymm11);  /* b  srcB(7..0, y+3),  srcB(7..0, y+3),  srcB(7..0, y+2),  srcB(7..0, y+2) */

            vpmaddubsw(Xbyak::Ymm(a1), ymm15); /* w  srcA(7,y+1)+srcA(6,y+1)..srcA(1,y+1)+srcA(0,y+1),  srcA(7,y+1)-srcA(6,y+1)..srcA(1,y+1)-srcA(0,y+1),  srcA(7,y+0)+srcA(6,y+0)..srcA(1,y+0)+srcA(0,y+0),  srcA(7,y+0)-srcA(6,y+0)..srcA(1,y+0)-srcA(0,y+0),  */
            vpmaddubsw(Xbyak::Ymm(a2), ymm15); /* w  srcB(7,y+1)+srcB(6,y+1)..srcB(1,y+1)+srcB(0,y+1),  srcB(7,y+1)-srcB(6,y+1)..srcB(1,y+1)-srcB(0,y+1),  srcB(7,y+0)+srcB(6,y+0)..srcB(1,y+0)+srcB(0,y+0),  srcB(7,y+0)-srcB(6,y+0)..srcB(1,y+0)-srcB(0,y+0), */
            vpmaddubsw(ymm10, ymm15); /* w  srcA(7,y+3)+srcA(6,y+3)..srcA(1,y+3)+srcA(0,y+3),  srcA(7,y+3)-srcA(6,y+3)..srcA(1,y+3)-srcA(0,y+3),  srcA(7,y+2)+srcA(6,y+2)..srcA(1,y+2)+srcA(0,y+2),  srcA(7,y+2)-srcA(6,y+2)..srcA(1,y+2)-srcA(0,y+2), */
            vpmaddubsw(ymm11, ymm15); /* w  srcB(7,y+3)+srcB(6,y+3)..srcB(1,y+3)+srcB(0,y+3),  srcB(7,y+3)-srcB(6,y+3)..srcB(1,y+3)-srcB(0,y+3),  srcB(7,y+2)+srcB(6,y+2)..srcB(1,y+2)+srcB(0,y+2),  srcB(7,y+2)-srcB(6,y+2)..srcB(1,y+2)-srcB(0,y+2), */

            vpsubw(Xbyak::Ymm(a1), Xbyak::Ymm(a2)); /* w src(7,y+1)+src(6,y+1)..src(1,y+1)+src(0,y+1), src(7,y+1)-src(6,y+1)..src(1,y+1)-src(0,y+1), src(7,y+0)+src(6,y+0)..src(1,y+0)+src(0,y+0), src(7,y+0)-src(6,y+0)..src(1,y+0)-src(0,y+0),  */
            vpsubw(ymm10, ymm11); /* w src(7,y+3)+src(6,y+3)..src(1,y+3)+src(0,y+3), src(7,y+3)-src(6,y+3)..src(1,y+3)-src(0,y+3), src(7,y+2)+src(6,y+2)..src(1,y+2)+src(0,y+2), src(7,y+2)-src(6,y+2)..src(1,y+2)-src(0,y+2),  */
        }
        else
        {
            vmovdqu(Xbyak::Xmm(a1), ptr[r0]);  /*srcA(7..0, y) */
            vmovdqu(Xbyak::Xmm(a2), ptr[r2]);  /* srcB(7..0, y) */
            vmovdqu(xmm10, ptr[r0 + r1 * 2]);  /* srcA(7..0, y+2) */
            vmovdqu(xmm11, ptr[r2 + r3 * 2]);  /* srcB(7..0, y+2) */

            lea(r0, ptr[r0 + r1*(2 * 2)]);
            lea(r2, ptr[r2 + r3*(2 * 2)]);

            vinserti128(Xbyak::Ymm(a1), Xbyak::Ymm(a1), ptr[r0], 1);
            vinserti128(Xbyak::Ymm(a2), Xbyak::Ymm(a2), ptr[r2], 1);
            vinserti128(ymm10, ymm10, ptr[r0 + r1 * 2], 1);
            vinserti128(ymm11, ymm11, ptr[r2 + r3 * 2], 1);

            vpsubw(Xbyak::Ymm(a1), Xbyak::Ymm(a2));
            vpsubw(ymm10, ymm11);

            // 76543210 76543210
            // -+-+-+-+-+ ++++++++

            vphsubw(ymm7, Xbyak::Ymm(a1), ymm10); // 7-6 5-4 3-2 1-0
            vphaddw(ymm6, Xbyak::Ymm(a1), ymm10); // 7+6 5+4 3+2 1+0

            vpunpcklqdq(Xbyak::Ymm(a1), ymm6, ymm7); // 7+6 5+4 3+2 1+0 7-6 5-4 3-2 1-0
            vpunpckhqdq(ymm10, ymm6, ymm7); // 7+6 5+4 3+2 1+0 7-6 5-4 3-2 1-0

            // a1; /* w src(7,y+1)+src(6,y+1)..src(1,y+1)+src(0,y+1), src(7,y+1)-src(6,y+1)..src(1,y+1)-src(0,y+1), src(7,y+0)+src(6,y+0)..src(1,y+0)+src(0,y+0), src(7,y+0)-src(6,y+0)..src(1,y+0)-src(0,y+0),  */
            // ymm10; /* w src(7,y+3)+src(6,y+3)..src(1,y+3)+src(0,y+3), src(7,y+3)-src(6,y+3)..src(1,y+3)-src(0,y+3), src(7,y+2)+src(6,y+2)..src(1,y+2)+src(0,y+2), src(7,y+2)-src(6,y+2)..src(1,y+2)-src(0,y+2),  */
        }

        vphsubw(ymm12, Xbyak::Ymm(a1), ymm10);
        vphaddw(ymm13, Xbyak::Ymm(a1), ymm10);

        vpunpckldq(ymm4, ymm12, ymm13);
        vpunpckhdq(ymm5, ymm12, ymm13);

        vphsubw(ymm6, ymm4, ymm4);
        vphaddw(ymm7, ymm4, ymm4);
        vphsubw(ymm8, ymm5, ymm5);
        vphaddw(ymm9, ymm5, ymm5);

        vpunpcklwd(Xbyak::Ymm(a1), ymm6, ymm7);
        vpunpcklwd(Xbyak::Ymm(a2), ymm8, ymm9);
    }

    void assemble()
    {
        auto &m0 = ymm0;
        auto &m1 = ymm1;
        auto &m2 = ymm2;
        auto &m3 = ymm3;
        auto &m4 = ymm4;
        auto &m5 = ymm5;
        auto &m6 = ymm6;
        auto &m7 = ymm7;
        auto &m8 = ymm8;
        auto &m9 = ymm9;
        auto &m10 = ymm10;
        auto &m11 = ymm11;
        auto &m12 = ymm12;
        auto &m13 = ymm13;
        auto &m14 = ymm14;
        auto &m15 = ymm15;

        regXmm(15);

        auto &r0 = arg64(0);
        auto &r1 = arg64(1);
        auto &r2 = arg64(2);
        auto &r3 = arg64(3);

        if (sizeof(Sample) == 1)
            vmovdqa(m15, ptr[rip + constant_010101010101010101ff01ff01ff01ff010101010101010101ff01ff01ff01ff]);

        hadamardHorizontal8(0, 1);

        lea(r0, ptr[r0 + r1 * (2 * sizeof(Sample))]);
        lea(r2, ptr[r2 + r3 * (2 * sizeof(Sample))]);
        hadamardHorizontal8(2, 3);

        // horizontal transform now done - output order is incorrect but same in each row so OK for SATD

        // (x, 0) in m0 low
        // (x, 1) in m1 low
        // (x, 2) in m0 high
        // (x, 3) in m1 high
        // (x, 4) in m2 low
        // (x, 5) in m3 low
        // (x, 6) in m2 high
        // (x, 7) in m3 high

// 8-bit input means 12-bit values at this point (including sign bit)
// 10-bit input means 14-bit values at this point (including sign bit)

// m0  H = 00+00000 L = +0000000
// m1  H = 000+0000 L = 0+000000
// m2  H = 000000+0 L = 0000+000
// m3  H = 0000000+ L = 00000+00

        vpaddw(m4, m0, m2); // H = 00+000+0 L = +000+000
        vpaddw(m5, m1, m3); // H = 000+000+ L = 0+000+00
        vpsubw(m6, m0, m2); // H = 00+000-0 L = +000-000
        vpsubw(m7, m1, m3); // H = 000+000- L = 0+000-00

        // 8-bit input means 13-bit values at this point (including sign bit)
        // 10-bit input means 15-bit values at this point (including sign bit)

        vpaddw(m0, m4, m5); // H = 00++00++L = ++00++00
        vpsubw(m1, m4, m5); // H = 00+-00+-L = +-00+-00
        vpaddw(m2, m6, m7); // H = 00++00--L = ++00--00
        vpsubw(m3, m6, m7); // H = 00+-00-+L = +-00-+00

        // 8-bit input means 14-bit values at this point (including sign bit)
        // 10-bit input means 16-bit values at this point (including sign bit)

        PUNPCKHDQQQ(4, 0, 1); // H = 00+-00+- L = 00++00++
        PUNPCKLDQQQ(5, 0, 1); // H = +-00+-00 L = ++00++00
        PUNPCKHDQQQ(6, 2, 3); // H = 00+-00-+ L = 00++00--
        PUNPCKLDQQQ(7, 2, 3); // L = +-00-+00 L = ++00--00

        if (sizeof(Sample) == 2)
        {
            vextracti128(xm8, m4, 1);
            vextracti128(xm9, m5, 1);
            vextracti128(xm10, m6, 1);
            vextracti128(xm11, m7, 1);

            vpmovsxwd(m4, m4);
            vpmovsxwd(m5, m5);
            vpmovsxwd(m6, m6);
            vpmovsxwd(m7, m7);
            vpmovsxwd(m8, m8);
            vpmovsxwd(m9, m9);
            vpmovsxwd(m10, m10);
            vpmovsxwd(m11, m11);

            vpaddd(m0, m4, m5); // H = +-+-+-+- L = ++++++++
            vpsubd(m1, m5, m4); // H = +--++--+ L = ++--++--
            vpaddd(m2, m6, m7); // H = +-+--+-+ L = ++++----
            vpsubd(m3, m7, m6); // H = +--+-++- L = ++----++
            vpaddd(m4, m8, m9); // H = +-+-+-+- L = ++++++++
            vpsubd(m5, m9, m8); // H = +--++--+ L = ++--++--
            vpaddd(m6, m10, m11); // H = +-+--+-+ L = ++++----
            vpsubd(m7, m11, m10); // H = +--+-++- L = ++----++
            // 10-bit input means 17-bit values at this point (including sign bit)

            // vertical transform now done too.

            vpabsd(m0, m0);
            vpabsd(m1, m1);
            vpabsd(m2, m2);
            vpabsd(m3, m3);
            vpabsd(m4, m4);
            vpabsd(m5, m5);
            vpabsd(m6, m6);
            vpabsd(m7, m7);
            // 10-bit input means 16-bit unsigned values at this point

            vpaddd(m0, m1);
            vpaddd(m2, m3);
            vpaddd(m4, m5);
            vpaddd(m6, m7);
            // 8-bit input means 17-bit unsigned values at this point

            vpaddd(m0, m2);
            vpaddd(m4, m6);
            // 8-bit input means 18-bit unsigned values at this point

            vpaddd(m0, m4);
            // 8-bit input means 19-bit unsigned values at this point
        }
        else
        {
            vpaddw(m0, m4, m5); // H = +-+-+-+- L = ++++++++
            vpsubw(m1, m5, m4); // H = +--++--+ L = ++--++--
            vpaddw(m2, m6, m7); // H = +-+--+-+ L = ++++----
            vpsubw(m3, m7, m6); // H = +--+-++- L = ++----++
            // 8-bit input means 15-bit values at this point (including sign bit)

                                        // vertical transform now done too.

            vpabsw(m0, m0);
            vpabsw(m1, m1);
            vpabsw(m2, m2);
            vpabsw(m3, m3);

            // 8-bit input means 14-bit unsigned values at this point

                                        // transformed absolute differences now computed, just need to sum them and return

            vpaddw(m0, m1);
            vpaddw(m2, m3);
            // 8-bit input means 15-bit unsigned values at this point

            vpaddw(m0, m2);
            // 8-bit input means 16-bit unsigned values at this point

            vpmaddwd(m0, m0, ptr[rip + constant_times_16_dw_1]);
            // 8-bit input means 17-bit unsigned values at this point
        }

        vextracti128(xmm1, m0, 1);

        vzeroupper();

        paddd(xmm0, xmm1);
        movhlps(xmm1, xmm0);
        paddd(xmm0, xmm1);
        pshuflw(xmm1, xmm0, 0xe);
        paddd(xmm0, xmm1);
        movd(eax, xmm0);
        add(eax, 2);
        shr(eax, 2 + (sizeof(Sample) == 2 ? 2 : 0));
    }
};


template <typename Sample>
void havoc_populate_hadamard_satd(havoc_table_hadamard_satd<Sample> *table, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    *havoc_get_hadamard_satd(table, 1) = 0;
    *havoc_get_hadamard_satd(table, 2) = 0;
    *havoc_get_hadamard_satd(table, 3) = 0;

    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
    {
        *havoc_get_hadamard_satd(table, 1) = compute_satd_c_ref<2, Sample>;
        *havoc_get_hadamard_satd(table, 2) = compute_satd_c_ref<4, Sample>;
        *havoc_get_hadamard_satd(table, 3) = compute_satd_c_ref<8, Sample>;
    }

#if USE_HM_DERIVED
    if (buffer.isa & HAVOC_C_OPT)
    {
        *havoc_get_hadamard_satd(table, 1) = compute_satd_c_opt_2x2<Sample>;
        *havoc_get_hadamard_satd(table, 2) = compute_satd_c_opt_4x4<Sample>;
        *havoc_get_hadamard_satd(table, 3) = compute_satd_c_opt_8x8<Sample>;
    }
#endif

#ifdef HAVOC_X64
    if (buffer.isa & HAVOC_SSE2)
    {
        Satd4<Sample> satd4(&buffer);
        *havoc_get_hadamard_satd(table, 2) = satd4;
    }
    if (buffer.isa & HAVOC_AVX2)
    {
        Satd8<Sample> satd8(&buffer);
        *havoc_get_hadamard_satd(table, 3) = satd8;
    }
#endif
}


struct BoundHadamardSatdBase
{
    int log2TrafoSize;
    int satd;
    int bits;
};

template <typename Sample>
struct bound_hadamard_satd
    :
    BoundHadamardSatdBase
{
    Sample *srcA, *srcB;
    havoc_hadamard_satd<Sample> *f;

    int init(havoc_code code)
    {
        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

        auto s = this;

        havoc_table_hadamard_satd<Sample> table;

        havoc_populate_hadamard_satd(&table, code);

        s->f = *havoc_get_hadamard_satd(&table, s->log2TrafoSize);

        if (buffer.isa == HAVOC_C_REF)
        {
            const int nCbS = 1 << s->log2TrafoSize;
            printf("\t%d bits %dx%d : ", s->bits, nCbS, nCbS);
        }

        return !!s->f;
    }

    void invoke(int n)
    {
        auto s = this;

        const int nCbS = 1 << s->log2TrafoSize;
        while (n--)
        {
            s->satd = s->f(s->srcA, 2 * nCbS, s->srcB, 2 * nCbS);
        }
    }
};

int init_hadamard_satd(void *p, havoc_code code)
{
    BoundHadamardSatdBase *b = (BoundHadamardSatdBase *)p;
    if (b->bits == 8)
        return static_cast<bound_hadamard_satd<uint8_t> *>(b)->init(code);
    else
        return static_cast<bound_hadamard_satd<uint16_t> *>(b)->init(code);
}


void invoke_hadamard_satd(void *p, int n)
{
    BoundHadamardSatdBase *b = (BoundHadamardSatdBase *)p;
    if (b->bits == 8)
        return static_cast<bound_hadamard_satd<uint8_t> *>(b)->invoke(n);
    else
        return static_cast<bound_hadamard_satd<uint16_t> *>(b)->invoke(n);
}


int mismatch_hadamard_satd(void *boundRef, void *boundTest)
{
    BoundHadamardSatdBase *ref = (BoundHadamardSatdBase *)boundRef;
    BoundHadamardSatdBase *test = (BoundHadamardSatdBase *)boundTest;

    return ref->satd != test->satd;
}


template <typename Sample>
void testHadamardSatd(int *error_count, havoc_instruction_set mask)
{
    auto constexpr bits = sizeof(Sample) == 2 ? 10 : 8;

    Sample srcA[16 * 8];
    Sample srcB[16 * 8];

    for (int i = 0; i < 16 * 8; ++i)
    {
        srcA[i] = (1 << bits) - 1 - (rand() & 7);
        srcB[i] = rand() & 7;
    }

    bound_hadamard_satd<Sample> b[2];

    b[0].srcA = srcA;
    b[0].srcB = srcB;
    b[0].bits = bits;

    for (b[0].log2TrafoSize = 3; b[0].log2TrafoSize >= 1; --b[0].log2TrafoSize)
    {
        b[1] = b[0];
        *error_count += havoc_test(&b[0], &b[1], init_hadamard_satd, invoke_hadamard_satd, mismatch_hadamard_satd, mask, 10);
    }
}


void havoc_test_hadamard_satd(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_hadamard_satd - Hadamard Sum of Absolute Transformed Differences\n");
    testHadamardSatd<uint16_t>(error_count, mask);
    testHadamardSatd<uint8_t>(error_count, mask);
}

