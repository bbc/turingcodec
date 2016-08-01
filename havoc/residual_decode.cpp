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
//
// Portions marked "f265" derived from the f265 project also
// governed by a BSD-style license included in the COPYING file.

#include "residual_decode.h"
#include "havoc_test.h"
#include "Jit.h"
#include <stdlib.h>
#include <string.h>


#ifdef WIN32
#define FASTCALL __fastcall
#else
#define FASTCALL
#endif


static int Clip3(int min, int max, int x)
{
    if (x > max) return max;
    if (x < min) return min;
    return x;
}


static void havoc_inverse_partial_butterfly_4x4_dst_c_opt(int16_t dst[4 * 4], const int16_t src[4 * 4], int shift)
{
    const int add = 1 << (shift - 1);
    const int src_stride = 4;
    const int dst_stride = 4;

    for (int i = 0; i < 4; i++)
    {
        int c[4];
        c[0] = src[i] + src[2 * src_stride + i];
        c[1] = src[2 * src_stride + i] + src[3 * src_stride + i];
        c[2] = src[i] - src[3 * src_stride + i];
        c[3] = 74 * src[1 * src_stride + i];

        dst[dst_stride*i + 0] = Clip3(-32768, 32767, (29 * c[0] + 55 * c[1] + c[3] + add) >> shift);
        dst[dst_stride*i + 1] = Clip3(-32768, 32767, (55 * c[2] - 29 * c[1] + c[3] + add) >> shift);
        dst[dst_stride*i + 2] = Clip3(-32768, 32767, (74 * (src[i] - src[2 * src_stride + i] + src[3 * src_stride + i]) + add) >> shift);
        dst[dst_stride*i + 3] = Clip3(-32768, 32767, (55 * c[0] + 29 * c[2] - c[3] + add) >> shift);
    }
}


template <int nTbS>
void havoc_inverse_partial_butterfly_c_opt(int16_t dst[], int16_t const src[], int shift);


template <>
void havoc_inverse_partial_butterfly_c_opt<4>(int16_t dst[], const int16_t src[], int shift)
{
    const int add = 1 << (shift - 1);
    const int src_stride = 4;
    const int dst_stride = 4;

    for (int j = 0; j<4; ++j)
    {
        static const int16_t table[4][4] =
        {
                { 64, 64, 64, 64 },
                { 83, 36, -36, -83 },
                { 64, -64, -64, 64 },
                { 36, -83, 83, -36 }
        };

        int E[2], O[2];
        O[0] = table[1][0] * src[1 * src_stride] + table[3][0] * src[3 * src_stride];
        O[1] = table[1][1] * src[1 * src_stride] + table[3][1] * src[3 * src_stride];
        E[0] = table[0][0] * src[0 * src_stride] + table[2][0] * src[2 * src_stride];
        E[1] = table[0][1] * src[0 * src_stride] + table[2][1] * src[2 * src_stride];

        dst[0] = Clip3(-32768, 32767, (E[0] + O[0] + add) >> shift);
        dst[1] = Clip3(-32768, 32767, (E[1] + O[1] + add) >> shift);
        dst[2] = Clip3(-32768, 32767, (E[1] - O[1] + add) >> shift);
        dst[3] = Clip3(-32768, 32767, (E[0] - O[0] + add) >> shift);

        ++src;
        dst += dst_stride;
    }
}


template <>
void havoc_inverse_partial_butterfly_c_opt<8>(int16_t dst[], const int16_t src[], int shift)
{
    const int add = 1 << (shift - 1);
    const int src_stride = 8;
    const int dst_stride = 8;

    for (int j = 0; j<8; ++j)
    {
        static const int16_t table[8][8] =
        {
                { 64, 64, 64, 64, 64, 64, 64, 64 },
                { 89, 75, 50, 18, -18, -50, -75, -89 },
                { 83, 36, -36, -83, -83, -36, 36, 83 },
                { 75, -18, -89, -50, 50, 89, 18, -75 },
                { 64, -64, -64, 64, 64, -64, -64, 64 },
                { 50, -89, 18, 75, -75, -18, 89, -50 },
                { 36, -83, 83, -36, -36, 83, -83, 36 },
                { 18, -50, 75, -89, 89, -75, 50, -18 }
        };

        int O[4];
        for (int k = 0; k<4; ++k)
        {
            O[k] = table[1][k] * src[1 * src_stride] + table[3][k] * src[3 * src_stride] + table[5][k] * src[5 * src_stride] + table[7][k] * src[7 * src_stride];
        }

        int EE[2], EO[2];
        EO[0] = table[2][0] * src[2 * src_stride] + table[6][0] * src[6 * src_stride];
        EO[1] = table[2][1] * src[2 * src_stride] + table[6][1] * src[6 * src_stride];
        EE[0] = table[0][0] * src[0 * src_stride] + table[4][0] * src[4 * src_stride];
        EE[1] = table[0][1] * src[0 * src_stride] + table[4][1] * src[4 * src_stride];

        int E[4];
        E[0] = EE[0] + EO[0];
        E[3] = EE[0] - EO[0];
        E[1] = EE[1] + EO[1];
        E[2] = EE[1] - EO[1];

        for (int k = 0; k<4; ++k)
        {
            dst[k] = Clip3(-32768, 32767, (E[k] + O[k] + add) >> shift);
            dst[k + 4] = Clip3(-32768, 32767, (E[3 - k] - O[3 - k] + add) >> shift);
        }

        src++;
        dst += dst_stride;
    }
}


template <>
void havoc_inverse_partial_butterfly_c_opt<16>(int16_t dst[], const int16_t src[], int shift)
{
    const int add = 1 << (shift - 1);
    const int src_stride = 16;
    const int dst_stride = 16;

    for (int j = 0; j<16; ++j)
    {
        static const int16_t table[16][16] =
        {
                { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
                { 90, 87, 80, 70, 57, 43, 25, 9, -9, -25, -43, -57, -70, -80, -87, -90 },
                { 89, 75, 50, 18, -18, -50, -75, -89, -89, -75, -50, -18, 18, 50, 75, 89 },
                { 87, 57, 9, -43, -80, -90, -70, -25, 25, 70, 90, 80, 43, -9, -57, -87 },
                { 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83 },
                { 80, 9, -70, -87, -25, 57, 90, 43, -43, -90, -57, 25, 87, 70, -9, -80 },
                { 75, -18, -89, -50, 50, 89, 18, -75, -75, 18, 89, 50, -50, -89, -18, 75 },
                { 70, -43, -87, 9, 90, 25, -80, -57, 57, 80, -25, -90, -9, 87, 43, -70 },
                { 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64 },
                { 57, -80, -25, 90, -9, -87, 43, 70, -70, -43, 87, 9, -90, 25, 80, -57 },
                { 50, -89, 18, 75, -75, -18, 89, -50, -50, 89, -18, -75, 75, 18, -89, 50 },
                { 43, -90, 57, 25, -87, 70, 9, -80, 80, -9, -70, 87, -25, -57, 90, -43 },
                { 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36 },
                { 25, -70, 90, -80, 43, 9, -57, 87, -87, 57, -9, -43, 80, -90, 70, -25 },
                { 18, -50, 75, -89, 89, -75, 50, -18, -18, 50, -75, 89, -89, 75, -50, 18 },
                { 9, -25, 43, -57, 70, -80, 87, -90, 90, -87, 80, -70, 57, -43, 25, -9 }
        };

        int O[8];
        for (int k = 0; k<8; ++k)
        {
            O[k] = table[1][k] * src[src_stride] + table[3][k] * src[3 * src_stride] + table[5][k] * src[5 * src_stride] + table[7][k] * src[7 * src_stride] +
                    table[9][k] * src[9 * src_stride] + table[11][k] * src[11 * src_stride] + table[13][k] * src[13 * src_stride] + table[15][k] * src[15 * src_stride];
        }

        int EO[4];
        for (int k = 0; k<4; ++k)
        {
            EO[k] = table[2][k] * src[2 * src_stride] + table[6][k] * src[6 * src_stride] + table[10][k] * src[10 * src_stride] + table[14][k] * src[14 * src_stride];
        }

        int EEE[2], EEO[2];
        EEO[0] = table[4][0] * src[4 * src_stride] + table[12][0] * src[12 * src_stride];
        EEE[0] = table[0][0] * src[0 * src_stride] + table[8][0] * src[8 * src_stride];
        EEO[1] = table[4][1] * src[4 * src_stride] + table[12][1] * src[12 * src_stride];
        EEE[1] = table[0][1] * src[0 * src_stride] + table[8][1] * src[8 * src_stride];

        int EE[4];
        for (int k = 0; k<2; ++k)
        {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 2] = EEE[1 - k] - EEO[1 - k];
        }

        int E[8];
        for (int k = 0; k<4; ++k)
        {
            E[k] = EE[k] + EO[k];
            E[k + 4] = EE[3 - k] - EO[3 - k];
        }

        for (int k = 0; k<8; ++k)
        {
            dst[k] = Clip3(-32768, 32767, (E[k] + O[k] + add) >> shift);
            dst[k + 8] = Clip3(-32768, 32767, (E[7 - k] - O[7 - k] + add) >> shift);
        }
        src++;
        dst += dst_stride;
    }
}


template <>
void havoc_inverse_partial_butterfly_c_opt<32>(int16_t dst[], const int16_t src[], int shift)
{
    const int add = 1 << (shift - 1);
    const int src_stride = 32;
    const int dst_stride = 32;

    for (int j = 0; j<32; j++)
    {
        static const int16_t table[32][32] =
        {
                { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
                { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13, 4, -4, -13, -22, -31, -38, -46, -54, -61, -67, -73, -78, -82, -85, -88, -90, -90 },
                { 90, 87, 80, 70, 57, 43, 25, 9, -9, -25, -43, -57, -70, -80, -87, -90, -90, -87, -80, -70, -57, -43, -25, -9, 9, 25, 43, 57, 70, 80, 87, 90 },
                { 90, 82, 67, 46, 22, -4, -31, -54, -73, -85, -90, -88, -78, -61, -38, -13, 13, 38, 61, 78, 88, 90, 85, 73, 54, 31, 4, -22, -46, -67, -82, -90 },
                { 89, 75, 50, 18, -18, -50, -75, -89, -89, -75, -50, -18, 18, 50, 75, 89, 89, 75, 50, 18, -18, -50, -75, -89, -89, -75, -50, -18, 18, 50, 75, 89 },
                { 88, 67, 31, -13, -54, -82, -90, -78, -46, -4, 38, 73, 90, 85, 61, 22, -22, -61, -85, -90, -73, -38, 4, 46, 78, 90, 82, 54, 13, -31, -67, -88 },
                { 87, 57, 9, -43, -80, -90, -70, -25, 25, 70, 90, 80, 43, -9, -57, -87, -87, -57, -9, 43, 80, 90, 70, 25, -25, -70, -90, -80, -43, 9, 57, 87 },
                { 85, 46, -13, -67, -90, -73, -22, 38, 82, 88, 54, -4, -61, -90, -78, -31, 31, 78, 90, 61, 4, -54, -88, -82, -38, 22, 73, 90, 67, 13, -46, -85 },
                { 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83 },
                { 82, 22, -54, -90, -61, 13, 78, 85, 31, -46, -90, -67, 4, 73, 88, 38, -38, -88, -73, -4, 67, 90, 46, -31, -85, -78, -13, 61, 90, 54, -22, -82 },
                { 80, 9, -70, -87, -25, 57, 90, 43, -43, -90, -57, 25, 87, 70, -9, -80, -80, -9, 70, 87, 25, -57, -90, -43, 43, 90, 57, -25, -87, -70, 9, 80 },
                { 78, -4, -82, -73, 13, 85, 67, -22, -88, -61, 31, 90, 54, -38, -90, -46, 46, 90, 38, -54, -90, -31, 61, 88, 22, -67, -85, -13, 73, 82, 4, -78 },
                { 75, -18, -89, -50, 50, 89, 18, -75, -75, 18, 89, 50, -50, -89, -18, 75, 75, -18, -89, -50, 50, 89, 18, -75, -75, 18, 89, 50, -50, -89, -18, 75 },
                { 73, -31, -90, -22, 78, 67, -38, -90, -13, 82, 61, -46, -88, -4, 85, 54, -54, -85, 4, 88, 46, -61, -82, 13, 90, 38, -67, -78, 22, 90, 31, -73 },
                { 70, -43, -87, 9, 90, 25, -80, -57, 57, 80, -25, -90, -9, 87, 43, -70, -70, 43, 87, -9, -90, -25, 80, 57, -57, -80, 25, 90, 9, -87, -43, 70 },
                { 67, -54, -78, 38, 85, -22, -90, 4, 90, 13, -88, -31, 82, 46, -73, -61, 61, 73, -46, -82, 31, 88, -13, -90, -4, 90, 22, -85, -38, 78, 54, -67 },
                { 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64 },
                { 61, -73, -46, 82, 31, -88, -13, 90, -4, -90, 22, 85, -38, -78, 54, 67, -67, -54, 78, 38, -85, -22, 90, 4, -90, 13, 88, -31, -82, 46, 73, -61 },
                { 57, -80, -25, 90, -9, -87, 43, 70, -70, -43, 87, 9, -90, 25, 80, -57, -57, 80, 25, -90, 9, 87, -43, -70, 70, 43, -87, -9, 90, -25, -80, 57 },
                { 54, -85, -4, 88, -46, -61, 82, 13, -90, 38, 67, -78, -22, 90, -31, -73, 73, 31, -90, 22, 78, -67, -38, 90, -13, -82, 61, 46, -88, 4, 85, -54 },
                { 50, -89, 18, 75, -75, -18, 89, -50, -50, 89, -18, -75, 75, 18, -89, 50, 50, -89, 18, 75, -75, -18, 89, -50, -50, 89, -18, -75, 75, 18, -89, 50 },
                { 46, -90, 38, 54, -90, 31, 61, -88, 22, 67, -85, 13, 73, -82, 4, 78, -78, -4, 82, -73, -13, 85, -67, -22, 88, -61, -31, 90, -54, -38, 90, -46 },
                { 43, -90, 57, 25, -87, 70, 9, -80, 80, -9, -70, 87, -25, -57, 90, -43, -43, 90, -57, -25, 87, -70, -9, 80, -80, 9, 70, -87, 25, 57, -90, 43 },
                { 38, -88, 73, -4, -67, 90, -46, -31, 85, -78, 13, 61, -90, 54, 22, -82, 82, -22, -54, 90, -61, -13, 78, -85, 31, 46, -90, 67, 4, -73, 88, -38 },
                { 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36 },
                { 31, -78, 90, -61, 4, 54, -88, 82, -38, -22, 73, -90, 67, -13, -46, 85, -85, 46, 13, -67, 90, -73, 22, 38, -82, 88, -54, -4, 61, -90, 78, -31 },
                { 25, -70, 90, -80, 43, 9, -57, 87, -87, 57, -9, -43, 80, -90, 70, -25, -25, 70, -90, 80, -43, -9, 57, -87, 87, -57, 9, 43, -80, 90, -70, 25 },
                { 22, -61, 85, -90, 73, -38, -4, 46, -78, 90, -82, 54, -13, -31, 67, -88, 88, -67, 31, 13, -54, 82, -90, 78, -46, 4, 38, -73, 90, -85, 61, -22 },
                { 18, -50, 75, -89, 89, -75, 50, -18, -18, 50, -75, 89, -89, 75, -50, 18, 18, -50, 75, -89, 89, -75, 50, -18, -18, 50, -75, 89, -89, 75, -50, 18 },
                { 13, -38, 61, -78, 88, -90, 85, -73, 54, -31, 4, 22, -46, 67, -82, 90, -90, 82, -67, 46, -22, -4, 31, -54, 73, -85, 90, -88, 78, -61, 38, -13 },
                { 9, -25, 43, -57, 70, -80, 87, -90, 90, -87, 80, -70, 57, -43, 25, -9, -9, 25, -43, 57, -70, 80, -87, 90, -90, 87, -80, 70, -57, 43, -25, 9 },
                { 4, -13, 22, -31, 38, -46, 54, -61, 67, -73, 78, -82, 85, -88, 90, -90, 90, -90, 88, -85, 82, -78, 73, -67, 61, -54, 46, -38, 31, -22, 13, -4 }
        };

        int O[16];
        for (int k = 0; k<16; k++)
        {
            O[k] = table[1][k] * src[src_stride] + table[3][k] * src[3 * src_stride] + table[5][k] * src[5 * src_stride] + table[7][k] * src[7 * src_stride] +
                    table[9][k] * src[9 * src_stride] + table[11][k] * src[11 * src_stride] + table[13][k] * src[13 * src_stride] + table[15][k] * src[15 * src_stride] +
                    table[17][k] * src[17 * src_stride] + table[19][k] * src[19 * src_stride] + table[21][k] * src[21 * src_stride] + table[23][k] * src[23 * src_stride] +
                    table[25][k] * src[25 * src_stride] + table[27][k] * src[27 * src_stride] + table[29][k] * src[29 * src_stride] + table[31][k] * src[31 * src_stride];
        }

        int EO[8];
        for (int k = 0; k<8; k++)
        {
            EO[k] = table[2][k] * src[2 * src_stride] + table[6][k] * src[6 * src_stride] + table[10][k] * src[10 * src_stride] + table[14][k] * src[14 * src_stride] +
                    table[18][k] * src[18 * src_stride] + table[22][k] * src[22 * src_stride] + table[26][k] * src[26 * src_stride] + table[30][k] * src[30 * src_stride];
        }

        int EEO[4];
        for (int k = 0; k<4; k++)
        {
            EEO[k] = table[4][k] * src[4 * src_stride] + table[12][k] * src[12 * src_stride] + table[20][k] * src[20 * src_stride] + table[28][k] * src[28 * src_stride];
        }

        int EEEO[2];
        EEEO[0] = table[8][0] * src[8 * src_stride] + table[24][0] * src[24 * src_stride];
        EEEO[1] = table[8][1] * src[8 * src_stride] + table[24][1] * src[24 * src_stride];

        int EEEE[2];
        EEEE[0] = table[0][0] * src[0] + table[16][0] * src[16 * src_stride];
        EEEE[1] = table[0][1] * src[0] + table[16][1] * src[16 * src_stride];

        int EEE[4];
        EEE[0] = EEEE[0] + EEEO[0];
        EEE[3] = EEEE[0] - EEEO[0];
        EEE[1] = EEEE[1] + EEEO[1];
        EEE[2] = EEEE[1] - EEEO[1];

        int EE[8];
        for (int k = 0; k<4; k++)
        {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 4] = EEE[3 - k] - EEO[3 - k];
        }

        int E[16];
        for (int k = 0; k<8; k++)
        {
            E[k] = EE[k] + EO[k];
            E[k + 8] = EE[7 - k] - EO[7 - k];
        }
        for (int k = 0; k<16; k++)
        {
            dst[k] = Clip3(-32768, 32767, (E[k] + O[k] + add) >> shift);
            dst[k + 16] = Clip3(-32768, 32767, (E[15 - k] - O[15 - k] + add) >> shift);
        }
        src++;
        dst += dst_stride;
    }
}


template <int nTbS>
void havoc_idst_c_opt(int16_t dst[], int16_t const coeffs[], int bitDepth)
{
    static_assert(nTbS == 4, "");
    int16_t temp[4 * 4];
    havoc_inverse_partial_butterfly_4x4_dst_c_opt(temp, coeffs, 7);
    havoc_inverse_partial_butterfly_4x4_dst_c_opt(dst, temp, 20 - bitDepth);
}


template <int nTbS>
void havoc_idct_c_opt(int16_t dst[], int16_t const coeffs[], int bitDepth)
{
    int16_t temp[nTbS * nTbS];
    havoc_inverse_partial_butterfly_c_opt<nTbS>(temp, coeffs, 7);
    havoc_inverse_partial_butterfly_c_opt<nTbS>(dst, temp, 20 - bitDepth);
}


void havoc_idct_4x4_c_opt(uint8_t *dst, intptr_t stride_dst, const uint8_t *pred, intptr_t stride_pred, const int16_t coeffs[4 * 4], int bitDepth)
{
    int16_t temp[2][4 * 4];
    havoc_inverse_partial_butterfly_c_opt<4>(temp[0], coeffs, 7);
    havoc_inverse_partial_butterfly_c_opt<4>(temp[1], temp[0], 20 - bitDepth);
    havoc_add_residual(4, dst, stride_dst, pred, stride_pred, temp[1], 8);
}


void havoc_idct_4x4_16_c_opt(uint16_t *dst, intptr_t stride_dst, const uint16_t *pred, intptr_t stride_pred, const int16_t coeffs[4 * 4], int bitDepth)
{
    int16_t temp[2][4 * 4];
    havoc_inverse_partial_butterfly_c_opt<4>(temp[0], coeffs, 7);
    havoc_inverse_partial_butterfly_c_opt<4>(temp[1], temp[0], 20 - bitDepth);
    havoc_add_residual(4, dst, stride_dst, pred, stride_pred, temp[1], bitDepth);
}


void havoc_idst_4x4_c_opt(uint8_t *dst, intptr_t stride_dst, const uint8_t *pred, intptr_t stride_pred, const int16_t coeffs[4 * 4], int bitDepth)
{
    int16_t temp[2][4 * 4];
    havoc_inverse_partial_butterfly_4x4_dst_c_opt(temp[0], coeffs, 7);
    havoc_inverse_partial_butterfly_4x4_dst_c_opt(temp[1], temp[0], 20 - bitDepth);
    havoc_add_residual(4, dst, stride_dst, pred, stride_pred, temp[1], 8);
}


void havoc_idst_4x4_16_c_opt(uint16_t *dst, intptr_t stride_dst, const uint16_t *pred, intptr_t stride_pred, const int16_t coeffs[4 * 4], int bitDepth)
{
    int16_t temp[2][4 * 4];
    havoc_inverse_partial_butterfly_4x4_dst_c_opt(temp[0], coeffs, 7);
    havoc_inverse_partial_butterfly_4x4_dst_c_opt(temp[1], temp[0], 20 - bitDepth);
    havoc_add_residual(4, dst, stride_dst, pred, stride_pred, temp[1], bitDepth);
}


template <int nTbS, typename Sample>
void havoc_idct_add_c_opt(Sample *dst, intptr_t stride_dst, Sample const* pred, intptr_t stride_pred, int16_t const coeffs[8 * 8], int bitDepth)
{
    int16_t temp[2][nTbS * nTbS];
    havoc_inverse_partial_butterfly_c_opt<nTbS>(temp[0], coeffs, 7);
    havoc_inverse_partial_butterfly_c_opt<nTbS>(temp[1], temp[0], 20 - bitDepth);
    havoc_add_residual(nTbS, dst, stride_dst, pred, stride_pred, temp[1], bitDepth);
}


struct InverseTransformAdd :
        Jit::Function
        {
            InverseTransformAdd(Jit::Buffer *buffer, int trType, int log2TrafoSize)
                :
                    Jit::Function(buffer, Jit::CountArguments<havoc_inverse_transform_add<uint8_t>>::value),
                    trType(trType),
                    log2TrafoSize(log2TrafoSize)
            {
                if (log2TrafoSize == 4 && (this->isa() & HAVOC_AVX2))
                    this->buildSinglePass(4, 16, 512 + 512 + 512);
                else if (log2TrafoSize == 5 && (this->isa() & HAVOC_AVX2))
                    this->buildSinglePass(4, 16, 256 * 16);
                else
                    this->build();
            }

            Xbyak::Label cosine_inverse_4;
            Xbyak::Label cosine_inverse_4_h;
            Xbyak::Label cosine_inverse_8;
            Xbyak::Label dd_0040;
            Xbyak::Label dd_0800;
            Xbyak::Label shuffle_018945cd018945cd;
            Xbyak::Label shuffle_2367abef2367abef;
            Xbyak::Label cosine_inverse_8_h;
            Xbyak::Label cosine_inverse_8_h2;
            Xbyak::Label cosine_inverse_16;
            Xbyak::Label shuffle_2367236723672367;
            Xbyak::Label shuffle_abefabefabefabef;
            Xbyak::Label cosine_inverse_16_h;
            Xbyak::Label shuffle_45cd45cd45cd45cd;
            Xbyak::Label shuffle_0189018901890189;
            Xbyak::Label shuffle_efcdab8967452301;

            Xbyak::Label pat_idct4_pass1;
            Xbyak::Label pat_idct4_pass2;
            Xbyak::Label pat_idct4_sign;
            Xbyak::Label pat_idct4_shuf1;
            Xbyak::Label pat_dw_64;
            Xbyak::Label pat_dw_2048;
            Xbyak::Label pat_idct16_shuf1;
            Xbyak::Label pat_idst_pass1;
            Xbyak::Label pat_idst_shuf;
            Xbyak::Label pat_idst_pass2;
            Xbyak::Label pat_idct8_pass1;
            Xbyak::Label pat_idct8_pass2;
            Xbyak::Label pat_idct16_pass1;
            Xbyak::Label pat_idct16_pass2;
            Xbyak::Label pat_idct16_shuf3;
            Xbyak::Label pat_idct16_sign;
            Xbyak::Label pat_idct16_shuf2;
            Xbyak::Label pass1_process;
            Xbyak::Label pat_idct32_pass1;
            Xbyak::Label pat_idct32_sign;
            Xbyak::Label loop_idct32_pass1;
            Xbyak::Label loop_idct32_pass1_odd;
            Xbyak::Label loop_idct32_pass1_combine;
            Xbyak::Label pat_idct32_pass2;
            Xbyak::Label pat_idct32_shuf1;
            Xbyak::Label pat_idct32_shuf2;
            Xbyak::Label loop_idct32_pass2_combine2;
            Xbyak::Label loop_idct32_pass2_odd;

#if USE_F265_DERIVED
            void dataAvx2()
            {
                // ---------------------- DCT/IDCT 32 macros ---------------------
                // Load 4 rows of 32-bit data from consecutive memory locations.
                // a1: base address, a2: data location index.
#define LOAD_DATA(a1, a2) \
        vmovdqu(y0, ptr[g5 + a1 + a2]); \
        vmovdqu(y1, ptr[g5 + a1 + a2 + 32]); \
        vmovdqu(y2, ptr[g5 + a1 + a2 + 64]); \
        vmovdqu(y3, ptr[g5 + a1 + a2 + 96]); \

                // Store 4 rows of 32-bit data to consecutive memory locations.
                // a1: base address, a2 : data location index.
#define STORE_DATA(a1, a2) \
        vmovdqu(ptr[g5 + a1 + a2], y0); \
        vmovdqu(ptr[g5 + a1 + a2 + 32], y1); \
        vmovdqu(ptr[g5 + a1 + a2 + 64], y2); \
        vmovdqu(ptr[g5 + a1 + a2 + 96], y3); \

                // Multiply 2 rows of data and add to the destination.
                // a1: destination, a2: input 0, a3: input 1.
#define MULT_ADD(a1, a2, a3) \
        vpmaddwd(y4, a2, y6); \
        vpmaddwd(y5, a3, y7); \
        vpaddd(a1, a1, y4); \
        vpaddd(a1, a1, y5); \

                // Multiply 2 rows of data and subtract from the destination.
                // a1: destination, a2: input 0, a3: input 1.
#define MULT_SUB(a1, a2, a3) \
        vpmaddwd(y4, a2, y6); \
        vpmaddwd(y5, a3, y7); \
        vpsubd(a1, a1, y4); \
        vpsubd(a1, a1, y5); \

                // Multiply 2 rows of data and perform a horizontal addition of terms.
                // a1: result, a2: tmp, a3: DCT/IDCT factor 0, a4: DCT/IDCT factor 1.
#define MULT_HADD(a1, a2, a3, a4) \
        vpmaddwd(a1, a3, y0); \
        vpmaddwd(a2, a4, y0); \
        vphaddd(a1, a1, a2); \



                align(32);

                // from f265
                L(pat_idst_pass1);
                dw({ 29, 74, 84, 55, 55, 74, -29, -84, 74, 0, -74, 74, 84, -74, 55, -29 });

                L(pat_idct4_shuf1);
                db({ 0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15 });

                L(pat_idst_shuf);
                dd({ 0, 4, 1, 5, 2, 6, 3, 7 });

                L(pat_idst_pass2);
                dw({ 29, 74, 29, 74, 29, 74, 29, 74, 84, 55, 84, 55, 84, 55, 84, 55 });
                dw({ 74, 0, 74, 0, 74, 0, 74, 0, -74, 74, -74, 74, -74, 74, -74, 74 });
                dw({ 55, 74, 55, 74, 55, 74, 55, 74, -29, -84, -29, -84, -29, -84, -29, -84 });
                dw({ 84, -74, 84, -74, 84, -74, 84, -74, 55, -29, 55, -29, 55, -29, 55, -29 });

                // For IDCT4 Pass 1, the multiplication by 64 is done by a left shift operation. Hence, we need only the coefficients
                // for 'add_1' and 'sub_1' (the notation is from the C code). The first 8 values are for calculating 'add_1' and the
                // remaining 8 for 'sub_1'.

                L(pat_idct4_pass1);
                dw({ 83, 36, 83, 36, 83, 36, 83, 36, 36, -83, 36, -83, 36, -83, 36, -83 });

                L(pat_idct4_pass2);
                dw({ 64, 83, 64, 36, 64, 36, -64, -83, 64, -36, -64, 83, 64, -83, 64, -36 });

                L(pat_idct4_sign);
                dd({ 1, 1, 1, 1, -1, -1, -1, -1 });

                L(pat_idct8_pass1);
                dw({ 64, 64, 83, 36, 64, -64, 36, -83, 64, -64, -36, 83, 64, 64, -83, -36 });
                dw({ 89, 50, 75, 18, 75, -89, -18, -50, 50, 18, -89, 75, 18, 75, -50, -89 });

                L(pat_idct8_pass2);
                dw({ 64, 89, 83, 75, 64, 50, 36, 18, 64, 75, 36, -18, -64, -89, -83, -50 });
                dw({ 64, 50, -36, -89, -64, 18, 83, 75, 64, 18, -83, -50, 64, 75, -36, -89 });
                dw({ 64, -18, -83, 50, 64, -75, -36, 89, 64, -50, -36, 89, -64, -18, 83, -75 });
                dw({ 64, -75, 36, 18, -64, 89, -83, 50, 64, -89, 83, -75, 64, -50, 36, -18 });

                L(pat_idct16_shuf1);
                db({ 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15 });

                align(32);
                L(pat_idct16_pass1);
                db({ 64, 89, 83, 75, 64, 50, 36, 18, 64, 75, 36, -18, -64, -89, -83, -50 });
                db({ 64, 50, -36, -89, -64, 18, 83, 75, 64, 18, -83, -50, 64, 75, -36, -89 });
                db({ 64, -18, -83, 50, 64, -75, -36, 89, 64, -50, -36, 89, -64, -18, 83, -75 });
                db({ 64, -75, 36, 18, -64, 89, -83, 50, 64, -89, 83, -75, 64, -50, 36, -18 });
                db({ 90, 87, 80, 70, 57, 43, 25, 9, 87, 57, 9, -43, -80, -90, -70, -25 });
                db({ 80, 9, -70, -87, -25, 57, 90, 43, 70, -43, -87, 9, 90, 25, -80, -57 });
                db({ 57, -80, -25, 90, -9, -87, 43, 70, 43, -90, 57, 25, -87, 70, 9, -80 });
                db({ 25, -70, 90, -80, 43, 9, -57, 87, 9, -25, 43, -57, 70, -80, 87, -90 });

                L(pat_idct16_pass2);
                db({ 64, 89, 83, 75, 90, 87, 80, 70, 64, 50, 36, 18, 57, 43, 25, 9, });
                db({ 64, 75, 36, -18, 87, 57, 9, -43, -64, -89, -83, -50, -80, -90, -70, -25, });
                db({ 64, 50, -36, -89, 80, 9, -70, -87, -64, 18, 83, 75, -25, 57, 90, 43, });
                db({ 64, 18, -83, -50, 70, -43, -87, 9, 64, 75, -36, -89, 90, 25, -80, -57, });
                db({ 64, -18, -83, 50, 57, -80, -25, 90, 64, -75, -36, 89, -9, -87, 43, 70, });
                db({ 64, -50, -36, 89, 43, -90, 57, 25, -64, -18, 83, -75, -87, 70, 9, -80, });
                db({ 64, -75, 36, 18, 25, -70, 90, -80, -64, 89, -83, 50, 43, 9, -57, 87, });
                db({ 64, -89, 83, -75, 9, -25, 43, -57, 64, -50, 36, -18, 70, -80, 87, -90, });

                L(pat_idct16_shuf3);
                db({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 });
                db({ 2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13 });

                align(32);
                L(pat_idct32_pass1);
                // Since we process 4 rows at a time, 4 consecutive values are required to obtain one output row.
                // Each row below can be divided into sets of 4 consecutive values and are reused by the inputs.
                // - Note: The notation in the explanation below is from the C code.
                //
                // The first 8 rows are used to calculate the odd terms (viz, s0, s1, .. s15).
                // - The first 4 rows (i.e. 16 sets) are used sequentially by input rows 1, 3, 5, 7 (i.e. src[32],
                // src[96], src[160] and src[224]) to calculate s0, s1, .. s15 (e.g. <90, 90, 88, 85> is used to
                // calculate s0; <90, 82, 67, 46> for s2 ...).
                // - The values in the first 4 rows are alternately negated and used by input rows 31, 29, 28, 27
                // (i.e. src[992], src[928], src[864] and src[800]) to obtain s15, s14, ... s0 respectively
                // (e.g. <-90, 90, -88, 85> is used to obtain s15; <90, -82, 67, -46> is used for s14). Note the
                // descending order of input rows, as well as the descending order of outputs (i.e. <s15, s14 .. s0>
                // instead of <s0, s1, ... s15>).
                // - The second 4 rows are used by input rows 9, 11, 13, 15 (i.e. src[288], src[352], src[416] and
                // src[480]) to calculate s0, s1 .. s15 (e.g. <82, 78, 73, 67> is used to obtain s0;
                // <22, -4, -31, -54> is used to obtain s1 ..).
                // - The second 4 row values are alternately negated and used by input rows 23, 21, 19, 17 (i.e.
                // src[736], src[672], src[608] and src[544]) to obtain s15, s14 .. s0 (e.g. <-82, 78, -73, 67>
                // is used to calculate s15; <22, 4, -31, 54> for s14). Note the descending order of input rows,
                // as well as the descending order of outputs.
                dw({ 90, 90, 88, 85, 90, 82, 67, 46, 88, 67, 31, -13, 85, 46, -13, -67 });
                dw({ 82, 22, -54, -90, 78, -4, -82, -73, 73, -31, -90, -22, 67, -54, -78, 38 });
                dw({ 61, -73, -46, 82, 54, -85, -4, 88, 46, -90, 38, 54, 38, -88, 73, -4 });
                dw({ 31, -78, 90, -61, 22, -61, 85, -90, 13, -38, 61, -78, 4, -13, 22, -31 });
                dw({ 82, 78, 73, 67, 22, -4, -31, -54, -54, -82, -90, -78, -90, -73, -22, 38 });
                dw({ -61, 13, 78, 85, 13, 85, 67, -22, 78, 67, -38, -90, 85, -22, -90, 4 });
                dw({ 31, -88, -13, 90, -46, -61, 82, 13, -90, 31, 61, -88, -67, 90, -46, -31 });
                dw({ 4, 54, -88, 82, 73, -38, -4, 46, 88, -90, 85, -73, 38, -46, 54, -61 });
                // The following rows are used to calculate the even terms (viz, a0, a1, .. a15).
                // The following first two rows are used to calculate e_0, e_1 .. e_7.
                // - The 8 sets of values in the following 2 rows are used sequentially by input rows 2, 6, 10, 14
                // (i.e. src[64], src[192], src[320] and src[448]) to obtain e_0, e_1, .. e_7 (e.g. <90, 87, 80, 70>
                //  is used to obtain e_0; <87, 57, 9, -43> for e_1 ..).
                // - The values are alternately negated and used by input rows 30, 26, 22, 18 (i.e. src[960],
                // src[832], src[704] and src[576]) to calculate e_7, e_6, .. e_0 (e.g. <-90, 87, -80, 70> is used
                // for e_7; <87, -57, 9, 43> for e_6). Note the descending order of input rows, and the outputs.
                dw({ 90, 87, 80, 70, 87, 57, 9, -43, 80, 9, -70, -87, 70, -43, -87, 9 });
                dw({ 57, -80, -25, 90, 43, -90, 57, 25, 25, -70, 90, -80, 9, -25, 43, -57 });
                // The following row is used by input rows 4, 12, 20, 28 (i.e. src[128], src[384], src[640] and
                // src[896]) to calculate a_4, a_5, a_6 and a_7 (e.g. <89, 75, 50, 18> is used to obtain a_4).
                dw({ 89, 75, 50, 18, 75, -18, -89, -50, 50, -89, 18, 75, 18, -50, 75, -89 });
                // The following row is used to get a_0, a_1, a_2 and a_3 by input rows 0, 16, 8, 24 (i.e. src[0],
                // src[512], src[256] and src[768]). Note the order of input rows. <64, 64, 83, 36> is used to
                // obtain a_0 and a_3, and <64, -64, 36, -83> for a_1 and a_2.
                dw({ 64, 64, 83, 36, 64, -64, 36, -83 });

                L(pat_idct32_pass2); // The first 8 rows are used to calculate the even terms (viz, a0, a1 .. a15). The symmetry of the
                // even terms helps reduce multiplications, and thereby reduce the number of IDCT factors required.
                // While combining the odd and even terms in IDCT Pass 1, the terms used to calculate a_0 .. a_7
                // (i.e. columns 0, 4, 8, 12, 16, 20, 24 and 28) in IDCT Pass 2 are grouped together and the terms
                // for e_0, ..e_7 (i.e. columns 2, 6, 10, 14, 18, 22, 26 and 30) are grouped together. In IDCT
                // Pass 2, just one set of multiplications gives 2 sets of outputs (e.g. a0 = a_0 + e_0 + x,
                // a15 = a_0 - e_0 + x). Likewise, IDCT factors in the following 8 rows are also grouped together
                // (e.g. <64, 89, 83, 75, 64, 50, 36, 18> is used for a_0; <90, 87, 80, 70, 57, 43, 25, 9> for e_0).
                db({ 64, 89, 83, 75, 64, 50, 36, 18, 90, 87, 80, 70, 57, 43, 25, 9 });
                db({ 64, 75, 36, -18, -64, -89, -83, -50, 87, 57, 9, -43, -80, -90, -70, -25 });
                db({ 64, 50, -36, -89, -64, 18, 83, 75, 80, 9, -70, -87, -25, 57, 90, 43 });
                db({ 64, 18, -83, -50, 64, 75, -36, -89, 70, -43, -87, 9, 90, 25, -80, -57 });
                db({ 64, -18, -83, 50, 64, -75, -36, 89, 57, -80, -25, 90, -9, -87, 43, 70 });
                db({ 64, -50, -36, 89, -64, -18, 83, -75, 43, -90, 57, 25, -87, 70, 9, -80 });
                db({ 64, -75, 36, 18, -64, 89, -83, 50, 25, -70, 90, -80, 43, 9, -57, 87 });
                db({ 64, -89, 83, -75, 64, -50, 36, -18, 9, -25, 43, -57, 70, -80, 87, -90 });
                // The following 16 rows are used to calculate the 16 odd terms (s0, s1 .. s15).
                // Each row below gives a term (e.g. row 0 => s0, row 1 => s1, ... row 15 => s15).
                db({ 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13, 4 });
                db({ 90, 82, 67, 46, 22, -4, -31, -54, -73, -85, -90, -88, -78, -61, -38, -13 });
                db({ 88, 67, 31, -13, -54, -82, -90, -78, -46, -4, 38, 73, 90, 85, 61, 22 });
                db({ 85, 46, -13, -67, -90, -73, -22, 38, 82, 88, 54, -4, -61, -90, -78, -31 });
                db({ 82, 22, -54, -90, -61, 13, 78, 85, 31, -46, -90, -67, 4, 73, 88, 38 });
                db({ 78, -4, -82, -73, 13, 85, 67, -22, -88, -61, 31, 90, 54, -38, -90, -46 });
                db({ 73, -31, -90, -22, 78, 67, -38, -90, -13, 82, 61, -46, -88, -4, 85, 54 });
                db({ 67, -54, -78, 38, 85, -22, -90, 4, 90, 13, -88, -31, 82, 46, -73, -61 });
                db({ 61, -73, -46, 82, 31, -88, -13, 90, -4, -90, 22, 85, -38, -78, 54, 67 });
                db({ 54, -85, -4, 88, -46, -61, 82, 13, -90, 38, 67, -78, -22, 90, -31, -73 });
                db({ 46, -90, 38, 54, -90, 31, 61, -88, 22, 67, -85, 13, 73, -82, 4, 78 });
                db({ 38, -88, 73, -4, -67, 90, -46, -31, 85, -78, 13, 61, -90, 54, 22, -82 });
                db({ 31, -78, 90, -61, 4, 54, -88, 82, -38, -22, 73, -90, 67, -13, -46, 85 });
                db({ 22, -61, 85, -90, 73, -38, -4, 46, -78, 90, -82, 54, -13, -31, 67, -88 });
                db({ 13, -38, 61, -78, 88, -90, 85, -73, 54, -31, 4, 22, -46, 67, -82, 90 });
                db({ 4, -13, 22, -31, 38, -46, 54, -61, 67, -73, 78, -82, 85, -88, 90, -90 });

                L(pat_idct32_shuf1);
                dd({ 7, 6, 5, 4, 3, 2, 1, 0 });

                L(pat_idct32_shuf2);
                db({ 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9 });

                L(pat_dw_64);
                dd({ 64 });

                L(pat_dw_2048);
                dd({ 2048 });

                L(pat_idct16_sign);
                dd({ 1, -1 });

                L(pat_idct16_shuf2);
                db({ 0, 4, 2, 6, 7, 3, 5, 1 });

                L(pat_idct32_sign);
                dw({ -1, -1 });


                if (this->log2TrafoSize == 4)
                {
                    auto &y0 = ymm0; auto &x0 = regXmm(0);
                    auto &y1 = ymm1; auto &x1 = regXmm(1);
                    auto &y2 = ymm2; auto &x2 = regXmm(2);
                    auto &y3 = ymm3; auto &x3 = regXmm(3);
                    auto &y7 = ymm7; auto &x7 = regXmm(7);
                    auto &y8 = ymm8; auto &x8 = regXmm(8);
                    auto &y9 = ymm9; auto &x9 = regXmm(9);
                    auto &y10 = ymm10; auto &x10 = regXmm(10);
                    auto &y11 = ymm11; auto &x11 = regXmm(11);
                    auto &y12 = ymm12; auto &x12 = regXmm(12);
                    auto &y13 = ymm13; auto &x13 = regXmm(13);
                    auto &y14 = ymm14; auto &x14 = regXmm(14);
                    auto &g8 = reg64(8);

                    // Process all rows in pass 1.
                    L(pass1_process);
                    // Process a row.

                    // ; %1-2: out, %3-4: in, %5: DCT factors.
#define PROCESS_ROW(a1, a2, a3, a4, a5) \
        vpbroadcastd(a2, ptr[a5]); \
        vpmaddwd(a1, a3, a2); \
        vpmaddwd(a2, a4, a2); \

                    PROCESS_ROW(y0, y1, y14, y13, g8);
                    PROCESS_ROW(y2, y3, y12, y11, g8 + 4);
                    vpaddd(y0, y0, y2);
                    vpaddd(y1, y1, y3);

                    PROCESS_ROW(y2, y3, y10, y9, g8 + 8);
                    vpaddd(y0, y0, y2);
                    vpaddd(y1, y1, y3);

                    PROCESS_ROW(y2, y3, y8, y7, g8 + 12);
                    vpaddd(y0, y0, y2);
                    vpaddd(y1, y1, y3);

                    add(g8, 16);
                    ret();

#undef PROCESS_ROW
                }

                if (this->log2TrafoSize == 5)
                {
                    // IDCT32 Helper functions
                    auto &y0 = ymm0;
                    auto &y1 = ymm1;
                    auto &y2 = ymm2;
                    auto &y3 = ymm3;
                    auto &y4 = ymm4;
                    auto &y5 = ymm5;
                    auto &y6 = ymm6;
                    auto &y7 = ymm7;
                    auto &y8 = ymm8;
                    auto &y9 = ymm9;
                    auto &y10 = ymm10;
                    auto &y11 = ymm11;
                    auto &y12 = ymm12;
                    auto &y13 = ymm13;
                    auto &y14 = ymm14;
                    auto &y15 = ymm15;

                    auto &g0 = reg64(0);
                    auto &g1 = reg64(1);
                    auto &g2 = reg64(2);
                    auto &g3 = reg64(3);
                    auto &g4 = reg64(4);
                    auto &g5 = reg64(5);
                    auto &g6 = reg64(6);
                    auto &g7 = reg64(7);
                    auto &g8 = reg64(8);

                    // Perform IDCT32 Pass 1 related operations.
                    L(loop_idct32_pass1);
                    {
                        vpbroadcastd(y0, ptr[g7]);                // Load the DCT multiplication factors.
                        vpbroadcastd(y1, ptr[g7 + 4]);

                        // a1: in 0, a2: in 1, a3: store location 1.
#define PROCESS_ROW(a1, a2, a3) \
        vpmaddwd(y2, a1, y0); \
        vpmaddwd(y3, a2, y1); \
        vpaddd(y2, y2, y3); \
        vmovdqu(ptr[g5 + g6 + a3], y2); \

                        PROCESS_ROW(y15, y11, 0);
                        PROCESS_ROW(y14, y10, 32);
                        PROCESS_ROW(y13, y9, 64);
                        PROCESS_ROW(y12, y8, 96);

#undef PROCESS_ROW

                        add(g6, 128);
                        add(g7, 8);
                        cmp(g6, g8);
                        jne(loop_idct32_pass1);
                        ret();

                        // Perform IDCT32 Pass 1 related operations for odd rows alone.
                        L(loop_idct32_pass1_odd);
                        LOAD_DATA(g6, 0);
                        vpbroadcastd(y6, ptr[g7]);
                        vpbroadcastd(y7, ptr[g7 + 4]);
                        MULT_ADD(y0, y15, y11);
                        MULT_ADD(y1, y14, y10);
                        MULT_ADD(y2, y13, y9);
                        MULT_ADD(y3, y12, y8);
                        STORE_DATA(g6, 0);
                        sub(g7, 8);
                        vpbroadcastd(y6, ptr[g7]);
                        vpbroadcastd(y7, ptr[g7 + 4]);
                        LOAD_DATA(g6, 128);
                        MULT_SUB(y0, y15, y11);
                        MULT_SUB(y1, y14, y10);
                        MULT_SUB(y2, y13, y9);
                        MULT_SUB(y3, y12, y8);
                        STORE_DATA(g6, 128);
                        add(g6, 256);
                        sub(g7, 8);
                    }
                    cmp(g6, 256 * 8);
                    jne(loop_idct32_pass1_odd);
                    ret();

                    // Combine odd and even terms for IDCT32 Pass 1.
                    L(loop_idct32_pass1_combine);
                    {
                        // a1: in/out 0, a2: in 1, a3: out 1/in 2.
#define ADD_SUB(a1, a2, a3) \
        vpsubd(a1, a2, a3); \
        vpaddd(a3, a2, a3); \
        vpsrad(a1, a1, 7); \
        vpsrad(a3, a3, 7); \

                        // a1-2: tmp, a3: in/out 0, a4: in/out 1.
#define COMBINE(a1, a2, a3, a4); \
        ADD_SUB(y4, a1, a3); \
        ADD_SUB(y5, a2, a4); \
        vpackssdw(y4, y4, y5); \
        vpackssdw(a3, a3, a4); \
        vpshufb(a3, a3, y6); /* Rearrange all the odd & even terms together : odd, even | odd, even.*/\
        vpshufb(a4, y4, y6); \
        vpermq(a3, a3, 0xD8); \
        vpermq(a4, a4, 0xD8); \

                        // a1-2: in 0-1, a3: memory index to store the separated even and odd terms.
#define SHUFFLE(a1, a2, a3) \
        vperm2i128(y4, a1, a2, 0x20);        /* All even terms. */ \
        vperm2i128(y5, a1, a2, 0x31);        /* All odd terms. */ \
        vpshufb(y4, y4, y6); \
        vmovdqu(ptr[g5 + g6 + a3], y4); \
        vmovdqu(ptr[g5 + g6 + a3 + 32], y5); \

        LOAD_DATA(g6, 0);
        vmovdqu(y15, ptr[g5 + g7 + 0]);
        vmovdqu(y14, ptr[g5 + g7 + 32]);
        vmovdqu(y13, ptr[g5 + g7 + 64]);
        vmovdqu(y12, ptr[g5 + g7 + 96]);
        COMBINE(y15, y14, y0, y1);
        COMBINE(y13, y12, y2, y3);
        // Shuffle the even terms to reduce multiplications for the next pass (i.e. Pass 2). Notation below is from C code.
        // Group srcptr[0], srcptr[128], srcptr[256], srcptr[384], srcptr[512], srcptr[640], srcptr[768], srcptr[896] in the lower lane and
        // srcptr[64], srcptr[192], srcptr[320], srcptr[448], srcptr[576], srcptr[704], srcptr[832], srcptr[960] in the higher lane.
        SHUFFLE(y0, y2, 0);              // Separate even and odd terms of IDCT32 Pass1 row "n".
        SHUFFLE(y1, y3, 64);             // Separate even and odd terms of IDCT32 Pass1 row "31 - n".
#undef SHUFFLE
#undef ADD_SUB
#undef COMBINE
        add(g6, 128);
        add(g7, g4);
                    }
                    cmp(g6, g8);
                    jne(loop_idct32_pass1_combine);
                    ret();

                    // Perform IDCT32 Pass 2 related operations for odd columns alone.
                    L(loop_idct32_pass2_odd);
                    {
                        vmovdqu(y0, ptr[g5 + g6 + 32]);
                        MULT_HADD(y7, y6, y15, y14);
                        MULT_HADD(y5, y4, y13, y12);
                        MULT_HADD(y6, y3, y11, y10);
                        MULT_HADD(y4, y2, y9, y8);
                        vphaddd(y7, y7, y5);
                        vphaddd(y6, y6, y4);
                        vperm2i128(y3, y7, y6, 0x20);
                        vperm2i128(y2, y7, y6, 0x31);
                        vpaddd(y0, y3, y2);
                        vmovdqu(ptr[g8 + g6], y0);
                        add(g6, 64);
                    }
                    cmp(g6, 64 * 32);
                    jne(loop_idct32_pass2_odd);
                    ret();

                    // Combine odd and even terms for IDCT32 Pass 2.
                    L(loop_idct32_pass2_combine2);
                    {

                        // a1-2: out 0-1 (to hold the sum and difference of the odd and even terms
                        // respectively), a3-4: memory index of the even and odd terms respectively.
#define COMBINE(a1, a2, a3, a4) \
        vmovdqu(y0, ptr[g5 + g4 + 2048 + a3]); \
        vmovdqu(y2, ptr[g5 + g4 + a4]);      /* Odd terms. */ \
        vpaddd(a1, y0, y2); \
        vpsubd(a2, y0, y2); \
        vpsrad(a1, a1, 12); \
        vpsrad(a2, a2, 12);

                        COMBINE(y4, y6, 0, 0);
                        COMBINE(y5, y7, 32, 32);

#undef COMBINE

                        vpackssdw(y0, y4, y5);
                        vpackssdw(y1, y6, y7);
                        vpermq(y0, y0, 0xD8);
                        vpshufb(y1, y1, y10);
                        vpermq(y1, y1, 0x27);
                        // Reconstruction.
                        vpmovzxbw(y2, ptr[g2]);
                        vpmovzxbw(y3, ptr[g2 + 16]);
                        add(g2, g3);
                        vpaddw(y0, y0, y2);
                        vpaddw(y1, y1, y3);
                        vpackuswb(y0, y0, y1);
                        vpermq(y0, y0, 0xd8);
                        vmovdqu(ptr[g0], y0);
                        add(g0, g1);
                        add(g4, g6);
                    }
                    cmp(g4, g7);
                    jne(loop_idct32_pass2_combine2);
                    ret();
                }
            }
#endif

            void data()
            {
#if USE_F265_DERIVED
                if (this->isa() & HAVOC_AVX2)
                {
                    this->dataAvx2();
                    return;
                }
#endif

                align();

                L(cosine_inverse_4);
                dw({ 64, 64 }, 4);
                dw({ 64, -64 }, 4);
                dw({ 83, 36 }, 4);
                dw({ 36, -83 }, 4);

                L(cosine_inverse_4_h);
                dw({ 64, -64, 36, -83, 64, 64, 83, 36 });

                L(cosine_inverse_8);
                dw({ 89, 75 }, 4);
                dw({ 50, 18 }, 4);
                dw({ 75, -18 }, 4);
                dw({ -89, -50 }, 4);
                dw({ 50, -89 }, 4);
                dw({ 18, 75 }, 4);
                dw({ 18, -50 }, 4);
                dw({ 75, -89 }, 4);

                L(dd_0040);
                dd({ 0x0040 }, 4);

                L(dd_0800);
                dd({ 0x0800 }, 4);

                L(shuffle_018945cd018945cd);
                db({ 0, 1, 8, 9, 4, 5, 12, 13, 0, 1, 8, 9, 4, 5, 12, 13 });

                L(shuffle_2367abef2367abef);
                db({ 2, 3, 6, 7, 10, 11, 14, 15, 2, 3, 6, 7, 10, 11, 14, 15 });

                L(cosine_inverse_8_h2);
                H_TABLE_ENTRY_8(89, 75, 50, 18, 75, -18, -89, -50);
                H_TABLE_ENTRY_8(50, -89, 18, 75, 18, -50, 75, -89);

                L(cosine_inverse_8_h);
                dw({ 89, 75, 50, 18 });
                dw({ 75, -18, -89, -50 });
                dw({ 50, -89, 18, 75 });
                dw({ 18, -50, 75, -89 });

                L(cosine_inverse_16);
                V_TABLE_ENTRY_4({ 90, 87, 80, 70, 57, 43, 25, 9 });
                V_TABLE_ENTRY_4({ 87, 57, 9, -43, -80, -90, -70, -25 });
                V_TABLE_ENTRY_4({ 80, 9, -70, -87, -25, 57, 90, 43 });
                V_TABLE_ENTRY_4({ 70, -43, -87, 9, 90, 25, -80, -57 });
                V_TABLE_ENTRY_4({ 57, -80, -25, 90, -9, -87, 43, 70 });
                V_TABLE_ENTRY_4({ 43, -90, 57, 25, -87, 70, 9, -80 });
                V_TABLE_ENTRY_4({ 25, -70, 90, -80, 43, 9, -57, 87 });
                V_TABLE_ENTRY_4({ 9, -25, 43, -57, 70, -80, 87, -90 });

                L(shuffle_2367236723672367);
                db({ 2, 3, 6, 7 }, 4);

                L(shuffle_abefabefabefabef);
                db({ 10, 11, 14, 15 }, 4);

                L(cosine_inverse_16_h);
                H_TABLE_ENTRY({ 90, 87, 80, 70, 57, 43, 25, 9,  87, 57, 9, -43, -80, -90, -70, -25 });
                H_TABLE_ENTRY({ 80, 9, -70, -87, -25, 57, 90, 43,  70, -43, -87, 9, 90, 25, -80, -57 });
                H_TABLE_ENTRY({ 57, -80, -25, 90, -9, -87, 43, 70,  43, -90, 57, 25, -87, 70, 9, -80 });
                H_TABLE_ENTRY({ 25, -70, 90, -80, 43, 9, -57, 87,  9, -25, 43, -57, 70, -80, 87, -90 });

                L(shuffle_45cd45cd45cd45cd);
                db({ 4, 5, 12, 13 }, 4);

                L(shuffle_0189018901890189);
                db({ 0, 1, 8, 9 }, 4);

                L(shuffle_efcdab8967452301);
                db({ 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1 });
            }

            void V_TABLE_ENTRY_4(std::initializer_list<int> const &args)
            {
                auto i = args.begin();
                for (int j = 0; j < 4; ++j)
                {
                    auto a = *i++;
                    auto b = *i++;
                    for (int i = 0; i < 4; ++i)
                    {
                        dw({ a, b });
                    }
                }
            }

            struct Index
            {
                const std::initializer_list<int>& list;
                Index(const std::initializer_list<int>& list) : list(list) {}
                int const &operator[](size_t i) { return *(this->list.begin() + i); }
            };

            void H_TABLE_ENTRY(std::initializer_list<int> const &args)
            {
                Index i(args);
                dw({ i[0], i[8], i[1], i[9], i[2], i[10], i[3], i[11] });
                dw({ i[4], i[12], i[5], i[13], i[6], i[14], i[7], i[15] });
            }

            void H_TABLE_ENTRY_8(int c1, int c2, int c3, int c4, int c5, int c6, int c7, int c8)
            {
                dw({ c1, c5 , c2, c6, c3, c7, c4 , c8 });
            }

            void havoc_partial_butterfly_inverse_8v_ssse3()
            {
                //  this function potentially be combined with the horizontal feature with no
                //  need for temporary memory buffer in between.

                //  void transform_partial_butterfly_inverse_8v_ssse3(int16_t *dst, const int16_t *src, int shift);

                auto &r0 = arg64(0);
                auto &r1 = arg64(1);
                auto &r2 = arg64(2);
                auto &r3 = arg64(3);
                auto &r4 = arg64(4);
                auto &r5 = reg64(5);
                auto &r6 = reg64(6);
                auto &r7 = reg64(7);

                auto &m0 = regXmm(0);
                auto &m1 = regXmm(1);
                auto &m2 = regXmm(2);
                auto &m3 = regXmm(3);
                auto &m4 = regXmm(4);
                auto &m5 = regXmm(5);
                auto &m6 = regXmm(6);
                auto &m7 = regXmm(7);
                auto &m8 = regXmm(8);
                auto &m9 = regXmm(9);
                auto &m10 = regXmm(10);
                auto &m11 = regXmm(11);
                auto &m12 = regXmm(12);
                auto &m13 = regXmm(13);

                lea(r3, ptr[rip + cosine_inverse_4]);

                movdqa(m8, ptr[r1 + 0 * 2 * 8]);
                // m8 = src(7, 0), src(6, 0), src(5, 0), src(4, 0), src(3, 0), src(2, 0), src(1, 0), src(0, 0)

                movdqa(m1, ptr[r1 + 4 * 2 * 8]);
                // m1 = src(7, 4), src(6, 4), src(5, 4), src(4, 4), src(3, 4), src(2, 4), src(1, 4), src(0, 4)

                movdqa(m0, m8);
                punpcklwd(m0, m1);
                // m0 = src(3, 4), src(3, 0), src(2, 4), src(2, 0), src(1, 4), src(1, 0), src(0, 4), src(0, 0)

                punpckhwd(m8, m1);
                // m8 = src(7, 4), src(7, 0), src(6, 4), src(6, 0), src(5, 4), src(5, 0), src(4, 4), src(4, 0)

                movdqa(m1, m0);
                pmaddwd(m1, ptr[r3]);
                // m1 = EE[0] : 3, 2, 1, 0

                movdqa(m9, m8);
                pmaddwd(m9, ptr[r3]);
                // m9 = EE[0] : 7, 6, 5, 4

                pmaddwd(m0, ptr[r3 + 16]);
                // m0 = EE[1] : 3, 2, 1, 0

                pmaddwd(m8, ptr[r3 + 16]);
                // m8 = EE[1] : 7, 6, 5, 4

                movdqa(m10, ptr[r1 + 2 * 2 * 8]);
                // m10 = src(7, 2), src(6, 2), src(5, 2), src(4, 2), src(3, 2), src(2, 2), src(1, 2), src(0, 2)

                movdqa(m3, ptr[r1 + 6 * 2 * 8]);
                // m3 = src(7, 6), src(6, 6), src(5, 6), src(4, 6), src(3, 6), src(2, 6), src(1, 6), src(0, 6)

                movdqa(m2, m10);
                punpcklwd(m2, m3);
                //; m2 = src(3, 6), src(3, 2), src(2, 6), src(2, 2), src(1, 6), src(1, 2), src(0, 6), src(0, 2)

                punpckhwd(m10, m3);
                // m10 = src(7, 6), src(7, 2), src(6, 6), src(6, 2), src(5, 6), src(5, 2), src(4, 6), src(4, 2)

                movdqa(m3, m2);
                pmaddwd(m3, ptr[r3 + 32]);
                // m3 = EO[0] : 3, 2, 1, 0

                movdqa(m11, m10);
                pmaddwd(m11, ptr[r3 + 32]);
                // m11 = EO[0] : 7, 6, 5, 4

                pmaddwd(m2, ptr[r3 + 48]);
                // m2 = EO[1] : 3, 2, 1, 0

                pmaddwd(m10, ptr[r3 + 48]);
                // m10 = EO[1] : 7, 6, 5, 4

                movdqa(m4, m1);
                paddd(m4, m3);
                // m4 = E[0] : 3, 2, 1, 0

                movdqa(m12, m9);
                paddd(m12, m11);
                // m12 = E[0] : 7, 6, 5, 4

                psubd(m1, m3);
                // m1 = E[3] : 3, 2, 1, 0

                psubd(m9, m11);
                // m9 = E[3] : 7, 6, 5, 4

                movdqa(m3, m0);
                paddd(m3, m2);
                // m3 = E[1] : 3, 2, 1, 0

                movdqa(m13, m8);
                paddd(m13, m10);
                // m13 = E[1] : 7, 6, 5, 4

                psubd(m0, m2);
                // m0 = E[2] : 3, 2, 1, 0

                psubd(m8, m10);
                // m8 = E[2] : 7, 6, 5, 4

                // m4, m3, m0, m1 contain E[] for words 3, 2, 1, 0
                // m12, m13, m8, m9 contain E[] for words 7, 6, 5, 4
                // could add rounding here to E[] instead of later

                movdqa(m5, ptr[r1 + 1 * 2 * 8]);
                // m5 = src(7, 1), src(6, 1), src(5, 1), src(4, 1), src(3, 1), src(2, 1), src(1, 1), src(0, 1)

                movdqa(m2, ptr[r1 + 3 * 2 * 8]);
                // m2 = src(7, 3), src(6, 3), src(5, 3), src(4, 3), src(3, 3), src(2, 3), src(1, 3), src(0, 3)

                movdqa(m10, m5);
                punpckhwd(m10, m2);
                // m10 = src(7, 3), src(7, 1), src(6, 3), src(6, 1), src(5, 3), src(5, 1), src(4, 3), src(4, 1)

                punpcklwd(m5, m2);
                // m5 = src(3, 3), src(3, 1), src(2, 3), src(2, 1), src(1, 3), src(1, 1), src(0, 3), src(0, 1)

                movdqa(m6, ptr[r1 + 5 * 2 * 8]);
                // m6 = src(7, 5), src(6, 5), src(5, 5), src(4, 5), src(3, 5), src(2, 5), src(1, 5), src(0, 5)

                movdqa(m7, ptr[r1 + 7 * 2 * 8]);
                // m7 = src(7, 7), src(6, 7), src(5, 7), src(4, 7), src(3, 7), src(2, 7), src(1, 7), src(0, 7)

                movdqa(m11, m6);
                punpckhwd(m11, m7);
                // m11 = src(7, 7), src(7, 5), src(6, 7), src(6, 5), src(5, 7), src(5, 5), src(4, 7), src(4, 5)

                punpcklwd(m6, m7);
                // m6 = src(3, 7), src(3, 5), src(2, 7), src(2, 5), src(1, 7), src(1, 5), src(0, 7), src(0, 5)

                // m4, m3, m0, m1 contain E[] for words 3, 2, 1, 0
                // m12, m13, m8, m9 contain E[] for words 7, 6, 5, 4
                // m5 contains source rows 1 and 3 for words 3, 2, 1, 0
                // m6 contains source rows 5 and 7 for words 3, 2, 1, 0
                // m10 contains source rows 1 and 3 for words 7, 6, 5, 4
                // m11 contains source rows 5 and 7 for words 7, 6, 5, 4

                lea(r3, ptr[rip + cosine_inverse_8]);
                lea(r4, ptr[rip + dd_0040]);

                OUTPUTROWS(0, m4, m12);
                OUTPUTROWS(1, m3, m13);
                OUTPUTROWS(2, m0, m8);
                OUTPUTROWS(3, m1, m9);
            }

            void OUTPUTROWS(int n, Xbyak::Xmm const &p2, Xbyak::Xmm const &p3)
            {
                auto &r0 = arg64(0);
                auto &r3 = arg64(3);
                auto &r4 = arg64(4);

                auto &m0 = regXmm(0);
                auto &m1 = regXmm(1);
                auto &m2 = regXmm(2);
                auto &m3 = regXmm(3);
                auto &m4 = regXmm(4);
                auto &m5 = regXmm(5);
                auto &m6 = regXmm(6);
                auto &m7 = regXmm(7);
                auto &m8 = regXmm(8);
                auto &m9 = regXmm(9);
                auto &m10 = regXmm(10);
                auto &m11 = regXmm(11);
                auto &m12 = regXmm(12);
                auto &m13 = regXmm(13);
                auto &m14 = regXmm(14);
                auto &m15 = regXmm(15);

                movdqa(m2, m5);
                pmaddwd(m2, ptr[r3 + 32 * n + 0]);
                movdqa(m14, m10);
                pmaddwd(m14, ptr[r3 + 32 * n + 0]);
                movdqa(m7, m6);
                pmaddwd(m7, ptr[r3 + 32 * n + 16]);
                movdqa(m15, m11);
                pmaddwd(m15, ptr[r3 + 32 * n + 16]);
                paddd(m2, m7);
                paddd(m14, m15);
                // m2 = O[n] : 3, 2, 1, 0
                // m14 = O[n] : 7, 6, 5, 4

                movdqu(m7, ptr[r4]);
                movdqu(m15, ptr[r4]);
                paddd(m7, m2);
                paddd(m15, m14);
                paddd(m7, p2);
                paddd(m15, p3);
                psrad(m7, 7);
                psrad(m15, 7);
                // m7 = (E[n] + O[n] + add) >> shift : 3, 2, 1, 0
                // m15 = (E[n] + O[n] + add) >> shift : 7, 6, 5, 4

                packssdw(m7, m15);
                movdqa(ptr[r0 + 16 * n], m7);

                paddd(p2, ptr[r4]);
                paddd(p3, ptr[r4]);
                psubd(p2, m2);
                psubd(p3, m14);
                psrad(p2, 7);
                psrad(p3, 7);
                // p2 = (E[n] - O[n] + add) >> shift : 3, 2, 1, 0
                // p3 = (E[n] - O[n] + add) >> shift : 7, 6, 5, 4

                packssdw(p2, p3);
                movdqa(ptr[r0 + 16 * (7 - n)], p2);
            }

            void havoc_partial_butterfly_inverse_8h_ssse3()
            {
                // void transform_partial_butterfly_inverse_8h_ssse3(int16_t *dst, const int16_t *src, int shift);
                auto r0 = reg64(0);
                auto r1 = reg64(1);
                auto r4d = Xbyak::Reg32(reg64(4).getIdx());

                auto &m0 = regXmm(0);
                auto &m1 = regXmm(1);
                auto &m2 = regXmm(2);
                auto &m3 = regXmm(3);

                mov(r4d, 8);
                L("looph");
                {
                    movdqa(m2, ptr[r1]);
                    movdqa(m0, m2);
                    pshufb(m0, ptr[rip + shuffle_018945cd018945cd]);
                    // m0 = src[6], src[2], src[4], src[0], src[6], src[2], src[4], src[0]

                    pmaddwd(m0, ptr[rip + cosine_inverse_4_h]);
                    // m0 = EO[0], EE[0], EO[1], EE[1]

                    phsubd(m1, m0);
                    // m1 = E[3], E[2], ?, ?

                    pshufd(m0, m0, ORDER(1, 0, 3, 2));
                    // m0 = EO[1], EE[1], EO[0], EE[0]

                    phaddd(m0, m0);
                    // m0 = E[1], E[2], ?, ?

                    punpckhqdq(m0, m1);
                    // m0 = E[3], E[2], E[1], E[0]

                    paddd(m0, ptr[rip + dd_0800]);
                    // m0 = E[3]+add, E[2]+add, E[1]+add, E[0]+add

                    pshufb(m2, ptr[rip + shuffle_2367abef2367abef]);

                    movdqa(m3, m2);
                    pmaddwd(m3, ptr[rip + cosine_inverse_8_h]);
                    pmaddwd(m2, ptr[rip + cosine_inverse_8_h + 16]);
                    phaddd(m3, m2);
                    // m3 = O[3], O[2], O[1], O[0]
                    // m0 = E[3], E[2], E[1], E[0]

                    movdqa(m1, m0);
                    paddd(m1, m3);
                    // m1 = (E[k] + O[k] + add) >> shift : 3, 2, 1, 0

                    psrad(m1, 12);

                    psubd(m0, m3);
                    psrad(m0, 12);
                    // m0 = (E[k] - O[k] + add) >> shift : 3, 2, 1, 0

                    pshufd(m0, m0, ORDER(0, 1, 2, 3));
                    // m0 = (E[k] - O[k] + add) >> shift : 0, 1, 2, 3

                    packssdw(m1, m0);

                    movdqa(ptr[r0], m1);

                    add(r0, 16);
                    add(r1, 16);
                }
                dec(r4d);
                jg("looph");
            }

            void INTERLEAVE(int rowA, Xbyak::Xmm const &mA, int rowB, Xbyak::Xmm const &mB, Xbyak::Xmm const &mTemp)
            {
                auto &r1 = reg64(1);

                movdqa(mA, ptr[r1 + rowA * 2 * 16]);
                // mA = src(7, rowA), src(6, rowA), src(5, rowA), src(4, rowA), src(3, rowA), src(2, rowA), src(1, rowA), src(0, rowA)

                movdqa(mTemp, ptr[r1 + rowB * 2 * 16]);
                // mTemp = src(7, rowB), src(6, rowB), src(5, rowB), src(4, rowB), src(3, rowB), src(2, rowB), src(1, rowB), src(0, rowB)

                movdqa(mB, mA);
                punpckhwd(mB, mTemp);
                // %4 = src(7, rowB), src(7, rowA), src(6, rowB), src(6, rowA), src(5, rowB), src(5, rowA), src(4, rowB), src(4, rowA)

                punpcklwd(mA, mTemp);
                // mA = src(3, rowB), src(3, rowA), src(2, rowB), src(2, rowA), src(1, rowB), src(1, rowA), src(0, rowB), src(0, rowA)
            }

            void havoc_partial_butterfly_inverse_16v_ssse3()
            {
                auto &r0 = reg64(0);
                auto &r1 = reg64(1);
                auto &r2 = reg64(2);
                auto &r3 = reg64(3);
                auto &r4 = reg64(4);
                auto r4d = Xbyak::Reg32(r4.getIdx());
                auto r5d = Xbyak::Reg32(reg64(5).getIdx());

                auto &m0 = regXmm(0);
                auto &m1 = regXmm(1);
                auto &m2 = regXmm(2);
                auto &m3 = regXmm(3);
                auto &m4 = regXmm(4);
                auto &m5 = regXmm(5);
                auto &m6 = regXmm(6);
                auto &m7 = regXmm(7);
                auto &m8 = regXmm(8);
                auto &m9 = regXmm(9);
                auto &m10 = regXmm(10);
                auto &m11 = regXmm(11);
                auto &m12 = regXmm(12);
                auto &m13 = regXmm(13);

                sub(rsp, 8 * 2 * 16);

                mov(r5d, 2);
                L("loop_horizontal");
                {
                    INTERLEAVE(1, m0, 3, m1, m9);
                    INTERLEAVE(5, m2, 7, m3, m9);
                    INTERLEAVE(9, m4, 11, m5, m9);
                    INTERLEAVE(13, m6, 15, m7, m9);

                    // m0 to m7 contain interleaved odd source rows

                    lea(r3, ptr[rip + cosine_inverse_16]);
                    mov(r4d, 8);

                    L("loop_compute_Ok");
                    {
                        movdqa(m8, m0);
                        pmaddwd(m8, ptr[r3 + 0 * 16]);
                        movdqa(m9, m1);
                        pmaddwd(m9, ptr[r3 + 0 * 16]);

                        movdqa(m10, m2);
                        pmaddwd(m10, ptr[r3 + 1 * 16]);
                        movdqa(m11, m3);
                        pmaddwd(m11, ptr[r3 + 1 * 16]);
                        paddd(m8, m10);
                        paddd(m9, m11);

                        movdqa(m10, m4);
                        pmaddwd(m10, ptr[r3 + 2 * 16]);
                        movdqa(m11, m5);
                        pmaddwd(m11, ptr[r3 + 2 * 16]);
                        paddd(m8, m10);
                        paddd(m9, m11);

                        movdqa(m10, m6);
                        pmaddwd(m10, ptr[r3 + 3 * 16]);
                        movdqa(m11, m7);
                        pmaddwd(m11, ptr[r3 + 3 * 16]);
                        paddd(m8, m10);
                        paddd(m9, m11);

                        movdqa(ptr[rsp], m8);
                        // review: unusual usage of stack pointer
                        movdqa(ptr[rsp + 16], m9);

                        lea(r3, ptr[r3 + 4 * 16]);
                        add(rsp, 2 * 16);
                    }
                    dec(r4d);
                    jg("loop_compute_Ok");

                    sub(rsp, 8 * 2 * 16);

                    // now at esp is O[k] for 8 vertical positions (32-bit signed values)

                    // the following is almost identical to 8x8 inverse transform

                    lea(r3, ptr[rip + cosine_inverse_4]);
                    lea(r4, ptr[rip + dd_0040]);

                    movdqa(m8, ptr[r1 + 0 * 2 * 16]);
                    movdqa(m1, ptr[r1 + 8 * 2 * 16]);
                    movdqa(m0, m8);
                    punpcklwd(m0, m1);
                    punpckhwd(m8, m1);

                    movdqa(m1, m0);
                    pmaddwd(m1, ptr[r3]);
                    paddd(m1, ptr[r4]);
                    // m1 = EEE[0] + add: 3, 2, 1, 0

                    movdqa(m9, m8);
                    pmaddwd(m9, ptr[r3]);
                    paddd(m9, ptr[r4]);
                    // m9 = EEE[0] + add : 7, 6, 5, 4

                    pmaddwd(m0, ptr[r3 + 16]);
                    paddd(m0, ptr[r4]);
                    // m0 = EEE[1] + add : 3, 2, 1, 0

                    pmaddwd(m8, ptr[r3 + 16]);
                    paddd(m8, ptr[r4]);
                    // m8 = EEE[1] + add : 7, 6, 5, 4

                    movdqa(m10, ptr[r1 + 4 * 2 * 16]);
                    movdqa(m3, ptr[r1 + 12 * 2 * 16]);
                    movdqa(m2, m10);
                    punpcklwd(m2, m3);
                    punpckhwd(m10, m3);
                    movdqa(m3, m2);
                    pmaddwd(m3, ptr[r3 + 32]);
                    // m3 = EEO[0] : 3, 2, 1, 0

                    movdqa(m11, m10);
                    pmaddwd(m11, ptr[r3 + 32]);
                    // m11 = EEO[0] : 7, 6, 5, 4

                    pmaddwd(m2, ptr[r3 + 48]);
                    // m2 = EEO[1] : 3, 2, 1, 0

                    pmaddwd(m10, ptr[r3 + 48]);
                    // m10 = EEO[1] : 7, 6, 5, 4

                    movdqa(m4, m1);
                    paddd(m4, m3);
                    // m4 = EE[0] : 3, 2, 1, 0

                    movdqa(m12, m9);
                    paddd(m12, m11);
                    // m12 = EE[0] : 7, 6, 5, 4

                    psubd(m1, m3);
                    // m1 = EE[3] : 3, 2, 1, 0

                    psubd(m9, m11);
                    // m9 = EE[3] : 7, 6, 5, 4

                    movdqa(m3, m0);
                    paddd(m3, m2);
                    // m3 = E[1] : 3, 2, 1, 0

                    movdqa(m13, m8);
                    paddd(m13, m10);
                    // m13 = E[1] : 7, 6, 5, 4

                    psubd(m0, m2);
                    // m0 = E[2] : 3, 2, 1, 0

                    psubd(m8, m10);
                    // m8 = E[2] : 7, 6, 5, 4

                    // m4, m3, m0, m1 contain EE[] for words 3, 2, 1, 0
                    // m12, m13, m8, m9 contain EE[] for words 7, 6, 5, 4

                    movdqa(m5, ptr[r1 + 2 * 2 * 16]);
                    movdqa(m2, ptr[r1 + 6 * 2 * 16]);
                    movdqa(m10, m5);
                    punpckhwd(m10, m2);
                    punpcklwd(m5, m2);

                    movdqa(m6, ptr[r1 + 10 * 2 * 16]);
                    movdqa(m7, ptr[r1 + 14 * 2 * 16]);
                    movdqa(m11, m6);
                    punpckhwd(m11, m7);
                    punpcklwd(m6, m7);

                    // m4, m3, m0, m1 contain EE[] for words 3, 2, 1, 0
                    // m12, m13, m8, m9 contain EE[] for words 7, 6, 5, 4
                    // m5 contains source rows 2 and 6 for words 3, 2, 1, 0
                    // m6 contains source rows 10 and 14 for words 3, 2, 1, 0
                    // m10 contains source rows 2 and 6 for words 7, 6, 5, 4
                    // m11 contains source rows 10 and 14 for words 7, 6, 5, 4

                    lea(r3, ptr[rip + cosine_inverse_8]);

                    OUTPUTROWS_B(0, m4, m12);
                    OUTPUTROWS_B(1, m3, m13);
                    OUTPUTROWS_B(2, m0, m8);
                    OUTPUTROWS_B(3, m1, m9);

                    add(r0, 16);
                    add(r1, 16);
                }
                dec(r5d);
                jg("loop_horizontal");

                add(rsp, 8 * 2 * 16);
            }

            void OUTPUTROWS_B(int k, Xbyak::Xmm EEk_l, Xbyak::Xmm EEk_h)
            {
                // k = k
                // EE[k] is in registers EEk_l (3:0) and EEk_h (7:4)
                // O[k] is at [esp + 8 * 4 * k]
                // EO[k] = 4-tap FIR
                // E[k] = EE[k] + EO[k]
                // dst[k] = E[k] + O[k]
                // dst[15 - k] = E[k] - O[k]
                // E[7-k] = EE[k] - EO[k]
                // dst[7-k] = E[7-k] + O[7-k]
                // dst[k+8] = E[7-k] - O[7-k]

                auto &r0 = reg64(0);
                auto &r3 = reg64(3);
                auto &m2 = regXmm(2);
                auto &m5 = regXmm(5);
                auto &m6 = regXmm(6);
                auto &m7 = regXmm(7);
                auto &m10 = regXmm(10);
                auto &m11 = regXmm(11);
                auto &m14 = regXmm(14);
                auto &m15 = regXmm(15);

                movdqa(m2, m5);
                pmaddwd(m2, ptr[r3+32*k+0]);
                movdqa(m14, m10);
                pmaddwd(m14, ptr[r3+32*k+0]);
                movdqa(m7, m6);
                pmaddwd(m7, ptr[r3+32*k+16]);
                movdqa(m15, m11);
                pmaddwd(m15, ptr[r3+32*k+16]);
                paddd(m2, m7);
                paddd(m14, m15);
                // m2 = EO[k] : 3, 2, 1, 0
                // m14 = EO[k] : 7, 6, 5, 4

                movdqa(m7, m2 );
                movdqa(m15, m14 );
                paddd(m7, EEk_l);
                paddd(m15, EEk_h);
                // m7 = (E[k] + add) = (EE[k] + EO[k] + add) : 3, 2, 1, 0
                // m15 = (E[k] + add) = (EE[k] + EO[k] + add) : 7, 6, 5, 4

                // at this point, all registers are in use

                psubd(EEk_l, m2);
                psubd(EEk_h, m14);
                // EEk_l = (E[7-k] + add) = (EE[k] - EO[k] + add) : 3, 2, 1, 0
                // EEk_h = (E[7-k] + add) = (EE[k] - EO[k] + add) : 7, 6, 5, 4

                movdqa(m2, m7);
                paddd(m2, ptr[rsp+8*4*k]);
                // m2 = (E[k] + O[k] + add) : 3, 2, 1, 0

                movdqa(m14, m15);
                paddd(m14, ptr[rsp+8*4*k + 16]);
                // m14 = (E[k] + O[k] + add) : 7, 6, 5, 4

                psrad(m2, 7);
                psrad(m14, 7);
                packssdw(m2, m14);
                movdqa(ptr[r0+16*2*k], m2);

                psubd(m7, ptr[rsp+8*4*k]);
                // m7 = (E[k] - O[k] + add) : 3, 2, 1, 0

                psubd(m15, ptr[rsp+8*4*k + 16]);
                // m15 = (E[k] - O[k] + add) : 7, 6, 5, 4

                psrad(m7, 7);
                psrad(m15, 7);
                packssdw(m7, m15);
                movdqa(ptr[r0+16*2*(15-k)], m7);

                movdqa(m2, EEk_l);
                paddd(m2, ptr[rsp+8*4*(7-k)]);
                // m2 = (E[7-k] + O[7-k] + add) : 3, 2, 1, 0

                movdqa(m14, EEk_h);
                paddd(m14, ptr[rsp+8*4*(7-k) + 16]);
                // m14 = (E[7-k] + O[7-k] + add) : 7, 6, 5, 4

                psrad(m2, 7);
                psrad(m14, 7);
                packssdw(m2, m14);
                movdqa(ptr[r0+16*2*(7-k)], m2);

                psubd(EEk_l, ptr[rsp+8*4*(7-k)]);
                // EEk_l = (E[k] - O[k] + add) : 3, 2, 1, 0

                psubd(EEk_h, ptr[rsp+8*4*(7-k) + 16]);
                // EEk_h = (E[k] - O[k] + add) : 7, 6, 5, 4

                psrad(EEk_l, 7);
                psrad(EEk_h, 7);
                packssdw(EEk_l, EEk_h);
                movdqa(ptr[r0+16*2*(8+k)], EEk_l);
            }

            void havoc_partial_butterfly_inverse_16h_ssse3()
            {
                //// void transform_partial_butterfly_inverse_16h_ssse3(int16_t *dst, const int16_t *src, int shift)//
                //INIT_XMM ssse3
                //cglobal partial_butterfly_inverse_16h, 3, 5, 16
                auto &r0 = reg64(0);
                auto &r1 = reg64(1);
                auto r4d = Xbyak::Reg32(reg64(4).getIdx());
                auto &m0 = regXmm(0);
                auto &m1 = regXmm(1);
                auto &m2 = regXmm(2);
                auto &m3 = regXmm(3);
                auto &m4 = regXmm(4);
                auto &m5 = regXmm(5);
                auto &m6 = regXmm(6);
                auto &m7 = regXmm(7);
                auto &m8 = regXmm(8);
                auto &m9 = regXmm(9);
                auto &m10 = regXmm(10);
                auto &m11 = regXmm(11);
                auto &m12 = regXmm(12);

                mov(r4d, 16);
                L("loop16h");
                {
                    movdqa(m0, ptr[r1]);
                    // m0 = src[7:0]

                    movdqa(m1, ptr[r1 + 16]);
                    // m1 = src[15:8]

                    add(r1, 32);

                    // 8 x 8-tap horizontal filters on odd positions

                    movdqa(m6, m0);
                    pshufb(m6, ptr[rip + shuffle_2367236723672367]); // 3, 1
                    movdqa(m7, m0);
                    pshufb(m7, ptr[rip + shuffle_abefabefabefabef]); // 7, 5
                    movdqa(m8, m1);
                    pshufb(m8, ptr[rip + shuffle_2367236723672367]);  // 11, 9
                    movdqa(m9, m1);
                    pshufb(m9, ptr[rip + shuffle_abefabefabefabef]); // 15, 13

                    movdqa(m2, m6);
                    pmaddwd(m2, ptr[rip + cosine_inverse_16_h + 0 * 16]);
                    pmaddwd(m6, ptr[rip + cosine_inverse_16_h + 1 * 16]);

                    movdqa(m3, m7);
                    pmaddwd(m3, ptr[rip + cosine_inverse_16_h + 2 * 16]);
                    pmaddwd(m7, ptr[rip + cosine_inverse_16_h + 3 * 16]);

                    movdqa(m4, m8);
                    pmaddwd(m4, ptr[rip + cosine_inverse_16_h + 4 * 16]);
                    pmaddwd(m8, ptr[rip + cosine_inverse_16_h + 5 * 16]);

                    movdqa(m5, m9);
                    pmaddwd(m5, ptr[rip + cosine_inverse_16_h + 6 * 16]);
                    pmaddwd(m9, ptr[rip + cosine_inverse_16_h + 7 * 16]);

                    paddd(m2, m3);
                    paddd(m4, m5);
                    paddd(m2, m4);
                    // m2 = O[3:0]

                    paddd(m6, m7);
                    paddd(m8, m9);
                    paddd(m6, m8);
                    // m6 = O[7:4]

                    movdqa(m3, m0);
                    pshufb(m3, ptr[rip + shuffle_45cd45cd45cd45cd]); // 6, 2
                    movdqa(m4, m1);
                    pshufb(m4, ptr[rip + shuffle_45cd45cd45cd45cd]);  // 14, 10

                    pmaddwd(m3, ptr[rip + cosine_inverse_8_h2 + 0 * 16]);
                    pmaddwd(m4, ptr[rip + cosine_inverse_8_h2 + 1 * 16]);

                    paddd(m3, m4);
                    // m3 = EO[3:0]

                    movdqa(m5, m0);
                    pshufb(m5, ptr[rip + shuffle_0189018901890189]);
                    // m5 = src[4, 0, 4, 0, 4, 0, 4, 0]

                    movdqa(m7, m1);
                    pshufb(m7, ptr[rip + shuffle_0189018901890189]);
                    // m7 = src[12, 8, 12, 8, 12, 8, 12, 8]

                    punpcklwd(m5, m7);
                    // m5 = src[12, 4, 8, 0, 12, 4, 8, 0]

                    pmaddwd(m5, ptr[rip + cosine_inverse_4_h]);
                    // m5 = EEO[0], EEE[0], EEO[1], EEE[1]

                    phsubd(m4, m5);
                    // m4 = EE[3], EE[2], ?, ?

                    pshufd(m5, m5, ORDER(1, 0, 3, 2));
                    // m5 = EEO[1], EEE[1], EEO[0], EEE[0]

                    phaddd(m5, m5);
                    // m5 = EE[1], EE[2], EE[1], EE[2]

                    punpckhqdq(m5, m4);
                    // m5 = EE[3:0]

                    paddd(m5, ptr[rip + dd_0800]);
                    // m5 = EE[3:0] + add

                    movdqa(m0, m3);
                    paddd(m0, m5);
                    // m0 = E[3:0] + add

                    psubd(m5, m3);
                    // m5 = E[4:7] + add

                    movdqa(m10, m5);
                    pshufd(m10, m5, ORDER(0, 1, 2, 3));
                    // m11 = E[7:4] + add

                    movdqa(m11, m10);
                    paddd(m11, m6);
                    psrad(m11, 12);
                    // m11 = dst[7:4]

                    movdqa(m12, m0);
                    paddd(m12, m2);
                    psrad(m12, 12);
                    // m12 = dst[3:0]

                    packssdw(m12, m11);
                    // m12 = dst[7:0]

                    movdqa(ptr[r0], m12);

                    movdqa(m11, m10);
                    psubd(m11, m6);
                    psrad(m11, 12);
                    // m11 = dst[8:11]

                    movdqa(m12, m0);
                    psubd(m12, m2);
                    psrad(m12, 12);
                    // m12 = dst[12:15]

                    packssdw(m12, m11);
                    // m12 = dst[8:15]

                    pshufb(m12, ptr[rip + shuffle_efcdab8967452301]);

                    movdqa(ptr[r0 + 16], m12);

                    add(r0, 32);
                }
                dec(r4d);
                jg("loop16h");
            }

            void assemble() override
                    {
#if USE_F265_DERIVED
                if (this->isa() & HAVOC_AVX2)
                {
                    if (this->log2TrafoSize == 2) this->assemble4x4Avx2();
                    if (this->log2TrafoSize == 3) this->assemble8x8Avx2();
                    if (this->log2TrafoSize == 4) this->assemble16x16Avx2();
                    if (this->log2TrafoSize == 5) this->assemble32x32Avx2();
                    return;
                }
#endif
                // void havoc_inverse_transform_add(uint8_t *dst, intptr_t stride_dst, const uint8_t *pred, intptr_t stride_pred, const int16_t *coeffs);

                auto &r0 = arg64(0);
                auto &r1 = arg64(1);
                auto &r2 = arg64(2);
                auto &r3 = arg64(3);
                auto &r4 = arg64(4);

                const auto intermediateBufferSize = 2 * (1 << (2 * this->log2TrafoSize));

                this->stackSize = 2 * intermediateBufferSize + 32;

                // save arguments
                mov(ptr[rsp], r0); // dst
                mov(ptr[rsp + 8], r1); // stride_dst
                mov(ptr[rsp + 16], r2); // pred
                mov(ptr[rsp + 24], r3); // stride_pred

                lea(r0, ptr[rsp + 32]); // temp[0]
                mov(r1, r4); // coeffs

                if (this->log2TrafoSize == 3)
                    this->havoc_partial_butterfly_inverse_8v_ssse3();
                else
                {
                    this->havoc_partial_butterfly_inverse_16v_ssse3();
                }

                lea(r0, ptr[rsp + 32 + intermediateBufferSize]); // temp[1]
                lea(r1, ptr[rsp + 32]); // temp[0]

                if (this->log2TrafoSize == 3)
                    this->havoc_partial_butterfly_inverse_8h_ssse3();
                else
                    this->havoc_partial_butterfly_inverse_16h_ssse3();

                // add residual
                mov(r0, ptr[rsp]); // dst
                mov(r1, ptr[rsp + 8]); // stride_dst
                mov(r2, ptr[rsp + 16]); // pred
                mov(r3, ptr[rsp + 24]); // stride_pred
                lea(r4, ptr[rsp + 32 + intermediateBufferSize]); // temp[1] (residual)

                auto &m0 = regXmm(0);
                auto &m1 = regXmm(1);
                auto &m2 = regXmm(2);
                pxor(m2, m2);

                if (this->log2TrafoSize == 3)
                    for (int y = 0; y < 8; ++y)
                    {
                        movdqu(m0, ptr[r4 + 16 * y]);
                        movq(m1, ptr[r2]);
                        punpcklbw(m1, m2);
                        paddsw(m0, m1);
                        packuswb(m0, m0);
                        movq(ptr[r0], m0);
                        lea(r0, ptr[r0 + r1]);
                        lea(r2, ptr[r2 + r3]);
                    }

                if (this->log2TrafoSize == 4)
                    for (int y = 0; y < 16; ++y)
                    {
                        auto &m3 = regXmm(3);
                        auto &m4 = regXmm(4);

                        movdqa(m0, ptr[r4 + 32 * y]);
                        movdqa(m1, ptr[r4 + 32 * y + 16]);
                        movdqa(m3, ptr[r2]);
                        movdqa(m4, m3);
                        punpcklbw(m3, m2);
                        punpckhbw(m4, m2);
                        paddsw(m0, m3);
                        paddsw(m1, m4);
                        packuswb(m0, m1);
                        movdqa(ptr[r0], m0);
                        lea(r0, ptr[r0 + r1]);
                        lea(r2, ptr[r2 + r3]);
                    }


                //	havoc_add_residual(8, dst, stride_dst, pred, stride_pred, temp[1]);
                //this->havoc_add_residual();
                    }

#if USE_F265_DERIVED

            // Interleave two registers.
            // a1: lower output, a2 : higher output, a3 : input 0, a4 : input 1.
#define INTERLEAVE_COL(a1, a2, a3, a4) \
        vpunpcklwd(a1, a3, a4); \
        vpunpckhwd(a2, a3, a4); \


            // Expand the DCT factors from bytes to words and store them in spill memory. We
            // could store them expanded in the data segment, but that would impact the
            // cache. We could also expand them on the fly, but that kills performance. This
            // operation has a limited impact on performance, so it's probably for the best.
            // a1: destination, a2: source, a3: loop counter, a4: tmp, a5: factor size.
#define EXPAND_FACTOR(a1, a2, a3, a4, a5) \
        { \
            xor(a3, a3); \
            Xbyak::Label loop_expand; \
            L(loop_expand); \
            vpmovsxbw(a4, ptr[a2 + a3]);           /* Load and expand the factors.*/ \
        vmovdqu(ptr[a1 + a3 * 2], a4);         /* Store. */ \
        add(a3, 16); \
        cmp(a3, a5); \
        jne(loop_expand); \
        } \

            void assemble32x32Avx2()
        {
            // from f265
            // IDCT 32x32.
            //
            // Input parameters:
            // - g0:     destination.
            // - g1:     destination stride.
            // - g2:     prediction.
            // - g3:     prediction stride.
            // - g4:     coefficients.
            // - g5:     spill buffer.
            //DEFFUN f265_lbd_idct_32_avx2, ia=6, at=848488, fa=0, ti=3, tv=16, ym=1

            auto &g0 = arg64(0);
            auto &g1 = arg64(1);
            auto &g2 = arg64(2);
            auto &g3 = arg64(3);
            auto &g4 = arg64(4);
            auto &g5 = reg64(5);
            auto &g6 = reg64(6);
            auto &g7 = reg64(7);
            auto &g8 = reg64(8);

            this->stackSize = 128 * 16;
            mov(g5, rsp);

            auto &y0 = ymm0; auto &x0 = regXmm(0);
            auto &y1 = ymm1; auto &x1 = regXmm(1);
            auto &y2 = ymm2; auto &x2 = regXmm(2);
            auto &y3 = ymm3; auto &x3 = regXmm(3);
            auto &y4 = ymm4; auto &x4 = regXmm(4);
            auto &y5 = ymm5; auto &x5 = regXmm(5);
            auto &y6 = ymm6; auto &x6 = regXmm(6);
            auto &y7 = ymm7; auto &x7 = regXmm(7);
            auto &y8 = ymm8; auto &x8 = regXmm(8);
            auto &y9 = ymm9; auto &x9 = regXmm(9);
            auto &y10 = ymm10; auto &x10 = regXmm(10);
            auto &y11 = ymm11; auto &x11 = regXmm(11);
            auto &y12 = ymm12; auto &x12 = regXmm(12);
            auto &y13 = ymm13; auto &x13 = regXmm(13);
            auto &y14 = ymm14; auto &x14 = regXmm(14);
            auto &y15 = ymm15; auto &x15 = regXmm(15);

            // Note : To make it easier, the notation used in the C
            // code is used to describe the output terms in each code section.

            // Load 4 consecutive, ascending rows of input.
            // a1: row offset, a2: row spacing, a3: g7 value, a4: g6 value.
#define LOAD_COL(a1, a2, a3, a4) \
        vmovdqu(y0, ptr[g4 + (0 + a1)*32]); \
        vmovdqu(y1, ptr[g4 + (1 + a1)*32]); \
        vmovdqu(y2, ptr[g4 + (0 + a2 + a1)*32]); \
        vmovdqu(y3, ptr[g4 + (1 + a2 + a1)*32]); \
        lea(g7, ptr[rip + pat_idct32_pass1 + a3]); \
        mov(g6, a4); \
        INTERLEAVE_COL(y15, y14, y0, y2) \
        INTERLEAVE_COL(y13, y12, y1, y3) \
        vmovdqu(y4, ptr[g4 + (0 + 2*a2 + a1)*32]); \
        vmovdqu(y5, ptr[g4 + (1 + 2*a2 + a1)*32]); \
        vmovdqu(y6, ptr[g4 + (0 + 3*a2 + a1)*32]); \
        vmovdqu(y7, ptr[g4 + (1 + 3*a2 + a1)*32]); \
        INTERLEAVE_COL(y11, y10, y4, y6) \
        INTERLEAVE_COL(y9, y8, y5, y7) \


            // Load 4 consecutive, descending rows of input, and negate alternate rows. The alternating is done to include
            // the sign of the IDCT factors. This helps reuse the IDCT factors.
            // a1: row offset, a2: row spacing, a3: g7 value, a4: g6 value.
#define LOAD_COL_REV(a1, a2, a3, a4) \
        vpbroadcastd(y4, ptr[rip + pat_idct32_sign]); \
        vmovdqu(y0, ptr[g4 + (1 + 3 * a2 + a1) * 32]); \
        vmovdqu(y1, ptr[g4 + (0 + 3 * a2 + a1) * 32]); \
        vmovdqu(y2, ptr[g4 + (1 + 2 * a2 + a1) * 32]); \
        vmovdqu(y3, ptr[g4 + (0 + 2 * a2 + a1) * 32]); \
        vpsignw(y2, y2, y4); \
        vpsignw(y3, y3, y4); \
        INTERLEAVE_COL(y15, y14, y1, y3) \
        INTERLEAVE_COL(y13, y12, y0, y2) \
        vmovdqu(y0, ptr[g4 + (1 + a2 + a1)*32]); \
        vmovdqu(y1, ptr[g4 + (0 + a2 + a1)*32]); \
        vmovdqu(y2, ptr[g4 + (1 + a1)*32]); \
        vmovdqu(y3, ptr[g4 + (0 + a1)*32]); \
        vpsignw(y2, y2, y4); \
        vpsignw(y3, y3, y4); \
        INTERLEAVE_COL(y11, y10, y1, y3) \
        INTERLEAVE_COL(y9, y8, y0, y2) \
        lea(g7, ptr[rip + pat_idct32_pass1 + a3]); \
        mov(g6, a4); \


            // IDCT Pass 1.
            // Each row of input requires 2 ymm registers => 8 ymm registers for 4 rows of inputs.

            // First, process the odd rows to obtain s0, s1 .. s15.
            // Since we can process only 4 rows at a time, and we have 16 odd rows, we do the processing 4 times.
            // 1. Process rows 1, 3, 5 and 7 (i.e. srcptr[32], srcptr[96], srcptr[160] and srcptr[224]).
            LOAD_COL(2, 4, 0, 0);
            mov(g8, 128 * 16);

            call(loop_idct32_pass1);

            // 2. Process rows 9, 11, 13, 15 (i.e. srcptr[288], srcptr[352], srcptr[416] and srcptr[480]).
            LOAD_COL(18, 4, 128, 0);

            align(16);
            // Added for speed optimization: Aligns the loops '.loop_pass1_odd',
            // '.loop_pass1_even1a' and the subroutines '.loop_idct32_pass1' and
            // '.loop_idct32_pass2_odd' at 16-byte memory boundaries.

            L("loop_pass1_odd");
            {
                LOAD_DATA(g6, 0);
                vpbroadcastd(y6, ptr[g7]);
                vpbroadcastd(y7, ptr[g7 + 4]);
                MULT_ADD(y0, y15, y11);
                MULT_ADD(y1, y14, y10);
                MULT_ADD(y2, y13, y9);
                MULT_ADD(y3, y12, y8);
                STORE_DATA(g6, 0);
                add(g6, 128);
                add(g7, 8);
            }
            cmp(g6, 128 * 16);
            jne("loop_pass1_odd");

            // 3. Process rows 23, 21, 19, 17 (i.e. srcptr[736], srcptr[672], srcptr[608] and srcptr[544]).
            // This seems to loop over too much data
            LOAD_COL_REV(34, 4, 248, 0);
            call(loop_idct32_pass1_odd);

            // 4. Process rows 31, 29, 27, 25 (i.e. srcptr[992], srcptr[928], srcptr[864] and srcptr[800]).
            // This seems to loop over too much data
            LOAD_COL_REV(50, 4, 120, 0);
            call(loop_idct32_pass1_odd);

            // Now process the even rows. Split this into 4.
            // 1. Process rows 0, 16, 8, 24 (srcptr[0], srcptr[512], srcptr[256] and srcptr[768]) to calculate a_0, a_1, a_2, a_3.
            vmovdqu(y0, ptr[g4]);                // Row 0.
            vmovdqu(y1, ptr[g4 + 32]);
            vmovdqu(y2, ptr[g4 + 1024]);         // Row 16.
            vmovdqu(y3, ptr[g4 + 1056]);
            xor (g6, g6);
            mov(g8, 384);
            INTERLEAVE_COL(y15, y14, y0, y2);
            INTERLEAVE_COL(y13, y12, y1, y3);
            vmovdqu(y4, ptr[g4 + 512]);          // Row 8.
            vmovdqu(y5, ptr[g4 + 544]);
            vmovdqu(y6, ptr[g4 + 1536]);         // Row 24.
            vmovdqu(y7, ptr[g4 + 1568]);
            lea(g7, ptr[rip + pat_idct32_pass1 + 352]);
            INTERLEAVE_COL(y11, y10, y4, y6);
            INTERLEAVE_COL(y9, y8, y5, y7);
            vpbroadcastd(y6, ptr[rip + pat_dw_64]);        // Bias.

            L("loop_pass1_even1a");
            {
                vpbroadcastd(y0, ptr[g7]);                // Load the DCT multiplication factors.
                vpbroadcastd(y1, ptr[g7 + 4]);

                // a1: input 0, a2: input 1, a3: store location's index.
#define PROCESS_ROW(a1, a2, a3) \
        vpmaddwd(y2, a1, y0); \
        vpmaddwd(y3, a2, y1); \
        vpaddd(y4, y2, y3); \
        vpsubd(y5, y2, y3); \
        vpaddd(y4, y4, y6); /* Add the bias. */ \
        vpaddd(y5, y5, y6); \
        vmovdqu(ptr[g5 + g6 + 2048 + a3], y4); \
        vmovdqu(ptr[g5 + g8 + 2048 + a3], y5);

                PROCESS_ROW(y15, y11, 0);
                PROCESS_ROW(y14, y10, 32);
                PROCESS_ROW(y13, y9, 64);
                PROCESS_ROW(y12, y8, 96);
#undef PROCESS_ROW

                add(g6, 128);
                sub(g8, 128);
                add(g7, 8);
            }
            cmp(g6, 128 * 2);
            jne("loop_pass1_even1a");

            // 2. Process rows 4, 12, 20, 28 (srcptr[128], srcptr[384], srcptr[640] and srcptr[896]) to get a_4, a_5, a_6 and a_7.
            // Also calculate (a_0 +/- a_4), (a_1 +/- a_5), (a_2 +/- a_6), (a_3 +/- a_7).
            LOAD_COL(8, 16, 320, 0);
            mov(g8, 896);

            L("loop_pass1_even1b");
            {
                vpbroadcastd(y6, ptr[g7]);
                vpbroadcastd(y7, ptr[g7 + 4]);
                // a1: input 0, a2: input 1, a3: store location's index.
#define PROCESS_ROW(a1, a2, a3) \
        vpmaddwd(y2, a1, y6); \
        vpmaddwd(y3, a2, y7); \
        vpaddd(y4, y2, y3); \
        vmovdqu(y5, ptr[g5 + 2048 + g6 + a3]); \
        vpaddd(y2, y4, y5); \
        vpsubd(y4, y5, y4); \
        vmovdqu(ptr[g5 + 2048 + g6 + a3], y2); \
        vmovdqu(ptr[g5 + 2048 + g8 + a3], y4);

                PROCESS_ROW(y15, y11, 0);
                PROCESS_ROW(y14, y10, 32);
                PROCESS_ROW(y13, y9, 64);
                PROCESS_ROW(y12, y8, 96);
#undef PROCESS_ROW

                add(g6, 128);
                sub(g8, 128);
                add(g7, 8);
            }
            cmp(g6, 128 * 4);
            jne("loop_pass1_even1b");

            // 3. Process rows 2, 6, 10, 14 to calculate the 4 terms of e_0, e_1, e_2, e_3, e_4, e_5, e_6, and e_7.
            LOAD_COL(4, 8, 256, 2048 + 1024);
            mov(g8, 2048 + 1024 + 128 * 8);
            call(loop_idct32_pass1);

            // 4. Process rows 30, 26, 22, 18 to get e_0, e_1, e_2, e_3, e_4, e_5, e_6, and e_7. Also obtain a0 .. a15.
            LOAD_COL_REV(36, 8, 312, 2048);
            //xor g8, g8
            L("loop_pass1_even1d");
            {
                // a1: input 0, a2: input 1, a3: memory index for storing and loading.
#define LOAD_ADD_STORE(a1, a2, a3) \
        vmovdqu(y4, ptr[g5 + g6 + a3]); \
        vmovdqu(y5, ptr[g5 + g6 + a3 + 32]); \
        vpaddd(y6, a1, y4); \
        vpaddd(y7, a2, y5); \
        vmovdqu(ptr[g5 + g6 + a3], y6); \
        vmovdqu(ptr[g5 + g6 + a3 + 32], y7); \
        vpsubd(y4, y4, a1); \
        vpsubd(y5, y5, a2); \
        vmovdqu(ptr[g5 + g6 + 1024 + a3], y4); \
        vmovdqu(ptr[g5 + g6 + 1024 + a3 + 32], y5);

                LOAD_DATA(g6, 1024);                // To obtain e_0, e_2, e_4, e_6.
                vpbroadcastd(y6, ptr[g7]);
                vpbroadcastd(y7, ptr[g7 + 4]);
                MULT_ADD(y0, y15, y11);
                MULT_ADD(y1, y14, y10);
                MULT_ADD(y2, y13, y9);
                MULT_ADD(y3, y12, y8);
                LOAD_ADD_STORE(y0, y1, 0);              // Calculate a0 & a15, a2 & a13 ...
                LOAD_ADD_STORE(y2, y3, 64);
                add(g6, 128);
                sub(g7, 8);
                vpbroadcastd(y6, ptr[g7]);
                vpbroadcastd(y7, ptr[g7 + 4]);
                LOAD_DATA(g6, 1024);               // To obtain e_1, e_3, e_5, e_7.
                MULT_SUB(y0, y15, y11);
                MULT_SUB(y1, y14, y10);
                MULT_SUB(y2, y13, y9);
                MULT_SUB(y3, y12, y8);
                LOAD_ADD_STORE(y0, y1, 0);            // Calculate a1 & a14, a3 & a12 ...
                LOAD_ADD_STORE(y2, y3, 64);

#undef  LOAD_ADD_STORE
                add(g6, 128);
                sub(g7, 8);
            }
            cmp(g6, 2048 + 128 * 8);
            jne("loop_pass1_even1d");

            // Combine the odd (s0, s1 .. s15) and even (a0, a1, .. a15) terms.
            // a1: g4 value, a2: g6 value, a3: g7 value, a4: g8 value.
#define SET_GPR(a1, a2, a3, a4) \
        mov(g4, a1); \
        mov(g6, a2); \
        mov(g7, a3); \
        mov(g8, a4);

            vbroadcasti128(y6, ptr[rip + pat_idct16_shuf1]);
            SET_GPR(128, 0, 2048, 1024);
            call(loop_idct32_pass1_combine);
            SET_GPR(-128, 1024, 3968, 1024 + 1024);
            call(loop_idct32_pass1_combine);
#undef SET_GPR

            // Now, perform IDCT Pass 2.
            // In the previous pass (IDCT32 Pass 1), after combining, the even terms were written at
            // even locations (i.e. ptr[g5 + 0], ptr[g5 + 64], .. ptr[g5 + 960]), while the odd terms were written at odd locations (i.e.
            // ptr[g5 + 32], ptr[g5 + 96], .. ptr[g5 + 992]). This division helps reduce the number of multiplications.
            // Further note that the even terms were split so that all terms that enable calculation of a_x were grouped together
            // and those that enable calculation of e_x (x = 0, 1 .. 7) were grouped together.

            // First process the even terms.
            // a1: g7 value, a2: g8 value.
#define LOAD_IDCT_FACT(a1, a2) \
        xor(g6, g6); \
        lea(g7, ptr[rip + pat_idct32_pass2+ a1]); \
        lea(g8, ptr[g5 + a2]); \
        vpmovsxbw(y15, ptr[g7]); \
        vpmovsxbw(y14, ptr[g7 + 16]); \
        vpmovsxbw(y13, ptr[g7 + 32]); \
        vpmovsxbw(y12, ptr[g7 + 48]); \
        vpmovsxbw(y11, ptr[g7 + 64]); \
        vpmovsxbw(y10, ptr[g7 + 80]); \
        vpmovsxbw(y9, ptr[g7 + 96]); \
        vpmovsxbw(y8, ptr[g7 + 112]); \

            LOAD_IDCT_FACT(0, 2048);
            vmovdqu(y2, ptr[rip + pat_idct32_shuf1]);
            vpbroadcastd(y1, ptr[rip + pat_dw_2048]);       // Bias.

            L("loop_pass2_even");
            {
                vpermq(y0, ptr[g5 + g6], 0xD8);
                MULT_HADD(y7, y6, y15, y14);
                MULT_HADD(y5, y4, y13, y12);
                MULT_HADD(y6, y4, y11, y10);
                MULT_HADD(y4, y3, y9, y8);
                vphaddd(y7, y7, y5);
                vphaddd(y6, y6, y4);
                vperm2i128(y5, y7, y6, 0x20);
                vperm2i128(y4, y7, y6, 0x31);
                vpaddd(y0, y5, y4);              // Calculate a0, a1, a2, ... a7.
                vpsubd(y3, y5, y4);              // Calculate a15, a14, a13 ... a8.
                vpermd(y3, y2, y3);              // Rearrange to give: a8, a9, ... a15.
                vpaddd(y0, y0, y1);              // Add the bias.
                vpaddd(y3, y3, y1);
                vmovdqu(ptr[g8 + g6], y0);
                vmovdqu(ptr[g8 + g6 + 32], y3);
                add(g6, 64);
            }
            cmp(g6, 64 * 32);
            jne("loop_pass2_even");


            // Calculate the odd terms.
            LOAD_IDCT_FACT(128, 0);
            call(loop_idct32_pass2_odd);
            LOAD_IDCT_FACT(256, 32);
            call(loop_idct32_pass2_odd);

            // Combine the even and odd terms.

            // a1: g4 value, a2: g6 value, a3: g7 value.
#define SET_GPR(a1, a2, a3) \
        mov(g4, a1); \
        mov(g6, a2); \
        mov(g7, a3); \

            vbroadcasti128(y10, ptr[rip + pat_idct32_shuf2]);
            SET_GPR(0, 128, 128 * 16);
            call(loop_idct32_pass2_combine2);
            SET_GPR(1984, -128, -64);
            call(loop_idct32_pass2_combine2);


#undef SET_GPR

#undef LOAD_COL
        }

        void assemble16x16Avx2()
        {
            // from f265

            // IDCT 16x16.
            // DEFFUN(f265_lbd_idct_16_avx2, ia=6, at=848488, fa=0, ti=3, tv=16, ym=1

            auto &g0 = arg64(0);
            auto &g1 = arg64(1);
            auto &g2 = arg64(2);
            auto &g3 = arg64(3);
            auto &g4 = arg64(4);
            auto &g5 = reg64(5);
            auto &g6 = reg64(6);
            auto &g7 = reg64(7);
            auto &g8 = reg64(8);

            // Spill buffer layout:
            // - 512 bytes for pass 1 even rows, and pass 2 results.
            // - 512 bytes for pass 1 results.
            // - 256 bytes for factor expansion.
            this->stackSize = 512 + 512 + 512;
            mov(g5, rsp);

            auto &y0 = ymm0; auto &x0 = regXmm(0);
            auto &y1 = ymm1; auto &x1 = regXmm(1);
            auto &y2 = ymm2; auto &x2 = regXmm(2);
            auto &y3 = ymm3; auto &x3 = regXmm(3);
            auto &y4 = ymm4; auto &x4 = regXmm(4);
            auto &y5 = ymm5; auto &x5 = regXmm(5);
            auto &y6 = ymm6; auto &x6 = regXmm(6);
            auto &y7 = ymm7; auto &x7 = regXmm(7);
            auto &y8 = ymm8; auto &x8 = regXmm(8);
            auto &y9 = ymm9; auto &x9 = regXmm(9);
            auto &y10 = ymm10; auto &x10 = regXmm(10);
            auto &y11 = ymm11; auto &x11 = regXmm(11);
            auto &y12 = ymm12; auto &x12 = regXmm(12);
            auto &y13 = ymm13; auto &x13 = regXmm(13);
            auto &y14 = ymm14; auto &x14 = regXmm(14);
            auto &y15 = ymm15; auto &x15 = regXmm(15);

            // DCT pass 1.

            // Expand the factors for pass 1.
            lea(g8, ptr[g5 + 1024]);         // Factor destination.
            lea(g7, ptr[rip + pat_idct16_pass1]);

            // Factor source.
            EXPAND_FACTOR(g8, g7, g6, y0, 16 * 16);

#if 1

            // Load half of the rows in 8 registers.
            // a1: row offset.
#define LOAD_COL(a1) \
        vmovdqu(y0, ptr[g4 + 0*64 + a1]); \
        vmovdqu(y1, ptr[g4 + 1*64 + a1]); \
        vmovdqu(y2, ptr[g4 + 2*64 + a1]); \
        vmovdqu(y3, ptr[g4 + 3*64 + a1]); \
        vmovdqu(y4, ptr[g4 + 4*64 + a1]); \
        vmovdqu(y5, ptr[g4 + 5*64 + a1]); \
        INTERLEAVE_COL(y14, y13, y0, y1); \
        INTERLEAVE_COL(y12, y11, y2, y3); \
        INTERLEAVE_COL(y10, y9, y4, y5); \
        vmovdqu(y0, ptr[g4 + 6*64 + a1]); \
        vmovdqu(y1, ptr[g4 + 7*64 + a1]); \
        INTERLEAVE_COL(y8, y7, y0, y1); \

            // Process the even rows.
            LOAD_COL(0)
            xor (g6, g6);                  // Loop counter.
            L("loop_pass1a");
            {
                call(pass1_process);          // Process all the rows.
                vmovdqu(ptr[g5 + g6], y0);           // Store the temporary sums.
                vmovdqu(ptr[g5 + g6 + 32], y1);
                add(g6, 64);
            }
            cmp(g6, 8 * 64);
            jne("loop_pass1a");

            // Process the odd rows.
            LOAD_COL(32);
            xor (g6, g6);                  // Loop counter.
            mov(g7, 480);                 // Column location 15-n.
            vpbroadcastd(y6, ptr[rip + pat_dw_64]);         // Bias.
            L("loop_pass1b");
            {
                call(pass1_process);          // Process all the rows.

                vmovdqu(y4, ptr[g5 + g6 * 2]);         // Load the temporary sums.
                vmovdqu(y5, ptr[g5 + g6 * 2 + 32]);
                vpaddd(y4, y4, y6);              // Add the bias.
                vpaddd(y5, y5, y6);

                vpsubd(y2, y4, y0);              // Subtract the odd terms.
                vpsubd(y3, y5, y1);
                vpaddd(y0, y0, y4);              // Add the odd terms.
                vpaddd(y1, y1, y5);

                vpsrad(y0, y0, 7);               // Shift.
                vpsrad(y1, y1, 7);
                vpsrad(y2, y2, 7);
                vpsrad(y3, y3, 7);

                vpackssdw(y0, y0, y1);              // Pack.
                vpackssdw(y1, y2, y3);

                vmovdqu(ptr[g5 + g6 + 512], y0);     // Column "n".
                vmovdqu(ptr[g5 + g7 + 512], y1);     // Column "15 - n".

                add(g6, 32);
                sub(g7, 32);
            }
            cmp(g6, 8 * 32);
            jne("loop_pass1b");

            // DCT pass 2.
            lea(g6, ptr[rip + pat_idct16_pass2]);  // Load the multiplication factors.
            vpmovsxbw(y15, ptr[g6 + 0 * 16]);
            vpmovsxbw(y14, ptr[g6 + 1 * 16]);
            vpmovsxbw(y13, ptr[g6 + 2 * 16]);
            vpmovsxbw(y12, ptr[g6 + 3 * 16]);
            vpmovsxbw(y11, ptr[g6 + 4 * 16]);
            vpmovsxbw(y10, ptr[g6 + 5 * 16]);
            vpmovsxbw(y9, ptr[g6 + 6 * 16]);
            vpmovsxbw(y8, ptr[g6 + 7 * 16]);
            vbroadcasti128(y7, ptr[rip + pat_idct16_shuf1]);  // Shuffle.
            vpbroadcastq(y6, ptr[rip + pat_idct16_sign]);   // Sign.
            vpbroadcastd(y5, ptr[rip + pat_dw_2048]);       // Bias.

            xor (g6, g6);                  // Loop counter.
            L("loop_pass2");
            {
                vmovdqu(y0, ptr[g5 + g6 + 512]);     // Load the column.
                vpshufb(y0, y0, y7);              // Put the even and odd terms together.

                // Multiply and add.
                // a1: result, a2: tmp, a3: DCT factors 0, a4: DCT factors 1.
#define MULTIPLY_ADD(a1, a2, a3, a4) \
        vpmaddwd(a2, a3, y0); \
        vpmaddwd(a1, a4, y0); \
        vphaddd(a1, a2, a1); \

                // Process the column.
                // a1: result, a2-3: tmp, a4-7: DCT factors.
#define PROCESS_COLUMN(a1, a2, a3, a4, a5, a6, a7) \
        MULTIPLY_ADD(a1, a3, a4, a5); \
        MULTIPLY_ADD(a2, a3, a6, a7); \
        vperm2i128(a3, a1, a2, 0x31); /* Align the lanes and add up the terms. */ \
        vperm2i128(a2, a1, a2, 0x20); \
        vpaddd(a1, a3, a2); \

                // Finish adding up the terms.
        // a1: result, a2: tmp.
#define FINISH(a1, a2) \
        vpsignd(a2, a1, y6);              /* Change the sign of the odd terms. */ \
        vphaddd(a1, a1, a2);             /* Finish summing up the terms. */ \
        vpaddd(a1, a1, y5);              /* Add the bias and shift. */ \
        vpsrad(a1, a1, 12);

        PROCESS_COLUMN(y1, y3, y4, y15, y14, y13, y12);
        PROCESS_COLUMN(y2, y3, y4, y11, y10, y9, y8);
        FINISH(y1, y3);
        FINISH(y2, y3);

#undef PROCESS_COLUMN

        vpackssdw(y0, y1, y2);              // Pack and store without reordering.
        vmovdqu(ptr[g5 + g6], y0);

        add(g6, 32);
            }
            cmp(g6, 512);
            jne("loop_pass2");

            // Reconstruct.
            vpmovzxbd(y2, ptr[rip + pat_idct16_shuf2]);  // Load the shuffle patterns.
            vmovdqu(y3, ptr[rip + pat_idct16_shuf3]);
            xor (g6, g6);                  // Loop counter.
            L("loop_rec");
            {
                vpermd(y0, y2, ptr[g5 + g6]);      // Load and reorder the residual.
                vpshufb(y0, y0, y3);
                vpmovzxbw(y1, ptr[g2]);                // Load the prediction.
                add(g2, g3);
                vpaddw(y0, y0, y1);              // Add the residual and pack.
                vpackuswb(y0, y0, y0);
                vpermq(y0, y0, 0xd8);            // Combine and store.
                vmovdqu(ptr[g0], x0);
                add(g0, g1);
                add(g6, 32);
            }
            cmp(g6, 512);
            jne("loop_rec");

#endif // 0
        }

        void assemble8x8Avx2()
        {
            // from f265

            // IDCT 8x8.
            //
            // Input parameters:
            // - g0:     destination.
            // - g1:     destination stride.
            // - g2:     prediction.
            // - g3:     prediction stride.
            // - g4:     coefficients.
            // - g5:     spill buffer.
            //DEFFUN f265_lbd_idct_8_avx2, ia=6, at=848488, fa=0, ti=1, tv=13, ym=1

            auto &g0 = arg64(0);
            auto &g1 = arg64(1);
            auto &g2 = arg64(2);
            auto &g3 = arg64(3);
            auto &g4 = arg64(4);
            auto &g5 = reg64(5);
            auto &g6 = reg64(6);

            this->stackSize = 128;
            mov(g5, rsp);

            auto &y0 = ymm0; auto &x0 = regXmm(0);
            auto &y1 = ymm1; auto &x1 = regXmm(1);
            auto &y2 = ymm2; auto &x2 = regXmm(2);
            auto &y3 = ymm3; auto &x3 = regXmm(3);
            auto &y4 = ymm4; auto &x4 = regXmm(4);
            auto &y5 = ymm5; auto &x5 = regXmm(5);
            auto &y6 = ymm6; auto &x6 = regXmm(6);
            auto &y7 = ymm7; auto &x7 = regXmm(7);
            auto &y8 = ymm8; auto &x8 = regXmm(8);
            auto &y9 = ymm9; auto &x9 = regXmm(9);
            auto &y10 = ymm10; auto &x10 = regXmm(10);
            auto &y11 = ymm11; auto &x11 = regXmm(11);
            auto &y12 = ymm12; auto &x12 = regXmm(12);

            // I. Load the 16-bit input data.

            // Load and interleave adjacent rows.
            //   Load: [ H1 G1 F1 E1 D1 C1 B1 A1 | H0 G0 F0 E0 D0 C0 B0 A0 ].
            //   Perm: [ H1 G1 F1 E1 H0 G0 F0 E0 | D1 C1 B1 A1 D0 C0 B0 A0 ].
            vpermq(y4, ptr[g4], 0xd8);          // rows 0 & 1.
            vpermq(y5, ptr[g4 + 64], 0xd8); // rows 4 & 5.
            vpermq(y6, ptr[g4 + 32], 0xd8); // rows 2 & 3.
            vpermq(y7, ptr[g4 + 96], 0xd8); // rows 6 & 7.

            // Prepare the linear combination factors.
            //   Perm: [ H4 H0 G4 G0 F0 F4 E4 E0 | D4 D0 C4 C0 B4 B0 A4 A0 ].
            vpunpcklwd(y0, y4, y5); // rows 0 & 4.
            vpunpcklwd(y2, y6, y7); // rows 2 & 6.
            vpunpckhwd(y1, y4, y5); // rows 1 & 5.
            vpunpckhwd(y3, y6, y7); // rows 3 & 7.

            lea(g6, ptr[rip + pat_idct8_pass1]); // Table address.
            xor(g4, g4); // Loop counter.
            vpbroadcastd(y8, ptr[rip + pat_dw_64]); // Load the bias (64).

            // II. Perform the 1st pass of IDCT8x8.
            L("loop_idct8_pass1");
            {

                // a1: output register, a2: 1st input, a3: 2nd input, a4: IDCT factor index.
#define PROCESS_COLUMN(a1, a2, a3, a4) \
        vpbroadcastd(y7, ptr[g6 + a4]); \
        vpbroadcastd(y6, ptr[g6 + a4 + 4]); \
        vpmaddwd(y7, y7, a2); \
        vpmaddwd(a1, y6, a3); \
        vpaddd(a1, y7, a1); \

                PROCESS_COLUMN(y4, y0, y2, 0); // Combine rows 0, 2, 4, 6.
                PROCESS_COLUMN(y5, y1, y3, 32);  // Combine rows 1, 3, 5, 7.

#undef PROCESS_COLUMN

                vpaddd(y6, y4, y5); // Do two linear combinations per itereration: 0-7, 1-6, 2-5, 3-4.
                vpsubd(y7, y4, y5);

                vpaddd(y6, y8, y6); // Add constant factor and shift.
                vpaddd(y7, y8, y7);
                vpsrad(y6, y6, 7);
                vpsrad(y7, y7, 7);

                vpackssdw(y5, y6, y7); // Pack.
                vpermq(y5, y5, 0xd8); // Put back the column values together.
                vmovdqu(ptr[g5 + g4], y5); // Store.

                add(g4, 32);
                add(g6, 8);
            }
            cmp(g4, 128);
            jne("loop_idct8_pass1");

            // III. Rearrange the rows of data.
            vmovdqu(y0, ptr[g5]); // Columns 0 & 7.
            vmovdqu(y1, ptr[g5 + 32]); // Columns 1 & 6.
            vmovdqu(y2, ptr[g5 + 64]); // Columns 2 & 5.
            // Columns 3 & 4 are already in y5.

            vinserti128(y3, y0, x1, 1); // Columns 0 & 1.
            vinserti128(y4, y2, x5, 1); // Columns 2 & 3.
            vperm2i128(y5, y5, y2, 0x31); // Columns 4 & 5.
            vperm2i128(y0, y0, y1, 0x13); // Columns 6 & 7.

            vmovdqu(ptr[g5], y3); // Store in order.
            vmovdqu(ptr[g5 + 32], y4);
            vmovdqu(ptr[g5 + 64], y5);
            vmovdqu(ptr[g5 + 96], y0);

            // IV. Perform the 2nd pass of IDCT8x8.
            // FIXME: might want to do like DCT 32x32 here.
            lea(g4, ptr[rip + pat_idct8_pass2]);   // Load the multiplication factors.
            vbroadcasti128(y12, ptr[g4]);
            vbroadcasti128(y11, ptr[g4 + 16]);
            vbroadcasti128(y10, ptr[g4 + 32]);
            vbroadcasti128(y9, ptr[g4 + 48]);
            vbroadcasti128(y8, ptr[g4 + 64]);
            vbroadcasti128(y7, ptr[g4 + 80]);
            vbroadcasti128(y6, ptr[g4 + 96]);
            vbroadcasti128(y5, ptr[g4 + 112]);

            vpbroadcastd(y4, ptr[rip + pat_dw_2048]); // Load the bias (2048).

            xor(g4, g4);
            L("loop_idct8_pass2");
            {
                vmovdqu(y0, ptr[g5 + g4]);// Load the next column

                // a1: output register, a2: IDCT factors 0, a3: IDCT factors 1.
#define PROCESS_COLUMN(a1, a2, a3) \
        vpmaddwd(y2, a3, y0); \
        vpmaddwd(a1, a2, y0); /* Multiply by the factors and add the terms once. */ \
        vphaddd(a1, a1, y2); \

                PROCESS_COLUMN(y1, y12, y11); // Do the 8 linear combinations & add the terms.
        PROCESS_COLUMN(y3, y10, y9);
        vphaddd(y1, y1, y3);

        PROCESS_COLUMN(y3, y8, y7);
        PROCESS_COLUMN(y0, y6, y5);
        vphaddd(y0, y3, y0);

#undef PROCESS_COLUMN

        vpaddd(y1, y4, y1); // Add the bias and shift.
        vpaddd(y0, y4, y0);
        vpsrad(y1, y1, 12);
        vpsrad(y0, y0, 12);

        vpackssdw(y1, y1, y0); // Pack.

        // Perform the reconstruction.
        vmovq(x0, ptr[g2]); // Load the prediction.
        vmovhps(x0, ptr[g2 + g3]);
        lea(g2, ptr[g2 + g3 * 2]);
        vpmovzxbw(y2, x0); // Convert the prediction to 16-bit.
        vpaddw(y1, y1, y2); // Add the residual.
        vpackuswb(y0, y1, y1); // Convert the reconstruction to 8-bit.
        vpermq(y0, y0, 0xd8); // Merge the data in the low lane.
        vmovq(ptr[g0], x0); // Store.
        vmovhps(ptr[g0 + g1], x0);
        lea(g0, ptr[g0 + g1 * 2]);
            }
            add(g4, 32);
            cmp(g4, 128);
            jne("loop_idct8_pass2");
        }

        void assemble4x4Avx2()
        {
            auto &g0 = arg64(0);
            auto &g1 = arg64(1);
            auto &g2 = arg64(2);
            auto &g3 = arg64(3);
            auto &g4 = arg64(4);
            auto &g5 = reg64(5);

            auto &y0 = ymm0; auto &x0 = regXmm(0);
            auto &y1 = ymm1; auto &x1 = regXmm(1);
            auto &y2 = ymm2; auto &x2 = regXmm(2);
            auto &y3 = ymm3; auto &x3 = regXmm(3);
            auto &y4 = ymm4; auto &x4 = regXmm(4);
            auto &y5 = ymm5; auto &x5 = regXmm(5);
            auto &y6 = ymm6; auto &x6 = regXmm(6);
            auto &y7 = ymm7; auto &x7 = regXmm(7);

            // from f265

            //; Input parameters:
            //; - g0:     destination.
            //; - g1:     destination stride.
            //; - g2:     prediction.
            //; - g3:     prediction stride.
            //; - g4:     coefficients.
            //; - g5:     spill buffer.

            if (this->trType)
            {
                // DEFFUN f265_lbd_idct_dst_avx2, ia=6, at=848488, fa=0, ti=0, tv=9, ym=1

                auto &y8 = ymm8; auto &x8 = regXmm(8);

                // Do IDST Pass 1.
                vmovdqu(y0, ptr[g4]);
                vbroadcasti128(y3, ptr[rip + pat_idct4_shuf1]);   // Load shuffle pattern.
                lea(g5, ptr[rip + pat_idst_pass1]);
                vmovdqu(y8, ptr[rip + pat_idst_shuf]);
                vpbroadcastd(y6, ptr[rip + pat_dw_64]);
                vpshufb(y0, y0, y3);
                LOAD_64BIT(y1, y2, g5);
                vpermd(y0, y8, y0);
                LOAD_64BIT(y5, y4, g5 + 16);
                MULT_4({ y1, y2, y5, y4 }, y0);
                lea(g5, ptr[rip + pat_idst_pass2]);
                vmovdqu(y0, ptr[g5]);
                vmovdqu(y7, ptr[g5 + 32]);
                vphaddd(y1, y1, y2);
                vphaddd(y5, y5, y4);
                vmovdqu(y2, ptr[g5 + 64]);
                vmovdqu(y4, ptr[g5 + 96]);
                vpaddd(y1, y1, y6);
                vpaddd(y5, y5, y6);
                vpsrad(y1, y1, 7);
                vpsrad(y5, y5, 7);
                vpackssdw(y1, y1, y5);              // Output: 15.11.14.10  13.9.12.8 | 7.3.6.2  5.1.4.0.
                lea(g4, ptr[g2 + g3 * 2]);

                // Do IDST Pass 2.
                MULT_4({ y0, y7, y2, y4 }, y1);
                vpbroadcastd(y6, ptr[rip + pat_dw_2048]);
                vmovd(x5, ptr[g2]);                // Load line 0 of prediction.

#define PERM_4(a1, a2, a3, a4, a5, a6) \
        /* %1: out 0, %2: out 1, %3: in 0/out 2, %4: in 1, %5: in 2/out 3, %6: in 3. */ \
        vperm2i128(a1, a3, a4, 0x20); \
        vperm2i128(a3, a3, a4, 0x31); \
        vperm2i128(a2, a5, a6, 0x20); \
        vperm2i128(a5, a5, a6, 0x31);

                PERM_4(y1, y7, y0, y7, y2, y4);
                vpaddd(y1, y1, y6);
                vpaddd(y7, y7, y6);
                vpunpckldq(x5, x5, ptr[g2 + g3]);       // Load line 1 of prediction.
                vmovd(x6, ptr[g4]);                // Load line 2 of prediction.
                vpaddd(y0, y1, y0);
                vpaddd(y2, y2, y7);
                vpunpckldq(x6, x6, ptr[g4 + g3]);       // Load line 3 of prediction.
                vpunpcklqdq(y1, y5, y6);
                vpsrad(y0, y0, 12);
                vpsrad(y2, y2, 12);
                vpmovzxbw(y1, x1);
                vpackssdw(y0, y0, y2);              // Pack output to 16-bit.
                vpshufb(y0, y0, y3);
                vpermd(y0, y8, y0);
            }
            else
            {
                //DEFFUN f265_lbd_idct_4_avx2, ia=6, at=848488, fa=0, ti=0, tv=8, ym=1

                // Do IDCT4 Pass 1.
                vmovdqu(y0, ptr[g4]);
                vmovdqu(y4, ptr[rip + pat_idct4_pass1]); // Load multiplication factors.
                vpermq(y1, y0, 0x00); // 3.2.1.0  3.2.1.0 | 3.2.1.0  3.2.1.0.
                vpermq(y2, y0, 0xAA); // 11.10.9.8  11.10.9.8 | 11.10.9.8  11.10.9.8.
                vpermq(y0, y0, 0xDD); // ; 15.14.13.12  7.6.5.4 | 15.14.13.12  7.6.5.4.
                vbroadcasti128(y3, ptr[rip + pat_idct4_shuf1]); // Load shuffle pattern.
                vpbroadcastd(y6, ptr[rip + pat_dw_64]); // Bias.
                vpmovsxwd(y1, x1);
                vpmovsxwd(y2, x2);
                vpsignd(y2, y2, ptr[rip + pat_idct4_sign]);
                vpshufb(y0, y0, y3);
                vpmaddwd(y0, y0, y4); // Obtain sub_1 | add_1.
                vpaddd(y1, y1, y2);
                vpslld(y1, y1, 6); // Obtain sub_0 | add_0.
                vpaddd(y1, y1, y6); // Add bias.
                vpaddd(y4, y1, y0); // sub_0 + sub_1 | add_0 + add_1.
                vpsubd(y5, y1, y0); // sub_0 - sub_1 | add_0 - add_1.
                vpsrad(y4, y4, 7); // Shift.
                vpsrad(y5, y5, 7);
                lea(g5, ptr[rip + pat_idct4_pass2]);
                vpackssdw(y0, y4, y5); // Convert to 16 - bit and pack.
                lea(g4, ptr[g2 + g3 * 2]);

                // Do IDCT4 Pass 2.
                LOAD_64BIT(y2, y3, g5);
                vpermq(y0, y0, 0x78); // Shuffle input : 15.11.7.3  14.10.6.2 | 13.9.5.1  12.8.4.0.
                LOAD_64BIT(y1, y4, g5 + 16);
                MULT_4({ y2, y3, y1, y4 }, y0);
                vpbroadcastd(y6, ptr[rip + pat_dw_2048]); // Bias.
                vmovd(x5, ptr[g2]); // Load line 0 of prediction.
                vphaddd(y3, y2, y3);
                vbroadcasti128(y7, ptr[rip + pat_idct16_shuf1]);
                vphaddd(y1, y1, y4);
                vpunpckldq(x5, x5, ptr[g2 + g3]); // Load line 1 of prediction.
                vmovd(x0, ptr[g4]); // Load line 2 of prediction.
                vpaddd(y3, y3, y6);
                vpaddd(y1, y1, y6);
                vpunpckldq(x0, x0, ptr[g4 + g3]); // Load line 3 of prediction.
                vpsrad(y1, y1, 12);
                vpunpcklqdq(y0, y5, y0);
                vpsrad(y3, y3, 12);
                vpmovzxbw(y0, x0);
                vpackssdw(y4, y3, y1);
                vpshufb(y1, y4, y7);
                //call recon_4x4               ; Reconstruction.
            }

            // IDCT4/IDST helper function - Reconstruction.
            // recon_4x4

            vpaddw(y0, y1, y0);
            vpackuswb(y1, y0, y0);
            vextracti128(x0, y1, 1);
            vpsrldq(y3, y1, 4);
            vmovd(ptr[g0], x1);
            vmovd(ptr[g0 + g1], x3);
            lea(g0, ptr[g0 + g1 * 2]);
            vpsrldq(y4, y0, 4);
            vmovd(ptr[g0], x0);
            vmovd(ptr[g0 + g1], x4);
        }

        template <class Address>
        void LOAD_64BIT(Xbyak::Ymm const &y1, Xbyak::Ymm const &y2, Address &&address)
        {
            vpbroadcastq(y1, ptr[address]);
            vpbroadcastq(y2, ptr[address + 8]);
        }

        void MULT_4(std::initializer_list<Xbyak::Ymm> regs, Xbyak::Ymm const &multiplier)
        {
            for (auto reg : regs) vpmaddwd(reg, multiplier, reg);
        }

#endif

        int trType;
        int log2TrafoSize;
    };



    static havoc_inverse_transform* get_inverse_transform(int trType, int log2TrafoSize, havoc_code code, int encoder)
    {
        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
        auto mask = buffer.isa;

        const int nCbS = 1 << log2TrafoSize;

        havoc_inverse_transform *f = 0;

        if (mask & (HAVOC_C_REF | HAVOC_C_OPT))
        {
            if (nCbS == 4) f = trType ? havoc_idst_c_opt<4> : havoc_idct_c_opt<4>;
            if (nCbS == 8) f = havoc_idct_c_opt<8>;
            if (nCbS == 16) f = havoc_idct_c_opt<16>;
            if (nCbS == 32) f = havoc_idct_c_opt<32>;
        }

        return f;
    }


    void havoc_populate_inverse_transform(havoc_table_inverse_transform* table, havoc_code code, int encoder)
    {
        *havoc_get_inverse_transform(table, 1, 2) = get_inverse_transform(1, 2, code, encoder);
        for (int log2TrafoSize = 2; log2TrafoSize <= 5; ++log2TrafoSize)
        {
            *havoc_get_inverse_transform(table, 0, log2TrafoSize) = get_inverse_transform(0, log2TrafoSize, code, encoder);
        }
    }


    template <typename Sample>
    havoc_inverse_transform_add<Sample>* get_inverse_transform_add(int trType, int log2TrafoSize, havoc_code code, int encoder);


    template<>
    havoc_inverse_transform_add<uint8_t>* get_inverse_transform_add<uint8_t>(int trType, int log2TrafoSize, havoc_code code, int encoder)
    {
        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
        auto mask = buffer.isa;
        const int nCbS = 1 << log2TrafoSize;

        havoc_inverse_transform_add<uint8_t> *f = 0;

        if (mask & (HAVOC_C_REF | HAVOC_C_OPT))
        {
            if (nCbS == 4) f = trType ? havoc_idst_4x4_c_opt : havoc_idct_4x4_c_opt;
            if (nCbS == 8) f = havoc_idct_add_c_opt<8, uint8_t>;
            if (nCbS == 16) f = havoc_idct_add_c_opt<16, uint8_t>;
            if (nCbS == 32) f = havoc_idct_add_c_opt<32, uint8_t>;
        }

        if (mask & HAVOC_SSSE3)
        {
            if (nCbS == 8)
            {
                InverseTransformAdd a(&buffer, trType, log2TrafoSize);
                f = a;
            }
            if (nCbS == 16)
            {
                InverseTransformAdd a(&buffer, trType, log2TrafoSize);
                f = a;
            }
        }

#if USE_F265_DERIVED
        if (mask & HAVOC_AVX2)
        {
            if (nCbS == 4 && !trType)
            {
                InverseTransformAdd a(&buffer, trType, log2TrafoSize);
                f = a;
            }
            if (nCbS == 4 && trType)
            {
                InverseTransformAdd a(&buffer, trType, log2TrafoSize);
                f = a;
            }
            if (nCbS == 8)
            {
                InverseTransformAdd a(&buffer, trType, log2TrafoSize);
                f = a;
            }
            if (nCbS == 16)
            {
                InverseTransformAdd a(&buffer, trType, log2TrafoSize);
                f = a;
            }
            if (nCbS == 32 && encoder)
            {
                // The 32x32 idct from f265 is non-conforming. Suspect arithmetic overflow.
                // This probem becomes apparent when decoding the DELTAQP_A conformance stream which has artificial, extreme coefficient values.
                InverseTransformAdd a(&buffer, trType, log2TrafoSize);
                f = a;
            }
        }
#endif

        return f;
    }


    template<>
    havoc_inverse_transform_add<uint16_t>* get_inverse_transform_add(int trType, int log2TrafoSize, havoc_code code, int encoder)
    {
        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
        auto mask = buffer.isa;

        const int nCbS = 1 << log2TrafoSize;

        havoc_inverse_transform_add<uint16_t> *f = 0;

        if (mask & (HAVOC_C_REF | HAVOC_C_OPT))
        {
            if (nCbS == 4) f = trType ? havoc_idst_4x4_16_c_opt : havoc_idct_4x4_16_c_opt;
            if (nCbS == 8) f = havoc_idct_add_c_opt<8, uint16_t>;
            if (nCbS == 16) f = havoc_idct_add_c_opt<16, uint16_t>;
            if (nCbS == 32) f = havoc_idct_add_c_opt<32, uint16_t>;
        }

        return f;
    }


    template <typename Sample>
    void havoc_populate_inverse_transform_add(havoc_table_inverse_transform_add<Sample> *table, havoc_code code, int encoder)
    {
        *havoc_get_inverse_transform_add(table, 1, 2) = get_inverse_transform_add<Sample>(1, 2, code, encoder);
        for (int log2TrafoSize = 2; log2TrafoSize <= 5; ++log2TrafoSize)
        {
            *havoc_get_inverse_transform_add(table, 0, log2TrafoSize) = get_inverse_transform_add<Sample>(0, log2TrafoSize, code, encoder);
        }
    }


    template void havoc_populate_inverse_transform_add<uint8_t>(havoc_table_inverse_transform_add<uint8_t> *table, havoc_code code, int encoder);
    template void havoc_populate_inverse_transform_add<uint16_t>(havoc_table_inverse_transform_add<uint16_t> *table, havoc_code code, int encoder);


    template <typename Sample>
    struct bind_inverse_transform_add
    {
        const int16_t *coefficients;
        const Sample *predicted;
        havoc_inverse_transform_add<Sample> *f;
        int log2TrafoSize;
        int trType;
        Sample dst[32 * 32];
    } ;


    template <typename Sample>
    int init_inverse_transform_add(void *p, havoc_code code)
    {
        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

        bind_inverse_transform_add<Sample> *s = (bind_inverse_transform_add<Sample> *)p;

        havoc_table_inverse_transform_add<Sample> table;

        havoc_populate_inverse_transform_add<Sample>(&table, code, 1);

        s->f = *havoc_get_inverse_transform_add<Sample>(&table, s->trType, s->log2TrafoSize);

        if (s->f) for (int i = 0; i < 32 * 32; ++i) s->dst[i] = 0xaa;

        if (s->f && buffer.isa == HAVOC_C_REF)
        {
            const int nCbS = 1 << s->log2TrafoSize;
            printf("\t%s %dx%d : ", s->trType ? "sine" : "cosine", nCbS, nCbS);
        }

        return !!s->f;
    }


    template <typename Sample>
    void invoke_inverse_transform_add(void *p, int n)
    {
        auto *s = (bind_inverse_transform_add<Sample> *)p;

        while (n--)
        {
            s->f(s->dst, (intptr_t)1 << s->log2TrafoSize, s->predicted, (intptr_t)1 << s->log2TrafoSize, s->coefficients, 8);
        }
    }


    template <typename Sample>
    int mismatch_transform_add(void *boundRef, void *boundTest)
    {
        bind_inverse_transform_add<Sample> *ref = (bind_inverse_transform_add<Sample> *)boundRef;
        bind_inverse_transform_add<Sample> *test = (bind_inverse_transform_add<Sample> *)boundTest;

        const int nCbS = 1 << ref->log2TrafoSize;

        return memcmp(ref->dst, test->dst, nCbS * nCbS * sizeof(Sample));
    }


    template <typename Sample>
    void havoc_test_inverse_transform_add(int *error_count, havoc_instruction_set mask)
    {
        printf("\ninverse_transform_add<uint%d_t> - Inverse Transform, then add to predicted\n", int(8 * sizeof(Sample)));

        HAVOC_ALIGN(32, int16_t, coefficients[32 * 32]);
        HAVOC_ALIGN(32, Sample, predicted[32 * 32]);

        for (int x = 0; x < 32 * 32; x++) coefficients[x] = (((rand() << 1) ^ rand()) & 0xff) - 128;
        for (int x = 0; x < 32 * 32; x++) predicted[x] = rand() & 0xff;

        bind_inverse_transform_add<Sample> b[2];
        b[0].coefficients = coefficients;
        b[0].predicted = predicted;

        for (int j = 2; j < 6; ++j)
        {
            b[0].trType = (j == 1) ? 1 : 0;
            b[0].log2TrafoSize = (j == 1) ? 2 : j;
            b[1] = b[0];

            *error_count += havoc_test(&b[0], &b[1], init_inverse_transform_add<Sample>, invoke_inverse_transform_add<Sample>, mismatch_transform_add<Sample>, mask, 100000);
        }
    }


    template void havoc_test_inverse_transform_add<uint8_t>(int *error_count, havoc_instruction_set mask);
    template void havoc_test_inverse_transform_add<uint16_t>(int *error_count, havoc_instruction_set mask);


    template <class Dst, class Src>
    Dst shiftRight(Src src, int shift)
    {
#if 1
        short temp = (src >> shift);
        src = temp;
#else
        src >>= shift;
#endif
        if (src > std::numeric_limits<Dst>::max()) return std::numeric_limits<Dst>::max();
        if (src < std::numeric_limits<Dst>::min()) return std::numeric_limits<Dst>::min();
        Dst dst = static_cast<Dst>(src);
        return dst;
    }


    template <typename Dst, typename Src>
    void havoc_partial_butterfly_4x4_dst_c_opt(Dst *dst, const Src *src, intptr_t src_stride, int shift)
    {
        const int add = 1 << (shift - 1);
        const int dst_stride = 4;

        for (int i = 0; i<4; i++)
        {
            int c[4];
            c[0] = src[src_stride*i + 0] + src[src_stride*i + 3];
            c[1] = src[src_stride*i + 1] + src[src_stride*i + 3];
            c[2] = src[src_stride*i + 0] - src[src_stride*i + 1];
            c[3] = 74 * src[src_stride*i + 2];

            dst[0 * dst_stride + i] = shiftRight<Dst>(29 * c[0] + 55 * c[1] + c[3] + add, shift);
            dst[1 * dst_stride + i] = shiftRight<Dst>(74 * (src[src_stride*i + 0] + src[src_stride*i + 1] - src[src_stride*i + 3]) + add, shift);
            dst[2 * dst_stride + i] = shiftRight<Dst>(29 * c[2] + 55 * c[0] - c[3] + add, shift);
            dst[3 * dst_stride + i] = shiftRight<Dst>(55 * c[2] - 29 * c[1] + c[3] + add, shift);
        }
    }


    template <typename Dst, typename Src>
    void havoc_partial_butterfly_4x4_c_opt(Dst *dst, const Src *src, intptr_t src_stride, int shift)
    {
        const int add = 1 << (shift - 1);
        const int dst_stride = 4;

        for (int j = 0; j<4; j++)
        {
            int E[2] = { src[0] + src[3], src[1] + src[2] };
            int O[2] = { src[0] - src[3], src[1] - src[2] };

            static const int16_t table[4][4] =
            {
                    { 64, 64, 64, 64 },
                    { 83, 36, -36, -83 },
                    { 64, -64, -64, 64 },
                    { 36, -83, 83, -36 }
            };

            dst[0 * dst_stride] = shiftRight<Dst>(table[0][0] * E[0] + table[0][1] * E[1] + add, shift);
            dst[2 * dst_stride] = shiftRight<Dst>(table[2][0] * E[0] + table[2][1] * E[1] + add, shift);
            dst[1 * dst_stride] = shiftRight<Dst>(table[1][0] * O[0] + table[1][1] * O[1] + add, shift);
            dst[3 * dst_stride] = shiftRight<Dst>(table[3][0] * O[0] + table[3][1] * O[1] + add, shift);

            src += src_stride;
            dst++;
        }
    }


    template <typename Dst, typename Src>
    void havoc_partial_butterfly_8x8_c_opt(Dst *dst, const Src *src, intptr_t src_stride, int shift)
    {
        const int add = 1 << (shift - 1);
        const int dst_stride = 8;

        for (int j = 0; j<8; j++)
        {
            int E[4], O[4];
            for (int k = 0; k<4; k++)
            {
                E[k] = src[k] + src[7 - k];
                O[k] = src[k] - src[7 - k];
            }

            int EE[2], EO[2];
            EE[0] = E[0] + E[3];
            EO[0] = E[0] - E[3];
            EE[1] = E[1] + E[2];
            EO[1] = E[1] - E[2];

            static const int16_t table[8][8] =
            {
                    { 64, 64, 64, 64, 64, 64, 64, 64 },
                    { 89, 75, 50, 18, -18, -50, -75, -89 },
                    { 83, 36, -36, -83, -83, -36, 36, 83 },
                    { 75, -18, -89, -50, 50, 89, 18, -75 },
                    { 64, -64, -64, 64, 64, -64, -64, 64 },
                    { 50, -89, 18, 75, -75, -18, 89, -50 },
                    { 36, -83, 83, -36, -36, 83, -83, 36 },
                    { 18, -50, 75, -89, 89, -75, 50, -18 }
            };

            dst[0 * dst_stride] = shiftRight<Dst>(table[0][0] * EE[0] + table[0][1] * EE[1] + add, shift);
            dst[4 * dst_stride] = shiftRight<Dst>(table[4][0] * EE[0] + table[4][1] * EE[1] + add, shift);
            dst[2 * dst_stride] = shiftRight<Dst>(table[2][0] * EO[0] + table[2][1] * EO[1] + add, shift);
            dst[6 * dst_stride] = shiftRight<Dst>(table[6][0] * EO[0] + table[6][1] * EO[1] + add, shift);

            dst[1 * dst_stride] = shiftRight<Dst>(table[1][0] * O[0] + table[1][1] * O[1] + table[1][2] * O[2] + table[1][3] * O[3] + add, shift);
            dst[3 * dst_stride] = shiftRight<Dst>(table[3][0] * O[0] + table[3][1] * O[1] + table[3][2] * O[2] + table[3][3] * O[3] + add, shift);
            dst[5 * dst_stride] = shiftRight<Dst>(table[5][0] * O[0] + table[5][1] * O[1] + table[5][2] * O[2] + table[5][3] * O[3] + add, shift);
            dst[7 * dst_stride] = shiftRight<Dst>(table[7][0] * O[0] + table[7][1] * O[1] + table[7][2] * O[2] + table[7][3] * O[3] + add, shift);

            src += src_stride;
            dst++;
        }
    }


    template <typename Dst, typename Src>
    void havoc_partial_butterfly_16x16_c_opt(Dst *dst, const Src *src, intptr_t src_stride, int shift)
    {
        const int add = 1 << (shift - 1);
        const intptr_t dst_stride = 16;

        for (int j = 0; j<16; ++j)
        {
            int E[8], O[8];
            for (int k = 0; k<8; ++k)
            {
                E[k] = src[k] + src[15 - k];
                O[k] = src[k] - src[15 - k];
            }

            int EE[4], EO[4];
            for (int k = 0; k<4; ++k)
            {
                EE[k] = E[k] + E[7 - k];
                EO[k] = E[k] - E[7 - k];
            }

            int EEE[2], EEO[2];
            EEE[0] = EE[0] + EE[3];
            EEO[0] = EE[0] - EE[3];
            EEE[1] = EE[1] + EE[2];
            EEO[1] = EE[1] - EE[2];

            static const int16_t table[16][16] =
            {
                    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
                    { 90, 87, 80, 70, 57, 43, 25, 9, -9, -25, -43, -57, -70, -80, -87, -90 },
                    { 89, 75, 50, 18, -18, -50, -75, -89, -89, -75, -50, -18, 18, 50, 75, 89 },
                    { 87, 57, 9, -43, -80, -90, -70, -25, 25, 70, 90, 80, 43, -9, -57, -87 },
                    { 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83 },
                    { 80, 9, -70, -87, -25, 57, 90, 43, -43, -90, -57, 25, 87, 70, -9, -80 },
                    { 75, -18, -89, -50, 50, 89, 18, -75, -75, 18, 89, 50, -50, -89, -18, 75 },
                    { 70, -43, -87, 9, 90, 25, -80, -57, 57, 80, -25, -90, -9, 87, 43, -70 },
                    { 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64 },
                    { 57, -80, -25, 90, -9, -87, 43, 70, -70, -43, 87, 9, -90, 25, 80, -57 },
                    { 50, -89, 18, 75, -75, -18, 89, -50, -50, 89, -18, -75, 75, 18, -89, 50 },
                    { 43, -90, 57, 25, -87, 70, 9, -80, 80, -9, -70, 87, -25, -57, 90, -43 },
                    { 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36 },
                    { 25, -70, 90, -80, 43, 9, -57, 87, -87, 57, -9, -43, 80, -90, 70, -25 },
                    { 18, -50, 75, -89, 89, -75, 50, -18, -18, 50, -75, 89, -89, 75, -50, 18 },
                    { 9, -25, 43, -57, 70, -80, 87, -90, 90, -87, 80, -70, 57, -43, 25, -9 }
            };

            dst[0 * dst_stride] = shiftRight<Dst>(table[0][0] * EEE[0] + table[0][1] * EEE[1] + add, shift);
            dst[8 * dst_stride] = shiftRight<Dst>(table[8][0] * EEE[0] + table[8][1] * EEE[1] + add, shift);
            dst[4 * dst_stride] = shiftRight<Dst>(table[4][0] * EEO[0] + table[4][1] * EEO[1] + add, shift);
            dst[12 * dst_stride] = shiftRight<Dst>(table[12][0] * EEO[0] + table[12][1] * EEO[1] + add, shift);

            for (int k = 2; k<16; k += 4)
            {
                dst[k*dst_stride] = shiftRight<Dst>(table[k][0] * EO[0] + table[k][1] * EO[1] + table[k][2] * EO[2] + table[k][3] * EO[3] + add, shift);
            }

            for (int k = 1; k<16; k += 2)
            {
                dst[k*dst_stride] = shiftRight<Dst>(table[k][0] * O[0] + table[k][1] * O[1] + table[k][2] * O[2] + table[k][3] * O[3] +
                                                    table[k][4] * O[4] + table[k][5] * O[5] + table[k][6] * O[6] + table[k][7] * O[7] + add, shift);
            }

            src += src_stride;
            ++dst;
        }
    }


    template <typename Dst, typename Src>
    void havoc_partial_butterfly_32x32_c_opt(Dst *dst, const Src *src, intptr_t src_stride, int shift)
    {
        const int add = 1 << (shift - 1);
        const int dst_stride = 32;

        for (int j = 0; j<32; j++)
        {
            int E[16], O[16];
            for (int k = 0; k<16; k++)
            {
                E[k] = src[k] + src[31 - k];
                O[k] = src[k] - src[31 - k];
            }

            int EE[8], EO[8];
            for (int k = 0; k<8; k++)
            {
                EE[k] = E[k] + E[15 - k];
                EO[k] = E[k] - E[15 - k];
            }

            int EEE[4], EEO[4];
            for (int k = 0; k<4; k++)
            {
                EEE[k] = EE[k] + EE[7 - k];
                EEO[k] = EE[k] - EE[7 - k];
            }

            int EEEE[2], EEEO[2];
            EEEE[0] = EEE[0] + EEE[3];
            EEEO[0] = EEE[0] - EEE[3];
            EEEE[1] = EEE[1] + EEE[2];
            EEEO[1] = EEE[1] - EEE[2];

            static const int16_t table[32][32] =
            {
                    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
                    { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13, 4, -4, -13, -22, -31, -38, -46, -54, -61, -67, -73, -78, -82, -85, -88, -90, -90 },
                    { 90, 87, 80, 70, 57, 43, 25, 9, -9, -25, -43, -57, -70, -80, -87, -90, -90, -87, -80, -70, -57, -43, -25, -9, 9, 25, 43, 57, 70, 80, 87, 90 },
                    { 90, 82, 67, 46, 22, -4, -31, -54, -73, -85, -90, -88, -78, -61, -38, -13, 13, 38, 61, 78, 88, 90, 85, 73, 54, 31, 4, -22, -46, -67, -82, -90 },
                    { 89, 75, 50, 18, -18, -50, -75, -89, -89, -75, -50, -18, 18, 50, 75, 89, 89, 75, 50, 18, -18, -50, -75, -89, -89, -75, -50, -18, 18, 50, 75, 89 },
                    { 88, 67, 31, -13, -54, -82, -90, -78, -46, -4, 38, 73, 90, 85, 61, 22, -22, -61, -85, -90, -73, -38, 4, 46, 78, 90, 82, 54, 13, -31, -67, -88 },
                    { 87, 57, 9, -43, -80, -90, -70, -25, 25, 70, 90, 80, 43, -9, -57, -87, -87, -57, -9, 43, 80, 90, 70, 25, -25, -70, -90, -80, -43, 9, 57, 87 },
                    { 85, 46, -13, -67, -90, -73, -22, 38, 82, 88, 54, -4, -61, -90, -78, -31, 31, 78, 90, 61, 4, -54, -88, -82, -38, 22, 73, 90, 67, 13, -46, -85 },
                    { 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83, 83, 36, -36, -83, -83, -36, 36, 83 },
                    { 82, 22, -54, -90, -61, 13, 78, 85, 31, -46, -90, -67, 4, 73, 88, 38, -38, -88, -73, -4, 67, 90, 46, -31, -85, -78, -13, 61, 90, 54, -22, -82 },
                    { 80, 9, -70, -87, -25, 57, 90, 43, -43, -90, -57, 25, 87, 70, -9, -80, -80, -9, 70, 87, 25, -57, -90, -43, 43, 90, 57, -25, -87, -70, 9, 80 },
                    { 78, -4, -82, -73, 13, 85, 67, -22, -88, -61, 31, 90, 54, -38, -90, -46, 46, 90, 38, -54, -90, -31, 61, 88, 22, -67, -85, -13, 73, 82, 4, -78 },
                    { 75, -18, -89, -50, 50, 89, 18, -75, -75, 18, 89, 50, -50, -89, -18, 75, 75, -18, -89, -50, 50, 89, 18, -75, -75, 18, 89, 50, -50, -89, -18, 75 },
                    { 73, -31, -90, -22, 78, 67, -38, -90, -13, 82, 61, -46, -88, -4, 85, 54, -54, -85, 4, 88, 46, -61, -82, 13, 90, 38, -67, -78, 22, 90, 31, -73 },
                    { 70, -43, -87, 9, 90, 25, -80, -57, 57, 80, -25, -90, -9, 87, 43, -70, -70, 43, 87, -9, -90, -25, 80, 57, -57, -80, 25, 90, 9, -87, -43, 70 },
                    { 67, -54, -78, 38, 85, -22, -90, 4, 90, 13, -88, -31, 82, 46, -73, -61, 61, 73, -46, -82, 31, 88, -13, -90, -4, 90, 22, -85, -38, 78, 54, -67 },
                    { 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64, 64, -64, -64, 64 },
                    { 61, -73, -46, 82, 31, -88, -13, 90, -4, -90, 22, 85, -38, -78, 54, 67, -67, -54, 78, 38, -85, -22, 90, 4, -90, 13, 88, -31, -82, 46, 73, -61 },
                    { 57, -80, -25, 90, -9, -87, 43, 70, -70, -43, 87, 9, -90, 25, 80, -57, -57, 80, 25, -90, 9, 87, -43, -70, 70, 43, -87, -9, 90, -25, -80, 57 },
                    { 54, -85, -4, 88, -46, -61, 82, 13, -90, 38, 67, -78, -22, 90, -31, -73, 73, 31, -90, 22, 78, -67, -38, 90, -13, -82, 61, 46, -88, 4, 85, -54 },
                    { 50, -89, 18, 75, -75, -18, 89, -50, -50, 89, -18, -75, 75, 18, -89, 50, 50, -89, 18, 75, -75, -18, 89, -50, -50, 89, -18, -75, 75, 18, -89, 50 },
                    { 46, -90, 38, 54, -90, 31, 61, -88, 22, 67, -85, 13, 73, -82, 4, 78, -78, -4, 82, -73, -13, 85, -67, -22, 88, -61, -31, 90, -54, -38, 90, -46 },
                    { 43, -90, 57, 25, -87, 70, 9, -80, 80, -9, -70, 87, -25, -57, 90, -43, -43, 90, -57, -25, 87, -70, -9, 80, -80, 9, 70, -87, 25, 57, -90, 43 },
                    { 38, -88, 73, -4, -67, 90, -46, -31, 85, -78, 13, 61, -90, 54, 22, -82, 82, -22, -54, 90, -61, -13, 78, -85, 31, 46, -90, 67, 4, -73, 88, -38 },
                    { 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36, 36, -83, 83, -36, -36, 83, -83, 36 },
                    { 31, -78, 90, -61, 4, 54, -88, 82, -38, -22, 73, -90, 67, -13, -46, 85, -85, 46, 13, -67, 90, -73, 22, 38, -82, 88, -54, -4, 61, -90, 78, -31 },
                    { 25, -70, 90, -80, 43, 9, -57, 87, -87, 57, -9, -43, 80, -90, 70, -25, -25, 70, -90, 80, -43, -9, 57, -87, 87, -57, 9, 43, -80, 90, -70, 25 },
                    { 22, -61, 85, -90, 73, -38, -4, 46, -78, 90, -82, 54, -13, -31, 67, -88, 88, -67, 31, 13, -54, 82, -90, 78, -46, 4, 38, -73, 90, -85, 61, -22 },
                    { 18, -50, 75, -89, 89, -75, 50, -18, -18, 50, -75, 89, -89, 75, -50, 18, 18, -50, 75, -89, 89, -75, 50, -18, -18, 50, -75, 89, -89, 75, -50, 18 },
                    { 13, -38, 61, -78, 88, -90, 85, -73, 54, -31, 4, 22, -46, 67, -82, 90, -90, 82, -67, 46, -22, -4, 31, -54, 73, -85, 90, -88, 78, -61, 38, -13 },
                    { 9, -25, 43, -57, 70, -80, 87, -90, 90, -87, 80, -70, 57, -43, 25, -9, -9, 25, -43, 57, -70, 80, -87, 90, -90, 87, -80, 70, -57, 43, -25, 9 },
                    { 4, -13, 22, -31, 38, -46, 54, -61, 67, -73, 78, -82, 85, -88, 90, -90, 90, -90, 88, -85, 82, -78, 73, -67, 61, -54, 46, -38, 31, -22, 13, -4 }
            };

            dst[0 * dst_stride] = shiftRight<Dst>(table[0][0] * EEEE[0] + table[0][1] * EEEE[1] + add, shift);
            dst[16 * dst_stride] = shiftRight<Dst>(table[16][0] * EEEE[0] + table[16][1] * EEEE[1] + add, shift);
            dst[8 * dst_stride] = shiftRight<Dst>(table[8][0] * EEEO[0] + table[8][1] * EEEO[1] + add, shift);
            dst[24 * dst_stride] = shiftRight<Dst>(table[24][0] * EEEO[0] + table[24][1] * EEEO[1] + add, shift);
            for (int k = 4; k<32; k += 8)
            {
                dst[k*dst_stride] = shiftRight<Dst>(table[k][0] * EEO[0] + table[k][1] * EEO[1] + table[k][2] * EEO[2] + table[k][3] * EEO[3] + add, shift);
            }
            for (int k = 2; k<32; k += 4)
            {
                dst[k*dst_stride] = shiftRight<Dst>(table[k][0] * EO[0] + table[k][1] * EO[1] + table[k][2] * EO[2] + table[k][3] * EO[3] +
                                                    table[k][4] * EO[4] + table[k][5] * EO[5] + table[k][6] * EO[6] + table[k][7] * EO[7] + add, shift);
            }
            for (int k = 1; k<32; k += 2)
            {
                dst[k*dst_stride] = shiftRight<Dst>(table[k][0] * O[0] + table[k][1] * O[1] + table[k][2] * O[2] + table[k][3] * O[3] +
                                                    table[k][4] * O[4] + table[k][5] * O[5] + table[k][6] * O[6] + table[k][7] * O[7] +
                                                    table[k][8] * O[8] + table[k][9] * O[9] + table[k][10] * O[10] + table[k][11] * O[11] +
                                                    table[k][12] * O[12] + table[k][13] * O[13] + table[k][14] * O[14] + table[k][15] * O[15] + add, shift);
            }
            src += src_stride;
            dst++;
        }
    }


    template <int bitDepth>
    void havoc_dst_4x4_c_opt(int16_t coeffs[4 * 4], const int16_t *src, intptr_t src_stride)
    {
        int16_t temp[4 * 4];
        havoc_partial_butterfly_4x4_dst_c_opt(temp, src, src_stride,1 + bitDepth - 8);
        havoc_partial_butterfly_4x4_dst_c_opt(coeffs, temp, 4, 8);
    }


    template <int bitDepth>
    void havoc_dct_4x4_c_opt(int16_t coeffs[4 * 4], const int16_t *src, intptr_t src_stride)
    {
        int16_t temp[4 * 4];
        havoc_partial_butterfly_4x4_c_opt(temp, src, src_stride, 1 + bitDepth - 8);
        havoc_partial_butterfly_4x4_c_opt(coeffs, temp, 4, 8);
    }


    template <int bitDepth>
    void havoc_dct_8x8_c_opt(int16_t coeffs[8 * 8], const int16_t *src, intptr_t src_stride)
    {
        int16_t temp[8 * 8];
        havoc_partial_butterfly_8x8_c_opt(temp, src, src_stride, 2 + bitDepth - 8);
        havoc_partial_butterfly_8x8_c_opt(coeffs, temp, 8, 9);
    }


    template <int bitDepth>
    void havoc_dct_16x16_c_opt(int16_t coeffs[16 * 16], const int16_t *src, intptr_t src_stride)
    {
        int16_t temp[16 * 16];
        havoc_partial_butterfly_16x16_c_opt(temp, src, src_stride, 3 + bitDepth - 8);
        havoc_partial_butterfly_16x16_c_opt(coeffs, temp, 16, 10);
    }


    template <int bitDepth>
    void havoc_dct_32x32_c_opt(int16_t coeffs[32 * 32], const int16_t *src, intptr_t src_stride)
    {
        int16_t temp[32 * 32];
        havoc_partial_butterfly_32x32_c_opt(temp, src, src_stride, 4 + bitDepth - 8);
        havoc_partial_butterfly_32x32_c_opt(coeffs, temp, 32, 11);
    }


    template <int bitDepth>
    struct ForwardDct16x16
            :
                    Jit::Function
                    {
                        ForwardDct16x16(Jit::Buffer *buffer)
                :
                    Jit::Function(buffer, 3)
                        {
                            this->buildSinglePass(1, 16, 32 + 2 * 2 * 16*16);
                        }

                        Xbyak::Label shuffle_efcdab8967452301;
                        Xbyak::Label dd_times_4_add;
                        Xbyak::Label cosine_8_h;
                        Xbyak::Label cosine_4_h;
                        Xbyak::Label cosine_1_h;
                        Xbyak::Label cosine_8_new;
                        Xbyak::Label cosine_4_new;
                        Xbyak::Label cosine_2_new;
                        Xbyak::Label const_00000200000002000000020000000200;
                        Xbyak::Label const_00000008000000080000000800000008;

                        void data()
                        {
                            align();

                            L(shuffle_efcdab8967452301);
                            db({ 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1 });

                            L(dd_times_4_add);
                            dd({ 4 << (bitDepth - 8) }, 4);

                            L(cosine_8_h);
                            dw({ 90, 87, 80, 70, 57, 43, 25, 9, 87, 57, 9, -43, -80, -90, -70, -25 });
                            dw({ 80, 9, -70, -87, -25, 57, 90, 43, 70, -43, -87, 9, 90, 25, -80, -57 });
                            dw({ 57, -80, -25, 90, -9, -87, 43, 70, 43, -90, 57, 25, -87, 70, 9, -80 });
                            dw({ 25, -70, 90, -80, 43, 9, -57, 87, 9, -25, 43, -57, 70, -80, 87, -90 });

                            L(cosine_4_h);
                            dw({ 89, 75, 50, 18,  50, 89, 18, -75 });
                            dw({ 50, -89, 18, 75,  89, -75, 50, -18 });

                            L(cosine_1_h);
                            dw({ 64, 64, -83, -36, -64, 64, -83, 36 });

                            L(cosine_8_new);
                            dd({ 90, 87, 80, 70, 57, 43, 25, 9, 87, 57, 9, -43, -80, -90, -70, -25 });
                            dd({ 80, 9, -70, -87, -25, 57, 90, 43, 70, -43, -87, 9, 90, 25, -80, -57 });
                            dd({ 57, -80, -25, 90, -9, -87, 43, 70, 43, -90, 57, 25, -87, 70, 9, -80 });
                            dd({ 25, -70, 90, -80, 43, 9, -57, 87, 9, -25, 43, -57, 70, -80, 87, -90 });

                            L(cosine_4_new);
                            dd({ 89, 75, 50, 18, 75, -18, -89, -50 });
                            dd({ 50, -89, 18, 75, 18, -50, 75, -89 });

                            L(cosine_2_new);
                            dd({ 83, 36, 36, -83 });

                            L(const_00000200000002000000020000000200);
                            dd({ 0x200 }, 4);

                            L(const_00000008000000080000000800000008);
                            dd({ 8 }, 4);

                        }

                        void assemble()
                        {
                            //void havoc_dct_16x16_ssse3(int16_t *coeffs, const int16_t *src, intptr_t src_stride)
                            //{
                            //	HAVOC_ALIGN(32, int16_t, temp[16 * 16]);
                            //	havoc_partial_butterfly_16h_ssse3(temp, src, src_stride, 3);
                            //	havoc_partial_butterfly_16v_ssse3(coeffs, temp, 10);
                            //}

                            auto &r0 = arg64(0);
                            auto &r1 = arg64(1);
                            auto &r2 = arg64(2);
                            auto &r3 = reg64(3);
                            Xbyak::Reg32 r3d{ r3.getIdx() };
                            auto &r4 = reg64(4);
                            Xbyak::Reg32 r4d{ r4.getIdx() };
                            auto &r5 = reg64(5);
                            Xbyak::Reg32 r5d{ r5.getIdx() };
                            //auto &r6 = reg64(6);
                            //auto &r7 = reg64(7);

                            auto &m0 = regXmm(0);
                            auto &m1 = regXmm(1);
                            auto &m2 = regXmm(2);
                            auto &m3 = regXmm(3);
                            auto &m4 = regXmm(4);
                            auto &m5 = regXmm(5);
                            auto &m6 = regXmm(6);
                            auto &m7 = regXmm(7);
                            auto &m8 = regXmm(8);
                            auto &m9 = regXmm(9);
                            auto &m10 = regXmm(10);
                            auto &m11 = regXmm(11);
                            auto &m12 = regXmm(12);
                            auto &m13 = regXmm(13);
                            auto &m14 = regXmm(14);
                            auto &m15 = regXmm(15);

                            mov(ptr[rsp], r0);
                            lea(r0, ptr[rsp + 32]);
                            //db({ 0xcc });
                            // void transform_partial_butterfly_16h_ssse3(int16_t *dst, const int16_t *src, std::intptr_t srcStride, int shift)//);
                            // shift parameter ignored (r3));
                            // INIT_XMM ssse3
                            // cglobal partial_butterfly_16h, 4, 7, 16
                            movdqa(m3, ptr[rip + shuffle_efcdab8967452301]);
                            movdqa(m4, ptr[rip + dd_times_4_add]);
                            mov(r4d, 16);

                            L("loop");
                            {
                                movdqu(m0, ptr[r1]);
                                movdqu(m1, ptr[r1 + 16]);
                                pshufb(m1, m3);

                                movdqa(m2, m0); paddw(m2, m1);
                                // m2 = Eptr[7:0]);

                                psubw(m0, m1);
                                // m0 = Optr[7:0]);

                                movdqa(m8, m0); pmaddwd(m8, ptr[rip + cosine_8_h + 0 * 16]);
                                movdqa(m9, m0); pmaddwd(m9, ptr[rip + cosine_8_h + 1 * 16]);
                                movdqa(m10, m0); pmaddwd(m10, ptr[rip + cosine_8_h + 2 * 16]);
                                movdqa(m11, m0); pmaddwd(m11, ptr[rip + cosine_8_h + 3 * 16]);
                                movdqa(m12, m0); pmaddwd(m12, ptr[rip + cosine_8_h + 4 * 16]);
                                movdqa(m13, m0); pmaddwd(m13, ptr[rip + cosine_8_h + 5 * 16]);
                                movdqa(m14, m0); pmaddwd(m14, ptr[rip + cosine_8_h + 6 * 16]);
                                movdqa(m15, m0); pmaddwd(m15, ptr[rip + cosine_8_h + 7 * 16]);

                                phaddd(m8, m9);
                                phaddd(m10, m11);
                                phaddd(m12, m13);
                                phaddd(m14, m15);

                                phaddd(m8, m10);
                                phaddd(m12, m14);

                                paddd(m8, m4);
                                psrad(m8, 3 + bitDepth - 8);
                                paddd(m12, m4);
                                psrad(m12, 3 + bitDepth - 8);

                                packssdw(m8, m12);
                                // m8 = dstptr[15,13,11, 9, 7, 5, 3, 1]);

                                movdqa(m5, m2); pshufb(m5, m3);
                                // m5 = Eptr[0:7]);

                                movdqa(m7, m2); psubw(m7, m5);
                                // m7 = -EOptr[0:3], EOptr[3:0]);

                                movdqa(m6, m7); pmaddwd(m6, ptr[rip + cosine_4_h]);
                                pmaddwd(m7, ptr[rip + cosine_4_h + 16]);
                                phaddd(m6, m7);
                                paddd(m6, m4);
                                psrad(m6, 3 + bitDepth - 8);
                                // dstptr[14, 10, 6, 2]);

                                pslld(m6, 16);
                                // m6 = dstptr[14], 0, dstptr[10], 0, dstptr[6], 0, dstptr[2], 0);

                                paddw(m2, m5);
                                // m12 = EEptr[0,1,2,3,3,2,1,0]);

                                pshufd(m11, m2, ORDER(3, 2, 3, 2));
                                // m11 = EEptr[0,1,2,3,0,1,2,3]);

                                movdqa(m10, m11); paddw(m10, m2);
                                // m10 = 2*EEptr[0,1,2,3], EEEptr[0,1,1,0]);

                                psubw(m11, m2);
                                // m11= 0,0,0,0,EEOptr[0,1],-EEOptr[1,0]);

                                punpckldq(m10, m11);
                                // m10 = EEOptr[0,1], EEEptr[0,1], -EEOptr[1,0], EEEptr[1,0]);

                                pmaddwd(m10, ptr[rip + cosine_1_h]);

                                paddd(m10, m4);
                                pslld(m10, 16 - (3 + bitDepth - 8));
                                // m10 = dstptr[12, 8, 4, 0] << 16);
                                // m10 = dstptr[12], 0, dstptr[8], 0, dstptr[4], 0, dstptr[0], 0);
                                psrld(m10, 16);

                                por(m10, m6);
                                // m10 = dstptr[14,12,10, 8, 6, 4, 2, 0]);

                                movdqa(m6, m10);
                                punpcklwd(m6, m8);
                                movdqu(ptr[r0], m6);

                                punpckhwd(m10, m8);
                                movdqu(ptr[r0 + 16], m10);

                                lea(r0, ptr[r0 + 2 * 16]);
                                lea(r1, ptr[r1 + r2 * 2]);
                                dec(r4d);
                            }
                            jg("loop");

                            // H transform done
                            //db({ 0xcc });

                            // V transform follows:
                            mov(r0, ptr[rsp]);
                            lea(r1, ptr[rsp + 32]);
                            mov(r2, 32);


                            mov(r5d, 4);
                            L("loop_left_right");
                            {
#define COMPUTE_4(a1, a2, a3) \
        movq(m##a1, ptr[r1 + a1 * 16 * 2]); \
        punpcklwd(m##a1, m##a1); \
        psrad(m##a1, 16); \
        /* m%1 = src(3:0, %1) */ \
        \
        movq(m8, ptr[r1 + (15 - a1) * 16 * 2]); \
        punpcklwd(m8, m8); \
        psrad(m8, 16); \
        /* m8 = src(3:0, 15-a1) */ \
        \
        a3(m##a1, m8); \
        /* m##a1 = Eptr[0](3:0) when a3 is paddd */ \
        /* m##a1 = Optr[0](3:0) when a3 is psubd */ \
        \
        movq(m##a2, ptr[r1 + (a2) * 16 * 2]); \
        punpcklwd(m##a2, m##a2); \
        psrad(m##a2, 16); \
        /* m##a2 = src(3:0, a2) */ \
        \
        movq(m9, ptr[r1 + (15 - a2) * 16 * 2]); \
        punpcklwd(m9, m9); \
        psrad(m9, 16); \
        /* m9 = src(3:0, 15-a2) */ \
        \
        a3(m##a2, m9); \
        /* m##a2 = Eptr[1](3:0) when a3 is paddd */ \
        /* m##a2 = Optr[1](3:0) when a3 is psubd */ \

                                COMPUTE_4(0, 1, psubd);
        COMPUTE_4(2, 3, psubd);
        COMPUTE_4(4, 5, psubd);
        COMPUTE_4(6, 7, psubd);
        // mx = Optr[x](3:0));

        lea(r0, ptr[r0 + 1 * 16 * 2]);
        lea(r4, ptr[rip + cosine_8_new]);
        mov(r3d, 8);
        L("loop_write_odd");
        {
            movdqa(m8, ptr[rip + const_00000200000002000000020000000200]);

            movdqa(m10, ptr[r4 + 0 * 16]);
            pshufd(m11, m10, ORDER(0, 0, 0, 0));
            pmulld(m11, m0);
            paddd(m8, m11);
            pshufd(m12, m10, ORDER(1, 1, 1, 1));
            pmulld(m12, m1);
            paddd(m8, m12);
            pshufd(m13, m10, ORDER(2, 2, 2, 2));
            pmulld(m13, m2);
            paddd(m8, m13);
            pshufd(m14, m10, ORDER(3, 3, 3, 3));
            pmulld(m14, m3);
            paddd(m8, m14);

            movdqa(m10, ptr[r4 + 1 * 16]);
            pshufd(m11, m10, ORDER(0, 0, 0, 0));
            pmulld(m11, m4);
            paddd(m8, m11);
            pshufd(m12, m10, ORDER(1, 1, 1, 1));
            pmulld(m12, m5);
            paddd(m8, m12);
            pshufd(m13, m10, ORDER(2, 2, 2, 2));
            pmulld(m13, m6);
            paddd(m8, m13);
            pshufd(m14, m10, ORDER(3, 3, 3, 3));
            pmulld(m14, m7);
            paddd(m8, m14);

            lea(r4, ptr[r4 + 2 * 16]);

            psrad(m8, 10);

            packssdw(m8, m8);
            movq(ptr[r0], m8);

            lea(r0, ptr[r0 + 2 * 16 * 2]);

            dec(r3d);
        }
        jg("loop_write_odd");

        lea(r0, ptr[r0 - (1 + 8 * 2) * 16 * 2]);

        COMPUTE_4(0, 1, paddd);
        COMPUTE_4(2, 3, paddd);
        COMPUTE_4(5, 4, paddd);
        COMPUTE_4(7, 6, paddd);
        // mx = Eptr[x](3:0));

        movdqa(m8, m0); psubd(m8, m7);
        movdqa(m9, m1); psubd(m9, m6);
        movdqa(m10, m2); psubd(m10, m5);
        movdqa(m11, m3); psubd(m11, m4);
        // m8:11 = EOptr[0:3] (3:0));

        paddd(m0, m7);
        paddd(m1, m6);
        paddd(m2, m5);
        paddd(m3, m4);
        // m0:3 = EEptr[0:3] (3:0));

        movdqa(m4, m0); psubd(m4, m3);
        movdqa(m5, m1); psubd(m5, m2);
        // m4:5 = EEOptr[0:1] (3:0));

        paddd(m0, m3);
        paddd(m1, m2);
        // m0:1 = EEEptr[0:1] (3:0));

        paddd(m0, ptr[rip + const_00000008000000080000000800000008]);

        lea(r0, ptr[r0 + 2 * 16 * 2]);
        lea(r4, ptr[rip + cosine_4_new]);
        mov(r3d, 4);
        L("loop_write_even_odd");
        {
            movdqa(m2, ptr[rip + const_00000200000002000000020000000200]);

            movdqa(m3, ptr[r4 + 0 * 16]);
            pshufd(m12, m3, ORDER(0, 0, 0, 0));
            pmulld(m12, m8);
            paddd(m2, m12);
            pshufd(m13, m3, ORDER(1, 1, 1, 1));
            pmulld(m13, m9);
            paddd(m2, m13);
            pshufd(m14, m3, ORDER(2, 2, 2, 2));
            pmulld(m14, m10);
            paddd(m2, m14);
            pshufd(m15, m3, ORDER(3, 3, 3, 3));
            pmulld(m15, m11);
            paddd(m2, m15);

            psrad(m2, 10);
            packssdw(m2, m2);
            movq(ptr[r0], m2);

            lea(r4, ptr[r4 + 1 * 16]);
            lea(r0, ptr[r0 + 4 * 16 * 2]);
            dec(r3d);
        }
        jg("loop_write_even_odd");
        lea(r0, ptr[r0 - (2 + 4 * 4) * 16 * 2]);

        movdqa(m3, ptr[rip + cosine_2_new]);

        movdqa(m2, ptr[rip + const_00000200000002000000020000000200]);
        pshufd(m12, m3, ORDER(0, 0, 0, 0));
        pmulld(m12, m4);
        paddd(m2, m12);
        pshufd(m13, m3, ORDER(1, 1, 1, 1));
        pmulld(m13, m5);
        paddd(m2, m13);
        psrad(m2, 10);
        packssdw(m2, m2);
        movq(ptr[r0 + (4) * 16 * 2], m2);

        movdqa(m2, ptr[rip + const_00000200000002000000020000000200]);
        pshufd(m12, m3, ORDER(2, 2, 2, 2));
        pmulld(m12, m4);
        paddd(m2, m12);
        pshufd(m13, m3, ORDER(3, 3, 3, 3));
        pmulld(m13, m5);
        paddd(m2, m13);
        psrad(m2, 10);
        packssdw(m2, m2);
        movq(ptr[r0 + (12) * 16 * 2], m2);

        // m0:1 = EEEptr[0:1] (3:0));
        movdqa(m2, m0); paddd(m2, m1);
        psubd(m0, m1);

        psrad(m2, 4);
        psrad(m0, 4);

        packssdw(m2, m2); // 50% utilization
        movq(ptr[r0 + (0) * 16 * 2], m2);

        packssdw(m0, m0); // 50% utilization
        movq(ptr[r0 + (8) * 16 * 2], m0);

        lea(r1, ptr[r1 + 8]);
        lea(r0, ptr[r0 + 8]);
        dec(r5d);
                            }
                            jg("loop_left_right");
                        }
                    };


                    template <int bitDepth>
                    havoc_transform* get_transform(int trType, int log2TrafoSize, havoc_code code)
                    {
                        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
                        const int nCbS = 1 << log2TrafoSize;

                        havoc_transform *f = 0;

                        if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT))
                        {
                            if (nCbS == 4) f = trType ? havoc_dst_4x4_c_opt<bitDepth> : havoc_dct_4x4_c_opt<bitDepth>;
                            if (nCbS == 8) f = havoc_dct_8x8_c_opt<bitDepth>;
                            if (nCbS == 16) f = havoc_dct_16x16_c_opt<bitDepth>;
                            if (nCbS == 32) f = havoc_dct_32x32_c_opt<bitDepth>;
                        }

                        if (buffer.isa & HAVOC_SSSE3)
                        {
                            if (nCbS == 16)
                            {
                                ForwardDct16x16<bitDepth> a(&buffer);
                                f = a;
                            }
                        }

                        return f;
                    }


                    template <int bitDepth>
                    void havoc_populate_transform(havoc_table_transform<bitDepth> *table, havoc_code code)
                    {
                        *havoc_get_transform(table, 1, 2) = get_transform<bitDepth>(1, 2, code);
                        for (int log2TrafoSize = 2; log2TrafoSize <= 5; ++log2TrafoSize)
                        {
                            *havoc_get_transform(table, 0, log2TrafoSize) = get_transform<bitDepth>(0, log2TrafoSize, code);
                        }
                    }


                    template void havoc_populate_transform<8>(havoc_table_transform<8> *table, havoc_code code);
                    template void havoc_populate_transform<10>(havoc_table_transform<10> *table, havoc_code code);


                    typedef struct
                    {
                        havoc_transform *f;
                        HAVOC_ALIGN(32, int16_t, dst[32 * 32]);
                        int16_t *src;
                        intptr_t src_stride;
                        int trType;
                        int log2TrafoSize;
                        int bitDepth;
                    }
                    bound_transform;


                    int init_transform(void *p, havoc_code code)
                    {
                        auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

                        bound_transform *s = (bound_transform *)p;

                        if (s->bitDepth == 8)
                        {
                            havoc_table_transform<8> table;
                            havoc_populate_transform(&table, code);
                            s->f = *havoc_get_transform(&table, s->trType, s->log2TrafoSize);
                        }
                        else
                        {
                            assert(s->bitDepth == 10);
                            havoc_table_transform<10> table;
                            havoc_populate_transform(&table, code);
                            s->f = *havoc_get_transform(&table, s->trType, s->log2TrafoSize);
                        }

                        if (s->f && buffer.isa == HAVOC_C_REF)
                        {
                            const int nCbS = 1 << s->log2TrafoSize;
                            printf("\t%d-bit %s %dx%d : ", s->bitDepth, s->trType ? "sine" : "cosine", nCbS, nCbS);
                        }

                        for (int x = 0; x < 32 * 32; x++) s->dst[x] = 0xab;

                        return !!s->f;
                    }


                    void invoke_transform(void *p, int n)
                    {
                        bound_transform *s = (bound_transform *)p;

                        while (n--)
                        {
                            s->f(s->dst, s->src, s->src_stride);
                        }
                    }


                    int mismatch_transform(void *boundRef, void *boundTest)
                    {
                        bound_transform *ref = (bound_transform *)boundRef;
                        bound_transform *test = (bound_transform *)boundTest;

                        const int nCbS = 1 << ref->log2TrafoSize;

                        return memcmp(ref->dst, test->dst, nCbS * nCbS * sizeof(int16_t));
                    }


                    void havoc_test_transform(int *error_count, havoc_instruction_set mask)
                    {
                        printf("\nhavoc_transform - Forward Transform\n");

                        HAVOC_ALIGN(32, int16_t, src[32 * 32]);
                        for (int x = 0; x < 32 * 32; x++) src[x] = (rand() & 0x1ff) - 0x100;

                        bound_transform b[2];
                        b[0].src = src;
                        b[0].src_stride = 32;

                        for (b[0].bitDepth = 8; b[0].bitDepth <= 10; b[0].bitDepth += 2)
                        {
                            for (int j = 2; j < 6; ++j)
                            {
                                b[0].trType = (j == 1) ? 1 : 0;
                                b[0].log2TrafoSize = (j == 1) ? 2 : j;

                                b[1] = b[0];

                                *error_count += havoc_test(&b[0], &b[1], init_transform, invoke_transform, mismatch_transform, mask, 100000);
                            }
                        }
                    }
