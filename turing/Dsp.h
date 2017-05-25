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

#ifndef INCLUDED_Dsp_h
#define INCLUDED_Dsp_h

#pragma once

#include "IntraReferenceSamples.h"
#include "GlobalState.h"
#include "Global.h"
#include "Picture.h"
#include "Memory.h"
#include "havoc/pred_inter.h"
#include <cstring>

template <typename Sample>
static Sample p(const Sample *neighbours, const Sample corner, int dx, int dy)
{
    if (dx < 0 && dy < 0)
    {
    }
    else if (dy < 0)
        {
        assert(dy == -1);
        return neighbours[64 + dx];
        }
    else if (dx < 0)
        {
        assert(dx == -1);
        return neighbours[63 - dy];
    }

    assert(dx == -1);
    assert(dy == -1);
    return corner;
}


static bool inline filterFlag(int cIdx, int predModeIntra, int nTbS)
{
    uint8_t const lookup[] =
    {
        0b111000, 0b000000, 0b111000, 0b110000, 0b110000,
        0b110000, 0b110000, 0b110000, 0b110000, 0b100000,
        0b000000, 0b100000, 0b110000, 0b110000, 0b110000,
        0b110000, 0b110000, 0b110000, 0b111000, 0b110000,
        0b110000, 0b110000, 0b110000, 0b110000, 0b110000,
        0b100000, 0b000000, 0b100000, 0b110000, 0b110000,
        0b110000, 0b110000, 0b110000, 0b110000, 0b111000,
    };

    return cIdx == 0 && !!(lookup[predModeIntra] & nTbS);
}


// Scaling process for transform coefficients
// Output of this process is the (nTbS)x(nTbS) array d of scaled transform coefficients with elements d[ x ][ y ].
static void inverseQuantize(Raster<int16_t> d, Raster<const int16_t> transCoeffLevel, Raster<const int> m, int log2TrafoSize, int qP, int bitDepth, bool check = false)
{
    const int bdShift = bitDepth + log2TrafoSize - 5;

    for (int y=0; y<(1<<log2TrafoSize); ++y)
    {
        for (int x=0; x<(1<<log2TrafoSize); ++x)
        {
            static const int levelScale[6] = {40, 45, 51, 57, 64, 72};

            auto v = (int16_t)Clip3<int64_t>(
                    -32768,
                    32767,
                    (( int64_t(transCoeffLevel(x, y)) * m(x, y) * levelScale[qP % 6]  << (qP/6) ) + (int64_t(1)<<(bdShift-1)) ) >> bdShift);

            if (check && v != d(x, y)) throw 0;

            d(x, y) = v;
        }
    }
}


template <class T>
static T clipCidx1(T x, int bitDepth)
{
    return Clip3<T>(0, (1<<bitDepth) - 1, x);
}


template <typename Sample, class H>
int16_t LumaSampleInterpolate(H &h, int xIntL, int yIntL, int xFracL, int yFracL, Raster<Sample> &refPicLXL, int bitDepthY)
{
    int predSampleLXL;
    int toffset = 3;
    int xA[8][8];
    int yA[8][8];
    for(int i=-3; i<=4; i++)
    {
        for(int j=-3; j<=4; j++)
        {
            xA[i+toffset][j+toffset] = Clip3(0, h[pic_width_in_luma_samples()]-1, xIntL+i);
            yA[i+toffset][j+toffset] = Clip3(0, h[pic_height_in_luma_samples()]-1, yIntL+j);
        }
    }
    int shift1 = bitDepthY-8;	//shift1 = 0;
    int shift2 = 6;				//shift2 = 6;
    int shift3 = 14-bitDepthY;	//shift3 = 6;

    if(xFracL==0 && yFracL==0)
    {
        int cStride = toffset;
        predSampleLXL = refPicLXL(xA[cStride][cStride],yA[cStride][cStride]) << shift3;
    }
    if(xFracL!=0 && yFracL==0) //Apply Filter in Horizontal direction to obtain positions a, b, c
    {
        int cStride = toffset;
        if(xFracL==1)
            predSampleLXL = (refPicLXL(xA[-3+cStride][cStride],yA[-3+cStride][cStride]) * -1 + refPicLXL(xA[-2+cStride][cStride],yA[-2+cStride][cStride]) * 4 + refPicLXL(xA[-1+cStride][cStride],yA[-1+cStride][cStride])* -10 + refPicLXL(xA[cStride][cStride],yA[cStride][cStride]) * 58 + refPicLXL(xA[1+cStride][cStride],yA[1+cStride][cStride]) * 17 + refPicLXL(xA[2+cStride][cStride],yA[2+cStride][cStride]) * -5 +  refPicLXL(xA[3+cStride][cStride],yA[3+cStride][cStride]) * 1)																	>> shift1;
        if(xFracL==2)
            predSampleLXL = (refPicLXL(xA[-3+cStride][cStride],yA[-3+cStride][cStride]) * -1 + refPicLXL(xA[-2+cStride][cStride],yA[-2+cStride][cStride]) * 4 + refPicLXL(xA[-1+cStride][cStride],yA[-1+cStride][cStride])* -11 + refPicLXL(xA[cStride][cStride],yA[cStride][cStride]) * 40 + refPicLXL(xA[1+cStride][cStride],yA[1+cStride][cStride]) * 40 + refPicLXL(xA[2+cStride][cStride],yA[2+cStride][cStride]) * -11 + refPicLXL(xA[3+cStride][cStride],yA[3+cStride][cStride]) * 4 + refPicLXL(xA[4+cStride][cStride],yA[4+cStride][cStride]) * -1) >> shift1;
        if(xFracL==3)
            predSampleLXL =							                                         (refPicLXL(xA[-2+cStride][cStride],yA[-2+cStride][cStride]) * 1 + refPicLXL(xA[-1+cStride][cStride],yA[-1+cStride][cStride]) * -5 + refPicLXL(xA[cStride][cStride],yA[cStride][cStride]) * 17 + refPicLXL(xA[1+cStride][cStride],yA[1+cStride][cStride]) * 58 + refPicLXL(xA[2+cStride][cStride],yA[2+cStride][cStride]) * -10 + refPicLXL(xA[3+cStride][cStride],yA[3+cStride][cStride]) * 4 + refPicLXL(xA[4+cStride][cStride],yA[4+cStride][cStride]) * -1) >> shift1;
    }
    if(xFracL==0 && yFracL!=0) // Apply filter in vertical direction to obtain positions d, h, n
    {
        int cStride = toffset;
        if(yFracL==1)
            predSampleLXL = (refPicLXL(xA[cStride][-3+cStride],yA[cStride][-3+cStride]) * -1 + refPicLXL(xA[cStride][-2+cStride],yA[cStride][-2+cStride]) * 4 + refPicLXL(xA[cStride][-1+cStride],yA[cStride][-1+cStride])* -10 + refPicLXL(xA[cStride][cStride],yA[cStride][cStride]) * 58 + refPicLXL(xA[cStride][1+cStride],yA[cStride][1+cStride]) * 17 + refPicLXL(xA[cStride][2+cStride],yA[cStride][2+cStride]) *  -5 + refPicLXL(xA[cStride][3+cStride],yA[cStride][3+cStride]) * 1)																	>> shift1;
        if(yFracL==2)
            predSampleLXL = (refPicLXL(xA[cStride][-3+cStride],yA[cStride][-3+cStride]) * -1 + refPicLXL(xA[cStride][-2+cStride],yA[cStride][-2+cStride]) * 4 + refPicLXL(xA[cStride][-1+cStride],yA[cStride][-1+cStride])* -11 + refPicLXL(xA[cStride][cStride],yA[cStride][cStride]) * 40 + refPicLXL(xA[cStride][1+cStride],yA[cStride][1+cStride]) * 40 + refPicLXL(xA[cStride][2+cStride],yA[cStride][2+cStride]) * -11 + refPicLXL(xA[cStride][3+cStride],yA[cStride][3+cStride]) * 4 + refPicLXL(xA[cStride][4+cStride],yA[cStride][4+cStride]) * -1)	>> shift1;
        if(yFracL==3)
            predSampleLXL =							                                         (refPicLXL(xA[cStride][-2+cStride],yA[cStride][-2+cStride]) * 1 + refPicLXL(xA[cStride][-1+cStride],yA[cStride][-1+cStride]) * -5 + refPicLXL(xA[cStride][cStride],yA[cStride][cStride]) * 17 + refPicLXL(xA[cStride][1+cStride],yA[cStride][1+cStride]) * 58 + refPicLXL(xA[cStride][2+cStride],yA[cStride][2+cStride]) * -10 + refPicLXL(xA[cStride][3+cStride],yA[cStride][3+cStride]) * 4 + refPicLXL(xA[cStride][4+cStride],yA[cStride][4+cStride]) * -1)	>> shift1;
    }
    if(xFracL!=0 && yFracL!=0) // Apply filter in both horizontal and vertical directions
    {
        int a[8];
        int b[8];
        int c[8];
        for(int i = -3; i<=4; i++) //Apply filter in horizontal direction to obtain positions ai,bi,ci for i=-3:4;
        {
            int cStride = toffset;
            int cStride2 = 3;
            a[i+cStride] = (refPicLXL(xA[-3+cStride][cStride2+i],yA[-3+cStride][cStride2+i]) * -1 + refPicLXL(xA[-2+cStride][cStride2+i],yA[-2+cStride][cStride2+i]) * 4 + refPicLXL(xA[-1+cStride][cStride2+i],yA[-1+cStride][cStride2+i]) * -10 + refPicLXL(xA[cStride][cStride2+i],yA[cStride][cStride2+i]) * 58 + refPicLXL(xA[1+cStride][cStride2+i],yA[1+cStride][cStride2+i]) * 17  + refPicLXL(xA[2+cStride][cStride2+i],yA[2+cStride][cStride2+i]) * -5  + refPicLXL(xA[3+cStride][cStride2+i],yA[3+cStride][cStride2+i]) * 1)																			>> shift1;
            b[i+cStride] = (refPicLXL(xA[-3+cStride][cStride2+i],yA[-3+cStride][cStride2+i]) * -1 + refPicLXL(xA[-2+cStride][cStride2+i],yA[-2+cStride][cStride2+i]) * 4 + refPicLXL(xA[-1+cStride][cStride2+i],yA[-1+cStride][cStride2+i]) * -11 + refPicLXL(xA[cStride][cStride2+i],yA[cStride][cStride2+i]) * 40 + refPicLXL(xA[1+cStride][cStride2+i],yA[1+cStride][cStride2+i]) * 40 + refPicLXL(xA[2+cStride][cStride2+i],yA[2+cStride][cStride2+i]) * -11 + refPicLXL(xA[3+cStride][cStride2+i],yA[3+cStride][cStride2+i]) * 4 + refPicLXL(xA[4+cStride][cStride2+i],yA[4+cStride][cStride2+i]) * -1)	>> shift1;
            c[i+cStride] = 									                                       (refPicLXL(xA[-2+cStride][cStride2+i],yA[-2+cStride][cStride2+i]) * 1 + refPicLXL(xA[-1+cStride][cStride2+i],yA[-1+cStride][cStride2+i]) * -5  + refPicLXL(xA[cStride][cStride2+i],yA[cStride][cStride2+i]) * 17 + refPicLXL(xA[1+cStride][cStride2+i],yA[1+cStride][cStride2+i]) * 58  + refPicLXL(xA[2+cStride][cStride2+i],yA[2+cStride][cStride2+i]) * -10 + refPicLXL(xA[3+cStride][cStride2+i],yA[3+cStride][cStride2+i]) * 4 + refPicLXL(xA[4+cStride][cStride2+i],yA[4+cStride][cStride2+i]) * -1)	>> shift1;
        }
        if(xFracL==1)			//Apply filter in vertical direction to obtain positions e, i, p;
        {
            int cStride = toffset;
            if(yFracL==1)
                predSampleLXL = (a[-3+cStride] * -1 + a[-2+cStride] * 4 + a[-1+cStride]* -10 + a[0+cStride] * 58 + a[1+cStride] * 17 + a[2+cStride] * -5 + a[3+cStride] * 1)						>> shift2;
            if(yFracL==2)
                predSampleLXL = (a[-3+cStride] * -1 + a[-2+cStride] * 4 + a[-1+cStride]* -11 + a[0+cStride] * 40 + a[1+cStride] * 40 + a[2+cStride] * -11 + a[3+cStride] * 4 + a[4+cStride] * -1)	>> shift2;
            if(yFracL==3)
                predSampleLXL =					    (a[-2+cStride] * 1 + a[-1+cStride] * -5 + a[0+cStride] * 17 + a[1+cStride] * 58 + a[2+cStride] * -10 + a[3+cStride] * 4 + a[4+cStride] * -1)	>> shift2;
        }
        if(xFracL==2)			//Apply filter in vertical direction to obtain positions f, j, q;
        {
            int cStride = toffset;
            if(yFracL==1)
                predSampleLXL = (b[-3+cStride] * -1 + b[-2+cStride] * 4 + b[-1+cStride]* -10 + b[0+cStride] * 58 + b[1+cStride] * 17 + b[2+cStride] * -5 + b[3+cStride] * 1)						>> shift2;
            if(yFracL==2)
                predSampleLXL = (b[-3+cStride] * -1 + b[-2+cStride] * 4 + b[-1+cStride]* -11 + b[0+cStride] * 40 + b[1+cStride] * 40 + b[2+cStride] * -11 + b[3+cStride] * 4 + b[4+cStride] * -1)	>> shift2;
            if(yFracL==3)
                predSampleLXL =					    (b[-2+cStride] * 1 + b[-1+cStride] * -5 + b[0+cStride] * 17 + b[1+cStride] * 58 + b[2+cStride] * -10 + b[3+cStride] * 4 + b[4+cStride] * -1)	>> shift2;
        }
        if(xFracL==3)			//Apply filter in vertical direction to obtain positions g, k, r;
        {
            int cStride = toffset;
            if(yFracL==1)
                predSampleLXL = (c[-3+cStride] * -1 + c[-2+cStride] * 4 + c[-1+cStride]* -10 + c[0+cStride] * 58 + c[1+cStride] * 17 + c[2+cStride] * -5 + c[3+cStride] * 1)						>> shift2;
            if(yFracL==2)
                predSampleLXL = (c[-3+cStride] * -1 + c[-2+cStride] * 4 + c[-1+cStride]* -11 + c[0+cStride] * 40 + c[1+cStride] * 40 + c[2+cStride] * -11 + c[3+cStride] * 4 + c[4+cStride] * -1)	>> shift2;
            if(yFracL==3)
                predSampleLXL =					    (c[-2+cStride] * 1 + c[-1+cStride] * -5 + c[0+cStride] * 17 + c[1+cStride] * 58 + c[2+cStride] * -10 + c[3+cStride] * 4 + c[4+cStride] * -1)	>> shift2;
        }
    }
    return predSampleLXL;
}


template <typename Sample, class H>
int16_t ChromaSampleInterpolate(H &h, int xIntC, int yIntC, int xFracC, int yFracC, Raster<Sample> &refPicLXC, int bitDepthC)
{
    int toffset = 1;
    int xB[4][4];
    int yB[4][4];
    for(int i=-1; i<=2; i++)
    {
        for(int j=-1; j<=2; j++)
        {
            xB[i+toffset][j+toffset] = Clip3(0,(h[pic_width_in_luma_samples()]/h[SubWidthC()])-1,xIntC+i);
            yB[i+toffset][j+toffset] = Clip3(0,(h[pic_height_in_luma_samples()]/h[SubHeightC()])-1,yIntC+j);
        }
    }
    int shift1 = bitDepthC-8;	//shift1 = 0;
    int shift2 = 6;				//shift2 = 6;
    int shift3 = 14-bitDepthC;	//shift3 = 6;

    int predSampleLXC;
    if(xFracC==0 && yFracC==0)
    {
        int cStride = toffset;
        predSampleLXC = refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) << shift3;
    }
    if(xFracC!=0 && yFracC==0) // Apply filter in horizontal direction to obtain positions ab, ac, ad, ae, af, ag, ah
    {
        int cStride = toffset;
        if(xFracC==1)
            predSampleLXC = (refPicLXC(xB[-1+cStride][cStride],yB[-1+cStride][cStride])* -2 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 58 + refPicLXC(xB[1+cStride][cStride],yB[1+cStride][cStride]) * 10 + refPicLXC(xB[2+cStride][cStride],yB[2+cStride][cStride]) * -2)  >> shift1;
        if(xFracC==2)
            predSampleLXC = (refPicLXC(xB[-1+cStride][cStride],yB[-1+cStride][cStride])* -4 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 54 + refPicLXC(xB[1+cStride][cStride],yB[1+cStride][cStride]) * 16 + refPicLXC(xB[2+cStride][cStride],yB[2+cStride][cStride]) * -2) >> shift1;
        if(xFracC==3)
            predSampleLXC =	(refPicLXC(xB[-1+cStride][cStride],yB[-1+cStride][cStride]) * -6 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 46 + refPicLXC(xB[1+cStride][cStride],yB[1+cStride][cStride]) * 28 + refPicLXC(xB[2+cStride][cStride],yB[2+cStride][cStride]) * -4) >> shift1;
        if(xFracC==4)
            predSampleLXC = (refPicLXC(xB[-1+cStride][cStride],yB[-1+cStride][cStride])* -4 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 36 + refPicLXC(xB[1+cStride][cStride],yB[1+cStride][cStride]) * 36 + refPicLXC(xB[2+cStride][cStride],yB[2+cStride][cStride]) * -4) >> shift1;
        if(xFracC==5)
            predSampleLXC = (refPicLXC(xB[-1+cStride][cStride],yB[-1+cStride][cStride])* -4 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 28 + refPicLXC(xB[1+cStride][cStride],yB[1+cStride][cStride]) * 46 + refPicLXC(xB[2+cStride][cStride],yB[2+cStride][cStride]) * -6) >> shift1;
        if(xFracC==6)
            predSampleLXC =	(refPicLXC(xB[-1+cStride][cStride],yB[-1+cStride][cStride]) * -2 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 16 + refPicLXC(xB[1+cStride][cStride],yB[1+cStride][cStride]) * 54 + refPicLXC(xB[2+cStride][cStride],yB[2+cStride][cStride]) * -4) >> shift1;
        if(xFracC==7)
            predSampleLXC =	(refPicLXC(xB[-1+cStride][cStride],yB[-1+cStride][cStride]) * -2 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 10 + refPicLXC(xB[1+cStride][cStride],yB[1+cStride][cStride]) * 58 + refPicLXC(xB[2+cStride][cStride],yB[2+cStride][cStride]) * -2) >> shift1;
    }
    if(xFracC==0 && yFracC!=0) //Apply filter in vertical direction to obtain positions ba, ca, da, ea, fa, ga, ha
    {
        int cStride = toffset;
        if(yFracC==1)
            predSampleLXC = (refPicLXC(xB[cStride][-1+cStride],yB[cStride][-1+cStride])* -2 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 58 + refPicLXC(xB[cStride][1+cStride],yB[cStride][1+cStride]) * 10 + refPicLXC(xB[cStride][2+cStride],yB[cStride][2+cStride]) * -2)  >> shift1;
        if(yFracC==2)
            predSampleLXC = (refPicLXC(xB[cStride][-1+cStride],yB[cStride][-1+cStride])* -4 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 54 + refPicLXC(xB[cStride][1+cStride],yB[cStride][1+cStride]) * 16 + refPicLXC(xB[cStride][2+cStride],yB[cStride][2+cStride]) * -2) >> shift1;
        if(yFracC==3)
            predSampleLXC =	(refPicLXC(xB[cStride][-1+cStride],yB[cStride][-1+cStride]) * -6 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 46 + refPicLXC(xB[cStride][1+cStride],yB[cStride][1+cStride]) * 28 + refPicLXC(xB[cStride][2+cStride],yB[cStride][2+cStride]) * -4) >> shift1;
        if(yFracC==4)
            predSampleLXC = (refPicLXC(xB[cStride][-1+cStride],yB[cStride][-1+cStride])* -4 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 36 + refPicLXC(xB[cStride][1+cStride],yB[cStride][1+cStride]) * 36 + refPicLXC(xB[cStride][2+cStride],yB[cStride][2+cStride]) * -4) >> shift1;
        if(yFracC==5)
            predSampleLXC = (refPicLXC(xB[cStride][-1+cStride],yB[cStride][-1+cStride])* -4 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 28 + refPicLXC(xB[cStride][1+cStride],yB[cStride][1+cStride]) * 46 + refPicLXC(xB[cStride][2+cStride],yB[cStride][2+cStride]) * -6) >> shift1;
        if(yFracC==6)
            predSampleLXC =	(refPicLXC(xB[cStride][-1+cStride],yB[cStride][-1+cStride]) * -2 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 16 + refPicLXC(xB[cStride][1+cStride],yB[cStride][1+cStride]) * 54 + refPicLXC(xB[cStride][2+cStride],yB[cStride][2+cStride]) * -4) >> shift1;
        if(yFracC==7)
            predSampleLXC =	(refPicLXC(xB[cStride][-1+cStride],yB[cStride][-1+cStride]) * -2 + refPicLXC(xB[cStride][cStride],yB[cStride][cStride]) * 10 + refPicLXC(xB[cStride][1+cStride],yB[cStride][1+cStride]) * 58 + refPicLXC(xB[cStride][2+cStride],yB[cStride][2+cStride]) * -2) >> shift1;
    }
    if(xFracC!=0 && yFracC!=0)
    {
        int ab[4];
        int ac[4];
        int ad[4];
        int ae[4];
        int af[4];
        int ag[4];
        int ah[4];

        for(int i = -1; i<=2; i++) //Apply filter in horizontal direction to obtain positions abi,aci,adi,aei,afi,agi,ahi for i=-1:2;
        {
            int cStride = toffset;
            int cStride2 = 1;//2*cStride;
            ab[i+cStride] =  (refPicLXC(xB[-1+cStride][cStride2+i],yB[-1+cStride][cStride2+i]) * -2 + refPicLXC(xB[cStride][cStride2+i],yB[cStride][cStride2+i]) * 58 + refPicLXC(xB[1+cStride][cStride2+i],yB[1+cStride][cStride2+i]) * 10  + refPicLXC(xB[2+cStride][cStride2+i],yB[2+cStride][cStride2+i]) * -2) >> shift1;
            ac[i+cStride] =  (refPicLXC(xB[-1+cStride][cStride2+i],yB[-1+cStride][cStride2+i]) * -4 + refPicLXC(xB[cStride][cStride2+i],yB[cStride][cStride2+i]) * 54 + refPicLXC(xB[1+cStride][cStride2+i],yB[1+cStride][cStride2+i]) * 16 + refPicLXC(xB[2+cStride][cStride2+i],yB[2+cStride][cStride2+i]) * -2 ) >> shift1;
            ad[i+cStride] =  (refPicLXC(xB[-1+cStride][cStride2+i],yB[-1+cStride][cStride2+i]) * -6  + refPicLXC(xB[cStride][cStride2+i],yB[cStride][cStride2+i]) *46 + refPicLXC(xB[1+cStride][cStride2+i],yB[1+cStride][cStride2+i]) * 28  + refPicLXC(xB[2+cStride][cStride2+i],yB[2+cStride][cStride2+i]) * -4) >> shift1;
            ae[i+cStride] =  (refPicLXC(xB[-1+cStride][cStride2+i],yB[-1+cStride][cStride2+i]) * -4 + refPicLXC(xB[cStride][cStride2+i],yB[cStride][cStride2+i]) * 36 + refPicLXC(xB[1+cStride][cStride2+i],yB[1+cStride][cStride2+i]) * 36  + refPicLXC(xB[2+cStride][cStride2+i],yB[2+cStride][cStride2+i]) * -4) >> shift1;
            af[i+cStride] =  (refPicLXC(xB[-1+cStride][cStride2+i],yB[-1+cStride][cStride2+i]) * -4 + refPicLXC(xB[cStride][cStride2+i],yB[cStride][cStride2+i]) * 28 + refPicLXC(xB[1+cStride][cStride2+i],yB[1+cStride][cStride2+i]) * 46 + refPicLXC(xB[2+cStride][cStride2+i],yB[2+cStride][cStride2+i]) * -6 ) >> shift1;
            ag[i+cStride] =  (refPicLXC(xB[-1+cStride][cStride2+i],yB[-1+cStride][cStride2+i]) * -2  + refPicLXC(xB[cStride][cStride2+i],yB[cStride][cStride2+i]) *16 + refPicLXC(xB[1+cStride][cStride2+i],yB[1+cStride][cStride2+i]) * 54  + refPicLXC(xB[2+cStride][cStride2+i],yB[2+cStride][cStride2+i]) * -4) >> shift1;
            ah[i+cStride] =  (refPicLXC(xB[-1+cStride][cStride2+i],yB[-1+cStride][cStride2+i]) * -2 + refPicLXC(xB[cStride][cStride2+i],yB[cStride][cStride2+i]) * 10 + refPicLXC(xB[1+cStride][cStride2+i],yB[1+cStride][cStride2+i]) * 58  + refPicLXC(xB[2+cStride][cStride2+i],yB[2+cStride][cStride2+i]) * -2) >> shift1;
        }
        if(xFracC==1)
        {
            int cStride = toffset;
            if(yFracC==1)
            {
                predSampleLXC = (ab[-1+cStride]* -2 + ab[0+cStride] * 58 + ab[1+cStride] * 10 + ab[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==2)
            {
                predSampleLXC = (ab[-1+cStride]* -4 + ab[0+cStride] * 54 + ab[1+cStride] * 16 + ab[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==3)
            {
                predSampleLXC =	(ab[-1+cStride] * -6 + ab[0+cStride] * 46 + ab[1+cStride] * 28 + ab[2+cStride] * -4)>> shift2;
            }
            if(yFracC==4)
            {
                predSampleLXC =	(ab[-1+cStride] * -4 + ab[0+cStride] * 36 + ab[1+cStride] * 36 + ab[2+cStride] * -4)>> shift2;
            }
            if(yFracC==5)
            {
                predSampleLXC =	(ab[-1+cStride] * -4 + ab[0+cStride] * 28 + ab[1+cStride] * 46 + ab[2+cStride] * -6)>> shift2;
            }
            if(yFracC==6)
            {
                predSampleLXC =	(ab[-1+cStride] * -2 + ab[0+cStride] * 16 + ab[1+cStride] * 54 + ab[2+cStride] * -4)>> shift2;
            }
            if(yFracC==7)
            {
                predSampleLXC =	(ab[-1+cStride] * -2 + ab[0+cStride] * 10 + ab[1+cStride] * 58 + ab[2+cStride] * -2)>> shift2;
            }
        }
        if(xFracC==2)
        {
            int cStride = toffset;
            if(yFracC==1)
            {
                predSampleLXC = (ac[-1+cStride]* -2 + ac[0+cStride] * 58 + ac[1+cStride] * 10 + ac[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==2)
            {
                predSampleLXC = (ac[-1+cStride]* -4 + ac[0+cStride] * 54 + ac[1+cStride] * 16 + ac[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==3)
            {
                predSampleLXC =	(ac[-1+cStride] * -6 + ac[0+cStride] * 46 + ac[1+cStride] * 28 + ac[2+cStride] * -4)>> shift2;
            }
            if(yFracC==4)
            {
                predSampleLXC =	(ac[-1+cStride] * -4 + ac[0+cStride] * 36 + ac[1+cStride] * 36 + ac[2+cStride] * -4)>> shift2;
            }
            if(yFracC==5)
            {
                predSampleLXC =	(ac[-1+cStride] * -4 + ac[0+cStride] * 28 + ac[1+cStride] * 46 + ac[2+cStride] * -6)>> shift2;
            }
            if(yFracC==6)
            {
                predSampleLXC =	(ac[-1+cStride] * -2 + ac[0+cStride] * 16 + ac[1+cStride] * 54 + ac[2+cStride] * -4)>> shift2;
            }
            if(yFracC==7)
            {
                predSampleLXC =	(ac[-1+cStride] * -2 + ac[0+cStride] * 10 + ac[1+cStride] * 58 + ac[2+cStride] * -2)>> shift2;
            }
        }
        if(xFracC==3)
        {
            int cStride = toffset;
            if(yFracC==1)
            {
                predSampleLXC = (ad[-1+cStride]* -2 + ad[0+cStride] * 58 + ad[1+cStride] * 10 + ad[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==2)
            {
                predSampleLXC = (ad[-1+cStride]* -4 + ad[0+cStride] * 54 + ad[1+cStride] * 16 + ad[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==3)
            {
                predSampleLXC =	(ad[-1+cStride] * -6 + ad[0+cStride] * 46 + ad[1+cStride] * 28 + ad[2+cStride] * -4)>> shift2;
            }
            if(yFracC==4)
            {
                predSampleLXC =	(ad[-1+cStride] * -4 + ad[0+cStride] * 36 + ad[1+cStride] * 36 + ad[2+cStride] * -4)>> shift2;
            }
            if(yFracC==5)
            {
                predSampleLXC =	(ad[-1+cStride] * -4 + ad[0+cStride] * 28 + ad[1+cStride] * 46 + ad[2+cStride] * -6)>> shift2;
            }
            if(yFracC==6)
            {
                predSampleLXC =	(ad[-1+cStride] * -2 + ad[0+cStride] * 16 + ad[1+cStride] * 54 + ad[2+cStride] * -4)>> shift2;
            }
            if(yFracC==7)
            {
                predSampleLXC =	(ad[-1+cStride] * -2 + ad[0+cStride] * 10 + ad[1+cStride] * 58 + ad[2+cStride] * -2)>> shift2;
            }
        }
        if(xFracC==4)
        {
            int cStride = toffset;
            if(yFracC==1)
            {
                predSampleLXC = (ae[-1+cStride]* -2 + ae[0+cStride] * 58 + ae[1+cStride] * 10 + ae[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==2)
            {
                predSampleLXC = (ae[-1+cStride]* -4 + ae[0+cStride] * 54 + ae[1+cStride] * 16 + ae[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==3)
            {
                predSampleLXC =	(ae[-1+cStride] * -6 + ae[0+cStride] * 46 + ae[1+cStride] * 28 + ae[2+cStride] * -4)>> shift2;
            }
            if(yFracC==4)
            {
                predSampleLXC =	(ae[-1+cStride] * -4 + ae[0+cStride] * 36 + ae[1+cStride] * 36 + ae[2+cStride] * -4)>> shift2;
            }
            if(yFracC==5)
            {
                predSampleLXC =	(ae[-1+cStride] * -4 + ae[0+cStride] * 28 + ae[1+cStride] * 46 + ae[2+cStride] * -6)>> shift2;
            }
            if(yFracC==6)
            {
                predSampleLXC =	(ae[-1+cStride] * -2 + ae[0+cStride] * 16 + ae[1+cStride] * 54 + ae[2+cStride] * -4)>> shift2;
            }
            if(yFracC==7)
            {
                predSampleLXC =	(ae[-1+cStride] * -2 + ae[0+cStride] * 10 + ae[1+cStride] * 58 + ae[2+cStride] * -2)>> shift2;
            }
        }
        if(xFracC==5)
        {
            int cStride = toffset;
            if(yFracC==1)
            {
                predSampleLXC = (af[-1+cStride]* -2 + af[0+cStride] * 58 + af[1+cStride] * 10 + af[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==2)
            {
                predSampleLXC = (af[-1+cStride]* -4 + af[0+cStride] * 54 + af[1+cStride] * 16 + af[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==3)
            {
                predSampleLXC =	(af[-1+cStride] * -6 + af[0+cStride] * 46 + af[1+cStride] * 28 + af[2+cStride] * -4)>> shift2;
            }
            if(yFracC==4)
            {
                predSampleLXC =	(af[-1+cStride] * -4 + af[0+cStride] * 36 + af[1+cStride] * 36 + af[2+cStride] * -4)>> shift2;
            }
            if(yFracC==5)
            {
                predSampleLXC =	(af[-1+cStride] * -4 + af[0+cStride] * 28 + af[1+cStride] * 46 + af[2+cStride] * -6)>> shift2;
            }
            if(yFracC==6)
            {
                predSampleLXC =	(af[-1+cStride] * -2 + af[0+cStride] * 16 + af[1+cStride] * 54 + af[2+cStride] * -4)>> shift2;
            }
            if(yFracC==7)
            {
                predSampleLXC =	(af[-1+cStride] * -2 + af[0+cStride] * 10 + af[1+cStride] * 58 + af[2+cStride] * -2)>> shift2;
            }		}
        if(xFracC==6)
        {
            int cStride = toffset;
            if(yFracC==1)
            {
                predSampleLXC = (ag[-1+cStride]* -2 + ag[0+cStride] * 58 + ag[1+cStride] * 10 + ag[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==2)
            {
                predSampleLXC = (ag[-1+cStride]* -4 + ag[0+cStride] * 54 + ag[1+cStride] * 16 + ag[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==3)
            {
                predSampleLXC =	(ag[-1+cStride] * -6 + ag[0+cStride] * 46 + ag[1+cStride] * 28 + ag[2+cStride] * -4)>> shift2;
            }
            if(yFracC==4)
            {
                predSampleLXC =	(ag[-1+cStride] * -4 + ag[0+cStride] * 36 + ag[1+cStride] * 36 + ag[2+cStride] * -4)>> shift2;
            }
            if(yFracC==5)
            {
                predSampleLXC =	(ag[-1+cStride] * -4 + ag[0+cStride] * 28 + ag[1+cStride] * 46 + ag[2+cStride] * -6)>> shift2;
            }
            if(yFracC==6)
            {
                predSampleLXC =	(ag[-1+cStride] * -2 + ag[0+cStride] * 16 + ag[1+cStride] * 54 + ag[2+cStride] * -4)>> shift2;
            }
            if(yFracC==7)
            {
                predSampleLXC =	(ag[-1+cStride] * -2 + ag[0+cStride] * 10 + ag[1+cStride] * 58 + ag[2+cStride] * -2)>> shift2;
            }
        }
        if(xFracC==7)
        {
            int cStride = toffset;
            if(yFracC==1)
            {
                predSampleLXC = (ah[-1+cStride]* -2 + ah[0+cStride] * 58 + ah[1+cStride] * 10 + ah[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==2)
            {
                predSampleLXC = (ah[-1+cStride]* -4 + ah[0+cStride] * 54 + ah[1+cStride] * 16 + ah[2+cStride] * -2)	>> shift2;
            }
            if(yFracC==3)
            {
                predSampleLXC =	(ah[-1+cStride] * -6 + ah[0+cStride] * 46 + ah[1+cStride] * 28 + ah[2+cStride] * -4)>> shift2;
            }
            if(yFracC==4)
            {
                predSampleLXC =	(ah[-1+cStride] * -4 + ah[0+cStride] * 36 + ah[1+cStride] * 36 + ah[2+cStride] * -4)>> shift2;
            }
            if(yFracC==5)
            {
                predSampleLXC =	(ah[-1+cStride] * -4 + ah[0+cStride] * 28 + ah[1+cStride] * 46 + ah[2+cStride] * -6)>> shift2;
            }
            if(yFracC==6)
            {
                predSampleLXC =	(ah[-1+cStride] * -2 + ah[0+cStride] * 16 + ah[1+cStride] * 54 + ah[2+cStride] * -4)>> shift2;
            }
            if(yFracC==7)
            {
                predSampleLXC =	(ah[-1+cStride] * -2 + ah[0+cStride] * 10 + ah[1+cStride] * 58 + ah[2+cStride] * -2)>> shift2;
            }
        }
    }
    return predSampleLXC;
}


/* brief Apply fractional Sample Interpolate Filter

xCb     x component of a luma location ( xCb, yCb ) specifying the top-left sample of the current luma coding block relative to the top-left luma sample of the current picture
yCb		y component of a luma location ( xCb, yCb ) specifying the top-left sample of the current luma coding block relative to the top-left luma sample of the current picture
xBl		x component of a luma location ( xBl, yBl ) specifying the top-left sample of the current luma prediction block relative to the top-left sample of the current luma coding block
yBl		y component of a luma location ( xBl, yBl ) specifying the top-left sample of the current luma prediction block relative to the top-left sample of the current luma coding block
nPbW	width of the luma prediction block
nPbH    height of the luma prediction block
mvLX	luma motion vector given in quarter-luma-sample units
mvCLX   chroma motion vector given in eighth-chroma-sample units
refPicLXL   the selected reference picture sample array (Luma)
refPicLXCb  the selected reference picture sample array (Chroma)
refPicLXCr  the selected reference picture sample array (Chroma)
predSampleLXL  an (nPbW)x(nPbH) array of prediction luma sample values
predSampleLXCb  an (nPbW / 2)x(nPbH / 2) array of prediction chroma sample values.
predSampleLXCr  an (nPbW / 2)x(nPbH / 2) array of prediction chroma sample values.
 */
template <typename Sample, class H>
void fractSampleFilter(
        H &h,
        int xCb, int yCb,
        int xBl, int yBl,
        int nPbW, int nPbH,
        int mvLX[2], int mvCLX[2],
        Raster<Sample> refPicLXL,
        Raster<Sample> refPicLXCb,
        Raster<Sample> refPicLXCr,
        Raster<int16_t> predSamplesLXL,
        Raster<int16_t> predSamplesLXCb,
        Raster<int16_t> predSamplesLXCr,
        int bitDepthY,
        int bitDepthC)
{
    int xPb = xCb + xBl;
    int yPb = yCb + yBl;

    // Luma samples
    int xBaseL = xPb + (mvLX[0] >> 2);
    int yBaseL = yPb + (mvLX[1] >> 2);

    if (xBaseL + nPbW + 4 < 0)
    {
        xBaseL = -nPbW - 4;
    }
    else if (xBaseL - 3 > h[pic_width_in_luma_samples()] - 1)
    {
        xBaseL = h[pic_width_in_luma_samples()] + 2;
    }

    if (yBaseL + nPbH + 4 < 0)
    {
        yBaseL = -nPbH - 4;
    }
    else if (yBaseL - 3 > h[pic_height_in_luma_samples()] - 1)
    {
        yBaseL = h[pic_height_in_luma_samples()] + 2;
    }

    const int xFracL = mvLX[0] & 0x3;
    const int yFracL = mvLX[1] & 0x3;

    for(int xL=0; xL<nPbW; xL++)
    {
        for(int yL=0; yL<nPbH; yL++)
        {
            const int xIntL = xBaseL + xL;
            const int yIntL = yBaseL + yL;
            predSamplesLXL(xL, yL) = LumaSampleInterpolate<Sample>(h, xIntL, yIntL, xFracL, yFracL, refPicLXL, bitDepthY);
        }
    }

    // Chroma samples
    const int xBaseC = (xPb/h[SubWidthC()]) + (mvCLX[0] >> (1 + h[SubWidthC()]));
    const int yBaseC = (yPb/ h[SubHeightC()]) + (mvCLX[1] >> (1 + h[SubHeightC()]));
    const int xFracC = mvLX[0] & 0x7;
    const int yFracC = mvLX[1] & 0x7;

    for(int xC = 0; xC<=nPbW / h[SubWidthC()] -1; xC++)
    {
        for(int yC = 0; yC<=nPbH / h[SubHeightC()] -1; yC++)
        {
            const int xIntC = xBaseC + xC;
            const int yIntC = yBaseC + yC;
            predSamplesLXCb(xC, yC) = ChromaSampleInterpolate<Sample>(h, xIntC, yIntC, xFracC, yFracC, refPicLXCb, bitDepthC);
            predSamplesLXCr(xC, yC) = ChromaSampleInterpolate<Sample>(h, xIntC, yIntC, xFracC, yFracC, refPicLXCr, bitDepthC);
        }

    }
}

template <typename Sample, class H>
void WeightedSamplePred(H &h, const PuData &puData, int nPbW, int nPbH, Raster<int16_t> predSamplesL0, Raster<int16_t> predSamplesL1, Raster<Sample> predSamples, int cIdx)
{
    int bitDepth = (!cIdx) ? h[BitDepthY()] : h[BitDepthC()];
    int weightedPredFlag = (h[slice_type()] == P) ? h[weighted_pred_flag()] : h[weighted_bipred_flag()];

    if(weightedPredFlag==0) //	Default Weighted Sample Prediction is invoked
    {
        int shift1 = 14 - bitDepth;
        int shift2 = 15 - bitDepth;
        int offset1 = (shift1 > 0) ? (1 << (shift1-1)) : 0;
        int offset2 = 1 << (shift2-1);
        if (puData.predFlag(L0) == 1 && puData.predFlag(L1) == 0)
        {
            for(int x=0; x<nPbW/(cIdx?h[SubWidthC()]:1); x++)
            {
                for(int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                {
                    predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, (predSamplesL0(x,y) + offset1) >> shift1);
                }
            }
        }
        else if (puData.predFlag(L0) == 0 && puData.predFlag(L1) == 1)
        {
            for (int x=0; x<nPbW/(cIdx? h[SubWidthC()] :1); x++)
            {
                for (int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                {
                    predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, (predSamplesL1(x,y) + offset1) >> shift1);
                }
            }
        }
        else //(predFlagL0==1 && predFlagL1==1)
        {
            for (int x=0; x<nPbW/(cIdx? h[SubWidthC()] :1); x++)
            {
                for (int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                {
                    predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, (predSamplesL0(x,y) + predSamplesL1(x,y) + offset2) >> shift2);
                }
            }
        }
    }
    else //Explicit Weighted Sample Prediction is invoked
    {
        int shift1 = 14 - bitDepth;
        int log2Wd, o0, o1, w0, w1;
        if(cIdx==0)
        {
            log2Wd = h[luma_log2_weight_denom()] + shift1;
            if (puData.predFlag(L0)) w0 = h[LumaWeightL0(puData.refIdx(L0))];
            if (puData.predFlag(L1)) w1 = h[LumaWeightL1(puData.refIdx(L1))];
            if (puData.predFlag(L0)) o0 = h[luma_offset_l0(puData.refIdx(L0))] * (1 << (bitDepth - 8));
            if (puData.predFlag(L1)) o1 = h[luma_offset_l1(puData.refIdx(L1))] * (1 << (bitDepth - 8));
        }
        else
        {
            log2Wd = h[ChromaLog2WeightDenom()] + shift1;
            if (puData.predFlag(L0)) w0 = h[ChromaWeightL0(puData.refIdx(L0), cIdx - 1)];
            if (puData.predFlag(L1)) w1 = h[ChromaWeightL1(puData.refIdx(L1), cIdx - 1)];
            if (puData.predFlag(L0)) o0 = h[ChromaOffsetL0(puData.refIdx(L0), cIdx - 1)] * (1 << (bitDepth - 8));
            if (puData.predFlag(L1)) o1 = h[ChromaOffsetL1(puData.refIdx(L1), cIdx - 1)] * (1 << (bitDepth - 8));
        }
        if (puData.predFlag(L0) == 1 && puData.predFlag(L1) == 0)
        {
            if(log2Wd >= 1)
            {
                for(int x=0; x<nPbW/(cIdx? h[SubWidthC()] :1); x++)
                {
                    for(int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                    {
                        predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, ((predSamplesL0(x,y) * w0 + (1 << (log2Wd-1))) >> log2Wd) + o0);
                    }
                }
            }
            else
            {
                for(int x=0; x<nPbW/(cIdx? h[SubWidthC()] :1); x++)
                {
                    for(int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                    {
                        predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, predSamplesL0(x,y) * w0 + o0);
                    }
                }

            }
        }
        else if (puData.predFlag(L0) == 0 && puData.predFlag(L1) == 1)
        {
            if(log2Wd >= 1)
            {
                for(int x=0; x<nPbW/(cIdx? h[SubHeightC()] :1); x++)
                {
                    for(int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                    {
                        predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, ((predSamplesL1(x,y) * w1 + (1 << (log2Wd-1))) >> log2Wd) + o1);
                    }
                }
            }
            else
            {
                for(int x=0; x<nPbW/(cIdx? h[SubWidthC()] :1); x++)
                {
                    for(int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                    {
                        predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, predSamplesL1(x,y) * w1 + o1);
                    }
                }

            }
        }
        else
        {
            assert(puData.predFlag(L0) == 1 && puData.predFlag(L1) == 1);

            for(int x=0; x<nPbW/(cIdx?h[SubWidthC()]:1); x++)
            {
                for(int y=0; y<nPbH/(cIdx? h[SubHeightC()] :1); y++)
                {
                    predSamples(x,y) = Clip3(0, (1 << bitDepth) - 1, (predSamplesL0(x,y) * w0 + predSamplesL1(x,y) * w1 + ((o0 + o1 + 1) << log2Wd)) >> (log2Wd + 1));
                }
            }
        }
    }
}

static int clipMvLumaComponent(int component, int nPbSize, int pictureSize)
{
    if (component + nPbSize + 4 < 0)
        return -nPbSize - 4;
    else if (component > pictureSize + 2)
        return pictureSize + 2;

    return component;
}

// copies a non-overlapping w*h rectangle of pixels from src to dst
template <class T>
void copyBlock(T *dst, std::intptr_t strideDst, const T *src, std::intptr_t strideSrc, size_t w, size_t h)
{
    while (h--)
    {
        memcpy(dst, src, w*sizeof(T));
        dst += strideDst;
        src += strideSrc;
    }
}


template <class T>
bool compareBlock(T *dst, std::intptr_t strideDst, const T *src, std::intptr_t strideSrc, size_t w, size_t h)
{
    bool different = false;
    while (h--)
    {
        different = different || !!memcmp(dst, src, w*sizeof(T));
        dst += strideDst;
        src += strideSrc;
    }
    return different;
}


template <class T, class U>
void copyBlock(Raster<T> dst, Raster<U> src, size_t w, size_t h)
{
    copyBlock(&dst(0, 0), dst.stride, &src(0, 0), src.stride, w, h);
}


template <class T, class U>
bool compareBlock(Raster<T> dst, Raster<U> src, size_t w, size_t h)
{
    return compareBlock(&dst(0, 0), dst.stride, &src(0, 0), src.stride, w, h);
}


template <class Sample>
static void predictUni(HavocTablePredUni <Sample> &table, Raster<Sample> pred0, Raster<Sample> pred1, Raster<Sample> pred2, Picture<Sample> &reference, MotionVector mv, int xPb, int yPb, int nPbW, int nPbH, int bitDepthY, int bitDepthC)
{
    auto const xFracL = mv[0] & 0x3;
    auto const yFracL = mv[1] & 0x3;

    auto const xBaseL = clipMvLumaComponent(xPb + (mv[0] >> 2), nPbW, reference[0].width);
    auto const yBaseL = clipMvLumaComponent(yPb + (mv[1] >> 2), nPbH, reference[0].height);

    Raster<Sample> ref0(reference[0], xBaseL, yBaseL);

    auto const index = (nPbW >> 2) - 1;

    {
        auto *f = *havocGetPredUni(&table, 8, nPbW, nPbH, xFracL, yFracL, bitDepthY);
        f(pred0.p, pred0.stride, ref0.p, ref0.stride, nPbW, nPbH, xFracL, yFracL, bitDepthY);
    }

    MotionVector mvCLX{
        MotionVector::ComponentType(mv[0] * 2 / reference.subWidthC),
                MotionVector::ComponentType(mv[1] * 2 / reference.subHeightC) };

    auto const xFracC = mvCLX[0] & 0x7;
    auto const yFracC = mvCLX[1] & 0x7;

    auto const xBaseC = xBaseL >> (reference.subWidthC - 1);
    auto const yBaseC = yBaseL >> (reference.subHeightC - 1);

    Raster<Sample> ref1(reference[1], xBaseC, yBaseC);
    Raster<Sample> ref2(reference[2], xBaseC, yBaseC);

    nPbW >>= (reference.subWidthC - 1);
    nPbH >>= (reference.subHeightC - 1);

    {
        auto *f = *havocGetPredUni(&table, 4, nPbW, nPbH, xFracC, yFracC, bitDepthC);
        f(pred1.p, pred1.stride, ref1.p, ref1.stride, nPbW, nPbH, xFracC, yFracC, bitDepthC);
        f(pred2.p, pred2.stride, ref2.p, ref2.stride, nPbW, nPbH, xFracC, yFracC, bitDepthC);
    }
}

template <class Sample>
static void predictBi(HavocTablePredBi<Sample> &table, Raster<Sample> pred0, Raster<Sample> pred1, Raster<Sample> pred2, Picture<Sample> &referenceA, MotionVector mvA, Picture<Sample> &referenceB, MotionVector mvB, int xPb, int yPb, int nPbW, int nPbH, int bitDepthY, int bitDepthC)
{
    const int xFracLA = mvA[0] & 0x3;
    const int yFracLA = mvA[1] & 0x3;
    const int xFracLB = mvB[0] & 0x3;
    const int yFracLB = mvB[1] & 0x3;

    const int xBaseLA = clipMvLumaComponent(xPb + (mvA[0] >> 2), nPbW, referenceA[0].width);
    const int yBaseLA = clipMvLumaComponent(yPb + (mvA[1] >> 2), nPbH, referenceA[0].height);
    const int xBaseLB = clipMvLumaComponent(xPb + (mvB[0] >> 2), nPbW, referenceB[0].width);
    const int yBaseLB = clipMvLumaComponent(yPb + (mvB[1] >> 2), nPbH, referenceB[0].height);

    Raster<Sample> refA0(referenceA[0], xBaseLA, yBaseLA);
    Raster<Sample> refB0(referenceB[0], xBaseLB, yBaseLB);

    auto f = *havocGetPredBi(&table, 8, nPbW, nPbH, xFracLA, yFracLA, xFracLB, yFracLB, bitDepthY);
    f(pred0.p, pred0.stride, refA0.p, refB0.p, refA0.stride, nPbW, nPbH, xFracLA, yFracLA, xFracLB, yFracLB, bitDepthY);

    MotionVector mvCLA{
        MotionVector::ComponentType(mvA[0] * 2 / referenceA.subWidthC),
                MotionVector::ComponentType(mvA[1] * 2 / referenceA.subHeightC) };

    MotionVector mvCLB{
        MotionVector::ComponentType(mvB[0] * 2 / referenceB.subWidthC),
                MotionVector::ComponentType(mvB[1] * 2 / referenceB.subHeightC) };

    auto const xFracCA = mvCLA[0] & 0x7;
    auto const yFracCA = mvCLA[1] & 0x7;
    auto const xFracCB = mvCLB[0] & 0x7;
    auto const yFracCB = mvCLB[1] & 0x7;

    auto const xBaseCA = xBaseLA >> (referenceA.subWidthC - 1);
    auto const yBaseCA = yBaseLA >> (referenceA.subHeightC - 1);
    auto const xBaseCB = xBaseLB >> (referenceB.subWidthC - 1);
    auto const yBaseCB = yBaseLB >> (referenceB.subHeightC - 1);

    Raster<Sample> refA1(referenceA[1], xBaseCA, yBaseCA);
    Raster<Sample> refA2(referenceA[2], xBaseCA, yBaseCA);
    Raster<Sample> refB1(referenceB[1], xBaseCB, yBaseCB);
    Raster<Sample> refB2(referenceB[2], xBaseCB, yBaseCB);

    nPbW >>= referenceA.subWidthC - 1;
    nPbH >>= referenceA.subHeightC - 1;

    auto fC = *havocGetPredBi(&table, 4, nPbW, nPbH, xFracCA, yFracCA, xFracCB, yFracCB, bitDepthC);
    fC(pred1.p, pred1.stride, refA1.p, refB1.p, refA1.stride, nPbW, nPbH, xFracCA, yFracCA, xFracCB, yFracCB, bitDepthC);
    fC(pred2.p, pred2.stride, refA2.p, refB2.p, refA2.stride, nPbW, nPbH, xFracCA, yFracCA, xFracCB, yFracCB, bitDepthC);
}


template <class H, class PredPicture>
void predictInter(PredPicture &predPicture, prediction_unit pu, H &h)
{
    typedef typename PredPicture::Type Sample;

    Snake<BlockData>::Cursor *cursor = h;
    const PuData &puData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);

    const auto xPb = pu.x0;
    const auto yPb = pu.y0;

    Picture<Sample> *referencePictureLX[2] = { 0, 0 };

    if (puData.getDpbIndex(L0) >= 0)
    {
        StateReconstructedPicture<Sample> *ref = static_cast<StateReconstructedPicture<Sample> *>(h[RefPicList(L0)][puData.refIdx(L0)].dp->reconstructedPicture.get());
        referencePictureLX[0] = ref->picture.get();
    }

    if (puData.getDpbIndex(L1) >= 0)
    {
        StateReconstructedPicture<Sample> *ref = static_cast<StateReconstructedPicture<Sample> *>(h[RefPicList(L1)][puData.refIdx(L1)].dp->reconstructedPicture.get());
        referencePictureLX[1] = ref->picture.get();
    }

    assert(referencePictureLX[0] || referencePictureLX[1]);

    Raster<Sample> pred0 = predPicture(xPb, yPb, 0);
    Raster<Sample> pred1 = predPicture(xPb, yPb, 1);
    Raster<Sample> pred2 = predPicture(xPb, yPb, 2);

    if (!h[weightedPredFlag()])
    {
        if (referencePictureLX[0] && !referencePictureLX[1])
        {
            HavocTablePredUni<Sample> *table = h;
            predictUni<Sample>(*table, pred0, pred1, pred2, *referencePictureLX[0], puData.mv(L0), xPb, yPb, pu.nPbW, pu.nPbH, h[BitDepthY()], h[BitDepthC()]);
        }
        else if (!referencePictureLX[0] && referencePictureLX[1])
        {
            HavocTablePredUni<Sample> *table = h;
            predictUni<Sample>(*table, pred0, pred1, pred2, *referencePictureLX[1], puData.mv(L1), xPb, yPb, pu.nPbW, pu.nPbH, h[BitDepthY()], h[BitDepthC()]);
        }
        else
        {
            HavocTablePredBi<Sample> *table = h;
            predictBi<Sample>(*table, pred0, pred1, pred2, *referencePictureLX[0], puData.mv(L0), *referencePictureLX[1], puData.mv(L1), xPb, yPb, pu.nPbW, pu.nPbH, h[BitDepthY()], h[BitDepthC()]);
        }
        return;
    }

    int16_t predSamplesLX[2][3][64*64];

    for (int refList = 0; refList < 2; ++refList)
    {
        auto x = refList;

        Picture<Sample> *reference = referencePictureLX[x];
        if (!reference) continue;

        const coding_quadtree *cqt = h;

        int mvLX[2] = { puData.mv(refList)[0], puData.mv(refList)[1] };

        auto &refPicLXL = (*reference)[0];
        auto &refPicLXCb = (*reference)[1];
        auto &refPicLXCr = (*reference)[2];

        Raster<int16_t> predSamplesLXL(predSamplesLX[x][0], 64);
        Raster<int16_t> predSamplesLXCb(predSamplesLX[x][1], 64);
        Raster<int16_t> predSamplesLXCr(predSamplesLX[x][2], 64);

        fractSampleFilter<Sample>(
                h,
                cqt->x0, cqt->y0,
                pu.x0-cqt->x0, pu.y0-cqt->y0,
                pu.nPbW, pu.nPbH,
                mvLX, mvLX,
                refPicLXL, refPicLXCb, refPicLXCr,
                predSamplesLXL, predSamplesLXCb, predSamplesLXCr,
                h[BitDepthY()], h[BitDepthC()]);
    }

    for (int cIdx=0; cIdx<3; ++cIdx)
    {
        Raster<int16_t> predSamplesL0(predSamplesLX[0][cIdx], 64);
        Raster<int16_t> predSamplesL1(predSamplesLX[1][cIdx], 64);
        Raster<Sample> predSamples = predPicture(xPb, yPb, cIdx);
        WeightedSamplePred<Sample>(h, puData, pu.nPbW, pu.nPbH, predSamplesL0, predSamplesL1, predSamples, cIdx);
    }
}

#endif
