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

// Rate and distortion measurement and estimation

#ifndef INCLUDED_Measure_h
#define INCLUDED_Measure_h

#pragma once

#include "Global.h"
#include "Dsp.h"
#include "Cost.h"
#include "Syntax.h"
#include "ContextModel.h"
#include "MotionVector.h"
#include "havoc/hadamard.h"
#ifdef WIN32
#include <intrin.h>
#endif

template <class F=void> struct Measure;
template <class> struct Reconstruct;


template <class H>
Lambda getReciprocalLambda(H &h)
{
    return static_cast<StateEncodePicture *>(h)->reciprocalLambda;
}


template <class H>
double getReciprocalSqrtLambda(H &h)
{
    return static_cast<StateEncodePicture *>(h)->reciprocalSqrtLambda;
}


template <class H>
double computeLambda(H &h)
{
    StateEncode *stateEncode = h;

    double qp = h[SliceQpY()];

    double qpFactor = static_cast<StateEncodePicture *>(h)->qpFactor;

    if (h[slice_type()] == I)
    {
        const double dLambda_scale = 1.0 - Clip3(0.0, 0.5, 0.05*(stateEncode->gopM - 1.0));
        qpFactor = 0.57*dLambda_scale;
    }
    double lambda = qpFactor * pow(2.0, (qp - 12.0) / 3.0 );

    if (h[PicOrderCntVal()] % stateEncode->gopM)
    {
        lambda *= Clip3(2.00, 4.00, ((qp - 12.0) / 6.0));
    }

    return lambda;
}


struct Result
{
    int32_t distortion[3];
    Cost rate;
    Cost value;
};

struct SumOfSquaredErrors;
struct Rectangle;

template <class H, class Metric, class Region>
struct Compute;


template <typename Sample>
int32_t measureSatd(havoc_table_hadamard_satd<Sample> *table, Raster<Sample> p1, Raster<Sample> p2, int width, int height)
{
    int32_t satd = 0;
    if ((width | height) & 0x3)
    {
        auto &f = *havoc_get_hadamard_satd<Sample>(table, 1);
        assert(!((width | height) & 0x1));
        for (int y = 0; y<height; y += 2)
        {
            for (int x = 0; x<width; x += 2)
            {
                satd += f(&p1(x, y), p1.stride, &p2(x, y), p2.stride);
            }
        }
    }
    else if ((width | height) & 0x7)
    {
        auto &f = *havoc_get_hadamard_satd<Sample>(table, 2);
        for (int y = 0; y<height; y += 4)
        {
            for (int x = 0; x<width; x += 4)
            {
                satd += f(&p1(x, y), p1.stride, &p2(x, y), p2.stride);
            }
        }
    }
    else
    {
        auto &f = *havoc_get_hadamard_satd<Sample>(table, 3);
        for (int y = 0; y<height; y += 8)
        {
            for (int x = 0; x<width; x += 8)
            {
                satd += f(&p1(x, y), p1.stride, &p2(x, y), p2.stride);
            }
        }
    }
    return satd;
}


template <class H>
struct Compute<H, struct SumOfAbsoluteTransformedDifferences, Turing::Rectangle>
{
    static int32_t go(H &h, int x0, int y0, int width, int height, int cIdx)
    {
        using Sample = typename SampleType<H>::Type;

        StateEncodePicture *stateEncodePicture = h;

        PictureWrapper &pictureWrapper = *stateEncodePicture->docket->picture;
        auto &pictureWrap = static_cast<PictureWrap<Sample> &>(pictureWrapper);
        auto &picture = static_cast<Picture<Sample> &>(pictureWrap);

        StateReconstructedPicture<Sample> *stateReconstructedPicture = h;
        Picture<Sample> &reconstructed = *stateReconstructedPicture->picture;

        if (cIdx)
        {
            width >>= 1;
            height >>= 1;
            if ((width | height) & 0x3) return 0;
        }

        havoc_table_hadamard_satd<Sample> *table = h;
        return measureSatd(table, picture(x0, y0, cIdx), reconstructed(x0, y0, cIdx), width, height);
    }
};


template <class H, class Metric>
struct Compute<H, Metric, prediction_unit>
{
    static int32_t go(H &h, prediction_unit pu, int cIdx)
    {
        return Compute<H, Metric, Turing::Rectangle>::go(h, pu.x0, pu.y0, pu.nPbW, pu.nPbH, cIdx);
    }
};


template <bool haveLzCnt>
static inline unsigned estimateRateOfMvdComponent(int d)
{
    unsigned u = abs(d);
    unsigned rate;

    if (haveLzCnt)
    {
#ifdef WIN32
        rate = 32 - __lzcnt(u);
#else
        rate = 32 - __builtin_clz(u);
#endif
    }
    else
    {
        rate = 0;
        for (int i = 0; i < 32; ++i)
            if (u & (1 << i))
                rate = i + 1;
    }

    constexpr bool selfCheck = false;
    if (selfCheck)
    {
        int rateCheck = 0;
        for (int i = 0; i < 32; ++i)
            if (u & (1 << i))
                rateCheck = i + 1;
        assert(rate == rateCheck);
    }

    return rate;
}


// estimate MVD rate for comparison of motion vector candidate: less accurate around zero and has constant offset
template <bool haveLzCnt>
static inline Cost rateOf(MotionVector const &mvd)
{
    return Cost::make(
            estimateRateOfMvdComponent<haveLzCnt>(mvd[0]) +
            estimateRateOfMvdComponent<haveLzCnt>(mvd[1]) + 1, -1);
}


template <class Metric, class H>
struct PuPredictionError;


template <class Metric, class H>
struct MvdCost
{
    static Cost measure(H &h, int refList, MotionVector mvd)
    {
        StateCodedData *stateCodedData = h;
        coding_quadtree const *cqt = h;
        prediction_unit const *pu = h;

        stateCodedData->codedPu.mvd(refList) = mvd;

        Snake<BlockData>::Cursor *cursor = h;
        BlockData &blockData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);
        computePuData(h, *pu, blockData);

        return PuPredictionError<Metric, H>::measure(h, *pu);
    }
};

#endif
