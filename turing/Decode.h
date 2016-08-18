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

// This file contains template specialisations for decoding of bitstreams to a YUV reconstructed output.

#ifndef INCLUDED_Decode_h
#define INCLUDED_Decode_h

#pragma once

#include "Global.h"
#include "Handlers.h"
#include "Picture.h"
#include "Mvp.h"
#include "havoc/residual_decode.h"
#include "havoc/quantize.h"


template <typename F>
struct Decode :
    Read<F>
    {
    };


template <> struct Decode<void> { };


template <> struct Decode<slice_segment_data>
{
    template <class H> static void go(slice_segment_data e, H &h)
    {
        if (h[BitDepthY()] > 8 || h[BitDepthC()] > 8)
            inner<uint16_t>(e, h);
        else
            inner<uint8_t>(e, h);
    }

    template <class Sample, class H> static void inner(slice_segment_data e, H &h)
    {
        StatePictures *statePictures = h;
        StatePicture *statePicture = h;

        if (h[first_slice_segment_in_pic_flag()])
        {
            auto *p = new StateReconstructedPicture<Sample>;
            const int pad = 96;// h[CtbSizeY()] + 16; // review: less padding will suffice
            p->picture.reset(new Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32));
            statePicture->reconstructedPicture.reset(p);
        }

        auto *p = static_cast<StateReconstructedPicture<Sample> *>(statePicture->reconstructedPicture.get());
        auto hPicture = h.extend(p);
        Read<slice_segment_data>::go(e, hPicture);
    }
};


template <> struct Decode<JustDecoded>
{
    template <class H> static void go(const JustDecoded &jd, H &h)
    {
        StatePicture *statePicture = h;
        StateReconstructedPictureBase *base = statePicture->reconstructedPicture.get();

        if (h[BitDepthY()] > 8 || h[BitDepthC()] > 8)
        {
            auto *derived = dynamic_cast<StateReconstructedPicture<uint16_t> *>(base);
            finishPicture<uint16_t>(h.extend(derived));
        }
        else
        {
            auto *derived = dynamic_cast<StateReconstructedPicture<uint8_t> *>(base);
            finishPicture<uint8_t>(h.extend(derived));
        }
    }
};


template <> struct Decode<OutputPicture>
{
    template <class H> static void go(const OutputPicture &op, H &h)
    {
        StatePictures *statePictures = h;
        StateDecode *stateDecode = h;

        auto &pic = *statePictures->getPicByPoc(op.poc);

        if (h[BitDepthY()] > 8 || h[BitDepthC()] > 8)
        {
            auto &dp = dynamic_cast<StateReconstructedPicture<uint16_t> &>(*pic.reconstructedPicture);
            stateDecode->deliver(dp, h[chroma_format_idc()], h[BitDepthY()], h[BitDepthC()]);
        }
        else
        {
            auto &dp = dynamic_cast<StateReconstructedPicture<uint8_t> &>(*pic.reconstructedPicture);
            stateDecode->deliver(dp, h[chroma_format_idc()], h[BitDepthY()], h[BitDepthC()]);
        }
    }
};


 template <>
struct Decode<coding_tree_unit>
{
    template <class H> static void go(const coding_tree_unit &ctu, H &h)
    {
        Read<coding_tree_unit>::go(ctu, h);
        StatePicture *statePicture = h;
        statePicture->loopFilterPicture->processCtu(h, ctu);
    }
};


template <>
struct Decode<coding_unit>
{
    template <class H> static void go(const coding_unit &cu, H &h)
    {
        Read<coding_unit>::go(cu, h);
        StatePicture *statePicture = h;
        statePicture->loopFilterPicture->processCu(h, cu);
    }
};


template <>
struct Decode<UnexpectedData> : Null<UnexpectedData> { };


template <>
struct Decode<sei_rbsp> : Null<sei_rbsp> { };


template <>
struct Decode<Process<prediction_unit>>
{
    template <class H> static void go(Process<prediction_unit>, H &h)
    {
        prediction_unit *pu = h;

        StatePicture *statePicture = h;
        statePicture->loopFilterPicture->processPu(h, *pu);

        using Sample = typename SampleType<H>::Type;
        StateReconstructedPicture<Sample> *stateReconstructedPicture = h;
        predictInter(*stateReconstructedPicture->picture, *pu, h);
    }
};



template <class V>
struct DecodePcmSample
{
    static bool constexpr isChroma = std::is_same<V, pcm_sample_chroma>::value;

    template <class H> static void go(Element<V, uv> fun, H &h)
    {
        Read<Element<V, uv>>::go(fun, h);

        coding_quadtree const *cqt = h;
        auto const nCbS = 1 << cqt->log2CbSize;
        auto const width = nCbS / (isChroma ? h[SubWidthC()] : 1);
        auto const height = nCbS / (isChroma ? h[SubHeightC()] : 1);

        auto const i = fun.v.i % width;
        auto const row = fun.v.i / width;
        auto const j = row % height;
        auto const cIdx = isChroma ? (1 + row / height) : 0;

        auto const shift = isChroma
            ? h[BitDepthC()] - h[PcmBitDepthC()]
            : h[BitDepthY()] - h[PcmBitDepthY()];

        auto const sample = h[fun.v];
        h[ReconstructedSamples(cqt->x0, cqt->y0, cIdx)](i, j) = sample << shift;
    }
};


template <> struct Decode<Element<pcm_sample_chroma, uv>> : DecodePcmSample<pcm_sample_chroma> {};
template <> struct Decode<Element<pcm_sample_luma, uv>> : DecodePcmSample<pcm_sample_luma> {};


template <>
struct Decode<transform_unit>
{
    template <class H> static void go(const transform_unit &tu, H &h)
    {
        StatePicture *statePictureBase = h;
        statePictureBase->loopFilterPicture->processTu(h, tu);

        Syntax<transform_unit>::go(tu, h);
    }
};


template <class V>
struct Decode<IfCbf<V, residual_coding>>
{
    template <class H> static void go(const IfCbf<V, residual_coding> &i, H &h)
    {
        const residual_coding &rc = i.f;

        if (h[current(CuPredMode(rc.x0, rc.y0))] == MODE_INTRA)
        {
            auto reconstructedSamples = h[ReconstructedSamples(rc.x0, rc.y0, rc.cIdx)];
            predictBlockIntra(h, reconstructedSamples, reconstructedSamples, rc);
        }

        Syntax<IfCbf<V, residual_coding>>::go(i, h);
    }
};


template <>
struct Decode<residual_coding>
{
    template <class H> static void go(const residual_coding &rc, H &h)
    {
        StatePicture *statePicture = h;
        statePicture->loopFilterPicture->processRc(h, rc);

        Read<residual_coding>::go(rc, h);

        const int nCbS = 1 << rc.log2TrafoSize;

    using Sample = typename SampleType<H>::Type;
    StateReconstructedPicture<Sample> *stateReconstructedPicture = h;
    auto recSamples = (*stateReconstructedPicture->picture)(rc.x0, rc.y0, rc.cIdx);
        auto predSamples = recSamples;

        auto const bitDepth = rc.cIdx ? h[BitDepthC()] : h[BitDepthY()];

        if (h[cu_transquant_bypass_flag()])
        {
            for (int y = 0; y<nCbS; ++y)
            {
                for (int x = 0; x<nCbS; ++x)
                {
                    recSamples(x, y) = static_cast<Sample>(clipCidx1<>(predSamples(x, y) + h[TransCoeffLevel(rc.x0, rc.y0, rc.cIdx, x, y)], bitDepth));
                }
            }
            return;
        }

        const int qpScaled = static_cast<QpState *>(h)->getQp(rc.cIdx);

        ALIGN(32, int16_t, coeffsQ[32 * 32]);
        ALIGN(32, int16_t, coeffs[32 * 32]);
        ALIGN(32, int16_t, res[32 * 32]);

        for (int yC = 0; yC<(1 << rc.log2TrafoSize); ++yC)
            for (int xC = 0; xC<(1 << rc.log2TrafoSize); ++xC)
                coeffsQ[(yC << rc.log2TrafoSize) + xC] = h[TransCoeffLevel(rc.x0, rc.y0, rc.cIdx, xC, yC)];

        Raster<int16_t> quantizedCoefficients(coeffsQ, 1ull << rc.log2TrafoSize);
        Raster<int16_t> coefficients(coeffs, 1ull << rc.log2TrafoSize);
        Raster<int16_t> resSamples(res, 1ull << rc.log2TrafoSize);

        const ScalingMatrices::Type matrix[32] = { 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 };
        Raster<const ScalingMatrices::Type> m(matrix, 0);

        Snake<BlockData>::Cursor *cursor = h;
        StateSpatial *stateSpatial = h;
        const coding_quadtree *cqt = h;
        const coding_unit cu(cqt->x0, cqt->y0, cqt->log2CbSize);

        if (h[scaling_list_enabled_flag()])
        {
            const int sizeId = rc.log2TrafoSize - 2;
            m = Raster<ScalingMatrices::Type>(static_cast<ScalingMatrices *>(h)->getMatrix(sizeId, ScalingMatrices::matrixId(sizeId, rc.cIdx, h[current(CuPredMode(rc.x0, rc.y0))])), 1ull << rc.log2TrafoSize);
            inverseQuantize(coefficients, quantizedCoefficients, m, rc.log2TrafoSize, qpScaled, bitDepth);
        }
        else
        {
            const int scale = static_cast<QpState *>(h)->getScale(rc.cIdx);
            const int shift = rc.log2TrafoSize - 1 + bitDepth - 8;
            const int n = 1 << 2 * rc.log2TrafoSize;

            havoc_quantize_inverse *inverseQuantize = *havoc_get_quantize_inverse(h, scale, shift);
            inverseQuantize(coefficients.p, quantizedCoefficients.p, scale, shift, n);
        }

        const int trType = (
                h[current(CuPredMode(rc.x0, rc.y0))] == MODE_INTRA &&
                nCbS == 4 &&
                rc.cIdx == 0) ? 1 : 0;

        if (h[transform_skip_flag()])
        {
            const int bdShift = 20 - bitDepth;
            for (int y = 0; y<nCbS; ++y)
            {
                for (int x = 0; x<nCbS; ++x)
                {
                    int r = coefficients(x, y) << 7;
                    r = (r + (1 << (bdShift - 1))) >> bdShift;
                    resSamples(x, y) = r;
                    recSamples(x, y) = static_cast<Sample>(clipCidx1<>(predSamples(x, y) + resSamples(x, y), bitDepth));
                }
            }
        }
        else
        {
            auto *inverseTransformAdd = *havoc_get_inverse_transform_add<Sample>(h, trType, rc.log2TrafoSize);
            inverseTransformAdd(recSamples.p, recSamples.stride, predSamples.p, predSamples.stride, coefficients.p, bitDepth);
        }
    }
};

#endif
