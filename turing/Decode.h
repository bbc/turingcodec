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
#include "IntraReferenceSamples.h"
#include "havoc/transform.h"
#include "havoc/quantize.h"


template <typename F>
struct Decode :
    Read<F>
{
};


template <> struct Decode<void> { };


template <> struct Decode<UnexpectedData> : Null<UnexpectedData> {};


template <> struct Decode<sei_rbsp> : Null<sei_rbsp> {};


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
        StatePicture *statePicture = h;

        auto *p = static_cast<StateReconstructedPicture<Sample> *>(statePicture->reconstructedPicture.get());
        auto hPicture = h.extend(p);
        Read<slice_segment_data>::go(e, hPicture);
    }
};



template <> struct Decode<PictureBegin>
{
    template <class H> static void go(PictureBegin pb, H &h)
    {
        Syntax<PictureBegin>::go(pb, h);
        if (h[BitDepthY()] > 8 || h[BitDepthC()] > 8)
            inner<uint16_t>(h);
        else
            inner<uint8_t>(h);
    }

    template <class Sample, class H> static void inner(H &h)
    {
        StatePicture *statePicture = h;

        auto *p = new StateReconstructedPicture<Sample>;

        const int pad = 96;// h[CtbSizeY()] + 16; // review: less padding will suffice
        p->picture.reset(new Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32));

        statePicture->reconstructedPicture.reset(p);

        setupStateReconstructedPicture(h[Concrete<StatePicture>()], h);
    }
};


// Invoked to deblock the current decoded picture
struct PictureDeblock {};


// Invoked to perform SAO filtering upon the current decoded picture
struct PictureSao {};


template <> struct Decode<PictureDeblock>
{
    template <class H> static void go(PictureDeblock, H &h)
    {
        if (!h[pps_deblocking_filter_disabled_flag()])
        {
            StatePicture *statePicture = h;
            Profiler::Timers *timers = h;

            auto *reconstructedPicture = &h[Concrete<StateReconstructedPictureBase>()];
            auto &picture = *reconstructedPicture->picture;

            auto &recPictureL = picture[0];
            auto &recPictureCb = picture[1];
            auto &recPictureCr = picture[2];

            timers->deblock.start();
            statePicture->loopFilterPicture->template deblock<EDGE_VER>(h, recPictureL, recPictureCb, recPictureCr, 0, 0, h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()]);
            statePicture->loopFilterPicture->template deblock<EDGE_HOR>(h, recPictureL, recPictureCb, recPictureCr, 0, 0, h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()]);
            timers->deblock.stop();
        }
    }
};


template <> struct Decode<PictureSao>
{
    template <class H> static void go(PictureSao, H &h)
    {
        using Sample = typename SampleType<H>::Type;

        if (h[sample_adaptive_offset_enabled_flag()])
        {
            StatePicture *statePicture = h;

            auto *reconstructedPicture = &h[Concrete<StateReconstructedPictureBase>()];
            auto &picture = *reconstructedPicture->picture;

            const int pad = 96;// h[CtbSizeY()] + 16; // less padding will suffice
            Picture<Sample> *pictureSao = new Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);

            statePicture->loopFilterPicture->applySao2<Sample>(h, *pictureSao, picture);

            reconstructedPicture->picture.reset(pictureSao);
        }

        StateReconstructedPicture<Sample> *reconstructedPicture = h;
        Padding::padPicture(*reconstructedPicture->picture);
    }
};


template <> struct Decode<PictureDone>
{
    template <class H> static void go(PictureDone pd, H &h)
    {
        if (h[BitDepthY()] > 8 || h[BitDepthC()] > 8)
            inner<uint16_t>(h);
        else
            inner<uint8_t>(h);

        Syntax<PictureDone>::go(pd, h);
    }

    template <class Sample, class H> static void inner(H &h)
    {
        StatePicture *statePicture = h;
        StateReconstructedPictureBase *base = statePicture->reconstructedPicture.get();

        auto *derived = static_cast<StateReconstructedPicture<Sample> *>(base);
        auto hPicture = h.extend(derived);

        hPicture(PictureDeblock());
        hPicture(PictureSao());
    }
};


template <typename Sample, class H>
void pictureOutput(PictureOutput &po, H &h)
{
    StatePictures *statePictures = h;
    StateDecode *stateDecode = h;

    Handler<typename H::Tag, StatePicture> hPicture;
    hPicture.state = po.statePicture.get();

    auto &dp = static_cast<StateReconstructedPicture<Sample> &>(*hPicture.state->reconstructedPicture);

    ThreePlanes<Sample> conformanceWindow = { *dp.picture,
        hPicture[SubWidthC()] * hPicture[conf_win_left_offset()],
        hPicture[SubHeightC()] * hPicture[conf_win_top_offset()],
        hPicture[SubWidthC()] * hPicture[conf_win_right_offset()],
        hPicture[SubHeightC()] * hPicture[conf_win_bottom_offset()] };


    // Conformance window parameters come from SPS of the decoded picture, not the current picture which may be in a subsequent CVS
    // (see RAP_B_Bossen_2)
    // Review: other parameters (below) may also be affected similarly: consider active parameter sets within StatePicture...
    int chromaFormat = hPicture[chroma_format_idc()];
    int bitDepthY = hPicture[BitDepthY()];
    int bitDepthC = hPicture[BitDepthC()];

    if (stateDecode->ofs)
    {
        // if requested and necessary, round picture to 8-bit precision
        if (stateDecode->vm["8-bit"].as<bool>() && std::is_same<Sample, uint16_t>::value)
        {
            std::shared_ptr<Picture<uint8_t>> picture8(new Picture<uint8_t>(
                conformanceWindow[0].width,
                conformanceWindow[0].height,
                chromaFormat, 0, 0, 32));

            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                auto const bitDepth = cIdx ? bitDepthC : bitDepthY;
                auto const add = 1 << bitDepth >> 8;
                auto const shift = bitDepth - 8;
                for (int y = 0; y < (*picture8)[cIdx].height; ++y)
                {
                    for (int x = 0; x < (*picture8)[cIdx].width; ++x)
                    {
                        int value = conformanceWindow[cIdx](x, y);
                        value += add;
                        value >>= shift;
                        if (value > 255)
                            value = 255;
                        (*picture8)[cIdx](x, y) = value;
                    }
                }
            }

            stateDecode->ofs << *picture8;
        }
        else
            stateDecode->ofs << conformanceWindow;
    }

    // review: remove decoder's md5 option (it's really an external functionality)
    if (stateDecode->vm.count("md5"))
    {
        std::ostringstream oss;
        oss << conformanceWindow;
        std::string const s = oss.str();
        md5_append(&stateDecode->md5Sum, &reinterpret_cast<md5_byte_t const&>(s.front()), static_cast<int>(s.size()));
    }

    stateDecode->progressReporter.progress(++stateDecode->n);
}


template <> struct Decode<PictureOutput>
{
    template <class H> static void go(PictureOutput &po, H &h)
    {
        Syntax<PictureOutput>::go(po, h);

        Handler<Decode<void>, StatePicture> hPicture;
        hPicture.state = po.statePicture.get();

        if (hPicture[BitDepthY()] > 8 || hPicture[BitDepthC()] > 8)
            pictureOutput<uint16_t>(po, h);
        else
            pictureOutput<uint8_t>(po, h);
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


template <> struct Decode<transform_unit>
{
    template <class H> static void go(const transform_unit &tu, H &h)
    {
        StatePicture *statePicture = h;
        statePicture->loopFilterPicture->processTu(h, tu);

        Read<transform_unit>::go(tu, h);
    }
};


struct IntraPrediction
{
    residual_coding rc;
};


template <> struct Decode<IntraPrediction>
{
    template <class H> static void go(IntraPrediction f, H &h)
    {
        auto const &rc = f.rc;

        auto samples = h[ReconstructedSamples(rc.x0, rc.y0, rc.cIdx)];
        using Sample = typename SampleType<H>::Type;
        IntraReferenceSamples<Sample> unfiltered;
        IntraReferenceSamples<Sample> filtered;

        unfiltered.substitute(h, samples, rc);

        const int predModeIntra = rc.cIdx ? h[IntraPredModeC(rc.x0, rc.y0)] : h[IntraPredModeY(rc.x0, rc.y0)];

        auto const ff = filterFlag(rc.cIdx, predModeIntra, 1 << rc.log2TrafoSize);

        if (ff)
            filtered.filter(unfiltered, h[strong_intra_smoothing_enabled_flag()], h[BitDepthY()], 1 << rc.log2TrafoSize);

        auto const &pF = ff ? filtered : unfiltered;

        auto const bitDepth = rc.cIdx ? h[BitDepthC()] : h[BitDepthY()];

        havoc::intra::Table<Sample> *table = h;
        auto *function = table->lookup(rc.cIdx, bitDepth, rc.log2TrafoSize, predModeIntra);

        (*function)(samples.p, samples.stride, &pF(0, -1), predModeIntra);
    }
};


template <class V> struct Decode<IfCbf<V, residual_coding>>
{
    template <class H> static void go(const IfCbf<V, residual_coding> &i, H &h)
    {
        const residual_coding &rc = i.f;

        if (h[current(CuPredMode(rc.x0, rc.y0))] == MODE_INTRA)
            h(IntraPrediction{ rc });

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

        QpState *qpState = h;

        Read<residual_coding>::go(rc, h);

        const int nCbS = 1 << rc.log2TrafoSize;

        using Sample = typename SampleType<H>::Type;
        StateReconstructedPicture<Sample> *stateReconstructedPicture = h;
        auto recSamples = (*stateReconstructedPicture->picture)(rc.x0, rc.y0, rc.cIdx);
        auto predSamples = recSamples;

        auto const bitDepth = rc.cIdx ? h[BitDepthC()] : h[BitDepthY()];

        if (h[cu_transquant_bypass_flag()])
        {
            for (int y = 0; y < nCbS; ++y)
                for (int x = 0; x < nCbS; ++x)
                    recSamples(x, y) = static_cast<Sample>(clipCidx1<>(predSamples(x, y) + h[TransCoeffLevel(rc.x0, rc.y0, rc.cIdx, x, y)], bitDepth));

            return;
        }

        ALIGN(32, int16_t, coeffsQ[32 * 32]);
        ALIGN(32, int16_t, coeffs[32 * 32]);
        ALIGN(32, int16_t, res[32 * 32]);

        for (int yC = 0; yC < (1 << rc.log2TrafoSize); ++yC)
            for (int xC = 0; xC < (1 << rc.log2TrafoSize); ++xC)
                coeffsQ[(yC << rc.log2TrafoSize) + xC] = h[TransCoeffLevel(rc.x0, rc.y0, rc.cIdx, xC, yC)];

        Raster<int16_t> quantizedCoefficients(coeffsQ, 1ull << rc.log2TrafoSize);
        Raster<int16_t> coefficients(coeffs, 1ull << rc.log2TrafoSize);
        Raster<int16_t> resSamples(res, 1ull << rc.log2TrafoSize);

        ScalingMatrices::Type const matrix[32] = { 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 };
        Raster<ScalingMatrices::Type const> m(matrix, 0);

        Snake<BlockData>::Cursor *cursor = h;
        StateSpatial *stateSpatial = h;
        const coding_quadtree *cqt = h;
        const coding_unit cu(cqt->x0, cqt->y0, cqt->log2CbSize);

        if (h[scaling_list_enabled_flag()])
        {
            ScalingMatrices *scalingMatrices = h;
            auto const sizeId = rc.log2TrafoSize - 2;
            auto const matrixId = ScalingMatrices::matrixId(sizeId, rc.cIdx, h[current(CuPredMode(rc.x0, rc.y0))]);
            m = { scalingMatrices->getMatrix(sizeId, matrixId), 1ull << rc.log2TrafoSize };
            inverseQuantize(coefficients, quantizedCoefficients, m, rc.log2TrafoSize, qpState->getQp(rc.cIdx), bitDepth);
        }
        else
        {
            const int scale = qpState->getScale(rc.cIdx);
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
            int const bdShift = 20 - bitDepth;
            for (int y = 0; y < nCbS; ++y)
                for (int x = 0; x < nCbS; ++x)
                {
                    int r = coefficients(x, y) << 7;
                    r = (r + (1 << (bdShift - 1))) >> bdShift;
                    resSamples(x, y) = r;
                    recSamples(x, y) = static_cast<Sample>(clipCidx1<>(predSamples(x, y) + resSamples(x, y), bitDepth));
                }
        }
        else
        {
            auto *inverseTransformAdd = *havoc::get_inverse_transform_add<Sample>(h, trType, rc.log2TrafoSize);
            inverseTransformAdd(recSamples.p, recSamples.stride, predSamples.p, predSamples.stride, coefficients.p, bitDepth);
        }
    }
};

#endif
