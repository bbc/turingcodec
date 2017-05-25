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

#include "SyntaxCtu.hpp"
#include "Reconstruct.h"
#include "Handlers.h"
#include "EstimateRate.h"
#include "Write.h"
#include "Syntax.h"
#include "Search.h"
#include "havoc/quantize.h"
#include "Rdoq.h"
#include <cstdlib> // temporary for rand()


// Base verb class for reconstruction / intra prediction actions
// Forces CBFs to 1, actual, encoded CBFs are ignored while this is in force.
// Review: is "Reconstruct" still the best name?
template <typename F> struct Reconstruct
    :
    Syntax<F>,
    ValueConst<rqt_root_cbf, 1>, // review: ValueConst<> only checked on Reconstruct<void> specialization
    ValueConst<cbf_luma, 1>,
    ValueConst<cbf_cb, 1>,
    ValueConst<cbf_cr, 1>
{
};


// do not write or measure any symbols when reconstructing
template <class V> struct Reconstruct<Element<V, ae>> : Null<Element<V, ae>> {};
template <class V> struct Reconstruct<Element<V, f>> : Null<Element<V, f>> {};
template <> struct Reconstruct<pcm_sample> : Null<pcm_sample> {};


// performs inter prediction, codes residual and computes SSD
template <typename F> struct ReconstructInter : Reconstruct<F> { };


template <>
struct ReconstructInter<transform_tree>
{
    template <class H> static void go(const transform_tree &tt, H &h)
    {
        using Sample = typename SampleType<H>::Type;

        StateEncode *stateEncode = h;
        StateCodedData *stateCodedData = h;
        coding_quadtree const *cqt = h;

        if (tt.trafoDepth && tt.blkIdx == 0)
        {
            assert(stateCodedData->transformTree.p == stateCodedData->transformTreeAncestry[tt.trafoDepth - 1].p);
            assert(!stateCodedData->transformTreeAncestry[tt.trafoDepth - 1].word0().cbfWord);
            // Advance to next Transform Tree in CodedData (have just split)
            ++stateCodedData->transformTree.p;
        }


        stateCodedData->transformTreeAncestry[tt.trafoDepth] = stateCodedData->transformTree;
        bool inferTransformDepth = !(tt.log2TrafoSize <= h[MaxTbLog2SizeY()] &&
            tt.log2TrafoSize > h[MinTbLog2SizeY()] &&
            tt.trafoDepth < h[MaxTrafoDepth()] && !(h[IntraSplitFlag()] && (tt.trafoDepth == 0)));

        Candidate<Sample> *candidate = h;
        stateCodedData->transformTree.init(cqt->log2CbSize - tt.log2TrafoSize, tt.blkIdx);
        stateCodedData->transformTreeAncestry[tt.trafoDepth] = stateCodedData->transformTree;
        h[split_transform_flag()] = inferTransformDepth ? infer(split_transform_flag(tt.x0, tt.y0, tt.trafoDepth), h) : ((tt.trafoDepth < candidate->rqtdepth) && !candidate->noresidual);
        stateCodedData->transformTree.word0().split_transform_flag = h[split_transform_flag()];

        Syntax<transform_tree>::go(tt, h);

        {
            auto &wordParent = stateCodedData->transformTreeAncestry[tt.trafoDepth - 1].word0().cbfWord;
            auto const &word = stateCodedData->transformTreeAncestry[tt.trafoDepth].word0().cbfWord;
            wordParent |= word;
        }
    }
};


template <>
struct ReconstructInter<transform_unit>
{
    template <class H> static void go(const transform_unit &tu, H &h)
    {
        StateCodedData *stateCodedData = h;
        coding_quadtree const *cqt = h;
        transform_tree const *tt = h;

        stateCodedData->transformTree.word0().split_transform_flag = 0;
        stateCodedData->residual = stateCodedData->transformTree.firstResidual();
        Syntax<transform_unit>::go(tu, h);

        // after each MODE_INTER transform_unit, update the coded data pointers
        stateCodedData->transformTree.p = stateCodedData->residual.p;
        stateCodedData->transformTreeChroma.p = stateCodedData->residual.p;
    }
};


template <>
struct Reconstruct<transform_tree>
{
    template <class H> static void go(const transform_tree &tt, H &h)
    {
        h[split_transform_flag()] = infer(split_transform_flag(tt.x0, tt.y0, tt.trafoDepth), h);
        Syntax<transform_tree>::go(tt, h);
    }
};




// performs luma intra prediction, codes residual and computes SSD
template <typename F> struct ReconstructIntraLuma : Reconstruct<F> { };

// performs chroma intra prediction, codes residual and computes SSD
template <typename F> struct ReconstructIntraChroma : Reconstruct<F> { };

template <>
struct ReconstructIntraChroma<transform_tree>
{
    template <class H>
    static void go(const transform_tree &tt, H &h)
    {
        bool const isIntraPartition = tt.trafoDepth == h[part_mode()];

        if (isIntraPartition)
        {
            Neighbourhood *neighbourhood = h;
            Snake<BlockData>::Cursor *cursor = h;
            cursor->relocate(neighbourhood->snake, tt, h[MinCbLog2SizeY()] - 1, true);
            cursor->value.intra.predModeY = static_cast<StateCodedData *>(h)->codedCu.IntraPredModeY(tt.blkIdx);
        }

        bool const split = (tt.trafoDepth == 0 && h[part_mode()]) || tt.log2TrafoSize > h[MaxTbLog2SizeY()];
        h[split_transform_flag()] = split;
        Syntax<transform_tree>::go(tt, h);
    }
};


int quantize_c(short *dst, const short *src, int quantiserScale, int quantiserShift, int quantiserOffset, int size)
{
    int cbf = 0;
    for (int coeffIdx = 0; coeffIdx < size; coeffIdx++)
    {
        int sign = src[coeffIdx] < 0 ? -1 : 1;
        int coeffAbs = abs(src[coeffIdx]) * quantiserScale;
        int level = (coeffAbs + quantiserOffset) >> quantiserShift;
        cbf |= level;
        dst[coeffIdx] = sign * level;
    }

    return cbf;
}


// review: has very similar code to PredictIntraLumaBlock
struct ReconstructIntraBlock
{
    template <class H> static void go(residual_coding const &rc, H h)
    {
        using Sample = typename SampleType<H>::Type;

        StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
        StateReconstructionCache<Sample> *stateReconstructionCache = h;
        StateCodedData *stateCodedData = h;
        NeighbourhoodEnc<Sample> *neighbourhood = h;
        Candidate<Sample> *candidate = h;
        transform_tree const *tt = h;
        coding_quadtree const *cqt = h;
        StateEncode *stateEncode = h;
        *static_cast<residual_coding*>(h) = rc;

        assert(h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA);
        assert(!h[cu_transquant_bypass_flag()]);

        int const nTbS = 1 << rc.log2TrafoSize;

        assert(candidate->zz[rc.cIdx] == (zPositionOf(rc) >> (rc.cIdx ? 2 : 0)));
        const auto z = candidate->zz[rc.cIdx];
        typename ReconstructionCache<Sample>::Piece piece = stateReconstructionCache->components[rc.cIdx].allocateBlock(rc.log2TrafoSize, z);
#ifdef DEBUG_PIECES
        piece.rc = rc;
#endif
        auto recSamples = stateReconstructionCache->components[rc.cIdx].get(piece);
        assert((int)recSamples.stride == 1 << rc.log2TrafoSize);
        auto predSamples = recSamples;

        candidate->appendPiece(rc.cIdx, rc, rc.log2TrafoSize, piece.i);

        const int rcOffsetX = (rc.x0 - cqt->x0) / (rc.cIdx ? 2 : 1);
        const int rcOffsetY = (rc.y0 - cqt->y0) / (rc.cIdx ? 2 : 1);
        const int stride = (1 << cqt->log2CbSize) / (rc.cIdx ? 2 : 1);

        int const bitDepth = rc.cIdx ? h[BitDepthC()] : h[BitDepthY()];

        if (rc.cIdx && h[IntraSplitFlag()])
        {
            Snake<BlockData>::Cursor *cursor = h;
            cursor->relocate(neighbourhood->snake, *cqt, h[MinCbLog2SizeY()] - 1, true);
            cursor->value.intra.predModeY = static_cast<StateCodedData *>(h)->codedCu.IntraPredModeY(0);
        }

        // intra prediction
        {
            int constexpr bitDepth = 2 * sizeof(Sample) + 6;

            auto &unfiltered = stateEncodeSubstream->unfiltered[rc.cIdx];
            auto &filtered = stateEncodeSubstream->filtered;
            const int predModeIntra = rc.cIdx ? h[IntraPredModeC(rc.x0, rc.y0)] : h[IntraPredModeY(rc.x0, rc.y0)];
            auto const ff = filterFlag(rc.cIdx, predModeIntra, nTbS);

            if (tt->trafoDepth && !h[part_mode()])
            {
                auto samples = neighbourhood->snakeIntraReferenceSamples[rc.cIdx].offset(rc.x0, rc.y0, rc.cIdx ? 1 : 0);
                unfiltered.substituteFast(h, samples, rc);

                if (ff)
                    filtered.filter(unfiltered, h[strong_intra_smoothing_enabled_flag()], bitDepth, nTbS);
            }

            auto const &pF = ff ? filtered : unfiltered;

            havoc::intra::Table<Sample> *table = h;
            havoc::intra::Function<Sample> *f = table->lookup(rc.cIdx, bitDepth, rc.log2TrafoSize, predModeIntra);
            f(&predSamples(0, 0), predSamples.stride, &pF(0, -1), predModeIntra);
        }

        // input picture
        auto &pictureInput = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
        auto sourceSamples = pictureInput(rc.x0, rc.y0, rc.cIdx);

        ALIGN(32, int16_t, resSamplesRawBuffer[32 * 32]);
        Raster<int16_t> resSamplesRaw(resSamplesRawBuffer, nTbS);

        // subtract prediction from input
        // review: SIMD optimisations, perhaps integrate into forward transform
        for (int y = 0; y < nTbS; ++y)
            for (int x = 0; x < nTbS; ++x)
                resSamplesRaw(x, y) = sourceSamples(x, y) - predSamples(x, y);

        ALIGN(32, int16_t, coefficientsBuffer[32 * 32]);
        Raster<int16_t> coefficients(coefficientsBuffer, nTbS, 0, 0);

        int const trType = nTbS == 4 && rc.cIdx == 0;
        bool const checkTSkip = nTbS == 4 && h[transform_skip_enabled_flag()];

        {
            // forward transform
            int constexpr bitDepth = 2 * sizeof(Sample) + 6;
            havoc::table_transform<bitDepth> *table = h;
            auto *transform = *havoc::get_transform<bitDepth>(table, trType, rc.log2TrafoSize);
            transform(coefficients.p, resSamplesRaw.p, resSamplesRaw.stride);
        }

        ALIGN(32, int16_t, quantizedCoefficientsBuffer[32 * 32]);
        Raster<int16_t> quantizedCoefficients(quantizedCoefficientsBuffer, nTbS, 0, 0);

        bool cbf;

        {
            // review: precompute most of this
            const int n = 1 << 2 * rc.log2TrafoSize;
            QpState *qpState = h;
            auto const &entry = qpState->lookup(rc.cIdx);
            auto const shiftQuantise = entry.quantizeShift - rc.log2TrafoSize;

            // quantization
            if (stateEncode->rdoq)
            {
                StateEncodePicture *stateEncodePicture = h;
                double lambda = stateEncodePicture->lambda;
                if(stateEncode->useRateControl)
                {
                    StateEncodePicture *stateEncodePicture = h;
                    int segmentPoc = stateEncodePicture->docket->segmentPoc;

                    stateEncode->rateControlParams->takeTokenOnRCLevel();
                    lambda = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                    stateEncode->rateControlParams->releaseTokenOnRCLevel();
                }

                Contexts *contexts = h;
                Rdoq rdoqEngine(lambda, contexts, entry.quantiseScale, entry.scale, rc.log2TrafoSize, bitDepth);
                cbf = !!rdoqEngine.runQuantisation(quantizedCoefficients.p, coefficients.p, entry.quantiseScale, shiftQuantise,
                    n, rc, h[scanIdx()], true, !!h[sign_data_hiding_enabled_flag()]);
            }
            else
            {
                havoc_quantize *quantize = *havoc_get_quantize(static_cast<havoc_table_quantize *>(h));
                cbf = !!quantize(quantizedCoefficients.p, coefficients.p, entry.quantiseScale, shiftQuantise, entry.offsetQuantiseShifted, n);
            }

            // inverse quantization
            auto const  shift = rc.log2TrafoSize - 1 + bitDepth - 8;
            havoc_quantize_inverse *inverseQuantize = *havoc_get_quantize_inverse(h, entry.scale, shift);
            inverseQuantize(coefficients.p, quantizedCoefficients.p, entry.scale, shift, n);
        }

        // review: would this be better as Sample instead of int16_t?
        ALIGN(32, int16_t, backupPredictionBuffer[32 * 32]);
        Raster<int16_t> backupPrediction(backupPredictionBuffer, 1 << rc.log2TrafoSize);

        ALIGN(32, int16_t, backupQuantCoeffsBuffer[32 * 32]);
        Raster<int16_t> backupQuantCoeffs(backupQuantCoeffsBuffer, 1 << rc.log2TrafoSize);

        // backup prediction and quantized coefficients in case we need to test transform skip
        if (checkTSkip && cbf)
        {
            for (int y = 0; y < nTbS; y++)
            {
                for (int x = 0; x < nTbS; x++)
                {
                    backupPrediction(x, y) = predSamples(x, y);
                    backupQuantCoeffs(x, y) = quantizedCoefficients(x, y);
                }
            }
        }

        {
            havoc::table_inverse_transform_add<Sample> *table = h;
            auto *inverseTransformAdd = *havoc::get_inverse_transform_add(table, trType, rc.log2TrafoSize);

            // inverse transform and add to predicted samples
            inverseTransformAdd(recSamples.p, recSamples.stride, predSamples.p, predSamples.stride, coefficients.p, bitDepth);
        }

        int32_t backupSSDNoTSkip; // save SSD without tskip
        {
            // compute distortion as SSD
            auto ssdFunction = *havoc_get_ssd<Sample>(h, rc.log2TrafoSize);
            StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
            backupSSDNoTSkip = stateEncodeSubstream->ssd[rc.cIdx] + ssdFunction(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, nTbS, nTbS);
            if (!(checkTSkip && cbf))
                stateEncodeSubstream->ssd[rc.cIdx] = backupSSDNoTSkip;
        }

        // backup CodedData pointers and words containing CBF flags before storing residual coded data (clean)
        auto backupStateClean = *stateCodedData;
        auto const backupCuWord0Clean = stateCodedData->codedCu.word0();
        auto const backupTuWord0Clean = stateCodedData->transformTree.word0();
        ContextsAndCost *contextsAndCost = h;
        auto backupContextsAndCostClean = *contextsAndCost;

        {
            // store residual information to CodedData without tskip
            StateCodedData *stateCodedData = h;
            auto &tu = rc.cIdx ? stateCodedData->transformTreeChroma : stateCodedData->transformTree;
            if (rc.cIdx != 2)
            {
                tu.init(cqt->log2CbSize - rc.log2TrafoSize - (rc.cIdx ? 1 : 0), tt->blkIdx);
                stateCodedData->residual = tu.firstResidual();
            }

            stateCodedData->residual.transformSkipFlag() = 0;

            assert(quantizedCoefficients.stride == nTbS);
            CodedData::storeResidual(stateCodedData->codedCu, stateCodedData->residual, quantizedCoefficients.p, rc.log2TrafoSize, h[scanIdx()], cbf, tu, rc.cIdx);
            if (rc.cIdx != 1)
                tu.p = stateCodedData->residual.p;
            stateCodedData->codedDataAfter = stateCodedData->residual.p;
        }

        if (checkTSkip && cbf)
        {
            // review: would this be better as Sample instead of int16_t?
            ALIGN(32, int16_t, recSamplesNoTSkipBuffer[32 * 32]);
            Raster<int16_t> recSamplesNoTSkip(recSamplesNoTSkipBuffer, 1 << rc.log2TrafoSize);

            // Backup reconstructed samples without tskip
            for (int y = 0; y < nTbS; y++)
            {
                for (int x = 0; x < nTbS; x++)
                {
                    recSamplesNoTSkip(x, y) = recSamples(x, y);
                }
            }

            // Return pointers to clean to estimate the bitrate without tskip
            *stateCodedData = backupStateClean;
            auto e = h.template change<EstimateRate<void>>();

            switch (rc.cIdx)
            {
                case 0:
                    e(IfCbf<cbf_luma, residual_coding>{cbf_luma(rc.x0, rc.y0, tt->trafoDepth), rc});
                    break;
                case 1:
                    e(IfCbf<cbf_cb, residual_coding>{cbf_cb(rc.x0, rc.y0, tt->trafoDepth), rc});
                    break;
                case 2:
                    e(IfCbf<cbf_cr, residual_coding>{cbf_cr(rc.x0, rc.y0, tt->trafoDepth), rc});
                    break;
                default: assert(0);
            }

            auto contextsCostNoTSkip = *static_cast<ContextsAndCost *>(candidate);

            // Return pointers to clean
            *contextsAndCost = backupContextsAndCostClean;
            *stateCodedData = backupStateClean;
            stateCodedData->codedCu.word0() = backupCuWord0Clean;
            stateCodedData->transformTree.word0() = backupTuWord0Clean;

            {
                // test transform skip
                {
                    const int bdShift = 13 - bitDepth;
                    for (int y = 0; y < nTbS; y++)
                        for (int x = 0; x < nTbS; x++)
                            coefficients(x, y) = resSamplesRaw(x, y) << bdShift;
                }

                {
                    const int n = 1 << 2 * rc.log2TrafoSize;
                    const int qpScaled = static_cast<QpState *>(h)->getQp(rc.cIdx);
                    const int scaleQuantise = static_cast<QpState *>(h)->getQuantiseScale(rc.cIdx);
                    const int shiftQuantise = 29 - bitDepth + qpScaled / 6 - rc.log2TrafoSize;
                    const int offsetQuantise = (h[slice_type()] == I ? 171 : 85) << (shiftQuantise - 9);
                    const int scale = static_cast<QpState *>(h)->getScale(rc.cIdx);
                    const int shift = rc.log2TrafoSize - 1 + bitDepth - 8;

                    auto const quantize = *havoc_get_quantize(h);

                    // quantization
                    if (stateEncode->rdoq)
                    {
                        double lambda;
                        if(stateEncode->useRateControl)
                        {
                            StateEncodePicture *stateEncodePicture = h;
                            int segmentPoc = stateEncodePicture->docket->segmentPoc;

                            stateEncode->rateControlParams->takeTokenOnRCLevel();
                            lambda = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                            stateEncode->rateControlParams->releaseTokenOnRCLevel();
                        }
                        else
                        {
                            lambda = static_cast<StateEncodePicture *>(h)->lambda;
                        }
                        Contexts *contexts = h;
                        bool isSdhEnabled = !!h[sign_data_hiding_enabled_flag()];
                        Rdoq rdoqEngine(lambda, contexts, scaleQuantise, scale, rc.log2TrafoSize, bitDepth);
                        cbf = !!rdoqEngine.runQuantisation(quantizedCoefficients.p, coefficients.p, scaleQuantise, shiftQuantise,
                                                            n, rc, h[scanIdx()], true, isSdhEnabled);
                    }
                    else
                    {
                        cbf = !!quantize(quantizedCoefficients.p, coefficients.p, scaleQuantise, shiftQuantise, offsetQuantise >> (shiftQuantise - 16), n);
                    }

                    auto *inverseQuantize = *havoc_get_quantize_inverse(h, scale, shift);

                    // inverse quantization
                    inverseQuantize(coefficients.p, quantizedCoefficients.p, scale, shift, n);
                }

                {
                    // load prediction from backup
                    for (int y = 0; y < nTbS; y++)
                        for (int x = 0; x < nTbS; x++)
                            predSamples(x, y) = static_cast<Sample>(backupPrediction(x, y));

                    // inverse transform skip and add to predicted samples
                    const int bdShift = 20 - bitDepth;
                    for (int y = 0; y < nTbS; y++)
                    {
                        for (int x = 0; x < nTbS; x++)
                        {
                            int r = coefficients(x, y) << 7;
                            r = (r + (1 << (bdShift - 1))) >> bdShift;
                            recSamples(x, y) = static_cast<Sample>(clipCidx1<>(predSamples(x, y) + r, bitDepth));
                        }
                    }
                }

                int32_t backupSSDTSkip;
                {
                    // compute distortion as SSD
                    auto ssdFunction = *havoc_get_ssd<Sample>(h, rc.log2TrafoSize);
                    StateEncodeSubstreamBase *stateEncodeSubstreamBase = h;
                    backupSSDTSkip = stateEncodeSubstreamBase->ssd[rc.cIdx] + ssdFunction(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, nTbS, nTbS);
                    stateEncodeSubstreamBase->ssd[rc.cIdx] = backupSSDTSkip;
                }

                {
                    // store residual information to CodedData with tskip
                    StateCodedData *stateCodedData = h;
                    auto &tu = rc.cIdx ? stateCodedData->transformTreeChroma : stateCodedData->transformTree;
                    if (rc.cIdx != 2)
                    {
                        tu.init(cqt->log2CbSize - rc.log2TrafoSize - (rc.cIdx ? 1 : 0), tt->blkIdx);
                        stateCodedData->residual = tu.firstResidual();
                    }

                    stateCodedData->residual.transformSkipFlag() = 1;

                    assert(quantizedCoefficients.stride == nTbS);
                    CodedData::storeResidual(stateCodedData->codedCu, stateCodedData->residual, quantizedCoefficients.p, rc.log2TrafoSize, h[scanIdx()], cbf, tu, rc.cIdx);
                    if (rc.cIdx != 1)
                    {
                        tu.p = stateCodedData->residual.p;
                    }
                    stateCodedData->codedDataAfter = stateCodedData->residual.p;
                }

                // Return pointers to clean
                *stateCodedData = backupStateClean;

                switch (rc.cIdx)
                {
                    case 0:
                        e(IfCbf<cbf_luma, residual_coding>{cbf_luma(rc.x0, rc.y0, tt->trafoDepth), rc});
                        break;
                    case 1:
                        e(IfCbf<cbf_cb, residual_coding>{cbf_cb(rc.x0, rc.y0, tt->trafoDepth), rc});
                        break;
                    case 2:
                        e(IfCbf<cbf_cr, residual_coding>{cbf_cr(rc.x0, rc.y0, tt->trafoDepth), rc});
                        break;
                    default: assert(0);
                }

                auto contextsCostTSkip = *static_cast<ContextsAndCost *>(candidate);
                *contextsAndCost = backupContextsAndCostClean;

                // Compute lambdaDistortion with and without tskip
                Lambda reciprocalLambda = getReciprocalLambda(h);
                if(stateEncode->useRateControl)
                {
                    StateEncodePicture *stateEncodePicture = h;
                    int segmentPoc = stateEncodePicture->docket->segmentPoc;

                    stateEncode->rateControlParams->takeTokenOnRCLevel();
                    double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                    stateEncode->rateControlParams->releaseTokenOnRCLevel();
                    reciprocalLambda.set(value);
                }
                contextsCostNoTSkip.lambdaDistortion += backupSSDNoTSkip * reciprocalLambda;
                contextsCostTSkip.lambdaDistortion += backupSSDTSkip   * reciprocalLambda;

                if (contextsCostNoTSkip.cost2() < contextsCostTSkip.cost2())
                {
                    // tskip is not cheaper, restore reconstructed samples and quantized coefficients without tskip
                    for (int y = 0; y < nTbS; y++)
                    {
                        for (int x = 0; x < nTbS; x++)
                        {
                            recSamples(x, y) = static_cast<Sample>(recSamplesNoTSkip(x, y));
                            quantizedCoefficients(x, y) = backupQuantCoeffs(x, y);
                        }
                    }

                    // restore SSD
                    auto const ssdFunction = *havoc_get_ssd<Sample>(h, rc.log2TrafoSize);
                    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
                    stateEncodeSubstream->ssd[rc.cIdx] = backupSSDNoTSkip;

                    // return pointers to clean
                    *stateCodedData = backupStateClean;
                    stateCodedData->codedCu.word0() = backupCuWord0Clean;
                    stateCodedData->transformTree.word0() = backupTuWord0Clean;
                    {
                        // store residual information to CodedData without tskip
                        StateCodedData *stateCodedData = h;
                        auto &tu = rc.cIdx ? stateCodedData->transformTreeChroma : stateCodedData->transformTree;
                        if (rc.cIdx != 2)
                        {
                            tu.init(cqt->log2CbSize - rc.log2TrafoSize - (rc.cIdx ? 1 : 0), tt->blkIdx);
                            stateCodedData->residual = tu.firstResidual();
                        }

                        stateCodedData->residual.transformSkipFlag() = 0;

                        assert(quantizedCoefficients.stride == nTbS);
                        CodedData::storeResidual(stateCodedData->codedCu, stateCodedData->residual, quantizedCoefficients.p, rc.log2TrafoSize, h[scanIdx()], true/*cbf*/, tu, rc.cIdx);
                        if (rc.cIdx != 1)
                        {
                            tu.p = stateCodedData->residual.p;
                        }
                        stateCodedData->codedDataAfter = stateCodedData->residual.p;
                    }
                }
            }
        } // if (checkTSkip && cbf)

        {
            // copy bottom and right edges of reconstructed block to inter prediction snake
            // review: may not always need to do this (could defer until later in some cases?)
            Raster<Sample> recPicture = recSamples.offset(-rc.x0 >> (rc.cIdx ? 1 : 0), -rc.y0 >> (rc.cIdx ? 1 : 0));
            int const size = 1 << rc.log2TrafoSize << (rc.cIdx ? 1 : 0);
            typename Snake<Sample>::Pointer &snake = neighbourhood->snakeIntraReferenceSamples[rc.cIdx];
            snake.copyFrom2D(Turing::Rectangle{ rc.x0, rc.y0, size, size }, recPicture, rc.cIdx ? 1 : 0);
        }


    }
};

template <> struct ReconstructIntraLuma<IfCbf<cbf_luma, residual_coding>> : ReconstructIntraBlock {};
template <> struct ReconstructIntraLuma<IfCbf<cbf_cb, residual_coding>> : Null<IfCbf<cbf_cb, residual_coding>> {};
template <> struct ReconstructIntraLuma<IfCbf<cbf_cr, residual_coding>> : Null<IfCbf<cbf_cr, residual_coding>> {};

template <> struct ReconstructIntraChroma<IfCbf<cbf_luma, residual_coding>> : Null<IfCbf<cbf_luma, residual_coding>> {};
template <> struct ReconstructIntraChroma<IfCbf<cbf_cb, residual_coding>> : ReconstructIntraBlock {};
template <> struct ReconstructIntraChroma<IfCbf<cbf_cr, residual_coding>> : ReconstructIntraBlock {};


struct PredictIntraLumaBlock
{
    template <class H> static void go(residual_coding const &rc, H h)
    {
        using Sample = typename SampleType<H>::Type;
        int constexpr bitDepth = 2 * sizeof(Sample) + 6;

        StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
        StateEncodePicture *stateEncodePicture = h;
        NeighbourhoodEnc<Sample> *neighbourhood = h;
        coding_quadtree const *cqt = h;
        transform_tree const *tt = h;

        assert(rc.cIdx == 0);
        assert(rc.log2TrafoSize >= 2);
        assert(rc.log2TrafoSize <= 5);

        int const nTbS = 1 << rc.log2TrafoSize;
        int const log2Resolution = rc.cIdx ? 1 : 0;

        // temporary, stack buffer to receive predicted samples:
        ALIGN(32, Sample, buffer[32 * 32]);
        Raster<Sample> predSamples{ buffer, 32 }; // review: smaller stride better?

        // perform intra prediction
        {
            auto &unfiltered = stateEncodeSubstream->unfiltered[rc.cIdx];
            auto &filtered = stateEncodeSubstream->filtered;
            const int predModeIntra = rc.cIdx ? h[IntraPredModeC(rc.x0, rc.y0)] : h[IntraPredModeY(rc.x0, rc.y0)];
            auto const ff = filterFlag(rc.cIdx, predModeIntra, nTbS);

            if (tt->trafoDepth && !h[part_mode()])
            {
                auto samples = neighbourhood->snakeIntraReferenceSamples[rc.cIdx].offset(rc.x0, rc.y0, rc.cIdx ? 1 : 0);
                unfiltered.substituteFast(h, samples, rc);

                if (ff)
                    filtered.filter(unfiltered, h[strong_intra_smoothing_enabled_flag()], bitDepth, nTbS);
            }

            auto const &pF = ff ? filtered : unfiltered;

            havoc::intra::Table<Sample> *table = h;
            havoc::intra::Function<Sample> *f = table->lookup(rc.cIdx, bitDepth, rc.log2TrafoSize, predModeIntra);
            f(&predSamples(0, 0), predSamples.stride, &pF(0, -1), predModeIntra);
        }

        using Sample = typename SampleType<H>::Type;

        // input picture
        auto &pictureInput = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
        auto source = pictureInput(rc.x0, rc.y0, 0);

        // compute SATD
        havoc_table_hadamard_satd<Sample> *table = h;
        if (rc.log2TrafoSize == 2)
        {
            auto const hadamard_satd = *havoc_get_hadamard_satd(table, 2);
            stateEncodeSubstream->satd += hadamard_satd(&source(0, 0), source.stride, &predSamples(0, 0), predSamples.stride);
        }
        else
        {
            auto const hadamard_satd = *havoc_get_hadamard_satd(table, 3);
            int const size = 1 << rc.log2TrafoSize;
            for (int dy = 0; dy < size; dy += 8)
            {
                for (int dx = 0; dx < size; dx += 8)
                {
                    stateEncodeSubstream->satd += hadamard_satd(&source(dx, dy), source.stride, &predSamples(dx, dy), predSamples.stride);
                }
            }
        }

        if (tt->trafoDepth && tt->blkIdx != 3)
        {
            // copy predicted samples to snake so that subsequent blocks in the transform tree can predict
            Raster<Sample> predPicture = predSamples.offset(-rc.x0 >> log2Resolution, -rc.y0 >> log2Resolution);
            typename Snake<Sample>::Pointer &snake = neighbourhood->snakeIntraReferenceSamples[rc.cIdx];
            int const size = 1 << (rc.log2TrafoSize + log2Resolution);
            snake.copyFrom2D(Turing::Rectangle{ rc.x0, rc.y0, size, size }, predPicture, log2Resolution);
        }
    }
};

// performs luma intra prediction with split_transform_flag = 0 and computes SATD
template <typename F> struct PredictIntraLuma : ReconstructIntraLuma<F> { };

template <> struct PredictIntraLuma<IfCbf<cbf_luma, residual_coding>> : PredictIntraLumaBlock {};
template <> struct PredictIntraLuma<IfCbf<cbf_cb, residual_coding>> : Null<IfCbf<cbf_cb, residual_coding>> {};
template <> struct PredictIntraLuma<IfCbf<cbf_cr, residual_coding>> : Null<IfCbf<cbf_cr, residual_coding>> { };

template <>
struct PredictIntraLuma<transform_tree>
{
    template <class H> static void go(const transform_tree &tt, H &h)
    {
        h[split_transform_flag()] = tt.log2TrafoSize > h[MaxTbLog2SizeY()];
        Syntax<transform_tree>::go(tt, h);
    }
};

template <class Cbf> struct ReconstructInterBlock
{
    template <class H> static void go(IfCbf<Cbf, residual_coding> const &ifcbf, H h)
    {
        using Sample = typename SampleType<H>::Type;

        residual_coding const &rc = ifcbf.f;
        StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
        StateCodedData *stateCodedData = h;
        Candidate<Sample> *candidate = h;
        coding_quadtree const *cqt = h;
        StateEncode *stateEncode = h;

        int const nTbS = 1 << rc.log2TrafoSize;

        assert(h[current(CuPredMode(cqt->x0, cqt->y0))] != MODE_INTRA);
        assert(!h[cu_transquant_bypass_flag()]);

        // input picture
        auto &pictureInput = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
        Raster<Sample> sourceSamples = pictureInput(rc.x0, rc.y0, rc.cIdx);

        Raster<int16_t> resSamplesRaw = stateEncodeSubstream->residual(rc.x0 - cqt->x0, rc.y0 - cqt->y0, rc.cIdx);//(resSamplesRawBuffer, 1 << rc.log2TrafoSize);

        // forward transform
        ALIGN(32, int16_t, coefficientsBuffer[32 * 32]);
        Raster<int16_t> coefficients(coefficientsBuffer, nTbS, 0, 0);
        int const trType = 0;
        bool const checkTSkip = nTbS == 4 && h[transform_skip_enabled_flag()];
        bool bUseTSkip = false;

        static int nn = 0;

        if (candidate->noresidual == 0)
        {
            int constexpr bitDepth = 2 * sizeof(Sample) + 6;
            havoc::table_transform<bitDepth> *table = h;
            auto *transform = *havoc::get_transform<bitDepth>(table, trType, rc.log2TrafoSize);
            transform(coefficients.p, resSamplesRaw.p, resSamplesRaw.stride);
        }

        // quantise
        ALIGN(32, int16_t, temp2[32 * 32]);
        Raster<int16_t> quantizedCoefficients(temp2, nTbS, 0, 0);

        // quantise
        ALIGN(32, int16_t, temp3[32 * 32]);
        Raster<int16_t> zeroCoefficients(temp3, nTbS, 0, 0);

        int const qpScaled = static_cast<QpState *>(h)->getQp(rc.cIdx);
        int const n = 1 << 2 * rc.log2TrafoSize;
        int const scaleQuantise = static_cast<QpState *>(h)->getQuantiseScale(rc.cIdx);
        int const shiftQuantise = 29 - (rc.cIdx ? h[BitDepthC()] : h[BitDepthY()]) + qpScaled / 6 - rc.log2TrafoSize;
        assert(h[slice_type()] != I);
        int const offsetQuantise = 85 << (shiftQuantise - 9);
        int const scale = static_cast<QpState *>(h)->getScale(rc.cIdx);
        int const bitDepth = rc.cIdx ? h[BitDepthC()] : h[BitDepthY()];
        const int shift = rc.log2TrafoSize - 1 + bitDepth - 8;
        auto const quantize = *havoc_get_quantize(static_cast<havoc_table_quantize *>(h));
        bool cbf;
        if (candidate->noresidual == 0)
        {
            if (stateEncode->rdoq)
            {
                double lambda;
                if(stateEncode->useRateControl)
                {
                    StateEncodePicture *stateEncodePicture = h;
                    int segmentPoc = stateEncodePicture->docket->segmentPoc;

                    stateEncode->rateControlParams->takeTokenOnRCLevel();
                    lambda = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                    stateEncode->rateControlParams->releaseTokenOnRCLevel();
                }
                else
                {
                    lambda = static_cast<StateEncodePicture*>(h)->lambda;
                }
                Contexts *contexts = h;
                bool isSdhEnabled = !!h[sign_data_hiding_enabled_flag()];
                Rdoq rdoqEngine(lambda, contexts, scaleQuantise, scale, rc.log2TrafoSize, rc.cIdx ? h[BitDepthC()] : h[BitDepthY()]);
                cbf = !!rdoqEngine.runQuantisation(quantizedCoefficients.p, coefficients.p, scaleQuantise, shiftQuantise,
                                                    n, rc, h[scanIdx()], false, isSdhEnabled);
            }
            else
            {
                cbf = !!quantize(quantizedCoefficients.p, coefficients.p, scaleQuantise, shiftQuantise, offsetQuantise >> (shiftQuantise - 16), n);
            }

            // inverse quantise
            {
                havoc_quantize_inverse *inverseQuantize = *havoc_get_quantize_inverse(static_cast<havoc_table_quantize_inverse *>(h), scale, shift);
                inverseQuantize(coefficients.p, quantizedCoefficients.p, scale, shift, n);
            }
        }
        else
        {
            cbf = 0;
        }
        auto recSamplesPiece = candidate->stateReconstructionCache->components[rc.cIdx].get(stateEncodeSubstream->interPieces[rc.cIdx][1 + candidate->rqtdepth]);
        auto recSamples = recSamplesPiece.offset((rc.x0 - cqt->x0) >> (rc.cIdx ? 1 : 0), (rc.y0 - cqt->y0) >> (rc.cIdx ? 1 : 0));
        auto recSamplesBackup = recSamples;

        auto predSamplesPiece = candidate->stateReconstructionCache->components[rc.cIdx].get(stateEncodeSubstream->interPieces[rc.cIdx][0]);
        auto predSamples = predSamplesPiece.offset((rc.x0 - cqt->x0) >> (rc.cIdx ? 1 : 0), (rc.y0 - cqt->y0) >> (rc.cIdx ? 1 : 0));
        if (candidate->noresidual == 0)
        {
            // add to prediction
            {
                havoc::table_inverse_transform_add<Sample> *table = h;
                auto *inverseTransformAdd = *havoc::get_inverse_transform_add<Sample>(table, trType, rc.log2TrafoSize);
                inverseTransformAdd(recSamples.p, recSamples.stride, predSamples.p, predSamples.stride, coefficients.p, rc.cIdx ? h[BitDepthC()] : h[BitDepthY()]);
            }
        }

        // measure SSD (reconstruction and prediction)
        int32_t backupSSDNoTSkip;
        {
            auto ssd = *havoc_get_ssd<Sample>(h, rc.log2TrafoSize);
            if (candidate->noresidual == 0)
            {
                backupSSDNoTSkip = stateEncodeSubstream->ssd[rc.cIdx] + ssd(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, nTbS, nTbS);
                if (!(checkTSkip && cbf))
                    stateEncodeSubstream->ssd[rc.cIdx] = backupSSDNoTSkip;
            }
            stateEncodeSubstream->ssdPrediction[rc.cIdx] += ssd(sourceSamples.p, sourceSamples.stride, predSamples.p, predSamples.stride, nTbS, nTbS);
        }


        if (checkTSkip && cbf)
        {

            ALIGN(32, int16_t, recSamplesNoTSkipBuffer[32 * 32]);
            Raster<int16_t> recSamplesNoTSkip(recSamplesNoTSkipBuffer, 1 << rc.log2TrafoSize);

            ALIGN(32, int16_t, backupQuantCoeffsBuffer[32 * 32]);
            Raster<int16_t> backupQuantCoeffs(backupQuantCoeffsBuffer, 1 << rc.log2TrafoSize);

            // backup reconstructed samples and quantized coefficients in case we need to test transform skip
            for (int y = 0; y < nTbS; y++)
            {
                for (int x = 0; x < nTbS; x++)
                {
                    recSamplesNoTSkip(x, y) = recSamples(x, y);
                    backupQuantCoeffs(x, y) = quantizedCoefficients(x, y);
                }
            }

            // populate coded data without tskip for rate estimation
            transform_tree const *tt = h;
            if (rc.cIdx == 0)
            {
                // luma
                stateCodedData->transformTree.init(cqt->log2CbSize - rc.log2TrafoSize - (rc.cIdx ? 1 : 0), tt->blkIdx);
                stateCodedData->residual = stateCodedData->transformTree.firstResidual();
            }

            // backup CodedData pointers and words containing CBF flags before storing residual coded data
            StateCodedData backupBefore = *stateCodedData;
            auto const backupCuWord0 = stateCodedData->codedCu.word0();
            auto const backupTuWord0 = stateCodedData->transformTree.word0();
            auto e = h.template change<EstimateRate<void>>();
            auto contextsCostBefore = *static_cast<ContextsAndCost *>(candidate);

            // store residual data to CodedData and make backup of state without tskip
            stateCodedData->residual.transformSkipFlag() = 0;
            assert(quantizedCoefficients.stride == nTbS);
            CodedData::storeResidual(stateCodedData->codedCu, stateCodedData->residual, quantizedCoefficients.p, rc.log2TrafoSize, h[scanIdx()], cbf, stateCodedData->transformTree, rc.cIdx);

            // prepare to measure rate of residual (rewind residual pointer)
            *stateCodedData = backupBefore;
            stateCodedData->transformTreeChroma = stateCodedData->transformTree;
            e(ifcbf);
            auto contextsCostNoTSkip = *static_cast<ContextsAndCost *>(candidate);

            // rewind pointers
            *stateCodedData = backupBefore;
            *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;
            stateCodedData->codedCu.word0() = backupCuWord0;
            stateCodedData->transformTree.word0() = backupTuWord0;


            {
                // test transform skip
                const int bdShift = 13 - bitDepth;
                for (int y = 0; y < nTbS; y++)
                    for (int x = 0; x < nTbS; x++)
                        coefficients(x, y) = resSamplesRaw(x, y) << bdShift;
            }

            {
                // quantise
                if (stateEncode->rdoq)
                {
                    double lambda;
                    if(stateEncode->useRateControl)
                    {
                        StateEncodePicture *stateEncodePicture = h;
                        int segmentPoc = stateEncodePicture->docket->segmentPoc;

                        stateEncode->rateControlParams->takeTokenOnRCLevel();
                        lambda = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                        stateEncode->rateControlParams->releaseTokenOnRCLevel();
                    }
                    else
                    {
                        lambda = static_cast<StateEncodePicture*>(h)->lambda;
                    }
                    Contexts *contexts = h;
                    bool const isSdhEnabled = !!h[sign_data_hiding_enabled_flag()];
                    Rdoq rdoqEngine(lambda, contexts, scaleQuantise, scale, rc.log2TrafoSize, rc.cIdx ? h[BitDepthC()] : h[BitDepthY()]);
                    cbf = !!rdoqEngine.runQuantisation(quantizedCoefficients.p, coefficients.p, scaleQuantise, shiftQuantise,
                                                        n, rc, h[scanIdx()], false, isSdhEnabled);
                }
                else
                {
                cbf = !!quantize(quantizedCoefficients.p, coefficients.p, scaleQuantise, shiftQuantise, offsetQuantise >> (shiftQuantise - 16), n);
                }

                // inverse quantise
                {
                    auto *inverseQuantize = *havoc_get_quantize_inverse(h, scale, shift);
                    inverseQuantize(coefficients.p, quantizedCoefficients.p, scale, shift, n);
                }
            }

            {
                // inverse transform skip and add to predicted samples
                const int bdShift = 20 - bitDepth;
                for (int y = 0; y < nTbS; y++)
                {
                    for (int x = 0; x < nTbS; x++)
                    {
                        int r = coefficients(x, y) << 7;
                        r = (r + (1 << (bdShift - 1))) >> bdShift;
                        recSamples(x, y) = static_cast<Sample>(clipCidx1<>(predSamples(x, y) + r, bitDepth));
                    }
                }
            }

            int32_t backupSSDTSkip;
            {
                // compute distortion as SSD
                auto const ssdFunction = *havoc_get_ssd<Sample>(h, rc.log2TrafoSize);
                backupSSDTSkip = stateEncodeSubstream->ssd[rc.cIdx] + ssdFunction(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, nTbS, nTbS);
                stateEncodeSubstream->ssd[rc.cIdx] = backupSSDTSkip;
            }

            // store residual data to CodedData and make backup of state without tskip
            stateCodedData->residual.transformSkipFlag() = 1;
            assert(quantizedCoefficients.stride == nTbS);
            CodedData::storeResidual(stateCodedData->codedCu, stateCodedData->residual, quantizedCoefficients.p, rc.log2TrafoSize, h[scanIdx()], cbf, stateCodedData->transformTree, rc.cIdx);

            // rewind pointers and estimate rate with tskip
            *stateCodedData = backupBefore;
            stateCodedData->transformTreeChroma = stateCodedData->transformTree;
            e(ifcbf);
            auto contextsCostTSkip = *static_cast<ContextsAndCost *>(candidate);

            // rewind pointers
            *stateCodedData = backupBefore;
            *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;
            stateCodedData->codedCu.word0() = backupCuWord0;
            stateCodedData->transformTree.word0() = backupTuWord0;

            // compute lambdaDistortion with and without tskip
            Lambda reciprocalLambda = getReciprocalLambda(h);
            if(stateEncode->useRateControl)
            {
                StateEncodePicture *stateEncodePicture = h;
                int segmentPoc = stateEncodePicture->docket->segmentPoc;

                stateEncode->rateControlParams->takeTokenOnRCLevel();
                double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                stateEncode->rateControlParams->releaseTokenOnRCLevel();
                reciprocalLambda.set(value);
            }
            contextsCostNoTSkip.lambdaDistortion += backupSSDNoTSkip * reciprocalLambda;
            contextsCostTSkip.lambdaDistortion += backupSSDTSkip   * reciprocalLambda;

            if (contextsCostNoTSkip.cost2() < contextsCostTSkip.cost2())
            {
                // tskip is not cheaper, restore reconstructed samples and quantized coefficients without tskip
                for (int y = 0; y < nTbS; y++)
                {
                    for (int x = 0; x < nTbS; x++)
                    {
                        recSamples(x, y) =static_cast<Sample>(recSamplesNoTSkip(x, y));
                        quantizedCoefficients(x, y) = backupQuantCoeffs(x, y);
                    }
                }

                // restore SSD
                stateEncodeSubstream->ssd[rc.cIdx] = backupSSDNoTSkip;
                cbf = true; // in case cbf became false when testing tskip
                bUseTSkip = false;
            }
            else
            {
                bUseTSkip = true;
            }
        }


        // populate coded data
        {
            transform_tree const *tt = h;

            if (rc.cIdx == 0)
            {
                // luma
                stateCodedData->transformTree.init(cqt->log2CbSize - rc.log2TrafoSize - (rc.cIdx ? 1 : 0), tt->blkIdx);
                stateCodedData->residual = stateCodedData->transformTree.firstResidual();
            }

            if (cbf)
            {
                // backup CodedData pointers and words containing CBF flags before storing residual coded data
                StateCodedData backupBefore = *stateCodedData;
                auto const backupCuWord0 = stateCodedData->codedCu.word0();
                auto const backupTuWord0 = stateCodedData->transformTree.word0();
                auto e = h.template change<EstimateRate<void>>();
                auto contextsCostBefore = *static_cast<ContextsAndCost *>(candidate);
                // This pointer needed for writing/measuring rate
                stateCodedData->transformTreeChroma = stateCodedData->transformTree;
                e(ifcbf.cbf, ae(v));

                StateCodedData backup0 = *stateCodedData;
                auto contextsCostCbf0 = *static_cast<ContextsAndCost *>(candidate);

                *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;
                *stateCodedData = backupBefore;
                // store residual data to CodedData and make backup of state with CBF=1
                stateCodedData->residual.transformSkipFlag() = bUseTSkip;
                assert(quantizedCoefficients.stride == nTbS);
                CodedData::storeResidual(stateCodedData->codedCu, stateCodedData->residual, quantizedCoefficients.p, rc.log2TrafoSize, h[scanIdx()], cbf, stateCodedData->transformTree, rc.cIdx);
                StateCodedData backup1 = *stateCodedData;
                // prepare to measure rate of residual (rewind residual pointer)
                stateCodedData->residual = backupBefore.residual;
                *stateCodedData = backupBefore;
                stateCodedData->transformTreeChroma = stateCodedData->transformTree;        // measure rate of residual_coding with cbf=1 and then restore contexts
                e(ifcbf);
                auto contextsCostCbf1 = *static_cast<ContextsAndCost *>(candidate);

                // compute distortion components for RDO
                // review: could measure distortion in transformed domain
                Lambda reciprocalLambda = getReciprocalLambda(h);
                if(stateEncode->useRateControl)
                {
                    StateEncodePicture *stateEncodePicture = h;
                    int segmentPoc = stateEncodePicture->docket->segmentPoc;

                    stateEncode->rateControlParams->takeTokenOnRCLevel();
                    double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                    stateEncode->rateControlParams->releaseTokenOnRCLevel();
                    reciprocalLambda.set(value);
                }
                contextsCostCbf1.lambdaDistortion += (stateEncodeSubstream->ssd[rc.cIdx]) * reciprocalLambda;
                contextsCostCbf0.lambdaDistortion += (stateEncodeSubstream->ssdPrediction[rc.cIdx]) * reciprocalLambda;

                if (contextsCostCbf0.cost2() < contextsCostCbf1.cost2())
                {
                    // cbf=0 case is cheaper
                    // restore CodedData as it was before storing residual data
                    *stateCodedData = backup0;
                    stateCodedData->codedCu.word0() = backupCuWord0;
                    stateCodedData->transformTree.word0() = backupTuWord0;

                    recSamples = recSamplesBackup;
                    for (int y = 0; y < nTbS; ++y)
                        for (int x = 0; x < nTbS; ++x)
                            recSamples(x, y) = predSamples(x, y);

                    stateEncodeSubstream->ssd[rc.cIdx] = stateEncodeSubstream->ssdPrediction[rc.cIdx];
                }
                else
                {
                    *stateCodedData = backup1;
                }
                *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;
            }
            else if (candidate->noresidual == 1)
            {
                stateCodedData->codedCu.word0().cbf[rc.cIdx] = 0;
                stateCodedData->transformTree.word0().cbf[rc.cIdx] = 0;
                stateEncodeSubstream->ssd[rc.cIdx] = stateEncodeSubstream->ssdPrediction[rc.cIdx];
                recSamples = recSamplesBackup;
                for (int y = 0; y < nTbS; ++y)
                    for (int x = 0; x < nTbS; ++x)
                        recSamples(x, y) = predSamples(x, y);
            }

            if (rc.cIdx == 2)
            {
                stateCodedData->transformTree.p = stateCodedData->residual.p;
            }
            stateCodedData->codedDataAfter = stateCodedData->residual.p;
        }
    }
};


template <class F> struct ReconstructInter<IfCbf<F, residual_coding>> : ReconstructInterBlock<F> {};

template <class H>
int32_t predictIntraLuma(transform_tree const &tt, H &h)
{
    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstreamBase *stateEncodeSubstream = h;
    stateEncodeSubstream->satd = 0;
    stateEncodeSubstream->ssd[0] = 0;
    stateEncodeSubstream->ssd[1] = 0;
    stateEncodeSubstream->ssd[2] = 0;

    {
        Neighbourhood *neighbourhood = h;
        Snake<BlockData>::Cursor *cursor = h;
        cursor->relocate(neighbourhood->snake, tt, h[MinCbLog2SizeY()] - 1, true);
        cursor->value.intra.predModeY = static_cast<StateCodedData *>(h)->codedCu.IntraPredModeY(tt.blkIdx);
    }

    auto r = h.template change<PredictIntraLuma<void>>();
    r(tt);

    return stateEncodeSubstream->satd;
}


template <class H>
int32_t reconstructIntraLuma(IntraPartition const &e, H &h)
{
    const int IntraSplitFlag = h[::IntraSplitFlag()];

    assert(e.split == IntraSplitFlag);

    h[MaxTrafoDepth()] = h[max_transform_hierarchy_depth_intra()] + IntraSplitFlag;

    transform_tree tt{ xPositionOf(e), yPositionOf(e), e.x0, e.y0, e.log2CbSize - e.split, IntraSplitFlag, e.blkIdx };

    {
        Neighbourhood *neighbourhood = h;
        Snake<BlockData>::Cursor *cursor = h;
        cursor->relocate(neighbourhood->snake, tt, h[MinCbLog2SizeY()] - 1, true);
        cursor->value.intra.predModeY = static_cast<StateCodedData *>(h)->codedCu.IntraPredModeY(tt.blkIdx);
    }

    if (tt.trafoDepth == 0)
    {
        StateCodedData *stateCodedData = h;

        if (h[current(CuPredMode(tt.x0, tt.y0))] == MODE_INTRA)
        {
            stateCodedData->transformTree = stateCodedData->codedCu.firstTransformTree();
        }
        else
        {
            stateCodedData->transformTree.p = stateCodedData->codedPu.p;
        }
        stateCodedData->codedDataAfter = stateCodedData->transformTree.p;
    }

    h[split_transform_flag(tt.x0, tt.y0)] = tt.log2TrafoSize > h[MaxTbLog2SizeY()];

    using Sample = typename SampleType<H>::Type;
    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    stateEncodeSubstream->ssd[0] = 0;
    stateEncodeSubstream->ssd[1] = 0;
    stateEncodeSubstream->ssd[2] = 0;
    auto r = h.template change<ReconstructIntraLuma<void>>();
    r(tt);
    return stateEncodeSubstream->ssd[0];
}


template <class H>
int32_t reconstructIntraChroma(transform_tree const &tt, H &h)
{
    if (tt.trafoDepth == 0)
    {
        StateCodedData *stateCodedData = h;

        if (h[current(CuPredMode(tt.x0, tt.y0))] == MODE_INTRA)
        {
            stateCodedData->transformTree = stateCodedData->codedCu.firstTransformTree();
        }
        else
        {
            stateCodedData->transformTree.p = stateCodedData->codedPu.p;
        }
        stateCodedData->codedDataAfter = stateCodedData->transformTree.p;
    }

    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    stateEncodeSubstream->ssd[0] = 0;
    stateEncodeSubstream->ssd[1] = 0;
    stateEncodeSubstream->ssd[2] = 0;
    auto r = h.template change<ReconstructIntraChroma<void>>();
    r(tt);
    return stateEncodeSubstream->ssd[1] + stateEncodeSubstream->ssd[2];
}


template <class H>
int32_t reconstructInter(transform_tree const &tt, H &h)
{
    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    StateEncode *stateEncode = h;

    assert(tt.trafoDepth == 0);
    assert(h[current(CuPredMode(tt.x0, tt.y0))] != MODE_INTRA);

    // subtract prediction from input to obtain residual
    for (int cIdx = 0; cIdx < 3; ++cIdx)
    {
        int const log2TrafoSize = tt.log2TrafoSize - (cIdx ? 1 : 0);
        int const nTbS = 1 << log2TrafoSize;

        using Sample = typename SampleType<H>::Type;
        
        auto &pictureInput = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
        auto sourceSamples = pictureInput(tt.x0, tt.y0, cIdx);

        Raster<Sample> predSamples = h[ReconstructedSamples{ tt.x0, tt.y0, cIdx }];

        auto residual = stateEncodeSubstream->residual(0, 0, cIdx);

        stateEncodeSubstream->ssdPrediction[cIdx] = 0;

        Candidate<Sample> *candidate = h;

        if (h[PartMode()] == PART_2Nx2N)
        {
            candidate->sadResidueQuad[0][0] = 0;
            candidate->sadResidueQuad[0][1] = 0;
            candidate->sadResidueQuad[1][0] = 0;
            candidate->sadResidueQuad[1][1] = 0;
        }

        for (int y = 0; y < nTbS; ++y)
        {
            for (int x = 0; x < nTbS; ++x)
            {
                int16_t const diff = sourceSamples(x, y) - predSamples(x, y);
                residual(x, y) = diff;
                if (h[PartMode()] == PART_2Nx2N)
                {
                    Candidate<Sample> *candidate = h;
                    candidate->sadResidueQuad[y >> (log2TrafoSize - 1)][x >> (log2TrafoSize - 1)] += abs(diff);
                }
            }
        }
    }

    StateCodedData *stateCodedData = h;
    stateCodedData->transformTree.p = stateCodedData->codedPu.p;
    stateCodedData->codedDataAfter = stateCodedData->transformTree.p;

    stateEncodeSubstream->ssd[0] = 0;
    stateEncodeSubstream->ssd[1] = 0;
    stateEncodeSubstream->ssd[2] = 0;

    auto r = h.template change<ReconstructInter<void>>();

    Candidate<Sample> *candidate = h;
    bool inferTransformDepth = !(tt.log2TrafoSize <= h[MaxTbLog2SizeY()] &&
        tt.log2TrafoSize > h[MinTbLog2SizeY()] &&
        tt.trafoDepth < h[MaxTrafoDepth()] && !(h[IntraSplitFlag()] && (tt.trafoDepth == 0)));
    
    bool checkRQT = !inferTransformDepth && !candidate->noresidual;

    if (!checkRQT)
    {
        candidate->rqtdepth = 0;
        r(tt);
    }
    else
    {
        StateCodedData *stateCodedData = h;
        Candidate<Sample> *candidate = h;
        CodedData::Type backupCodedData[3 * (8 * 8 * 8 + 64 * 64) / 2];
        int lengthBackup = 0, distscale = 4;
        int32_t backupssd[6];
        bool cbfZero = false;
        auto m = h.template change<EstimateRate<void>>();

        auto backupContextsAndCostBefore = *static_cast<ContextsAndCost *>(candidate);
        auto backupCuWord0Before = candidate->codedCu.word0().raw;

        //descend transform tree to perform decisions:
        candidate->rqtdepth = 1;
        r(tt);
        auto backupCuWord0One = candidate->codedCu.word0().raw;
        if (!r[rqt_root_cbf()])
            cbfZero = true;
        if (!cbfZero)
        {
            {
                CodedData::Type *src = stateCodedData->codedPu.p;
                CodedData::Type *srcend = stateCodedData->codedDataAfter;
                CodedData::Type *dst = backupCodedData;
                while (src != srcend)
                {
                    *dst++ = *src++;
                    lengthBackup++;
                }
            }
            stateCodedData->transformTree.p = stateCodedData->codedPu.p;
            stateCodedData->codedDataAfter = stateCodedData->transformTree.p;

            //descend transform tree to measure rate:
            if (m[rqt_root_cbf()]) 
                m(tt);

            auto contextsAndCostOne = *static_cast<ContextsAndCost *>(candidate);
            Lambda reciprocalLambda = getReciprocalLambda(h);
            if(stateEncode->useRateControl)
            {
                StateEncodePicture *stateEncodePicture = h;
                int segmentPoc = stateEncodePicture->docket->segmentPoc;

                stateEncode->rateControlParams->takeTokenOnRCLevel();
                double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                stateEncode->rateControlParams->releaseTokenOnRCLevel();
                reciprocalLambda.set(value);
            }
            contextsAndCostOne.lambdaDistortion += (stateEncodeSubstream->ssd[0] + distscale * stateEncodeSubstream->ssd[1] + distscale * stateEncodeSubstream->ssd[2]) * reciprocalLambda;
        
            // restore everything as it was before trying with depth = 1
            *static_cast<ContextsAndCost *>(candidate) = backupContextsAndCostBefore;
            candidate->codedCu.word0().raw = backupCuWord0Before;

            backupssd[0] = stateEncodeSubstream->ssd[0];
            backupssd[1] = stateEncodeSubstream->ssd[1];
            backupssd[2] = stateEncodeSubstream->ssd[2];
            backupssd[3] = stateEncodeSubstream->ssdPrediction[0];
            backupssd[4] = stateEncodeSubstream->ssdPrediction[1];
            backupssd[5] = stateEncodeSubstream->ssdPrediction[2];

            stateEncodeSubstream->ssd[0] = 0;
            stateEncodeSubstream->ssd[1] = 0;
            stateEncodeSubstream->ssd[2] = 0;
            stateEncodeSubstream->ssdPrediction[0] = 0;
            stateEncodeSubstream->ssdPrediction[1] = 0;
            stateEncodeSubstream->ssdPrediction[2] = 0;

            stateCodedData->transformTree.p = stateCodedData->codedPu.p;
            stateCodedData->codedDataAfter = stateCodedData->transformTree.p;

            candidate->rqtdepth = 0;
            r(tt);
            auto backupCuWord0Zero = candidate->codedCu.word0().raw;

            stateCodedData->transformTree.p = stateCodedData->codedPu.p;
            stateCodedData->codedDataAfter = stateCodedData->transformTree.p;

            if (m[rqt_root_cbf()])
                m(tt);
            candidate->lambdaDistortion += (stateEncodeSubstream->ssd[0] + distscale * stateEncodeSubstream->ssd[1] + distscale * stateEncodeSubstream->ssd[2]) * reciprocalLambda;
        
            if (cbfZero || candidate->cost2() < contextsAndCostOne.cost2())
            {
                candidate->codedCu.word0().raw = backupCuWord0Zero;
                candidate->rqtdepth = 0;
            }
            else
            {
                candidate->codedCu.word0().raw = backupCuWord0One;
                {
                    CodedData::Type *dst = stateCodedData->codedPu.p;
                    CodedData::Type *src = backupCodedData;
                    for (int i = 0; i < lengthBackup; i++)
                    {
                        *dst++ = *src++;
                    }
                }
                candidate->rqtdepth = 1;
                stateEncodeSubstream->ssd[0] = backupssd[0];
                stateEncodeSubstream->ssd[1] = backupssd[1];
                stateEncodeSubstream->ssd[2] = backupssd[2];
                stateEncodeSubstream->ssdPrediction[0] = backupssd[3];
                stateEncodeSubstream->ssdPrediction[1] = backupssd[4];
                stateEncodeSubstream->ssdPrediction[2] = backupssd[5];
            }
        }
        else
        {
            candidate->codedCu.word0().raw = backupCuWord0One;
            candidate->rqtdepth = 0;
        }
        *static_cast<ContextsAndCost *>(candidate) = backupContextsAndCostBefore;
    }
    return stateEncodeSubstream->ssd[0] + stateEncodeSubstream->ssd[1] + stateEncodeSubstream->ssd[2];
}


template <> struct ReconstructInter<Element<split_transform_flag, ae>>
{
    template <class H> static void go(Element<split_transform_flag, ae> e, H &h)
    {
        //StateCodedData *stateCodedData = h;
        //h[e.v] = stateCodedData->transformTree.word0().split_transform_flag;
        //Reconstruct<Element<split_transform_flag, ae>>::go(e, h);
    }
};


// explicit template instantiations - these were semi-automatically generated based on missing symbols reported by the linker
#define EXPLICIT_INSTANTIATIONS(mode, Sample) \
    template int32_t predictIntraLuma<Handler<Search<Mode<mode>>,Candidate<Sample>,StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> >(transform_tree const &,Handler<Search<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> &); \
    template int32_t reconstructIntraLuma<Handler<Search<Mode<mode>>,Candidate<Sample>,StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> >(IntraPartition const &,Handler<Search<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> &); \
    template int32_t reconstructIntraChroma<Handler<Search<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> >(transform_tree const &,Handler<Search<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> &); \
    template int32_t reconstructInter<Handler<Search<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> >(transform_tree const &,Handler<Search<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> &); \
    template int32_t reconstructInter<Handler<SearchMerge2Nx2N<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> >(transform_tree const &,Handler<SearchMerge2Nx2N<Mode<mode>>, Candidate<Sample>, StateEncodeSubstream<Sample>,StateEncodePicture2<Sample>,StateEncode> &);

EXPLICIT_INSTANTIATIONS(0, uint8_t)
EXPLICIT_INSTANTIATIONS(1, uint8_t)
EXPLICIT_INSTANTIATIONS(0, uint16_t)
EXPLICIT_INSTANTIATIONS(1, uint16_t)
