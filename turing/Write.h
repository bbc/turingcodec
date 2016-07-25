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

// Functions for writing bitstream

#ifndef INCLUDED_Write_h
#define INCLUDED_Write_H

#pragma once


#include "Global.h"
#include "GlobalState.h"
#include "HevcTypes.h"
#include "Cabac.h"
#include "Snake.h"
#include "Dsp.h"
#include "Mvp.h"
#include "Reconstruct.h"
#include "Picture.h"
#include "Binarization.h"
#include "Vanilla.h"
#include "Measure.h"
#include "Speed.h"
#include "Search.h"
#include "Cost.h"
#include "StateEncode.h"
#include "LoopFilter.h"
#include "havoc/ssd.h"
#include "SyntaxSei.h"
#include <cassert>
#include <iostream>
#include <cstdlib>


// write the specified element to bitstream
template <typename F> struct Write :
    Syntax<F>
    {
    };

template <> struct Write<void> { };


// encode the specified element
template <class F> struct Encode : Write<F> { };
template <> struct Encode<void> { };


// estimate the bitrate of the specified element
template <class F> struct EstimateRate;
template <> struct EstimateRate<void> { };


template <> struct Write<video_parameter_set_rbsp>
{
    template <class H> static void go(video_parameter_set_rbsp e, H &hhh)
    {
        std::shared_ptr<Vps> vps = hhh[Table<Vps>()][0];
        auto hh =  hhh.extend(&*vps);
        auto h =  hh.extend(&vps->ptl);

        Syntax<video_parameter_set_rbsp>::go(e, h);
    }
};


template <> struct Write<seq_parameter_set_rbsp>
{
    template <class H> static void go(seq_parameter_set_rbsp e, H &hhh)
    {
        std::shared_ptr<Sps> sps = hhh[Table<Sps>()][0];
        auto hh =  hhh.extend(&*sps);
        auto h =  hh.extend(&sps->ptl);

        Syntax<seq_parameter_set_rbsp>::go(e, h);
    }
};

template <>
struct Write<CabacRestart>
{
    template <class H> static void go(const CabacRestart &e, H &h)
    {
        CabacWriter &writer = *static_cast<CabacWriter *>(h);
        writer.restartCabac2();
    }
};

DEFINE_STRUCT_ARITY_1(FinishCabac, extraBit);

template <>
struct Write<FinishCabac>
{
    template <class H> static void go(FinishCabac &v, H &h)
    {
        CabacWriter &writer = *static_cast<CabacWriter *>(h);
        writer.finishCabac();
        if (v.extraBit)
        {
            writer.writeBit(1);
        }
    }
};

template <>
struct Write<EndOfSubStream>
{
    template <class H> static void go(EndOfSubStream &v, H &h)
    {
        assert(h[byte_aligned()]);

        StateEncode& stateEncode = *static_cast<StateEncode *>(h);

        // Set this flag so that we break out of the loop in Syntax<slice_segment_data, H>
        h[end_of_slice_segment_flag()] = 1;
    }
};

template <typename T, class H>
struct NextBits<T, H>
{
    static T go(H &h, int n, bool aligned)
    {
        // Writing bitstream but next_bits() is nonetheless invoked from Syntax.h.
        // Returning this value ensures that the logic in the syntax tables operates appropriately.
        return 0x000001;
    }
};

template <class V>
struct Write<Element<V, f>>
{
    template <class H> static void go(Element<V, f> fun, H &h)
    {
        static_cast<BitWriter *>(h)->writeBits(fun.m.n, Fixed<V>::value);
    }
};


template <class V>
struct Write<Element<V, u>>
{
    template <class H> static void go(Element<V, u> fun, H &h)
    {
        static_cast<BitWriter *>(h)->writeBits(fun.m.n, h[fun.v]);
    }
};


template <class V>
struct Write<Element<V, uv>>
{
    template <class H> static void go(Element<V, uv> e, H &h)
    {
        auto const n = NumberOfBitsUv<V>::get(e.v, h);
        static_cast<BitWriter *>(h)->writeBits(n, h[e.v]);
    }
};


template <class V>
struct Write<Element<V, b>>
{
    template <class H> static void go(Element<V, b> fun, H &h)
    {
        static_cast<BitWriter *>(h)->writeBits(fun.m.n, h[fun.v]);
    }
};


template <class V>
struct Write<Element<V, ue>>
{
    template <class H> static void go(Element<V, ue> fun, H &h)
    {
        int code = h[fun.v];
        assert(code >= 0);

        ++code;
        int len = 1;
        for (int temp = code; temp != 1; temp >>= 1)
        {
            len += 2;
        }

        BitWriter *bitWriter = h;
        bitWriter->writeBits(len / 2, 0);
        bitWriter->writeBits((len + 1) / 2, code);
    }
};


template <class V>
struct Write<Element<V, se>>
{
    template <class H> static void go(Element<V, se> fun, H &h)
    {
        int code = h[fun.v];
        code = ( code <= 0) ? -code<<1 : (code<<1)-1;

        ++code;
        int len = 1;
        for (int temp = code; temp != 1; temp >>= 1)
        {
            len += 2;
        }

        BitWriter *bitWriter = h;
        bitWriter->writeBits(len / 2, 0);
        bitWriter->writeBits((len + 1) / 2, code);
    }
};

template <>
struct Write<nal_unit>
{
    template <class H> static void go(const nal_unit &fun, H &h)
    {
        NalWriter *nalWriter = h;

        Syntax<nal_unit>::go(fun, h);

        auto startRbspPos = nalWriter->data->size();

        switch (h[nal_unit_type()])
        {

#define X(TYPE, NAME, DESCRIPTION, ELEMENT, VCL) \
case TYPE: h(ELEMENT()); break;

            NAL_UNIT_TYPES
#undef X

            default:
                assert(!"unrecognised nal_unit_type");
                break;
        }

        // Emulation prevention: insert some 0x03 bytes into the stream, as necessary.
        BitWriter::insertEp3Bytes(*nalWriter->data, startRbspPos);
    }
};

template <> struct Write<sei_message>
{
    template <class H> static void go(sei_message, H &h)
    {
        BitWriter *bitWriter = h;

        // write sei_payload to a temporary buffer first
        std::vector<uint8_t> seiData;
        auto *rbspData = bitWriter->data;
        bitWriter->data = &seiData;
        h(sei_payload(h[last_payload_type_byte()], 0));
        bitWriter->data = rbspData;

        // now we know the size of sei_payload() data and can write the sei_message() header
        h[last_payload_size_byte()] = static_cast<typename ValueType<last_payload_size_byte>::Type>(seiData.size());

        while (h[last_payload_type_byte()] > 255)
        {
            h(ff_byte() /* equal to 0xFF */, f(8));
            h[last_payload_type_byte()] -= 255;
        }
        h(last_payload_type_byte(), u(8));

        while (h[last_payload_size_byte()] > 255)
        {
            h(ff_byte() /* equal to 0xFF */, f(8));
            h[last_payload_size_byte()] -= 255;
        }
        h(last_payload_size_byte(), u(8));

        // append the sei_payload() data to the bitstream
        bitWriter->data->insert(bitWriter->data->end(), seiData.begin(), seiData.end());
    }
};

template <>
struct Write<access_unit_delimiter_rbsp>
{
    template <class H> static void go(access_unit_delimiter_rbsp, H &)
    {
        assert(!"not yet implemented - coherent state for pic_type needs to be arranged");
    }
};


template <class S>
struct Access<byte_aligned, S, typename std::enable_if<std::is_base_of<BitWriter, S>::value>::type>
{
    typedef bool Type;
    static Type get(byte_aligned, BitWriter &bitWriter)
    {
        return !bitWriter.shift;
    }
};


template <class S>
struct Access<more_data_in_payload, S, typename std::enable_if<std::is_base_of<BitWriter, S>::value>::type>
{
    typedef bool Type;
    static Type get(more_data_in_payload, BitWriter &s)
    {
        // There will be more data if stream is not aligned
        return !Access<byte_aligned, BitWriter>::get(byte_aligned(), s);
    }
};


struct NotImplemented
{
    template <class F, class H> static void go(F, H &)
    {
        throw std::runtime_error("not yet implemented");
    }
};

// The following are not yet used in encoding
template <> struct Write<scaling_list_data> : NotImplemented {};
template <> struct Write<pred_weight_table> : NotImplemented {};
template <> struct Write<hrd_parameters> : NotImplemented {};
template <> struct Write<Element<reserved_payload_extension_data, uv>> : NotImplemented {};


template <>
struct Write<vui_parameters>
{
    template <class H> static void go(const vui_parameters &fun, H &h)
    {
        Syntax<vui_parameters>::go(fun, h);
    }
};


// Suppress this call (this data has already been generated with emulation prevention - it will be appended later)
template <> struct Write<slice_segment_data> : Null<slice_segment_data> { };


// Suppress this call (rbsp_slice_segment_trailing_bits() is included in the last substream data that will be appended later)
template <> struct Write<rbsp_slice_segment_trailing_bits> : Null<rbsp_slice_segment_trailing_bits> { };


// Entropy estimation values from HM. These have radix position 15, i.e. rate of 1 bit is 0x8000.
static int32_t entropyBitsHm[] = {
        // Corrected table, most notably for last state
        0x07b23, 0x085f9, 0x074a0, 0x08cbc, 0x06ee4, 0x09354, 0x067f4, 0x09c1b, 0x060b0, 0x0a62a, 0x05a9c, 0x0af5b, 0x0548d, 0x0b955, 0x04f56, 0x0c2a9,
        0x04a87, 0x0cbf7, 0x045d6, 0x0d5c3, 0x04144, 0x0e01b, 0x03d88, 0x0e937, 0x039e0, 0x0f2cd, 0x03663, 0x0fc9e, 0x03347, 0x10600, 0x03050, 0x10f95,
        0x02d4d, 0x11a02, 0x02ad3, 0x12333, 0x0286e, 0x12cad, 0x02604, 0x136df, 0x02425, 0x13f48, 0x021f4, 0x149c4, 0x0203e, 0x1527b, 0x01e4d, 0x15d00,
        0x01c99, 0x166de, 0x01b18, 0x17017, 0x019a5, 0x17988, 0x01841, 0x18327, 0x016df, 0x18d50, 0x015d9, 0x19547, 0x0147c, 0x1a083, 0x0138e, 0x1a8a3,
        0x01251, 0x1b418, 0x01166, 0x1bd27, 0x01068, 0x1c77b, 0x00f7f, 0x1d18e, 0x00eda, 0x1d91a, 0x00e19, 0x1e254, 0x00d4f, 0x1ec9a, 0x00c90, 0x1f6e0,
        0x00c01, 0x1fef8, 0x00b5f, 0x208b1, 0x00ab6, 0x21362, 0x00a15, 0x21e46, 0x00988, 0x2285d, 0x00934, 0x22ea8, 0x008a8, 0x239b2, 0x0081d, 0x24577,
        0x007c9, 0x24ce6, 0x00763, 0x25663, 0x00710, 0x25e8f, 0x006a0, 0x26a26, 0x00672, 0x26f23, 0x005e8, 0x27ef8, 0x005ba, 0x284b5, 0x0055e, 0x29057,
        0x0050c, 0x29bab, 0x004c1, 0x2a674, 0x004a7, 0x2aa5e, 0x0046f, 0x2b32f, 0x0041f, 0x2c0ad, 0x003e7, 0x2ca8d, 0x003ba, 0x2d323, 0x0010c, 0x3bfbb};


template <class T>
struct EntropyEstimate
{
    T table[128];
    EntropyEstimate()
    {
        for (int i = 0; i < 128; ++i)
            table[i].set(entropyBitsHm[i], 15);
    }
};


static Cost measureEncodeDecision(ContextModel &contextModel, int binVal)
{
    static const EntropyEstimate<Cost> entropyEstimate;
    Cost const rate = entropyEstimate.table[contextModel.state ^ binVal];

    if (binVal != contextModel.getMps())
        contextModel.updateLps();
    else
        contextModel.updateMps();

    return rate;
}


template <class V>
struct Measure<EncodeDecision<V>>
{
    template <class H> static void go(EncodeDecision<V> v, H &h)
    {
        StateEstimateRate *stateEstimateRate = h;
        Contexts *contexts = h;
        ContextModel &contextModel = contexts->template get<typename Context<V>::Type>(v.ctxInc);

        auto const rate = measureEncodeDecision(contextModel, v.binVal);
        stateEstimateRate->rate += rate;
    }
};


template <class V>
struct Write<EncodeDecision<V>>
{
    template <class H> static void go(EncodeDecision<V> const& v, H &h)
    {
        CabacWriter *writer = h;
        Contexts *contexts = h;
        ContextModel &contextModel = contexts->template get<typename Context<V>::Type>(v.ctxInc);

        auto const pStateIdx = contextModel.getState();
        auto const valMPS = contextModel.getMps();
        auto const qRangeIdx = ( writer->ivlCurrRange  >>  6 ) & 3; // (9-45)
        auto const ivlLpsRange = rangeTabLPS(pStateIdx, qRangeIdx); // uiLPS

        writer->ivlCurrRange -= ivlLpsRange; // m_uiRange

        int numBits;
        if (v.binVal != valMPS)
        {
            static const int renormTable[] =
            {
                    6, 5, 4, 4, 3, 3, 3, 3,
                    2, 2, 2, 2, 2, 2, 2, 2,
                    1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1,
            };
            numBits = renormTable[ivlLpsRange >> 3];

            writer->ivlLow += writer->ivlCurrRange;
            writer->ivlCurrRange = ivlLpsRange;

            contextModel.updateLps();
        }
        else
        {
            numBits = writer->ivlCurrRange < 256 ? 1 : 0;

            contextModel.updateMps();
        }

        writer->ivlLow <<= numBits;
        writer->ivlCurrRange <<= numBits;
        writer->bitsLeft -= numBits;

        if (numBits)
            writer->testAndWriteOut();
    }
};


template <class V>
struct Measure<EncodeBypass<V>>
{
    template <class H> static void go(EncodeBypass<V> const&, H &h)
    {
        StateEstimateRate *stateEstimateRate = h;
        ++stateEstimateRate->rate;
    }
};


template <class V>
struct Write<EncodeBypass<V>>
{
    template <class H> static void go(EncodeBypass<V> const& encodeBypass, H &h)
    {
        CabacWriter *writer = h;
        writer->ivlLow <<= 1;
        if (encodeBypass.binVal)
            writer->ivlLow  += writer->ivlCurrRange;
        --writer->bitsLeft;
        writer->testAndWriteOut();
    }
};


template <class V>
struct Measure<EncodeTerminate<V>>
{
    template <class H> static void go(const EncodeTerminate<V> &, H &)
    {
        // should really count cost of pcm_flag here, but it is tiny
    }
};


template <class V>
struct Write<EncodeTerminate<V>>
{
    template <class H> static void go(const EncodeTerminate<V> &encodeTerminate, H &h)
    {
        CabacWriter *writer = h;

        writer->ivlCurrRange -= 2;

        int numBits;
        if( encodeTerminate.binVal )
        {
            writer->ivlLow  += writer->ivlCurrRange;
            writer->ivlCurrRange = 2;
            numBits = 7;
        }
        else
        {
            numBits = writer->ivlCurrRange < 256 ? 1 : 0;
        }

        if (numBits)
        {
            writer->ivlLow <<= numBits;
            writer->ivlCurrRange <<= numBits;
            writer->bitsLeft -= numBits;
            writer->testAndWriteOut();
        }
    }
};


template <>
struct Encode<coding_tree_unit>
{
    template <class H> static void go(const coding_tree_unit &ctu, H &h)
    {
        typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
        static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");


        //#ifdef SNAKE_DEBUG
        coding_quadtree cqt(
                (h[CtbAddrInRs()] % h[PicWidthInCtbsY()]) << h[CtbLog2SizeY()],
                (h[CtbAddrInRs()] / h[PicWidthInCtbsY()]) << h[CtbLog2SizeY()],
                h[CtbLog2SizeY()],
                0);
        //#endif

        preCtu(h);

        Syntax<coding_tree_unit>::go(ctu, h);

        postCtu(h);

#ifdef SNAKE_DEBUG
        bool lastCtuInRow = cqt.x0 + (1 << cqt.log2CbSize) >= h[pic_width_in_luma_samples()];

        if (lastCtuInRow)
        {
            cqt.x0 += h[CtbSizeY()];
            --cqt.log2CbSize;

            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                Picture<Sample> &reconstructed = h[ReconstructedPicture()];
                Raster<Sample> recSamples = reconstructed(-(1 << cqt.log2CbSize), 0, cIdx);

                NeighbourhoodEnc<Sample> *neighbourhood = h;
                neighbourhood->snakeIntraReferenceSamples[cIdx].copyFrom2D(cqt, recSamples, cIdx ? 1 : 0);
            }

            cqt.y0 += 1 << cqt.log2CbSize;

            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                Picture<Sample> &reconstructed = h[ReconstructedPicture()];
                Raster<Sample> recSamples = reconstructed(-(1 << cqt.log2CbSize), 0, cIdx);

                NeighbourhoodEnc<Sample> *neighbourhood = h;
                neighbourhood->snakeIntraReferenceSamples[cIdx].copyFrom2D(cqt, recSamples, cIdx ? 1 : 0);
            }
        }
#endif
    }
};


template <>
struct Write<coding_quadtree>
{
    template <class H> static void go(const coding_quadtree &cqt, H &h)
    {
        StateCodedData *stateCodedData = h;

        static_cast<Neighbourhood *>(h)->snake.checkPosition(before, cqt);

        auto const codedCtDepth = stateCodedData->codedCu.word0().CtDepth;

        h[split_cu_flag(cqt.x0, cqt.y0)] = cqt.cqtDepth < codedCtDepth;

        Syntax<coding_quadtree>::go(cqt, h);

        static_cast<Neighbourhood *>(h)->snake.checkPosition(after, cqt);
    }
};


template <>
struct Encode<coding_quadtree>
{
    // This function encodes a CTU. It is called only for depth=0.
    // Deeper levels recurse using Search<coding_quadtree>, for RDO search,
    // then Write<coding_quadtree> to write out the bitstream
    // Review: consistency - naming, access and usage of various Candidates
    template <class H> static void go(const coding_quadtree &cqt, H &h)
    {
        typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
        static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");

        StateReconstructionCache<Sample> *stateReconstructionCache = h;

        assert(cqt.cqtDepth == 0);

        NeighbourhoodEnc<Sample> *neighbourhood = h;
        neighbourhood->snake.checkPosition(before, cqt);

        stateReconstructionCache->reset();

        Candidate<Sample> **currentCandidate = h;
        Candidate<Sample> *originalCandidate = *currentCandidate;

        originalCandidate->stateReconstructionCache = stateReconstructionCache;
        originalCandidate->check = cartesianToZ(cqt.x0, cqt.y0);

        originalCandidate->resetZero();

        CandidateStash<Sample> candidate(*originalCandidate, cqt, *stateReconstructionCache);
        candidate.resetPieces();
        candidate.copy(*originalCandidate, before, cqt, true, true);
        candidate.copyIntraReferenceSamplesLuma(*originalCandidate, before, cqt);
        candidate.copyIntraReferenceSamplesChroma(*originalCandidate, before, cqt);
        candidate.ContextsAndCost::copy(*originalCandidate);

#ifdef _DEBUG
        // review: is it useful to surface these as a hidden command-line option?
        bool constexpr checkDistortion = false; // Review fix this for 10-bit content
        bool constexpr checkRate = true;
#else
        bool constexpr checkDistortion = false;
        bool constexpr checkRate = false;
#endif

        CandidateStash<Sample> candidateCheck(*originalCandidate, cqt, *stateReconstructionCache);
        candidateCheck.resetPieces();
        if (checkRate)
        {
            candidateCheck.copy(*originalCandidate, before, cqt, true, true);
            candidateCheck.copyIntraReferenceSamplesLuma(*originalCandidate, before, cqt);
            candidateCheck.copyIntraReferenceSamplesChroma(*originalCandidate, before, cqt);
            candidateCheck.ContextsAndCost::copy(*originalCandidate);
        }

#ifdef DEBUG_PIECES
        candidate.StatePieces<Sample>::z[0] = zPositionOf(cqt);
        candidate.StatePieces<Sample>::z[1] = zPositionOf(cqt) >> 2;
        candidate.StatePieces<Sample>::z[2] = zPositionOf(cqt) >> 2;
#endif

        *currentCandidate = &candidate;

        {
            Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->searchTotal);
            auto hSearch = h.template change<Search<void>>();

            if(h[cu_qp_delta_enabled_flag()])
            {
                static_cast<QpState *>(h)->setCanWrite(false);

                StateEncode *stateEncode = h;

                if(stateEncode->useRateControl)
                {
                    StateEncodePicture *stateEncodePicture = h;
                    int currentPictureLevel = stateEncodePicture->docket->sopLevel;
                    bool isIntraSlice = h[slice_type()] == I;
                    stateEncode->rateControlEngine->setValidityFlag(false, h[CtbAddrInRs()]);
                    double bpp = stateEncode->rateControlEngine->getCtuTargetBits(isIntraSlice, h[CtbAddrInRs()]);
                    double estLambda;
                    int qp;

                    if(isIntraSlice)
                    {
                        stateEncode->rateControlEngine->getCtuEstLambdaAndQp(bpp, h[SliceQpY()], h[CtbAddrInRs()], estLambda, qp);
                    }
                    else
                    {
                        estLambda = stateEncode->rateControlEngine->getCtuEstLambda(bpp, h[CtbAddrInRs()], currentPictureLevel);
                        qp = stateEncode->rateControlEngine->getCtuEstQp(h[CtbAddrInRs()], currentPictureLevel);
                    }

                    stateEncodePicture->lambda = estLambda;
                    stateEncodePicture->reciprocalLambda.set(1.0 / stateEncodePicture->lambda);
                    stateEncodePicture->reciprocalSqrtLambda = sqrt(1.0 / stateEncodePicture->lambda);
                    static_cast<QpState *>(h)->setQpInternal(0, 0, 64, qp, 0);
                }
                else
                    static_cast<QpState *>(h)->setQpInternal(0, 0, 64, h[SliceQpY()], 0);
            }
            else
                static_cast<QpState *>(h)->setCanWrite(true);

            StateFunctionTables *stateFunctionTables = h;

            // walks tree, finds best settings, stores then into CodedData sequential structure
            if (stateFunctionTables->instruction_set_support & HAVOC_LZCNT)
            {
                auto hSearch = h.template change<Search<Mode<1>>>();
                Search<coding_quadtree>::go(cqt, hSearch);
            }
            else
            {
                auto hSearch = h.template change<Search<Mode<0>>>();
                Search<coding_quadtree>::go(cqt, hSearch);
            }
        }

        candidate.snakeIntraReferenceSamples[0].checkPosition(after, cqt);
        candidate.snakeIntraReferenceSamples[1].checkPosition(after, cqt);
        candidate.snakeIntraReferenceSamples[2].checkPosition(after, cqt);

        Contexts contextsAfterSearch = candidate;

        if (checkRate)
        {
            // measure rate again to check against that obtained during searches

            static_cast<StateCodedData &>(candidateCheck) = static_cast<StateCodedData &>(candidate);
            candidateCheck.StateCodedData::startReading();

            candidateCheck.StateEstimateRate::rate.set(0, 0);

            *currentCandidate = &candidateCheck;
            auto w = h.template change<EstimateRate<void>>();
            w[split_cu_flag(cqt.x0, cqt.y0)] = cqt.cqtDepth < static_cast<StateCodedData *>(w)->codedCu.word0().CtDepth;

            // walk tree and estimate the rate of resultant bitstream
            Syntax<coding_quadtree>::go(cqt, w);

            candidateCheck.Contexts::checkSameAs(contextsAfterSearch);
            ASSERT(candidateCheck.StateEstimateRate::rate == candidate.StateEstimateRate::rate);// && "rate estimated during search doesn't match whole-CTU estimated rate");
        }

        static_cast<StateCodedData &>(*originalCandidate) = static_cast<StateCodedData &>(candidate);
        originalCandidate->StateCodedData::startReading();

        *currentCandidate = originalCandidate;

        {
            Candidate<Sample> *candidate = h;

            Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->write);

            // walk tree and writes bitstream
            // review: derive split_cu_flag directly from coded data
            h[split_cu_flag(cqt.x0, cqt.y0)] = cqt.cqtDepth < static_cast<StateCodedData *>(h)->codedCu.word0().CtDepth;

            if(h[cu_qp_delta_enabled_flag()])
            {
                static_cast<QpState *>(h)->setCanWrite(true);
                static_cast<QpState *>(h)->setQpInternal(0, 0, 64, h[SliceQpY()], 0);
            }

            BitWriter &writer = *static_cast<BitWriter *>(h);
            size_t bytesBeforeWriting = writer.pos();

            auto w = h.template change<Write<void>>();
            Syntax<coding_quadtree>::go(cqt, w);

            if(h[cu_qp_delta_enabled_flag()])
                static_cast<QpState *>(h)->swapInternalMemory();

            StateEncode *stateEncode = h;
            if(stateEncode->useRateControl)
            {
                size_t bytesAfterWriting = writer.pos();
                int codingBits = ((int)bytesAfterWriting - (int)bytesBeforeWriting) << 3;
                StateEncodePicture *stateEncodePicture = h;
                int currentPictureLevel = stateEncodePicture->docket->sopLevel;

                // Update the rate controller engine
                stateEncode->rateControlEngine->updateCtuController(codingBits, h[slice_type()] == I, h[CtbAddrInRs()], currentPictureLevel);
            }

            // review: test should pass even if this flag set
            if(!h[cu_qp_delta_enabled_flag()])
                candidate->checkSameAs(contextsAfterSearch);
        }

        // Copy reconstructed pieces into reconstructed picture buffer
        candidate.StatePieces<Sample>::commit(h[ReconstructedPicture()], cqt);

        if (checkDistortion)
        {
            // Measure distortion again, check that it is same as that measured during search
            typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
            static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");

            // input picture
            ThreePlanes<Sample> &pictureInput = dynamic_cast<ThreePlanes<Sample>&>(*static_cast<StateEncodePicture *>(h)->docket->picture);
            ThreePlanes<Sample> &pictureReconstructed = h[ReconstructedPicture()];

            int32_t ssd = 0;
            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                auto const log2TrafoSize = cqt.log2CbSize - (cIdx ? 1 : 0);
                auto const srcA = pictureInput(cqt.x0, cqt.y0, cIdx);
                auto const srcB = pictureReconstructed(cqt.x0, cqt.y0, cIdx);
                auto const width = std::min(1 << log2TrafoSize, (h[pic_width_in_luma_samples()] - cqt.x0) >> (cIdx ? 1 : 0));
                auto const height = std::min(1 << log2TrafoSize, (h[pic_height_in_luma_samples()] - cqt.y0) >> (cIdx ? 1 : 0));
                for (int y = 0; y < height; ++y)
                {
                    for (int x = 0; x < width; ++x)
                    {
                        int diff = srcA(x, y) - srcB(x, y);
                        ssd += diff * diff;
                    }
                }
            }

            Cost lambdaDistortion = getReciprocalLambda(h) * ssd;
            bool const ok = candidate.ContextsAndCost::lambdaDistortion == lambdaDistortion;
            if (!ok) std::cout << h[PicOrderCntVal()] << " " << cqt << " " << ssd << " " << lambdaDistortion << " " << candidate.ContextsAndCost::lambdaDistortion << "\n";
            ASSERT(ok && "distortion (or lambda) does not match that measured during searches");
        }

        // review: unnecessary copy
        originalCandidate->snakeIntraReferenceSamples[0].copyBlockFrom(candidate.snakeIntraReferenceSamples[0], after, cqt, 0);
        originalCandidate->snakeIntraReferenceSamples[1].copyBlockFrom(candidate.snakeIntraReferenceSamples[1], after, cqt, 1);
        originalCandidate->snakeIntraReferenceSamples[2].copyBlockFrom(candidate.snakeIntraReferenceSamples[2], after, cqt, 1);

        neighbourhood->snake.checkPosition(after, cqt);
    }
};


template <class Direction>
struct Write<Deleted<coding_quadtree, Direction>>
{
    template <class H> static void go(const Deleted<coding_quadtree, Direction> &cqt, H &h)
    {
        // Review: use Direction as a hint for optimisation
        Neighbourhood *neighbourhood = h;
        BlockData blockData = BlockData();
        neighbourhood->snake.commitRectangle(cqt, blockData, neighbourhood->MinCbLog2SizeYMinus1);
        neighbourhood->recordMerge(h, cqt, true);
    }
};


template <>
struct Write<coding_unit>
{
    template <class H> static void go(const coding_unit &cu, H &h)
    {
        StateEncodeSubstreamBase *stateEncodeSubstreamBase = h;
        StateCodedData *stateCodedData = h;
        coding_quadtree const *cqt = h;
        Neighbourhood *neighbourhood = h;
        Snake<BlockData>::Cursor *cursor = h;

        cursor->relocate(neighbourhood->snake, cu, neighbourhood->MinCbLog2SizeYMinus1);

        static_cast<QpState *>(h)->preCu(cu, h);

        stateEncodeSubstreamBase->cqt = *cqt;
        stateEncodeSubstreamBase->partIdx = 0;
        stateCodedData->partIdx = 0;

        auto const CuPredMode = h[current(::CuPredMode(cu.x0, cu.y0))];

        BlockData &blockData = cursor->current(cu.x0, cu.y0, neighbourhood->MinCbLog2SizeYMinus1);
        blockData.setup(cqt, CuPredMode);
        blockData.intra.pcm = stateCodedData->codedCu.word0().pcm_flag;

        cursor->relocate(neighbourhood->snake, cu, h[MinCbLog2SizeY()] - 1); // < review: remove?

        if (CuPredMode != MODE_INTRA) stateCodedData->firstPu();

        if(h[cu_qp_delta_enabled_flag()] && static_cast<QpState *>(h)->getCanWrite())
        {

            // Set up cu_qp_delta_abs and sign flag
            int rowQgModulo = (cu.y0 & (h[CtbSizeY()] - 1)) >> 3;
            int colQgModulo = (cu.x0 & (h[CtbSizeY()] - 1)) >> 3;

            StateEncode *stateEncode = h;
            int QpY;
            // Review for intra RC
            if(stateEncode->useRateControl)
                QpY = stateEncode->rateControlEngine->getCtuStoredQp(h[CtbAddrInRs()]);
            else
                QpY = h[SliceQpY()];

            if(stateEncode->useAq)
            {
                AdaptiveQuantisation &aqInfo = dynamic_cast<AdaptiveQuantisation&>(*static_cast<StateEncodePicture *>(h)->docket->aqInfo);
                int qpOffset = aqInfo.getAqOffset(cu.y0, cu.x0, cqt->cqtDepth);
                QpY = Clip3(0, 51, QpY + qpOffset);
            }

            int qpForPrediction;
            if (h[IsCuQpDeltaCoded()])
                qpForPrediction = static_cast<QpState *>(h)->getCodedQp();
            else
                qpForPrediction = static_cast<QpState *>(h)->getQpYPred(cu, h);
            int valueToWrite = QpY - qpForPrediction;
            valueToWrite = (valueToWrite + 78 + h[QpBdOffsetY()] + (h[QpBdOffsetY()]/2)) % (52 + h[QpBdOffsetY()]) - 26 - (h[QpBdOffsetY()]/2);
            static_cast<QpState *>(h)->setQpValue(QpY);

            // Set the internal QP for this CU to the predictor one (default one)
            int size = std::max<int>(1, 1 << ((cu.log2CbSize - 3) << 1));
            static_cast<QpState *>(h)->setQpInternal(rowQgModulo, colQgModulo, size, qpForPrediction);

            int signFlagToWrite = valueToWrite > 0 ? 0 : 1;
            h[cu_qp_delta_abs()] = abs(valueToWrite);
            h[cu_qp_delta_sign_flag()] = signFlagToWrite;
        }

        // walk syntax call tree to write constituent elements, PUs, TT, etc.
        Syntax<coding_unit>::go(cu, h);

        static_cast<QpState *>(h)->postCu(cu, h);

        neighbourhood->recordMerge(h, *cqt, CuPredMode == MODE_INTRA);

        if (std::is_same<typename H::Tag, Write<void>>::value)
        {

            static_cast<StatePicture *>(h)->loopFilterPicture->processCu(h, cu);

            if (CuPredMode == MODE_INTRA)
            {
                StatePicture *decodedPicture = h;
                decodedPicture->motion->fillRectangleIntra(cu);
            }
            StateEncode *stateEncode = h;
            if(stateEncode->useRateControl && h[slice_type()] != I)
            {
                stateEncode->rateControlEngine->updateValidityFlag(!(h[current(cu_skip_flag(cu.x0, cu.y0))]), h[CtbAddrInRs()]);
            }
        }

        // advance CodedData pointer to next CU
        stateCodedData->codedCu.p = stateCodedData->transformTreeChroma.p;
    }
};


template <> struct Write<pcm_sample>
{
    template <class H> static void go(pcm_sample, H &h)
    {
        Neighbourhood *neighbourhood = h;
        Snake<BlockData>::Cursor *cursor = h;
        coding_quadtree const *cqt = h;
        CabacWriter *cabacWriter = h;
        StateCodedData *stateCodedData = h;

        assert(h[PcmBitDepthY()] == 8);
        assert(h[PcmBitDepthC()] == 8);

        // Input picture
        typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
        static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");

        ThreePlanes<Sample> *picture = &dynamic_cast<ThreePlanes<Sample>&>(*static_cast<StateEncodePicture *>(h)->docket->picture);

        assert(!"this function writes PCM to bitstream but does not reconstruct - intra and inter prediction from here may fail");

        // Write IPCM samples to bitstream

        const int nCbS = 1 << cqt->log2CbSize;
        for (int j=0; j<nCbS; ++j)
        {
            for (int i=0; i<nCbS; ++i)
            {
                // pcm_sample_luma
                cabacWriter->writeBits(8, (*picture)[0](cqt->x0 + i, cqt->y0 + j));
            }
        }

        for (int j=0; j<nCbS/2; ++j)
        {
            for (int i=0; i<nCbS/2; ++i)
            {
                // pcm_sample_chroma
                cabacWriter->writeBits(8, (*picture)[1](cqt->x0 / 2 + i, cqt->y0 / 2 + j));
            }
        }

        for (int j=0; j<nCbS/2; ++j)
        {
            for (int i=0; i<nCbS/2; ++i)
            {
                // pcm_sample_chroma
                cabacWriter->writeBits(8, (*picture)[2](cqt->x0 / 2 + i, cqt->y0 / 2 + j));
            }
        }

        h(CabacRestart());

        // Update snake
        cursor->commit(*cqt, neighbourhood->MinCbLog2SizeYMinus1, true);

        // Update coded data read pointers appropriately
        stateCodedData->transformTreeChroma.p = stateCodedData->codedCu.firstTransformTree().p;
    }
};

template <class> struct MeasureLuma;
template <class> struct MeasureChroma;
template <class> struct EstimateRateLuma;
template <class> struct EstimateRateChroma;

template <>
struct Write<transform_tree>
{
    template <class H> static void go(const transform_tree &tt, H &h)
    {
        typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
        static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");

        coding_quadtree const *cqt = h;
        Snake<BlockData>::Cursor *cursor = h;
        Neighbourhood *neighbourhood = h;
        StateCodedData *stateCodedData = h;

        // true if transform_tree is same square block as an intra partition
        bool const intraPartition = h[current(CuPredMode(tt.x0, tt.y0))] == MODE_INTRA && (h[IntraSplitFlag()] == tt.trafoDepth);

        if (intraPartition)
        {
            cursor->relocate(neighbourhood->snake, tt, neighbourhood->MinCbLog2SizeYMinus1, true);
            cursor->value.intra.predModeY = static_cast<StateCodedData *>(h)->codedCu.IntraPredModeY(tt.blkIdx);
        }

        if (tt.trafoDepth == 0)
        {
            if (h[current(CuPredMode(tt.x0, tt.y0))] != MODE_INTRA && h[PartMode()] == PART_2Nx2N && h[merge_flag(tt.x0, tt.y0)])
            {
                stateCodedData->transformTreeAncestry[tt.trafoDepth].p = stateCodedData->codedPu.p;
                assert(h[cbf_luma(tt.x0, tt.y0)] || h[cbf_cr(tt.x0, tt.y0)] || h[cbf_cb(tt.x0, tt.y0)]);
            }

            if (h[current(CuPredMode(tt.x0, tt.y0))] == MODE_INTRA)
            {
                stateCodedData->transformTree = stateCodedData->codedCu.firstTransformTree();
                stateCodedData->transformTreeChroma = stateCodedData->codedCu.firstTransformTreeChroma();
            }
            else
            {
                stateCodedData->transformTree = { stateCodedData->codedPu.p };
                stateCodedData->transformTreeChroma = stateCodedData->transformTree;
                stateCodedData->transformTree.check(tt.trafoDepth, tt.blkIdx);
            }
        }
        else
        {
            if (h[current(CuPredMode(tt.x0, tt.y0))] == MODE_INTRA)
            {
            }
            else
            {
                if (tt.trafoDepth && tt.blkIdx == 0)
                {
                    // Advance to next Transform Tree in CodedData (have just split)
                    ++stateCodedData->transformTree.p;
                    stateCodedData->transformTreeChroma = stateCodedData->transformTree;
                    stateCodedData->transformTree.check(tt.trafoDepth, tt.blkIdx);
                }
                else if (tt.trafoDepth == 0)
                {
                    // very first TT
                    stateCodedData->transformTree = { stateCodedData->codedPu.p };
                    stateCodedData->transformTreeChroma = stateCodedData->transformTree;
                    stateCodedData->transformTree.check(tt.trafoDepth, tt.blkIdx);
                }
            }
        }

        stateCodedData->transformTreeAncestry[tt.trafoDepth] = stateCodedData->transformTree;

        if (h[current(CuPredMode(tt.x0, tt.y0))] != MODE_INTRA && tt.trafoDepth)
        {
            ASSERT(stateCodedData->transformTreeAncestry[tt.trafoDepth - 1].word0().split_transform_flag);
        }

        h[split_transform_flag()] = infer(split_transform_flag(tt.x0, tt.y0, tt.trafoDepth), h);

        //auto p = stateCodedData->transformTree.p;
        Syntax<transform_tree>::go(tt, h);

        //std::cout << "Write " << tt;
        //for (; p != stateCodedData->transformTree.p; ++p)
        //{
        //	std::cout << " " << std::hex << std::setw(4) << std::setfill('0') << *p;
        //}
        //std::cout << "\n";

        typedef std::is_same<typename H::Tag, MeasureLuma<void>> MeasuringLuma;
        typedef std::is_same<typename H::Tag, MeasureChroma<void>> MeasuringChroma;

        if (intraPartition &&
                !MeasuringLuma::value && !std::is_same<typename H::Tag, EstimateRateLuma<void>>::value)
        {
            cursor->value.intra.predModeY = static_cast<StateCodedData *>(h)->codedCu.IntraPredModeY(tt.blkIdx);
            cursor->commit(tt, neighbourhood->MinCbLog2SizeYMinus1, true);
        }
    }
};


template <>
struct Write<transform_unit>
{
    template <class H> static void go(const transform_unit &tu, H &h)
    {
        Syntax<transform_unit>::go(tu, h);

        if (h[current(CuPredMode(tu.x0, tu.y0))] != MODE_INTRA)
        {
            StateCodedData *stateCodedData = h;
            stateCodedData->transformTree.p = stateCodedData->residual.p;
            stateCodedData->transformTreeChroma.p = stateCodedData->residual.p;
        }

        // Review: replace if () with static techniques.
        if (std::is_same<typename H::Tag, Write<void>>::value)
        {
            static_cast<StatePicture *>(h)->loopFilterPicture->processTu(h, tu);
        }
    }
};

template <class Prefix, class Suffix, class H>
void setLastSignificantCoeff(H &h, residual_coding e, int n)
{
    static const std::array<int, 32> groupIdx = {0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9};
    static const std::array<int, 10> minInGroup = {0,1,2,3,4,6,8,12,16,24};

    h[Prefix()] = groupIdx[n];
    h[Suffix()] = n - minInGroup[groupIdx[n]];
}

template <int log2TrafoSize, class H>
int ctxIncSigCoeffFlag(H &h, int cIdx, int xC, int yC, int codedSubBlockFlagRight, int codedSubBlockFlagBelow)
{
    int sigCtx;
    if (log2TrafoSize == 2)
    {
        const int ctxIdxMap[16] =
        {
                0, 1, 4, 5,
                2, 3, 4, 5,
                6, 6, 8, 8,
                7, 7, 8, 8
        };
        sigCtx = ctxIdxMap[(yC << 2) + xC];
    }
    else if (xC + yC == 0)
    {
        sigCtx = 0;
    }
    else
    {
        const int xS = xC >> 2;
        const int yS = yC >> 2;
        int prevCsbf = 0;

        if (xS < (1 << (log2TrafoSize - 2)) - 1)
        {
            prevCsbf |= codedSubBlockFlagRight ? 1 : 0;
        }
        if (yS < (1 << (log2TrafoSize - 2)) - 1)
        {
            prevCsbf |= codedSubBlockFlagBelow ? 2 : 0;
        }

        const int xP = xC & 3;
        const int yP = yC & 3;

        if (prevCsbf == 0)
        {
            sigCtx = (xP + yP == 0) ? 2 : (xP + yP < 3) ? 1 : 0;
        }
        else if (prevCsbf == 1)
        {
            sigCtx = (yP == 0) ? 2 : (yP == 1) ? 1 : 0;
        }
        else if (prevCsbf == 2)
        {
            sigCtx = (xP == 0) ? 2 : (xP == 1) ? 1 : 0;
        }
        else
        {
            assert(prevCsbf == 3);
            sigCtx = 2;
        }

        if (cIdx == 0)
        {
            if (xS + yS > 0)
            {
                sigCtx += 3;
            }

            if (log2TrafoSize == 3)
            {
                sigCtx += (h[scanIdx()] == 0) ? 9 : 15;
            }
            else
            {
                sigCtx += 21;
            }
        }
        else
        {
            if (log2TrafoSize == 3)
            {
                sigCtx += 9;
            }
            else
            {
                sigCtx += 12;
            }
        }
    }

    assert(sigCtx >= 0);
    assert(sigCtx < 54);

    if (cIdx == 0)
    {
        return sigCtx;
    }
    else
    {
        return 27 + sigCtx;
    }
}


template <int log2TrafoSize>
struct OptimizedResidualCodingEncode
{
    template <class H> static void go(const residual_coding &rc, H &hParent)
    {
        typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
        static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");

        StateEncodeSubstreamBase *stateEncodeSubstreamBase = hParent;
        ResidualCodingState residualCodingState(hParent);
        auto h = hParent.extend(&residualCodingState);
        bool sdhEnabledFlag = !!h[sign_data_hiding_enabled_flag()];
        bool signHidden;

        StateReconstructionCache<Sample> *stateReconstructionCache = h;
        StateCodedData *stateCodedData = h;
        Candidate<Sample> *candidate = h;
        coding_quadtree const *cqt = h;

        CodedData::Residual &residual = stateCodedData->residual;
        CodedData::SubBlock subBlock = stateCodedData->residual.initialSubBlock(log2TrafoSize);

        // Find last subblock and last significant coefficient in that block
        int const lastScanPos = CodedData::Residual::firstSetBit(*subBlock.p);
        int const lastSubBlock = residual.lastSubBlock(rc.log2TrafoSize);
        {
            const int xS = ScanOrder(log2TrafoSize - 2, h[scanIdx()], lastSubBlock, 0);
            const int yS = ScanOrder(log2TrafoSize - 2, h[scanIdx()], lastSubBlock, 1);
            const int xC = (xS << 2) + ScanOrder(2, h[scanIdx()], lastScanPos, 0);
            const int yC = (yS << 2) + ScanOrder(2, h[scanIdx()], lastScanPos, 1);
            setLastSignificantCoeff<last_sig_coeff_x_prefix, last_sig_coeff_x_suffix>(h, rc, (h[scanIdx()] == 2) ? yC : xC);
            setLastSignificantCoeff<last_sig_coeff_y_prefix, last_sig_coeff_y_suffix>(h, rc, (h[scanIdx()] == 2) ? xC : yC);
        }

        //Encode transform_skip_flag
        if (h[transform_skip_enabled_flag()] && !h[cu_transquant_bypass_flag()] && (log2TrafoSize <= h[Log2MaxTransformSkipSize()]))
            h(transform_skip_flag(rc.x0, rc.y0, rc.cIdx), ae(v));

        h(last_sig_coeff_x_prefix(), ae(v));
        h(last_sig_coeff_y_prefix(), ae(v));
        if (h[last_sig_coeff_x_prefix()] > 3)
            h(last_sig_coeff_x_suffix(), ae(v));
        if (h[last_sig_coeff_y_prefix()] > 3)
            h(last_sig_coeff_y_suffix(), ae(v));

        int &greater1Ctx = stateEncodeSubstreamBase->greater1Ctx;
        int lastGreater1Flag = 1;
        int lastGreater1Ctx = -1;
        int previousGreater1Flag;

        // review: consider bitmask instead of array here
        uint8_t codedSubBlockFlagSnake[17] = { 0 };

        bool codedSubBlock0 = true;

        for (int i = lastSubBlock; i >= 0; --i)
        {
            int const xS = ScanOrder(log2TrafoSize - 2, h[scanIdx()], i, 0);
            int const yS = ScanOrder(log2TrafoSize - 2, h[scanIdx()], i, 1);

            bool codedSubBlockFlag = !!residual.codedSubBlockFlag(i);

            if (i == 0 && !codedSubBlockFlag)
            {
                codedSubBlock0 = false;
                codedSubBlockFlag = true;
                residual.p = subBlock.p;
                static CodedData::Type zero = 0;
                subBlock.p = &zero;
            }

            codedSubBlockFlagSnake[8 + xS - yS] = codedSubBlockFlag;

            int inferSbDcSigCoeffFlag = 0;
            if ((i < lastSubBlock) && (i > 0))
            {
                int const csbfCtx =
                        codedSubBlockFlagSnake[8 + (xS + 1) - yS] +
                        codedSubBlockFlagSnake[8 + xS - (yS + 1)];

                int const ctxInc = (rc.cIdx ? 2 : 0) + std::min(csbfCtx, 1);

                h(EncodeDecision<coded_sub_block_flag>(codedSubBlockFlag, ctxInc));

                inferSbDcSigCoeffFlag = 1;
            }

            if (codedSubBlockFlag)
            {
                // Sublock loop for sig_coeff_flag
                for (int n = (i == lastSubBlock) ? lastScanPos - 1 : 15; n >= 0; n--)
                {
                    if (n > 0 || !inferSbDcSigCoeffFlag)
                    {
                        auto const xC = (xS << 2) + ScanOrder(2, h[scanIdx()], n, 0);
                        auto const yC = (yS << 2) + ScanOrder(2, h[scanIdx()], n, 1);

                        auto const codedSubBlockFlagRight = codedSubBlockFlagSnake[8 + (xS + 1) - yS];
                        auto const codedSubBlockFlagBelow = codedSubBlockFlagSnake[8 + xS - (yS + 1)];
                        auto const ctxInc = ctxIncSigCoeffFlag<log2TrafoSize>(h, rc.cIdx, xC, yC, codedSubBlockFlagRight, codedSubBlockFlagBelow);

                        h(EncodeDecision<sig_coeff_flag>(subBlock.sigCoeffFlag(n), ctxInc));

                        if (subBlock.sigCoeffFlag(n))
                        {
                            inferSbDcSigCoeffFlag = 0;
                        }
                    }
                }

                // Sublock loop for coeff_abs_level_greater1_flag
                // Also derive coeff_abs_level_greater2_flag(n)
                int firstSigScanPos = 16;
                int lastSigScanPos = -1;
                int numGreater1Flag = 0;
                int lastGreater1ScanPos = -1;

                int& ctxSet = stateEncodeSubstreamBase->ctxSet;

                CodedData::Type *pAbsCoeff = subBlock.p + 3;

                for (int n = 15; n >= 0; n--)
                {
                    h[coeff_abs_level_greater2_flag(n)] = 0;

                    if (subBlock.sigCoeffFlag(n))
                    {
                        auto const xC = (xS << 2) + ScanOrder(2, h[scanIdx()], n, 0);
                        auto const yC = (yS << 2) + ScanOrder(2, h[scanIdx()], n, 1);

                        if (numGreater1Flag < 8)
                        {
                            h[coeff_abs_level_greater1_flag(n)] = subBlock.coeffGreater1Flag(n);

                            if (numGreater1Flag == 0)
                            {
                                ctxSet = (i == 0 || rc.cIdx > 0) ? 0 : 2;
                                if (lastGreater1Ctx < 0)
                                {
                                    lastGreater1Ctx = 1;
                                }
                                else
                                {
                                    lastGreater1Ctx = greater1Ctx;
                                    if (lastGreater1Ctx > 0)
                                    {
                                        lastGreater1Flag = previousGreater1Flag;
                                        lastGreater1Ctx = lastGreater1Flag ? 0 : lastGreater1Ctx + 1;
                                    }
                                }
                                if (lastGreater1Ctx == 0)
                                {
                                    ctxSet = ctxSet + 1;
                                    assert(ctxSet < 4);
                                }
                                greater1Ctx = 1;
                            }
                            else
                            {
                                if (greater1Ctx > 0)
                                {
                                    lastGreater1Flag = previousGreater1Flag;
                                }
                                greater1Ctx = lastGreater1Flag ? 0 : greater1Ctx + 1;
                            }

                            h(coeff_abs_level_greater1_flag(n), ae(v));
                            previousGreater1Flag = h[coeff_abs_level_greater1_flag(n)];
                            numGreater1Flag++;
                            if (h[coeff_abs_level_greater1_flag(n)] && lastGreater1ScanPos == -1)
                            {
                                lastGreater1ScanPos = n;
                                h[coeff_abs_level_greater2_flag(lastGreater1ScanPos)] = *pAbsCoeff > 2 ? 1 : 0;
                            }
                        }
                        else
                        {
                            h[coeff_abs_level_greater1_flag(n)] = 0;
                        }

                        if (lastSigScanPos == -1)
                            lastSigScanPos = n;

                        firstSigScanPos = n;

                        if (subBlock.coeffGreater1Flag(n))
                            ++pAbsCoeff;
                    }
                }

                if (lastGreater1ScanPos != -1)
                {
                    h(coeff_abs_level_greater2_flag(lastGreater1ScanPos), ae(v));
                }

                if(h[cu_transquant_bypass_flag()])
                    signHidden = false;
                else
                    signHidden = lastSigScanPos - firstSigScanPos > 3;

                // review - unnecessary copy
                for (int n = 15; n >= 0; n--)
                {
                    bool sdhCondition = ((n != firstSigScanPos) || !sdhEnabledFlag || !signHidden);
                    if (subBlock.sigCoeffFlag(n) && sdhCondition)
                    {
                        h[coeff_sign_flag(n)] = subBlock.coeffSignFlag(n);
                        h(coeff_sign_flag(n), ae(v));
                    }
                }

                // Sublock loop for coeff_abs_level_remaining
                pAbsCoeff = subBlock.p + 3;
                int numSigCoeff = 0;
                stateEncodeSubstreamBase->cLastAbsLevel = 0;
                stateEncodeSubstreamBase->cLastRiceParam = 0;
                for (int n = 15; n >= 0; n--)
                {
                    const int xC = (xS << 2) + ScanOrder(2, h[scanIdx()], n, 0);
                    const int yC = (yS << 2) + ScanOrder(2, h[scanIdx()], n, 1);
                    if (subBlock.sigCoeffFlag(n))
                    {
                        ASSERT(codedSubBlockFlag);
                        ASSERT(subBlock.sigCoeffFlag(n));

                        int16_t absCoeff = 1;
                        if (subBlock.coeffGreater1Flag(n))
                        {
                            absCoeff = *pAbsCoeff++;
                        }

                        int& baseLevel = stateEncodeSubstreamBase->baseLevel;
                        baseLevel = 1 + h[coeff_abs_level_greater1_flag(n)] + h[coeff_abs_level_greater2_flag(n)];
                        if (baseLevel == ((numSigCoeff < 8) ? ((n == lastGreater1ScanPos) ? 3 : 2) : 1))
                        {
                            h[coeff_abs_level_remaining(n)] = absCoeff - baseLevel;
                            h(coeff_abs_level_remaining(n), ae(v));
                        }
                        else
                        {
                            h[coeff_abs_level_remaining(n)] = 0;
                        }
                        numSigCoeff++;
                    }
                }

                if (subBlock.significant())
                {
                    subBlock.p = pAbsCoeff;
                }
            }
        }

        if (codedSubBlock0)
        {
            residual.p = subBlock.p;
        }
    }
};


template <class V>
struct Write<IfCbf<V, residual_coding>>
{
    template <class H> static void go(const IfCbf<V, residual_coding> &e, H &h)
    {
        residual_coding rc = e.f;
        coding_quadtree const *cqt = h;
        transform_tree const *tt = h;
        StateCodedData *stateCodedData = h;

        if (rc.cIdx == 0)
        {
            stateCodedData->residual = stateCodedData->transformTree.firstResidual();
        }
        else if (rc.cIdx == 1 && h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA)
        {
            stateCodedData->residual = stateCodedData->transformTreeChroma.firstResidual();
        }

        if (h[e.cbf])
        {
            CopyValueToState<residual_coding, H> copy{ h, rc };

            if (std::is_same<typename H::Tag, Write<void>>::value)
            {
                static_cast<StatePicture *>(h)->loopFilterPicture->processRc(h, rc);
            }

            switch (rc.log2TrafoSize)
            {
                case 2:
                    OptimizedResidualCodingEncode<2>::go(rc, h);
                    break;
                case 3:
                    OptimizedResidualCodingEncode<3>::go(rc, h);
                    break;
                case 4:
                    OptimizedResidualCodingEncode<4>::go(rc, h);
                    break;
                case 5:
                    OptimizedResidualCodingEncode<5>::go(rc, h);
                    break;
            }
        }

        if (rc.cIdx == 0)
        {
            if (h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA)
            {
                stateCodedData->transformTree.p = stateCodedData->residual.p;
            }
        }
        else if (rc.cIdx == 2)
        {
            if (h[current(CuPredMode(cqt->x0, cqt->y0))] != MODE_INTRA)
            {
                stateCodedData->transformTree.p = stateCodedData->residual.p;
            }
            stateCodedData->transformTreeChroma.p = stateCodedData->residual.p;
        }
        stateCodedData->codedDataAfter = stateCodedData->residual.p;
    }
};


template <>
struct Write<prediction_unit>
{
    template <class H> static void go(const prediction_unit &pu, H &h)
    {
        StateSubstream *stateSubstream = h;
        StateCodedData *stateCodedData = h;
        Neighbourhood *neighbourhood = h;
        Snake<BlockData>::Cursor *cursor = h;

        ASSERT(h[slice_type()] != I);

        cursor->relocate(neighbourhood->snake, pu, h[MinCbLog2SizeY()] - 1);
        BlockData &blockData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);
        processPredictionUnit(h, pu, blockData, stateSubstream->partIdx);

        Syntax<prediction_unit>::go(pu, h);

        if (std::is_same<typename H::Tag, Write<void>>::value)
        {
            // This is the CTU-final pass, writing bitstream.
            // Populate loop filter and temporal motion buffers accordingly.

            StatePicture *statePictureBase = h;
            statePictureBase->loopFilterPicture->processPu(h, pu);

            StatePicture *decodedPicture = h;
            decodedPicture->motion->fillRectangle(pu, blockData);
        }

        cursor->commit(pu, h[MinCbLog2SizeY()] - 1);

        neighbourhood->recordMerge(h, pu);

        if (!h[merge_flag(pu.x0, pu.y0)])
        {
            stateCodedData->codedPu.p = stateCodedData->codedPu.next();
        }

        if (h[PartMode()] != PART_2Nx2N)
        {
            ++stateSubstream->partIdx;
            ++stateCodedData->partIdx;
        }

        assert(stateCodedData->partIdx == stateSubstream->partIdx);

        stateCodedData->transformTree.p = stateCodedData->codedPu.p;
        stateCodedData->transformTreeChroma.p = stateCodedData->codedPu.p;
    }
};


template <>
struct Write<mvd_coding>
{
    template <class H> static void go(const mvd_coding &mvdc, H &h)
    {
        StateCodedData *stateCodedData = h;

        MotionVector const &mvd = stateCodedData->codedPu.mvd(mvdc.refList);

        // review: unncessary data wrangling, optimise using custom state(stores abs, sign) and accessors or inline element writing code here
        int const absMv[2] = { std::abs((int)mvd[0]), std::abs((int)mvd[1]) };
        h[abs_mvd_greater0_flag(0)] = absMv[0] > 0 ? 1 : 0;
        h[abs_mvd_greater0_flag(1)] = absMv[1] > 0 ? 1 : 0;
        h[abs_mvd_greater1_flag(0)] = absMv[0] > 1 ? 1 : 0;
        h[abs_mvd_greater1_flag(1)] = absMv[1] > 1 ? 1 : 0;
        h[abs_mvd_minus2(0)] = absMv[0] - 2;
        h[abs_mvd_minus2(1)] = absMv[1] - 2;
        h[mvd_sign_flag(0)] = mvd[0] < 0 ? 1 : 0;
        h[mvd_sign_flag(1)] = mvd[1] < 0 ? 1 : 0;

        Syntax<mvd_coding>::go(mvdc, h);
    }
};

template <>
struct Write<OutputPicture>
{
    template <class H> static void go(const OutputPicture &op, H &h)
    {
        std::shared_ptr<StatePicture> pic = static_cast<StatePictures *>(h)->getPicByPoc(op.poc);
        static_cast<StateEncode *>(h)->decodedPictures.push_back(pic);
    }
};


template <>
struct Write<DeletePicture>
{
    template <class H> static void go(const DeletePicture &dp, H &h)
    {
    }
};


template <>
struct Write<DpbClear>
{
    template <class H> static void go(const DpbClear &f, H &h)
    {
    }
};


template <> struct Write<pic_timing>
{
    template <class H> static void go(pic_timing e, H &h)
    {
        Hrd hrdDefault{};
        Hrd *hrd = getHrd(h);
        if (!hrd) hrd = &hrdDefault;
        auto h2 =  h.extend(&*hrd);
        Syntax<pic_timing>::go(e, h2);
    }
};

// Much SEI not yet implemented. Comment appropriate line here when writing of an SEI message is implemented.
template <> struct Write<reserved_sei_message> : NotImplemented { };
template <> struct Write<buffering_period> : NotImplemented { };
template <> struct Write<pan_scan_rect> : NotImplemented { };
template <> struct Write<filler_payload> : NotImplemented { };
template <> struct Write<user_data_registered_itu_t_t35> : NotImplemented { };
template <> struct Write<recovery_point> : NotImplemented { };
template <> struct Write<scene_info> : NotImplemented { };
template <> struct Write<picture_snapshot> : NotImplemented { };
template <> struct Write<progressive_refinement_segment_start> : NotImplemented { };
template <> struct Write<progressive_refinement_segment_end> : NotImplemented { };
template <> struct Write<film_grain_characteristics> : NotImplemented { };
template <> struct Write<post_filter_hint> : NotImplemented { };
template <> struct Write<tone_mapping_info> : NotImplemented { };
template <> struct Write<frame_packing_arrangement> : NotImplemented { };
template <> struct Write<display_orientation> : NotImplemented { };
template <> struct Write<structure_of_pictures_info> : NotImplemented { };
template <> struct Write<decoding_unit_info> : NotImplemented { };
template <> struct Write<temporal_sub_layer_zero_index> : NotImplemented { };
template <> struct Write<scalable_nesting> : NotImplemented { };
template <> struct Write<region_refresh_info> : NotImplemented { };
template <> struct Write<no_display> : NotImplemented { };
template <> struct Write<time_code> : NotImplemented { };
template <> struct Write<mastering_display_colour_volume> : NotImplemented { };
template <> struct Write<segmented_rect_frame_packing_arrangement> : NotImplemented { };
template <> struct Write<temporal_motion_constrained_tile_sets> : NotImplemented { };
template <> struct Write<chroma_resampling_filter_hint> : NotImplemented { };
template <> struct Write<knee_function_info> : NotImplemented { };
template <> struct Write<colour_remapping_info> : NotImplemented { };
template <> struct Write<deinterlaced_field_identification> : NotImplemented { };
template <> struct Write<content_light_level> : NotImplemented { };
template <> struct Write<layers_not_present> : NotImplemented { };
template <> struct Write<inter_layer_constrained_tile_sets> : NotImplemented { };
template <> struct Write<bsp_nesting> : NotImplemented { };
template <> struct Write<bsp_initial_arrival_time> : NotImplemented { };
template <> struct Write<sub_bitstream_property> : NotImplemented { };
template <> struct Write<alpha_channel_info> : NotImplemented { };
template <> struct Write<overlay_info> : NotImplemented { };
template <> struct Write<temporal_mv_prediction_constraints> : NotImplemented { };
template <> struct Write<frame_field_info> : NotImplemented { };
template <> struct Write<three_dimensional_reference_displays_info> : NotImplemented { };
template <> struct Write<depth_representation_info> : NotImplemented { };
template <> struct Write<multiview_scene_info> : NotImplemented { };
template <> struct Write<multiview_acquisition_info> : NotImplemented { };
template <> struct Write<multiview_view_position> : NotImplemented { };


// Review: deprecated, remove.  Newer code uses EstimateRate<> instead.
template <class F>
struct Measure :
    Write<F>
    {
    };

template <class F> struct MeasureLuma : Measure<F> { };
template <> struct MeasureLuma<Element<cbf_cb, ae>> : Null<Element<cbf_cb, ae>>{};
template <> struct MeasureLuma<Element<cbf_cr, ae>> : Null<Element<cbf_cr, ae>>{};
template <> struct MeasureLuma<Element<intra_chroma_pred_mode, ae>> : Null<Element<intra_chroma_pred_mode, ae>>{};


template <class F> struct MeasureChroma : Measure<F> { };
template <> struct MeasureChroma<Element<cbf_luma, ae>> : Null<Element<cbf_luma, ae>>{};
template <> struct MeasureChroma<Element<split_transform_flag, ae>> : Null<Element<split_transform_flag, ae>>{};


template <class H>
void writeHeaders(H &h)
{
    // zero_byte is required before parameter set NALUs
    h(zero_byte()  /* equal to 0x00 */, f(8));

    // Write the VPS NALU
    h[nal_unit_type()] = VPS_NUT;
    h(byte_stream_nal_unit(0));

    h(zero_byte()  /* equal to 0x00 */, f(8));

    // Write the SPS NALU
    h[nal_unit_type()] = SPS_NUT;
    h(byte_stream_nal_unit(0));

    h(zero_byte()  /* equal to 0x00 */, f(8));

    // Write the PPS NALU
    h[nal_unit_type()] = PPS_NUT;
    h(byte_stream_nal_unit(0));
}

#endif
