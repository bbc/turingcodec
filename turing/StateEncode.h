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

#ifndef INCLUDED_StateEncode_h
#define INCLUDED_StateEncode_h

#pragma once

// review: too many headers here - reorganise?
#include "StatePicture.h"
#include "StateWavefront.h"
#include "Global.h"
#include "Dsp.h"
#include "Search.h"
#include "Speed.h"
#include "Reconstruct.h"
#include "GlobalState.h"
#include "ReconstructionCache.h"
#include "ThreadPool.h"
#include "StateSpatial.h"
#include "CabacWriter.h"
#include "InputQueue.h"
#include "FixedPoint.h"
#include "Cost.h"
#include "Memory.h"
#include "CodedData.h"
#include "sei/alternative_transfer_characteristics.h"
#include "sei/decoded_picture_hash.h"
#include "sei/pic_timing.h"
#include "sei/active_parameter_sets.h"
#include "sei/user_data_unregistered.h"
#include "sei/mastering_display_colour_volume.h"
#include "AdaptiveQuantisation.h"
#include "RateControl.h"
#include <boost/program_options.hpp>
#include <condition_variable>
#include <deque>
#include <fstream>
#include <vector>


// State to track estimated bitrate.
struct StateEstimateRate
{
    // bitrate since start of CTU
    Cost rate;
};


struct ContextsAndCost :
    Contexts,
    StateEstimateRate
{
    // Lagrangian cost metric (R + lambda * D) accumulated since start of CTU
    Cost cost2()
    {
        return this->rate + this->lambdaDistortion;
    }

    // cumulative (lambda * D) since start of CTU
    Cost lambdaDistortion;

    // called at start of CTU encoding
    void resetZero()
    {
        this->rate.set(0, 0);
        this->lambdaDistortion.set(0, 0);
    }

    void setMax()
    {
        this->rate = std::numeric_limits<Cost>::max();
        this->lambdaDistortion.set(0, 0);
    }

    void copy(ContextsAndCost const &other)
    {
        *this = other;
    }
};


// Maintains lists of objects each providing access to cached reconstructed blocks for later assembly of reconstructed CTU.
template <typename Sample>
struct StatePieces
{
    StateReconstructionCache<Sample> *stateReconstructionCache; // review - better as function parameter than stored here?

    template <class Rectangular>
    void copyBefore(Rectangular const &rectangular, StatePieces &other, int mask = 0x7)
    {
        this->stateReconstructionCache = other.stateReconstructionCache;

        assert(this->pieces[0][before] == this->pieces[0][after]);
        this->parentPieces[0] = other.pieces[0][after];
        this->parentPieces[1] = other.pieces[1][after];
        this->parentPieces[2] = other.pieces[2][after];
        this->zz[0] = other.zz[0];
        this->zz[1] = other.zz[1];
        this->zz[2] = other.zz[2];
#ifdef DEBUG_PIECES
        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            typename ReconstructionCache<Sample>::Piece *p = this->pieces[cIdx][before];
            intptr_t d = this->pieces[cIdx][after] - p;
            assert(d >= 0 && d <= 256);
        }
#endif
    }

    // Review - the [before] entries may be duplication of CandidateStash<Sample> buffer pointer??
    template <class Rectangular>
    void copyAfter(Rectangular const &rectangular, const StatePieces &other, int mask = 0x7)
    {
        other.keep = true;
        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            if (!(mask & (1 << cIdx))) continue;

            assert(other.parentPieces[cIdx] >= this->pieces[cIdx][before]);
            assert(other.parentPieces[cIdx] <= this->pieces[cIdx][after]);

            while (this->pieces[cIdx][after] != other.parentPieces[cIdx])
            {
                --this->pieces[cIdx][after];
                this->StatePieces<Sample>::zz[cIdx] -= 1 << (2 * this->pieces[cIdx][after]->log2Size);
                this->stateReconstructionCache->components[cIdx].freeBlock(*this->pieces[cIdx][after]);
            }

            typename ReconstructionCache<Sample>::Piece const *src = other.pieces[cIdx][before];
            typename ReconstructionCache<Sample>::Piece const *srcEnd = other.pieces[cIdx][after];

            while (src != srcEnd)
            {
                this->appendPiece(
                    cIdx,
#ifdef DEBUG_PIECES
                    src->rc,
#else
                    residual_coding(),
#endif
                    src->log2Size,
                    src->i);

                ++src;
            }
        }
    }

    // Check that there are the expected topology of pieces in the respective buffer.
    template <class Rectangular>
    void checkPieces(int cIdx, When when, const Rectangular &rectangular)
    {
#ifdef DEBUG_PIECES
        int z0 = cartesianToZ(xPositionOf(rectangular), yPositionOf(rectangular));
        if (when == before)
        {
            if (this->pieces[cIdx][before] != this->pieces[cIdx][after])
            {
                auto previous = this->pieces[cIdx][after][-1];
                int z = cartesianToZ(previous.rc.x0, previous.rc.y0);
                const int size = 1 << (previous.rc.log2TrafoSize + (cIdx ? 1 : 0));
                z += size * size;
                if (z != z0)
                {
                    std::cout << "expected position before " << rectangular << " actually at (" << previous.rc.x0 + size << ", " << previous.rc.y0 + size << ")\n";
                    ASSERT(0);
                }
            }
        }
        else
        {
            ASSERT(this->pieces[cIdx][before] != this->pieces[cIdx][after]);
            z0 += widthOf(rectangular) * heightOf(rectangular);
            auto previous = this->pieces[cIdx][after][-1];
            int z = cartesianToZ(previous.rc.x0, previous.rc.y0);
            const int size = 1 << (previous.rc.log2TrafoSize + (cIdx ? 1 : 0));
            z += size * size;
            if (z != z0)
            {
                std::cout << "expected position after " << rectangular << " actually at (" << previous.rc.x0 + size << ", " << previous.rc.y0 + size << ")\n";
                ASSERT(0);
            }
        }
        assert(this->zz[cIdx] == (z0 >> (cIdx ? 2 : 0)));
#endif
    }

    void appendPiece(int cIdx, const residual_coding &rc, int log2Size, int16_t i)
    {
#ifdef DEBUG_PIECES
        int z = zPositionOf(rc);
        ASSERT(z == (this->zz[cIdx] << (cIdx ? 2 : 0)));
        if (this->pieces[cIdx][after] != this->pieces[cIdx][before])
        {
            auto &previous = this->pieces[cIdx][after][-1];
            int zPrev = cartesianToZ(previous.rc.x0, previous.rc.y0);
            ASSERT(z == zPrev + (1 << (2 * previous.rc.log2TrafoSize + (previous.rc.cIdx ? 2 : 0))));
        }
        this->pieces[cIdx][after]->rc = rc;
#endif
        this->zz[cIdx] += 1 << (2 * log2Size);
        this->pieces[cIdx][after]->log2Size = log2Size;
        this->pieces[cIdx][after]->i = i;
        ++this->pieces[cIdx][after];
    }

    // Copy from list of cached reconstructed pieces to reconstructed picture buffer
    void commit(ThreePlanes<Sample> &picture, coding_quadtree cqt)
    {
        {
            auto p = this->pieces[0][before];
#ifdef DEBUG_PIECES
            ASSERT(p->rc.x0 == cqt.x0);
            ASSERT(p->rc.y0 == cqt.y0);
#endif
            this->stateReconstructionCache->components[0].commit(picture[0].offset(cqt.x0, cqt.y0), cqt.log2CbSize, p);
            ASSERT(p == this->pieces[0][after]);
        }

        {
            auto p = this->pieces[1][before];
#ifdef DEBUG_PIECES
            ASSERT(p->rc.x0 == cqt.x0);
            ASSERT(p->rc.y0 == cqt.y0);
#endif
            this->stateReconstructionCache->components[1].commit(picture[1].offset(cqt.x0 / 2, cqt.y0 / 2), cqt.log2CbSize - 1, p);
            ASSERT(p == this->pieces[1][after]);
        }

        {
            auto p = this->pieces[2][before];
#ifdef DEBUG_PIECES
            ASSERT(p->rc.x0 == cqt.x0);
            ASSERT(p->rc.y0 == cqt.y0);
#endif
            this->stateReconstructionCache->components[2].commit(picture[2].offset(cqt.x0 / 2, cqt.y0 / 2), cqt.log2CbSize - 1, p);
            ASSERT(p == this->pieces[2][after]);
        }
    }

    void resetPieces()
    {
        this->pieces[0][before] = this->pieces[0][after] = 0;
        this->pieces[1][before] = this->pieces[1][after] = 0;
        this->pieces[2][before] = this->pieces[2][after] = 0;
    }

    // Returns a pointer to the first entry in the specified color component's array of entries
    typename ReconstructionCache<Sample>::Piece *firstPiece(int cIdx, int x, int y)
    {
#ifdef DEBUG_PIECES
        ASSERT(x == this->pieces[cIdx][before]->rc.x0);
        ASSERT(y == this->pieces[cIdx][before]->rc.y0);
#endif
        return this->pieces[cIdx][before];
    }

    // review: why mutable? So that keep can be modified on a const instance?
    bool mutable keep;

    StatePieces() : keep(false) { }

    ~StatePieces()
    {
        if (!this->keep)
            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                typename ReconstructionCache<Sample>::Piece *p = this->pieces[cIdx][before];
                intptr_t d = this->pieces[cIdx][after] - p;
                ASSERT(d >= 0 && d <= 256);
                while (p != this->pieces[cIdx][after])
                    this->stateReconstructionCache->components[cIdx].freeBlock(*p++);
            }
    }

    typename ReconstructionCache<Sample>::Piece *pieces[3/* cIdx */][2 /* after */];
    typename ReconstructionCache<Sample>::Piece *parentPieces[3/* cIdx */];

    // current Z-order address
    int zz[3/* cIdx */];
};

// Review: Candidate<Sample> deserves its own header file

// Object containing information that is discarded or kept depending on the outcome of a RDO decision
// Typically the encoder will maintain two or more Candidates pending a decision. It will keep or copy the
// Candidate<Sample> representing the lowest RDO cost and discard the other.
// Candidate<Sample> contains the following information:
//   * CABAC contexts' state
//   * Neighbourhood information (Snakes)
//   * An array of Piece objects providing access to cached reconstructed blocks for later assembly of reconstructed CTU.
template <typename Sample>
struct Candidate :
    ContextsAndCost,
    NeighbourhoodEnc<Sample>,
    StateCodedData,
    StatePieces<Sample>
    {
        template <class Rectangular>
        void copy(const Candidate<Sample> &other, When when, Rectangular const &rectangular, bool merge, bool blockData)
        {
            if (blockData)
            {
                this->snake.copyBlockFrom(other.Neighbourhood::snake, when, rectangular, this->MinCbLog2SizeYMinus1, 1, 1);
            }

            if (merge)
            {
                this->snakeMerge.copyBlockFrom(other.snakeMerge, when, rectangular, this->MinCbLog2SizeYMinus1, 1, 1);
            }

            static_cast<Snake<BlockData>::Cursor &>(*this) = static_cast<const Snake<BlockData>::Cursor &>(other);
        }

        template <class Rectangular>
        void copyIntraReferenceSamplesLuma(const Candidate<Sample> &other, When when, Rectangular const &rectangular)
        {
            int intraReferencePadding = std::min(widthOf(rectangular), 32); // review: check
            this->snakeIntraReferenceSamples[0].copyBlockFrom(other.snakeIntraReferenceSamples[0], when, rectangular, 0, intraReferencePadding, intraReferencePadding);
        }

        template <class Rectangular>
        void copyIntraReferenceSamplesChroma(const Candidate<Sample> &other, When when, Rectangular const &rectangular)
        {
            int intraReferencePadding = std::min(widthOf(rectangular) >> 1, 32); // review: check
            this->snakeIntraReferenceSamples[1].copyBlockFrom(other.snakeIntraReferenceSamples[1], when, rectangular, 1, intraReferencePadding, intraReferencePadding);
            this->snakeIntraReferenceSamples[2].copyBlockFrom(other.snakeIntraReferenceSamples[2], when, rectangular, 1, intraReferencePadding, intraReferencePadding);
        }

        Candidate() { }

public:
        int checkParent;
        int check;
        int rcudepthstatus;
        int noresidual;
        int rqtdepth;

        template <class T>
        Candidate<Sample> &operator<<(const T& t)
        {
#ifdef _DEBUG
            this->name << t;
#endif
            return *this;
        }
        int32_t sadResidueQuad[2 /* y */][2 /* x */];

private:
#ifdef _DEBUG
        std::ostringstream name;
#endif

private:
        Candidate<Sample> const &operator=(Candidate<Sample> const &) = delete;
        Candidate(Candidate<Sample> const &) = delete;
    };


// CandidateStash<Sample> is a Candidate<Sample> with its own local memories for neighbourhood information and reconstruction entries.
template <typename Sample>
struct CandidateStash :
    Candidate<Sample>
    {
        CandidateStash()
        {
            // prevent destructor causing havoc
            this->keep = true;
        }

        template <class Rectangular>
        CandidateStash(const Candidate<Sample> &other, Rectangular const & rectangular, StateReconstructionCache<Sample> &stateReconstructionCache)
        {
            this->stateReconstructionCache = &stateReconstructionCache;

            this->MinCbLog2SizeYMinus1 = other.MinCbLog2SizeYMinus1;

            this->snakeVector.resize(rectangular, this->MinCbLog2SizeYMinus1);
            this->snakeVectorMerge.resize(rectangular, this->MinCbLog2SizeYMinus1);
            this->snakeVectorIntraReferenceSamples[0].resize(rectangular, 0);
            this->snakeVectorIntraReferenceSamples[1].resize(rectangular, 1);
            this->snakeVectorIntraReferenceSamples[2].resize(rectangular, 1);

            this->snake = this->snakeVector;
            this->snakeMerge = this->snakeVectorMerge;
            this->snakeIntraReferenceSamples[0] = this->snakeVectorIntraReferenceSamples[0];
            this->snakeIntraReferenceSamples[1] = this->snakeVectorIntraReferenceSamples[1];
            this->snakeIntraReferenceSamples[2] = this->snakeVectorIntraReferenceSamples[2];

            static_cast<Snake<BlockData>::Pointer &>(*this) = this->snake;
            this->Snake<BlockData>::Cursor::p = this->snake.origin;
        }

        void resetPieces()
        {
            this->pieces[0][before] = this->pieces[0][after] = this->entryArray[0];
            this->pieces[1][before] = this->pieces[1][after] = this->entryArray[1];
            this->pieces[2][before] = this->pieces[2][after] = this->entryArray[2];

            this->StateCodedData::reset(this->codedDataArray);
        }

    private:
        // This lot ends up on the stack: review - consider ways of being more cache friendly.
        typename Snake<BlockData>::Array<16, 16, 1, 1> snakeVector;
        typename Snake<BlockData>::Array<16, 16, 1, 1> snakeVectorMerge;
        typename Snake<Sample>::template Array<64, 64, 32, 32> snakeVectorIntraReferenceSamples[3];
        typename ReconstructionCache<Sample>::Piece entryArray[3][256];
        CodedData::Type codedDataArray[(16*8*8 + 64*64) * 3 / 2];
    };


struct StateEncodeSubstreamBase :
    StateSubstream,
    CabacWriter,
    ValueConst<constrained_intra_pred_flag, 0>,
    ValueConst<cross_component_prediction_enabled_flag, 0>,
    ValueConst<ChromaArrayType, 1>,
    ValueCache<MaxTbLog2SizeY>,
    ValueCache<MinTbLog2SizeY>
{
    template <class H>
    StateEncodeSubstreamBase(H &h, std::vector<uint8_t> &buffer) :
        CabacWriter(buffer),
        StateSubstream(h),
        ValueConst<constrained_intra_pred_flag, 0>(h),
        ValueConst<cross_component_prediction_enabled_flag, 0>(h),
        ValueConst<::ChromaArrayType, 1>(h),
        ValueCache<MaxTbLog2SizeY>(h),
        ValueCache<MinTbLog2SizeY>(h),
        residual(64, 64, 1, 0, 0, 32)
    {
        this->mvPreviousInteger2Nx2N[0] = { 0, 0 };
        this->mvPreviousInteger2Nx2N[1] = { 0, 0 };
    }

    Picture<int16_t> residual;

    //coding_quadtree cqt;

    int ctxSet;
    int cRiceParam;

    // sum of square differences for each colour component
    // review: perhaps just need luma and chroma here
    // review: rename to "distortion" so this can also be used for SATD?
    int32_t ssd[3 /* cidx */];
    int32_t ssdPrediction[3 /* cIdx */];
    uint32_t satd;

    // best vector from most recent integer motion esarch
    MotionVector mvPreviousInteger2Nx2N[2 /* refList */];

    Cost costMvdZero[2 /* refList */][2 /* mvp_lX_flag */];
};


template <typename Sample>
struct StateEncodeSubstream :
    StateEncodeSubstreamBase,
    StateReconstructionCache<Sample>
    {
        template <class H>
        StateEncodeSubstream(H &h, std::vector<uint8_t> &buffer) :
        StateEncodeSubstreamBase(h, buffer),
        StateReconstructionCache<Sample>(h, 16)
        {
        }

        Raster<Sample> reconstructed0[3 /* cIdx */];

        // represents state before coding the current CU
        Candidate<Sample> *originalCandidate;

        typename ReconstructionCache<Sample>::Piece interPieces[3/* cIdx */][3 /* cbf */];

        IntraReferenceSamples<Sample> filtered;
        IntraReferenceSamples<Sample> unfiltered[3];
    };


struct StateEncodePictureSei :
    StateSei,
    ValueHolder<payload_extension_present>,
    Hrd,
    PicTiming,
    DecodedPictureHash,
    ActiveParameterSets2,
    AlternativeTransferCharacteristics,
    UserDataUnregistered,
    MasteringDisplayColourVolume
    {
    };


struct StateEncodePicture :
    AccessOperators<StateEncodePicture>,
    StatePicture,
    StateWavefront,
    Strps,
    NalWriter,
    StateEncodePictureSei
    {
        using AccessOperators<StateEncodePicture>::operator[];

        StateEncodePicture()
        {
        }

        StateEncodePicture(std::shared_ptr<InputQueue::Docket> docket)
        :
            docket(docket)
        {
        }

        template <typename Sample, class H>
        void resize(H &h)
        {
            this->StateWavefront::resize(h);

            Turing::Rectangle rectangle;
            rectangle.x0 = 0;
            rectangle.y0 = 0;
            rectangle.width = h[PicWidthInCtbsY()];
            rectangle.height = h[PicHeightInCtbsY()];
            rectangle.width <<= h[CtbLog2SizeY()];
            rectangle.height <<= h[CtbLog2SizeY()];

            if (std::is_same<Sample, uint16_t>::value)
            {
                this->StateSpatial::snakeVectorIntraReferenceSamples16[0].resize(rectangle, 0);
                this->StateSpatial::snakeVectorIntraReferenceSamples16[1].resize(rectangle, 1);
                this->StateSpatial::snakeVectorIntraReferenceSamples16[2].resize(rectangle, 1);

                auto *p = new NeighbourhoodEnc<uint16_t>;
                p->snakeIntraReferenceSamples[0] = this->snakeVectorIntraReferenceSamples16[0];
                p->snakeIntraReferenceSamples[1] = this->snakeVectorIntraReferenceSamples16[1];
                p->snakeIntraReferenceSamples[2] = this->snakeVectorIntraReferenceSamples16[2];
                this->neighbourhood2.reset(p);
            }
            else
            {
                this->StateSpatial::snakeVectorIntraReferenceSamples8[0].resize(rectangle, 0);
                this->StateSpatial::snakeVectorIntraReferenceSamples8[1].resize(rectangle, 1);
                this->StateSpatial::snakeVectorIntraReferenceSamples8[2].resize(rectangle, 1);

                auto *p = new NeighbourhoodEnc<uint8_t>;
                p->snakeIntraReferenceSamples[0] = this->snakeVectorIntraReferenceSamples8[0];
                p->snakeIntraReferenceSamples[1] = this->snakeVectorIntraReferenceSamples8[1];
                p->snakeIntraReferenceSamples[2] = this->snakeVectorIntraReferenceSamples8[2];
                this->neighbourhood2.reset(p);
            }

            this->StateSpatial::resize(h, *this->neighbourhood2);

            const int nSubstreams =
                    h[entropy_coding_sync_enabled_flag()]
                      ? h[PicHeightInCtbsY()]
                          : 1;

            this->substreams.resize(nSubstreams);
        }

        std::vector<std::vector<uint8_t>> substreams;
        std::shared_ptr<InputQueue::Docket> docket;

        double qpFactor;
        double lambda;
        Lambda reciprocalLambda;
        double reciprocalSqrtLambda;

        std::shared_ptr<Neighbourhood> neighbourhood2;

    };


template <class Sample>
struct StateEncodePicture2 :
    StateEncodePicture,
    StateReconstructedPicture<Sample>
    {
        using StateEncodePicture::StateEncodePicture;
    };


struct BackupPredictionInfo
{
    template <class Action, class V>
    void perform(Neighbourhood &neighbourhood, const V &v)
    {
        this->backup.perform<Action>(neighbourhood.snake, v, neighbourhood.MinCbLog2SizeYMinus1);
        this->backupMerge.perform<Action>(neighbourhood.snakeMerge, v, neighbourhood.MinCbLog2SizeYMinus1);
    }

    Snake<BlockData>::Backup<16, 16> backup;
    Snake<BlockData>::Backup<16, 16> backupMerge;
};


struct PsnrAnalysis
{
    struct PsnrFrame
    {
        double currentPsnr[3] = {0.0, 0.0, 0.0};
    };

    PsnrAnalysis(int bitDepth=0)
    {
        this->bitDepth = bitDepth;
        this->pictures = 0;
    }

    template <typename Sample>
    void analyse(int poc, Picture<Sample> &picture, Picture<Sample> &reference)
    {
        std::unique_lock<mutex> lockFrameDistortion(frameDistortionToken);
        assert(psnrOnFrameBasis.find(poc) == psnrOnFrameBasis.end());
        double pictureSse[4] = { 0.0, 0.0, 0.0, 0.0 };
        double pictureSamples[4] = { 0.0, 0.0, 0.0, 0.0 };
        PsnrFrame currentDistortion;

        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            for (int y = 0; y < picture[cIdx].height; ++y)
            {
                for (int x = 0; x < picture[cIdx].width; ++x)
                {
                    const int diff = picture[cIdx](x, y) - reference[cIdx](x, y);
                    pictureSse[cIdx] += diff * diff;
                }
            }

            pictureSamples[cIdx] = picture[cIdx].width * picture[cIdx].height;
            this->sse[cIdx] += pictureSse[cIdx];
            this->samples[cIdx] += pictureSamples[cIdx];

            this->sse[3] += pictureSse[cIdx];
            this->samples[3] += pictureSamples[cIdx];

            pictureSse[3] += pictureSse[cIdx];
            pictureSamples[3] += pictureSamples[cIdx];
        }

        for (int cIdx = 0; cIdx < 4; ++cIdx)
        {
            const double max = (1 << this->bitDepth) - 1;
            double mse = pictureSse[cIdx] / pictureSamples[cIdx];
            double psnr = 999.0;
            if (mse > 0.0) psnr = 10.0 * log10(max * max / mse);
            this->sumPsnr[cIdx] += psnr;
            if(cIdx < 3) currentDistortion.currentPsnr[cIdx] = psnr;
        }

        ++this->pictures;
        psnrOnFrameBasis[poc] = currentDistortion;
    }
    void report(std::ostream &os)
    {
        std::unique_lock<mutex> lockFrameDistortion(frameDistortionToken);
        std::ostringstream o;

        o << "PSNR report (dB)\n";

        o << std::setprecision(4);
        o << std::fixed;

        o << "PSNR overall      [Y Cb Cr all]: ";
        for (int cIdx = 0; cIdx < 4; ++cIdx)
        {
            const double max = (1 << this->bitDepth) - 1;
            double mse = this->sse[cIdx] / this->samples[cIdx];
            double psnr = 999.0;
            if (mse > 0.0) psnr = 10.0 * log10(max * max / mse);

            o << " " << psnr;
        }
        o << "\n";

        o << "PSNR picture mean [Y Cb Cr all]: ";
        for (int cIdx = 0; cIdx < 4; ++cIdx)
        {
            auto const psnr = this->sumPsnr[cIdx] / this->pictures;
            o << " " << psnr;
        }
        o << "\n";

        os << o.str();
    }
    void getPsnrFrameData(int poc, double &psnrY, double &psnrU, double &psnrV)
    {
        std::unique_lock<mutex> lockFrameDistortion(frameDistortionToken);
        auto currentFrameDistortion = psnrOnFrameBasis.find(poc);
        assert(currentFrameDistortion != psnrOnFrameBasis.end());
        psnrY = currentFrameDistortion->second.currentPsnr[0];
        psnrU = currentFrameDistortion->second.currentPsnr[1];
        psnrV = currentFrameDistortion->second.currentPsnr[2];
    }
    void removePsnrFrameData(int poc)
    {
        std::unique_lock<mutex> lockFrameDistortion(frameDistortionToken);
        auto currentFrameDistortion = psnrOnFrameBasis.find(poc);
        assert(currentFrameDistortion != psnrOnFrameBasis.end());
        psnrOnFrameBasis.erase(poc);
    }
    int bitDepth;
    int pictures;
    double sse[4]         = {0.0, 0.0, 0.0, 0.0};
    double samples[4]     = {0.0, 0.0, 0.0, 0.0};
    double sumPsnr[4]     = {0.0, 0.0, 0.0, 0.0};

    std::map<int, PsnrFrame> psnrOnFrameBasis;
    std::mutex frameDistortionToken;
};


struct StateEncode :
    Speed,
    ThreadPool,
    InputQueue,
    StateSequence,
    StateFunctionTables,
    StateWriteUserDataUnregistered,
    ValueHolder<rbsp_byte>,
    ValueHolder<more_data_in_byte_stream>,
    ValueHolder<more_rbsp_data>,
    ValueHolder<more_rbsp_trailing_data>
    {
        StateEncode(const boost::program_options::variables_map &vm) :
            InputQueue(vm["max-gop-n"].as<int>(), vm["max-gop-m"].as<int>(),
                       vm["field-coding"].as<bool>(),
                       vm["shot-change"].as<bool>(),
                       vm["segment"].as<int>(),
                       vm["qp"].as<int>()),
            ThreadPool((vm["no-parallel-processing"].as<bool>())? 1 : vm["threads"].as<int>()),
            StateFunctionTables(true,
#ifdef VALGRIND_FRIENDLY
            0
#else
            vm["asm"].as<int>()
#endif
            ? havoc_instruction_set_support() : havoc_instruction_set(HAVOC_C_OPT | HAVOC_C_REF)),
            userDataUnregSeiWritten(false)
        {
        }

        StateEncode(const StateEncode &);

        static void writePicture(std::ostream &o, StatePicture &statePicture, int bpp)
        {
            if (bpp == 8)
            {
                StateReconstructedPicture<uint8_t> &dp = static_cast<StateEncodePicture2<uint8_t> &>(statePicture);
                o << *dp.picture;
            }
            else if (bpp == 16)
            {
                StateReconstructedPicture<uint16_t> &dp = static_cast<StateEncodePicture2<uint16_t> &>(statePicture);
                o << *dp.picture;
            }
            else
            {
                assert(!"unexpected bpp");
            }
        }

        // purge reconstructed pictures, in output order
        void decodedPicturesFlush(int internalSampleSize, int externalSampleSize)
        {
            if (!fieldcoding)
            {
                while (!this->decodedPictures.empty() && this->decodedPictures.front()->reconstructed)
                {
                    if (this->fileOutYuvPictures)
                        writePicture(this->fileOutYuvPictures, *this->decodedPictures.front(), internalSampleSize);

                    if (this->fileOutYuvFrames)
                        writePicture(this->fileOutYuvFrames, *this->decodedPictures.front(), internalSampleSize);

                    this->decodedPictures.pop_front();
                }
            }
            else
            {
                if (!this->decodedPictures.empty() && (this->decodedPictures.size() % 2 == 0))
                {
                    auto hstate = *this->decodedPictures.front();
                    int size = static_cast<int>(this->decodedPictures.size());
                    if (internalSampleSize == 8 && externalSampleSize == 8)
                    {
                        for (int i = 0; i < size; i += 2)
                        {
                            StateReconstructedPicture<uint8_t> &dptop = static_cast<StateEncodePicture2<uint8_t> &>(*this->decodedPictures[0]);
                            StateReconstructedPicture<uint8_t> &dpbottom = static_cast<StateEncodePicture2<uint8_t> &>(*this->decodedPictures[1]);

                            if (this->decodedPictures[0]->reconstructed && this->decodedPictures[1]->reconstructed)
                            {

                                if (this->fileOutYuvFrames)
                                {
                                    ThreePlanes<uint8_t> conformanceWindowTop = { *dptop.picture,
                                        hstate[SubWidthC()] * hstate[conf_win_left_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_top_offset()],
                                        hstate[SubWidthC()] * hstate[conf_win_right_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_bottom_offset()] };
                                    ThreePlanes<uint8_t> conformanceWindowBottom = { *dpbottom.picture,
                                        hstate[SubWidthC()] * hstate[conf_win_left_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_top_offset()],
                                        hstate[SubWidthC()] * hstate[conf_win_right_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_bottom_offset()] };

                                    auto &top = conformanceWindowTop;
                                    auto &bottom = conformanceWindowBottom;

                                    writeFields(this->fileOutYuvFrames, top, bottom);
                                }

                                if (this->fileOutYuvPictures)
                                    writePicture(this->fileOutYuvPictures, *this->decodedPictures.front(), internalSampleSize);

                                this->decodedPictures.pop_front();

                                if (this->fileOutYuvPictures)
                                    writePicture(this->fileOutYuvPictures, *this->decodedPictures.front(), internalSampleSize);

                                this->decodedPictures.pop_front();
                            }
                        }
                    }
                    else
                    {
                        for (int i = 0; i < size; i += 2)
                        {
                            StateReconstructedPicture<uint16_t> &dptop = static_cast<StateEncodePicture2<uint16_t> &>(*this->decodedPictures[0]);
                            StateReconstructedPicture<uint16_t> &dpbottom = static_cast<StateEncodePicture2<uint16_t> &>(*this->decodedPictures[1]);

                            if (this->decodedPictures[0]->reconstructed && this->decodedPictures[1]->reconstructed)
                            {

                                if (this->fileOutYuvFrames)
                                {
                                    ThreePlanes<uint16_t> conformanceWindowTop = { *dptop.picture,
                                        hstate[SubWidthC()] * hstate[conf_win_left_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_top_offset()],
                                        hstate[SubWidthC()] * hstate[conf_win_right_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_bottom_offset()] };
                                    ThreePlanes<uint16_t> conformanceWindowBottom = { *dpbottom.picture,
                                        hstate[SubWidthC()] * hstate[conf_win_left_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_top_offset()],
                                        hstate[SubWidthC()] * hstate[conf_win_right_offset()],
                                        hstate[SubHeightC()] * hstate[conf_win_bottom_offset()] };

                                    auto &top = conformanceWindowTop;
                                    auto &bottom = conformanceWindowBottom;
                                    writeFields(this->fileOutYuvFrames, top, bottom);
                                }

                                if (this->fileOutYuvPictures)
                                    writePicture(this->fileOutYuvPictures, *this->decodedPictures.front(), internalSampleSize);

                                this->decodedPictures.pop_front();

                                if (this->fileOutYuvPictures)
                                    writePicture(this->fileOutYuvPictures, *this->decodedPictures.front(), internalSampleSize);

                                this->decodedPictures.pop_front();
                            }
                        }
                    }
                }
            }
        }
        bool scd;
        bool fieldcoding;
        bool sao;
        bool framedoubling;
        bool saoslow;
        bool amp;
        bool smp;
        bool nosmp;
        bool ecu;
        bool fdm;
        bool fdam;
        bool esd;
        bool cfm;
        bool met;
        bool rdoq;
        bool rcudepth;
        bool rqt;
        bool sdh;
        bool tskip;
        bool aps;
        bool wpp;
        bool decodedHashSei;
        int verbosity;
        int gopM;
        int baseQp;
        int maxnummergecand;
        size_t concurrentFrames;
        int internalbitdepth;
        int externalbitdepth;
        HashType hashType;
        bool useAq;
        bool useRateControl;
        int preferredTransferCharacteristics;
        bool userDataUnregSeiWritten;
        int userDataUnregMsgLen;
        bool repeatHeaders;
        bool masteringDisplayInfoPresent;

        DecodedPictureHash decodedPictureHash;
        std::ofstream fileOutYuvPictures;
        std::ofstream fileOutYuvFrames;

        //std::unique_ptr<SequenceController> rateControlEngine;
        std::shared_ptr<RateControlParameters> rateControlParams;
        map<int, SequenceController*> rateControlMap;
        std::unique_ptr<PsnrAnalysis> psnrAnalysis;
        bool enableProfiler;

        struct Response
        {
            bool eos;
            bool hungry;
            std::shared_ptr<StateEncodePicture> picture;
            bool done;
            bool keyframe; //< true if IDR, CDR or BLA picture (for FFmpeg etc.)
        };

        struct FrameHash
        {
            std::vector<int> hash;
        };

        struct MasteringDisplayInfo
        {
            uint16_t displayPrimariesX[3];
            uint16_t displayPrimariesY[3];
            uint16_t whitePointX;
            uint16_t whitePointY;
            uint32_t maxDisplayMasteringLuma;
            uint32_t minDisplayMasteringLuma;
        };

        // Structure to store the information about the mastering display colour volume as specified in ST.2086
        MasteringDisplayInfo masterDisplayInfo;

        // Map to store hashes on a frame basis, associated mutex and helper functions
        std::map<int, FrameHash> hashOnFrameBasis;
        std::mutex hashFrameToken;
        void addFrameHash(int poc, FrameHash &element)
        {
            std::unique_lock<std::mutex> lockFrameBasisMap(hashFrameToken);
            assert(hashOnFrameBasis.find(poc) == hashOnFrameBasis.end());
            hashOnFrameBasis[poc] = element;
        }
        void removeFrameHash(int poc)
        {
            std::unique_lock<std::mutex> lockFrameBasisMap(hashFrameToken);
            assert(hashOnFrameBasis.find(poc) != hashOnFrameBasis.end());
            hashOnFrameBasis.erase(poc);
        }
        FrameHash& getFrameHash(int poc)
        {
            std::unique_lock<std::mutex> lockFrameBasisMap(hashFrameToken);
            assert(hashOnFrameBasis.find(poc) != hashOnFrameBasis.end());
            auto currentFrameHash = hashOnFrameBasis.find(poc);
            return currentFrameHash->second;
        }

        // pictures currently being encoded (bitstream order)
        std::deque<Response> responses;
        std::condition_variable responsesAvailable;

        // reconstructed pictures ready for output (display order)
        std::deque<std::shared_ptr<StatePicture>> decodedPictures;
    };


template <typename Sample>
struct Access<ReconstructedSamples, StateEncodeSubstream<Sample>>
{
    typedef Raster<Sample> Type;
    static Type get(ReconstructedSamples e, StateEncodeSubstream<Sample> &s)
    {
        return s.reconstructed0[e.cIdx].offset((e.x) >> (e.cIdx ? 1 : 0), (e.y) >> (e.cIdx ? 1 : 0));
    }
};


template <class S>
struct Access<intra_chroma_pred_mode, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef ValueType<intra_chroma_pred_mode>::Type Type;
    static Type get(intra_chroma_pred_mode, StateCodedData &s)
    {
        return s.codedCu.word0().intra_chroma_pred_mode;
    }
};


template <class S>
struct Access<Mvd, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef ValueType<Mvd>::Type &Type;
    static Type get(Mvd mvd, StateCodedData &s)
    {
        return s.codedPu.mvd(mvd.refList);
    }
};


template <class S>
struct Access<inter_pred_idc, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef ValueType<inter_pred_idc>::Type Type;
    static Type get(inter_pred_idc, StateCodedData &s)
    {
        if (!s.codedPu.word0().metadata[L0].predFlag)
        {
            return PRED_L1;
        }
        else if (!s.codedPu.word0().metadata[L1].predFlag)
        {
            return PRED_L0;
        }
        else
        {
            return PRED_BI;
        }
    }
};


template <class S>
struct Access<ref_idx_l0, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef ValueType<ref_idx_l0>::Type Type;
    static Type get(ref_idx_l0 e, StateCodedData &s)
    {
        return s.codedPu.word0().metadata[L0].ref_idx_lX;
    }
};


template <class S>
struct Access<ref_idx_l1, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef ValueType<ref_idx_l1>::Type Type;
    static Type get(ref_idx_l1 e, StateCodedData &s)
    {
        return s.codedPu.word0().metadata[L1].ref_idx_lX;
    }
};


template <class S>
struct Access<mvp_l0_flag, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef ValueType<mvp_l0_flag>::Type Type;
    static Type get(mvp_l0_flag e, StateCodedData &s)
    {
        return s.codedPu.word0().metadata[L0].mvp_lX_flag;
    }
};


template <class S>
struct Access<mvp_l1_flag, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef ValueType<mvp_l1_flag>::Type Type;
    static Type get(mvp_l1_flag e, StateCodedData &s)
    {
        return s.codedPu.word0().metadata[L1].mvp_lX_flag;
    }
};


template <class S>
struct Access<part_mode, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef int Type;
    static Type get(part_mode v, StateCodedData &s)
    {
        return s.codedCu.word0().part_mode;
    }
};

#endif
