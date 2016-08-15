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

// Encoding top-level task.

#include "TaskEncodeSubstream.h"
#include "CabacWriter.h"
#include "StateEncode.h"
#include "Handlers.h"
#include "Write.h"
#include "SyntaxRbsp.hpp"
#include "SyntaxCtu.hpp"
#include "EstimateRate.h"
#include "RateControl.h"
#include "EstimateIntraComplexity.h"


template <typename Sample>
TaskEncodeSubstream<Sample>::TaskEncodeSubstream(PointerTuple<StateEncodePicture2<Sample>, StateEncode> pointers, std::shared_ptr<StateEncodePicture> stateEncodePicture, int begin, int end) :
    stateEncodeSubstream(pointers, static_cast<StateEncodePicture *>(pointers)->substreams[begin / pointers[PicWidthInCtbsY()]]),
    h(pointers.template change<Encode<void>>().extend(&this->stateEncodeSubstream).extend(&this->candidate)),
    begin(begin),
    end(end),
    isLastSubstream(end == h[PicSizeInCtbsY()]),
    stateEncodePicture(stateEncodePicture)
    {
        h[CtbAddrInRs()] = begin;
        h[CtbAddrInTs()] = h[CtbAddrRsToTs(begin)];
        this->encodingDone = false;

        static_cast<NeighbourhoodEnc<Sample>&>(this->candidate) = *static_cast<NeighbourhoodEnc<Sample>*>(stateEncodePicture->neighbourhood2.get());

        this->candidate.resetPieces();
    }


template <typename Sample>
bool TaskEncodeSubstream<Sample>::blocked()
{
    StateWavefront *stateWavefront = h;

    const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
    const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];

    if (rx != h[PicWidthInCtbsY()] - 1)
    {
        if (!stateWavefront->encoded(rx + 1, ry - 1))
        {
            // blocked by previous wavefront thread in this picture
            return true;
        }
    }

    // To encode the current CTU at (rx, ry), the CTU (rx+2, ry+1) of each reference picture must be ready for prediction.
    int depX = std::min(rx + 2, h[PicWidthInCtbsY()] - 1);
    int depY = std::min(ry + 1, h[PicHeightInCtbsY()] - 1);

    if (h[slice_type()] != I)
    {
        for (int rIdx = 0; rIdx <= h[num_ref_idx_l0_active_minus1()]; rIdx++)
        {
            StateEncodePicture *stateEncodePicture = dynamic_cast<StateEncodePicture *>(h[RefPicList(L0)][rIdx].dp.get());
            StateWavefront *stateWavefrontRef = stateEncodePicture;
            if (!stateWavefrontRef->deblocked(depX, depY)) return true;
        }
    }

    if (h[slice_type()] == B)
    {
        for (int rIdx = 0; rIdx <= h[num_ref_idx_l1_active_minus1()]; rIdx++)
        {
            StateEncodePicture *stateEncodePicture = dynamic_cast<StateEncodePicture *>(h[RefPicList(L1)][rIdx].dp.get());
            StateWavefront *stateWavefrontRef = stateEncodePicture;
            if (!stateWavefrontRef->deblocked(depX, depY)) return true;
        }
    }

    return false;
}


template <typename Sample>
bool TaskEncodeSubstream<Sample>::blockedLock()
{
    ThreadPool *threadPool = h;
    threadPool->lock();
    const bool blocked = this->blocked();
    threadPool->unlock();
    return blocked;
}


template <typename Sample>
bool TaskEncodeSubstream<Sample>::run()
{
    Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->encode);

    StateWavefront *stateWavefront = h;
    StateEncode *stateEncode = h;
    InputQueue *inputQueue = h;
    ThreadPool *threadPool = h;

    if(h[cu_qp_delta_enabled_flag()] && (h[CtbAddrInRs()] % h[PicWidthInCtbsY()] == 0))
    {
        static_cast<QpState *>(h)->initInternalMemory(h);
    }

    if(stateEncode->useRateControl && h[CtbAddrInRs()] == 0)
    {
        bool isShotChange = this->stateEncodePicture->docket->isShotChange;
        int currentPictureLevel = this->stateEncodePicture->docket->sopLevel;
        int sopSize = this->stateEncodePicture->docket->currentGopSize;
        if(sopSize > 1)
            stateEncode->rateControlEngine->setSopSize(sopSize);
        if(h[slice_type()] == I)
        {
            if(isShotChange)
            {
                std::cout<<"Resetting the rate control engine with POC: "<<this->stateEncodePicture->docket->poc<<"\n";
                stateEncode->rateControlEngine->reset();
            }
            EstimateIntraComplexity &icInfo = dynamic_cast<EstimateIntraComplexity&>(*static_cast<StateEncodePicture *>(h)->docket->icInfo);
            stateEncode->rateControlEngine->pictureRateAllocationIntra(icInfo);
            stateEncode->rateControlEngine->initNewSop();
        }
        else
        {
            if(currentPictureLevel == 1)
            {
                // New SOP starts, set the rate budget for GOP and this current picture
                stateEncode->rateControlEngine->initNewSop();
            }
            stateEncode->rateControlEngine->pictureRateAllocation(currentPictureLevel);
        }
        // Compute lambda
        this->stateEncodePicture->lambda = stateEncode->rateControlEngine->estimatePictureLambda(h[slice_type()] == I, currentPictureLevel);

        // Derive QP from lambda
        int currentQP = stateEncode->rateControlEngine->deriveQpFromLambda(this->stateEncodePicture->lambda, h[slice_type()] == I, currentPictureLevel);
        h[slice_qp_delta()] = currentQP - stateEncode->rateControlEngine->getBaseQp();
        this->stateEncodePicture->reciprocalLambda.set(1.0 / this->stateEncodePicture->lambda);
        this->stateEncodePicture->reciprocalSqrtLambda = sqrt(1.0 / this->stateEncodePicture->lambda);
        char data[100];
        sprintf(data, "| %06d | %10d | %9.2f | %4d |", this->stateEncodePicture->docket->poc, stateEncode->rateControlEngine->getPictureTargetBits(), this->stateEncodePicture->lambda, currentQP);
        stateEncode->rateControlEngine->writetoLogFile(data);
    }

    while (h[CtbAddrInRs()] != this->end)
    {
        const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
        const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];

        // Is this task blocked, waiting for previous substreams or other dependencies?
        if (this->blockedLock()) return true;

        if (!this->isLastSubstream && rx == 0)
        {
            assert(h[entropy_coding_sync_enabled_flag()]);

            // Create the next substream task.
            // Review, we could do this a little later (after iteration call and/or at CTU rx=1)
            threadPool->add(*new TaskEncodeSubstream(h, this->stateEncodePicture, this->end, this->end + h[PicWidthInCtbsY()]));
        }

        Syntax<slice_segment_data>::iteration(h);

        threadPool->lock();
        stateWavefront->encoded.set(rx, ry);
        threadPool->nudge();
        threadPool->unlock();
    }

    assert(h[CtbAddrInRs()] == this->end);

    // Perform end-of-substream tasks

    if (h[CtbAddrInRs()] == h[PicSizeInCtbsY()])
    {
        // End of picture (and slice segment)
        Syntax<rbsp_slice_segment_trailing_bits>::go(rbsp_slice_segment_trailing_bits(), h);

        this->stateEncodePicture->clearRpls();
    }

    BitWriter::insertEp3Bytes(*static_cast<CabacWriter *>(h)->data, 0);

    if(stateEncode->useRateControl && h[CtbAddrInRs()] == h[PicSizeInCtbsY()])
    {
        // Only set coding rate at the end of one frame
        BitWriter *bitWrite = h;
        size_t rate = (bitWrite->data->size()) << 3;
        stateEncode->rateControlEngine->setCodingBits(static_cast<int>(rate));
        char data[100];
        sprintf(data, " %10d |", (int)rate);
        stateEncode->rateControlEngine->writetoLogFile(data);
    }

    threadPool->lock();
    const int n = --stateWavefront->nSubstreamsUnfinished;
    threadPool->nudge();
    threadPool->unlock();

    delete this;
    return false;
}


// explict template instantiations
template struct TaskEncodeSubstream<uint8_t>;
template struct TaskEncodeSubstream<uint16_t>;
