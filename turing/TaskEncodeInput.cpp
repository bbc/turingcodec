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

// Encoder frame input top-level task.

#include "TaskEncodeInput.h"
#include "Write.h"
#include "TaskEncodeSubstream.h"
#include "TaskDeblock.h"
#include "TaskSao.h"
#include "SyntaxNal.hpp"
#include "SyntaxRbsp.hpp"


template <class H>
TaskEncodeInput<H>::TaskEncodeInput(H &h)
:
h(h)
{
}


template <class H>
bool TaskEncodeInput<H>::blocked()
{
    StateEncode *stateEncode = this->h;
    if (stateEncode->responses.size() < stateEncode->concurrentFrames)
    {
        // The number of frames presently being processed has not reached the maximum.
        return false;
    }


    return true;
}


template <class H>
void setupSliceHeader(H &h, const InputQueue::Docket *docket)
{
    StateEncode *stateEncode = h;
    h[nal_unit_type()] = docket->nut;
    h[slice_type()] = docket->sliceType;
    h[first_slice_segment_in_pic_flag()] = 1;
    h[slice_segment_address()] = 0;
    h[slice_pic_order_cnt_lsb()] = docket->poc & (h[MaxPicOrderCntLsb()] - 1);

    h[NumPositivePics(0)] = h[num_positive_pics()] = 0;
    h[NumNegativePics(0)] = h[num_negative_pics()] = 0;

    h[pic_output_flag()] = 1;

    h[slice_qp_delta()] = docket->qpOffset;
    h[five_minus_max_num_merge_cand()] = (5 - stateEncode->maxnummergecand);

    // review: this needs to be set with appropriate logic. E.g. set to 1 when mvd_l1_zero_flag is set on the first reference picture.
    h[collocated_from_l0_flag()] = 0;

    if (h[sample_adaptive_offset_enabled_flag()])
    {
        h[slice_sao_luma_flag()] = true;
        if (h[ChromaArrayType()] != 0)
            h[slice_sao_chroma_flag()] = true;
    }

    if (docket->poc)
    {
        // A std::set is sorted numerically but all values in docket->references.negative are negative so we need to iterate backwards in order to process closest references first
        int i = 0;
        for (auto delta = docket->references.negative.rbegin(); delta != docket->references.negative.rend(); ++delta)
        {
            auto const deltaPoc = *delta;
            h[DeltaPocS0(0, i)] = deltaPoc;
            h[delta_poc_s0_minus1(i)] = -(h[DeltaPocS0(0, i)] - (i ? h[DeltaPocS0(0, i - 1)] : 0)) - 1;
            ++i;
        }

        h[NumNegativePics(0)] = h[num_negative_pics()] = i;

        h[UsedByCurrPicS0(0, 0)] = h[used_by_curr_pic_s0_flag(0)] = (h[slice_type()] == I ? 0 : 1);
    }

    if (h[slice_type()] == B)
    {
        int i = 0;
        for (int deltaPoc : docket->references.positive)
        {
            h[DeltaPocS1(0, i)] = deltaPoc;
            h[delta_poc_s1_minus1(i)] = (h[DeltaPocS1(0, i)] - (i ? h[DeltaPocS1(0, i - 1)] : 0)) - 1;
            ++i;
        }

        h[NumPositivePics(0)] = h[num_positive_pics()] = i;

        if (h[NumPositivePics(0)] > 0)
        {
            h[UsedByCurrPicS1(0, 0)] = h[used_by_curr_pic_s1_flag(0)] = 1;
        }

        h[mvd_l1_zero_flag()] = 0;

        if (i == 0)
        {
            // There are no positive (future) references so copy negative (historic) references from L0 to L1
            h[DeltaPocS1(0, 0)] = h[DeltaPocS0(0, 0)];
            h[delta_poc_s1_minus1(0)] = h[delta_poc_s0_minus1(0)];
            h[mvd_l1_zero_flag()] = 1;
        }
    }

    h[slice_temporal_mvp_enabled_flag()] = 1;
}


template <class H>
template <typename Sample>
void TaskEncodeInput<H>::startPictureEncode(StateEncode::Response &response, std::shared_ptr<InputQueue::Docket> docket, H &hh)
{
    auto newPicture = new StateEncodePicture2<Sample>(docket);
    response.picture.reset(newPicture);
    StatePictures *statePictures = hh;

    auto h = hh.extend(&*newPicture);

    // Activate parameter sets
    h[Active<Vps>()] = h[Table<Vps>()][0];
    h[Active<Sps>()] = h[Table<Sps>()][0];
    h[Active<Pps>()] = h[Table<Pps>()][0];

    setupSliceHeader(h, response.picture->docket.get());
    // DPB state update - may bump pictures...
    statePictures->sliceHeaderDone(h);
    auto &picture = h[Concrete<StatePicture>()];
    StateEncode* stateEncode = h;
    setupStateReconstructedPicture(picture, h, stateEncode->saoslow);

    // Further DPB state update - may also bump pictures

    h(PictureDone());

    auto *p = new StateReconstructedPicture<Sample>;
    StateReconstructedPicture<Sample> *reconstructedPicture = h;
    p->picture = reconstructedPicture->picture;
    p->saoPicture = reconstructedPicture->saoPicture;
    p->deblockPicture = reconstructedPicture->deblockPicture;
    response.picture->reconstructedPicture.reset(p);

    // Allocate various picture memories
    response.picture->resize<Sample>(h);
    
    

    if (stateEncode->useRateControl)
    {
        stateEncode->rateControlParams->takeTokenOnRCLevel();
        bool isShotChange = response.picture->docket->isShotChange;
        int currentPictureLevel = response.picture->docket->sopLevel;
        int currentPoc = response.picture->docket->poc;
        int segmentPoc = response.picture->docket->segmentPoc;
        int sopSize = response.picture->docket->currentGopSize;
        int sopId = response.picture->docket->sopId;
        int pocInSop = response.picture->docket->pocInSop;

        if (h[slice_type()] == I)
        {
            if (docket->absolutePoc == docket->segmentPoc)
            {
                
                assert(stateEncode->rateControlMap.find(docket->segmentPoc) == stateEncode->rateControlMap.end());
                SequenceController *currentSequenceController = new SequenceController(stateEncode->rateControlParams);
                stateEncode->rateControlMap[docket->segmentPoc] = currentSequenceController;
            }

            stateEncode->rateControlMap.find(segmentPoc)->second->initNewIntraPeriod(response.picture->docket);
            stateEncode->rateControlMap.find(segmentPoc)->second->pictureRateAllocationIntra(response.picture->docket,stateEncode->rateControlParams->cpbInfo);
            stateEncode->rateControlMap.find(segmentPoc)->second->initNewSop(response.picture->docket);
        }
        else
        {
            //if (currentPictureLevel == 1)
            if(stateEncode->rateControlMap.find(segmentPoc)->second->bInitNewSop(response.picture->docket))
            {
                // New SOP starts, set the rate budget for GOP and this current picture
                stateEncode->rateControlMap.find(segmentPoc)->second->initNewSop(response.picture->docket);
            }
            stateEncode->rateControlMap.find(segmentPoc)->second->pictureRateAllocation(response.picture->docket,stateEncode->rateControlParams->cpbInfo);
        }

        // Compute lambda
        response.picture->lambda = stateEncode->rateControlMap.find(segmentPoc)->second->estimatePictureLambda(currentPoc);

        // Derive QP from lambda
        int currentQP = stateEncode->rateControlMap.find(segmentPoc)->second->deriveQpFromLambda(response.picture->lambda, h[slice_type()] == I, currentPictureLevel, h[PicOrderCntVal()]);

        response.picture->qpFactor = response.picture->docket->qpFactor;
        h[slice_qp_delta()] = currentQP - stateEncode->rateControlMap.find(segmentPoc)->second->getBaseQp();
        response.picture->reciprocalLambda.set(1.0 / response.picture->lambda);
        response.picture->reciprocalSqrtLambda = sqrt(1.0 / response.picture->lambda);
        stateEncode->rateControlParams->releaseTokenOnRCLevel();
    }
    else
    {
        // Lambda computations - review: better elsewhere?
        response.picture->qpFactor = response.picture->docket->qpFactor;
        response.picture->lambda = computeLambda(h);
        response.picture->reciprocalLambda.set(1.0 / response.picture->lambda);
        response.picture->reciprocalSqrtLambda = sqrt(1.0 / response.picture->lambda);
    }

    const int nCtusInFirstSubstream = h[entropy_coding_sync_enabled_flag()] ? h[PicWidthInCtbsY()] : h[PicSizeInCtbsY()];

    ThreadPool *threadPool = hh;

    {
        StateEncode *stateEncode = hh;

        std::unique_lock<std::mutex> lock(threadPool->mutex());
        stateEncode->responses.push_back(response);
        if (stateEncode->useRateControl)
        {
            int segmentPoc = response.picture->docket->segmentPoc;
            stateEncode->rateControlParams->takeTokenOnRCLevel();
            stateEncode->rateControlMap.find(segmentPoc)->second->decreaseNumLeftSameHierarchyLevel(response.picture->docket);
            stateEncode->rateControlParams->releaseTokenOnRCLevel();
        }
        stateEncode->responsesAvailable.notify_all();
    }

    // Enqueue first encoding task and deblocking tasks for theadpool execution
    threadPool->add(newTaskDeblock(h, response.picture, 0, nCtusInFirstSubstream));
    threadPool->add(newTaskSao(h, response.picture, 0, nCtusInFirstSubstream));
    threadPool->add(*new TaskEncodeSubstream<Sample>(h, response.picture, 0, nCtusInFirstSubstream));
}


template <class H>
bool TaskEncodeInput<H>::run()
{
    StateEncode::Response response;

    do
    {
        response.eos = false;
        response.hungry = false;
        response.done = false;
        response.picture.reset();

        ThreadPool *threadPool = this->h;
        StateEncode *stateEncode = this->h;
        InputQueue *inputQueue = this->h;

        threadPool->lock();
        inputQueue->preanalyse();
        if (this->blocked())
        {
            threadPool->unlock();
            return true;
        }

        auto docket = inputQueue->getDocket();
        const bool eos = inputQueue->eos();
        threadPool->unlock();

        //inputQueue->preanalyse();
        if (docket)
        {
            EstimateIntraComplexity *icInfo = new EstimateIntraComplexity();
            if (stateEncode->useRateControl)
            {
                if (docket->picture->sampleSize == 8)
                {
                    auto &pictureWrap = static_cast<PictureWrap<uint8_t> &>(*(docket->picture));
                    auto &orgPicture = static_cast<Picture<uint8_t> &>(pictureWrap);
                    int sliceType = docket->sliceType;

                    if (sliceType == 2) //INTRA
                    {
                        icInfo->allocate(orgPicture[0].height, orgPicture[0].width);
                        icInfo->preAnalysis<uint8_t>(docket->picture);
                    }
                }
                else
                {
                    auto &pictureWrap = static_cast<PictureWrap<uint16_t> &>(*(docket->picture));
                    auto &orgPicture = static_cast<Picture<uint16_t> &>(pictureWrap);
                    int sliceType = docket->sliceType;

                    if (sliceType == 2) //INTRA
                    {
                        icInfo->allocate(orgPicture[0].height, orgPicture[0].width);
                        icInfo->preAnalysis<uint16_t>(docket->picture);
                    }
                }
                docket->icInfo.reset(icInfo);
            }

            if (docket->picture->sampleSize == 16)
            {
                startPictureEncode<uint16_t>(response, docket, this->h);
            }
            else
            {
                startPictureEncode<uint8_t>(response, docket, this->h);
            }
        }
        else
        {
            if (eos)
            {
                response.eos = true;

                while (static_cast<StatePictures *>(h)->bumpingProcess(h, true));

                std::unique_lock<std::mutex> lock(threadPool->mutex());
                stateEncode->responses.push_back(response);
                stateEncode->responsesAvailable.notify_all();
            }
            else
            {
                response.hungry = true;
                response.done = true;

                std::unique_lock<std::mutex> lock(threadPool->mutex());
                stateEncode->responses.push_front(response);
                stateEncode->responsesAvailable.notify_all();
            }
        }

    } while (!response.eos);

    delete this;
    return false;
}



// Explict template instantiation
template struct TaskEncodeInput<struct Handler<Encode<void>, StateEncode>>;
