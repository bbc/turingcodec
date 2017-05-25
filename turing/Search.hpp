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

 // Encoder "Search" functions - all encoder decisions are made by functions here.
 // For example, CQT topology, CU prediction, mode, CU partitioning, intra prediction modes, inter modes, motion vectors, etc.


#include "Search.h"
#include "Global.h"
#include "Write.h"
#include "Syntax.h"
#include "SyntaxCtu.hpp"
#include "SyntaxElements.h"
#include "Reconstruct.h"
#include "Profiler.h"
#include "EstimateRate.h"
#include "havoc/sad.h"
#include "Aps.h"


// Finds the best value of IntraPredModeY for a single intra partition
template <class H>
void searchIntraPartition(H &h, IntraPartition intraPartition)
{
    using Sample = typename SampleType<H>::Type;
    int constexpr bitDepth = 6 + 2 * sizeof(Sample);

    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    coding_quadtree *cqt = h;

    Candidate<Sample> **activeCandidate = h;
    Candidate<Sample> *originalCandidate = *activeCandidate;

    auto const log2PartitionSize = intraPartition.log2CbSize - intraPartition.split;
    if (log2PartitionSize != 6)
    {
        auto samples = originalCandidate->snakeIntraReferenceSamples[0].offset(xPositionOf(intraPartition), yPositionOf(intraPartition), 0);
        auto &unfiltered = stateEncodeSubstream->unfiltered[0];
        auto &filtered = stateEncodeSubstream->filtered;

        unfiltered.substituteFast(h, samples, residual_coding(xPositionOf(intraPartition), yPositionOf(intraPartition), log2PartitionSize, 0));
        if (log2PartitionSize != 2)
            filtered.filter(unfiltered, h[strong_intra_smoothing_enabled_flag()], bitDepth, 1 << log2PartitionSize);
    }

    {
        Snake<BlockData>::Cursor *cursor = originalCandidate;
        cursor->relocate(originalCandidate->snake, intraPartition, h[MinCbLog2SizeY()] - 1);
    }

    CandModeList &candModeList = *static_cast<CandModeList *>(h);
    candModeList.populate(0, h, xPositionOf(intraPartition), yPositionOf(intraPartition));

    // There are only three possible different values for the rate of IntraPredModeY:
    // A. prev_intra_luma_pred_flag = 1 and mpm_idx = 0 - IntraPredModeY == CandModeList[0]
    Cost rateA;
    rateA.set(0, 0); // this case is always forced forward to refinement, forcing rate and distortion to zero acheives this.

    // B. prev_intra_luma_pred_flag = 1 and mpm_idx = 1 or 2  - IntraPredModeY == CandModeList[1 or 2]
    Cost rateB;
    rateB.set(2, 0);  // rate of mpm_idx
    rateB += estimateRateNoContextUpdate(h, EncodeDecision<prev_intra_luma_pred_flag>(1, 0));

    // C. prev_intra_luma_pred_flag = 0 - IntraPredModeY derived from rem_intra_luma_pred_mode
    Cost rateC;
    rateC.set(5, 0); // rate of rem_intra_luma_pred_mode
    rateC += estimateRateNoContextUpdate(h, EncodeDecision<prev_intra_luma_pred_flag>(0, 0));

    // Costs are offset by -costC so that we can zero-initialise here
    Cost costs[35] = { 0 };
    costs[candModeList[0]] = rateA - rateC;
    costs[candModeList[1]] = rateB - rateC;
    costs[candModeList[2]] = rateB - rateC;

    // Make initial RD cost estimates of all possible modes using Hadamard SATD
    {
        Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->searchIntraSatd);

        StateReconstructionCache<Sample> *stateReconstructionCache = h;
        CandidateStash<Sample> candidate(*originalCandidate, intraPartition, *stateReconstructionCache);
        *activeCandidate = &candidate;

        Lambda lambda;
        lambda.set(getReciprocalSqrtLambda(h));
        StateEncode *stateEncode = h;
        if(stateEncode->useRateControl)
        {
            StateEncodePicture *stateEncodePicture = h;
            int segmentPoc = stateEncodePicture->docket->segmentPoc;
            stateEncode->rateControlParams->takeTokenOnRCLevel();
            double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalSqrtLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
            stateEncode->rateControlParams->releaseTokenOnRCLevel();
            lambda.set(value);
        }

        for (int n = 0; n < 35; ++n)
        {
            candidate.resetPieces();
            candidate.copy(*originalCandidate, before, intraPartition, false, true);
            candidate.copyIntraReferenceSamplesLuma(*originalCandidate, before, intraPartition);
            candidate.StatePieces<Sample>::copyBefore(intraPartition, *originalCandidate);
            candidate.StateCodedData::copyBefore(intraPartition, *originalCandidate);

            candidate.codedCu.IntraPredModeY(intraPartition.blkIdx) = n;

            {
                Neighbourhood *neighbourhood = h;
                Snake<BlockData>::Cursor *cursor = h;
                BlockData &blockData = cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1);
                blockData.setup(cqt, MODE_INTRA);
            }

            int32_t distortion;

            {
                auto &e = intraPartition;

                const int i = (e.blkIdx & 1) << (e.log2CbSize - 1);
                const int j = (e.blkIdx >> 1) << (e.log2CbSize - 1);

                distortion = predictIntraLuma(transform_tree(e.x0 + i, e.y0 + j, e.x0, e.y0, e.log2CbSize - e.split, e.split, e.blkIdx), h);
            }

            costs[n] += lambda * distortion;
        }

        *activeCandidate = originalCandidate;
    }

    Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->searchIntraRd);

    StateReconstructionCache<Sample> *stateReconstructionCache = h;

    CandidateStash<Sample> candidateStash[2];
    CandidateStash<Sample> *champion = &candidateStash[0];
    CandidateStash<Sample> *challenger = &candidateStash[1];
    champion->ContextsAndCost::setMax();

    const auto max = static_cast<Speed *>(h)->nCandidatesIntraRefinement(intraPartition.log2CbSize >> intraPartition.split);

    // Refine using holistic rate and sum of square differences distortion
    auto nMpm = 0;
    for (auto j = 0; j < max + nMpm; ++j)
    {
        // Find champion candidate from previous search
        int IntraPredModeY = 0;
        Cost costBest = costs[0];
        for (int i = 1; i < 35; ++i)
        {
            if (costs[i] < costBest)
            {
                costBest = costs[i];
                IntraPredModeY = i;
            }
        }

        // Mark this candidate as done
        costs[IntraPredModeY] = std::numeric_limits<Cost>::max();

        if (j == max - 1)
        {
            // We have already refined (max) candidates: look at the MPM candidates
            for (int i = 0; i < candModeList.neighbourModes; ++i)
            {
                if (costs[candModeList[i]] != std::numeric_limits<Cost>::max())
                {
                    // This candidate hasn't been searched yet, mark it for refinement.
                    costs[candModeList[i]] = Cost{ 0 };
                    ++nMpm;
                }
            }
        }

        challenger->~CandidateStash();
        new (challenger) CandidateStash<Sample>(*originalCandidate, intraPartition, *stateReconstructionCache);
        challenger->resetPieces();

        challenger->copy(*originalCandidate, before, intraPartition, false, true);
        challenger->copyIntraReferenceSamplesLuma(*originalCandidate, before, intraPartition);
        challenger->StatePieces<Sample>::copyBefore(intraPartition, *originalCandidate);
        challenger->StateCodedData::copyBefore(intraPartition, *originalCandidate);
        challenger->ContextsAndCost::copy(*originalCandidate);

        challenger->codedCu.IntraPredModeY(intraPartition.blkIdx) = IntraPredModeY;
        challenger->transformTree.p = challenger->codedDataAfter;

        *activeCandidate = challenger;

        {
            {
                Neighbourhood *neighbourhood = h;
                Snake<BlockData>::Cursor *cursor = h;
                BlockData &blockData = cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1);
                blockData.setup(cqt, MODE_INTRA);
            }

            StateCodedData *stateCodedData = challenger;
            CodedData::TransformTree transformTree = stateCodedData->transformTree;

            // Reconstruct output partition and measure distortion using SSD
            auto ssd = reconstructIntraLuma(intraPartition, h);

            Lambda lambda = getReciprocalLambda(h);
            StateEncode *stateEncode = h;
            if (stateEncode->useRateControl)
            {
                StateEncodePicture *stateEncodePicture = h;
                int segmentPoc = stateEncodePicture->docket->segmentPoc;
                stateEncode->rateControlParams->takeTokenOnRCLevel();
                double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                stateEncode->rateControlParams->releaseTokenOnRCLevel();
                lambda.set(value);
            }
            (*activeCandidate)->ContextsAndCost::lambdaDistortion += lambda * ssd;

            stateCodedData->codedDataAfter = stateCodedData->residual.p;
            stateCodedData->codedCu.chromaOffset() = CodedData::Type(stateCodedData->codedDataAfter - stateCodedData->codedCu.firstTransformTree().p);
            ASSERT(stateCodedData->codedCu.firstTransformTreeChroma().p == stateCodedData->codedDataAfter);
            stateCodedData->codedCu.firstTransformTreeChroma().word0().raw = 0;
            stateCodedData->transformTree = transformTree;
            stateCodedData->transformTreeChroma = transformTree;
            // Measure bits used to send the partition
            // review - could reuse IntraPredModeY rate calculated in initial search (but would need to update context of prev_intra_luma_pred_flag)

            if (challenger->cost2() < champion->cost2())
            {
                auto m = h.template change<EstimateRateLuma<void>>();
                m(intraPartition);
            }
        }

        if (challenger->cost2() < champion->cost2())
        {
            std::swap(challenger, champion);
        }
    }

    Snake<BlockData>::Cursor *cursor = champion;
    cursor->current(0, 0, champion->MinCbLog2SizeYMinus1).intra.predModeY = champion->codedCu.IntraPredModeY(intraPartition.blkIdx);

    champion->snake.commitRectangle(intraPartition, cursor->value, champion->MinCbLog2SizeYMinus1);

    originalCandidate->copy(*champion, after, intraPartition, false, true);
    originalCandidate->copyIntraReferenceSamplesLuma(*champion, after, intraPartition);
    originalCandidate->ContextsAndCost::copy(*champion);
    originalCandidate->StatePieces<Sample>::copyAfter(intraPartition, *champion);
    originalCandidate->StateCodedData::copyAfter(intraPartition, *champion);

    *activeCandidate = originalCandidate;
}


template <class H>
void searchIntraChroma(H &h, const coding_quadtree &cqt)
{
    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    StateReconstructionCache<Sample> *stateReconstructionCache = h;
    Candidate<Sample> **activeCandidate = h;

    coding_unit const cu(cqt.x0, cqt.y0, cqt.log2CbSize);

    Candidate<Sample> *originalCandidate = *activeCandidate;

    CandidateStash<Sample> candidateStash[2];
    CandidateStash<Sample> *champion = &candidateStash[0];
    CandidateStash<Sample> *challenger = &candidateStash[1];

    champion->ContextsAndCost::setMax();

    for (int cIdx = 1; cIdx < 3; ++cIdx)
    {
        auto samples = originalCandidate->snakeIntraReferenceSamples[cIdx].offset(cqt.x0, cqt.y0, 1);
        auto &unfiltered = stateEncodeSubstream->unfiltered[cIdx];
        unfiltered.substituteFast(h, samples, residual_coding(cqt.x0, cqt.y0, cqt.log2CbSize - 1, cIdx));
    }

    for (int n = 0; n < 5; ++n)
    {
        challenger->~CandidateStash();
        new (challenger) CandidateStash<Sample>(*originalCandidate, cqt, *stateReconstructionCache);

        challenger->resetPieces();

        challenger->copy(*originalCandidate, before, cqt, false, false);
        challenger->copyIntraReferenceSamplesChroma(*originalCandidate, before, cqt);
        challenger->ContextsAndCost::copy(*originalCandidate);
        challenger->StatePieces<Sample>::copyBefore(cqt, *originalCandidate);

        challenger->StateCodedData::copyBefore(cqt, *originalCandidate, -1);
        challenger->transformTreeChroma = challenger->codedCu.firstTransformTreeChroma();
        challenger->codedCu.word0().intra_chroma_pred_mode = n;

        challenger->forcePosition(before, cu);

        *activeCandidate = challenger;

        //if (cu == coding_unit(64, 64, 6))
        //{
        //	std::cout << "searchIntraChroma " << cu << " challenger = " << &challenger << "\n";
        //	std::cout << "codedCu.p = " << challenger->codedCu.p << "\n";
        //}

        auto ssd = reconstructIntraChroma(transform_tree(cu.x0, cu.y0, cu.x0, cu.y0, cu.log2CbSize, 0, 0), h);

        challenger->transformTreeChroma = challenger->codedCu.firstTransformTreeChroma();

        {
            // measure bits used to send the chroma
            auto m = h.template change<EstimateRateChroma<void>>();

            m(intra_chroma_pred_mode(cu.x0, cu.y0), ae(v));
            if (m[rqt_root_cbf()])
            {
                m[MaxTrafoDepth()] = (m[current(CuPredMode(cu.x0, cu.y0))] == MODE_INTRA)
                    ? m[max_transform_hierarchy_depth_intra()] + m[IntraSplitFlag()]
                    : m[max_transform_hierarchy_depth_inter()];
                m(transform_tree(cu.x0, cu.y0, cu.x0, cu.y0, cu.log2CbSize, 0, 0));
            }
        }

        ASSERT(challenger->codedDataAfter == challenger->residual.p);

        Lambda reciprocalLambda = getReciprocalLambda(h);
        StateEncode *stateEncode = h;
        if(stateEncode->useRateControl)
        {
            StateEncodePicture *stateEncodePicture = h;
            int segmentPoc = stateEncodePicture->docket->segmentPoc;
            stateEncode->rateControlParams->takeTokenOnRCLevel();
            double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
            stateEncode->rateControlParams->releaseTokenOnRCLevel();
            reciprocalLambda.set(value);
        }
        challenger->ContextsAndCost::lambdaDistortion += reciprocalLambda * ssd;

        if (challenger->cost2() < champion->cost2())
        {
            std::swap(challenger, champion);
        }
    }

    originalCandidate->copy(*champion, after, cu, false, false);
    originalCandidate->copyIntraReferenceSamplesChroma(*champion, after, cu);
    originalCandidate->ContextsAndCost::copy(*champion);
    originalCandidate->StatePieces<Sample>::copyAfter(cu, *champion, 0x6 /* (chroma only) */);
    originalCandidate->StateCodedData::copyAfter(cu, *champion);

    *activeCandidate = originalCandidate;
}

#define FORCE_PCM(h) false


template <typename Sample, class H>
void searchIntraCu(H &h, CandidateStash<Sample> *&challenger, CandidateStash<Sample> *&champion)
{
    StateReconstructionCache<Sample> *stateReconstructionCache = h;
    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    coding_quadtree const *cqt = h;
    Candidate<Sample> **activeCandidate = h;
    Candidate<Sample> *originalCandidate = *activeCandidate;

    challenger->~CandidateStash<Sample>();
    new (challenger) CandidateStash<Sample>(*originalCandidate, *cqt, *stateReconstructionCache);

    challenger->resetPieces();
    challenger->copy(*originalCandidate, before, *cqt, true, true);
    challenger->copyIntraReferenceSamplesLuma(*originalCandidate, before, *cqt);
    challenger->copyIntraReferenceSamplesChroma(*originalCandidate, before, *cqt);
    challenger->ContextsAndCost::copy(*originalCandidate);
    challenger->StatePieces<Sample>::copyBefore(*cqt, *originalCandidate);
    challenger->StateCodedData::copyBefore(*cqt, *originalCandidate);

    *activeCandidate = challenger;

    Neighbourhood *neighbourhood = h;
    Snake<BlockData>::Cursor *cursor = h;

    BlockData &blockData = cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1);
    blockData.refIdxPlus1[0] = 0;
    blockData.refIdxPlus1[1] = 0;
    blockData.skip = false;

    {
        StateCodedData *stateCodedData = h;
        stateCodedData->startCu();
        stateCodedData->codedCu.init();
        stateCodedData->codedCu.word0().CtDepth = cqt->cqtDepth;
        stateCodedData->codedCu.word0().CuPredMode = MODE_INTRA;
        stateCodedData->codedCu.word0().intra_chroma_pred_mode = 4;
        stateCodedData->codedCu.word0().part_mode = 0;
        stateCodedData->codedCu.chromaOffset() = 0;
        stateCodedData->codedDataAfter = stateCodedData->codedCu.firstTransformTree().p;
    }

    // reset cached intra prediction source samples - mark as invalid so they're recreated at next CU - review
    for (int cIdx = 0; cIdx < 3; ++cIdx)
        new (&stateEncodeSubstream->unfiltered[cIdx]) IntraReferenceSamples<Sample>();
    new (&stateEncodeSubstream->filtered) IntraReferenceSamples<Sample>();

    RateDistortion::Metric cost2Nx2N;

    if (!static_cast<Speed *>(h)->doIntraSearch())
    {
        // We're running fast and don't have time for an intra search - always use planar mode
        assert(!"currently broken");
    }
    else
    {
        RateDistortion::Metric costNxN = RateDistortion::Metric::max();

        (*activeCandidate)->codedCu.word0().part_mode = 0;

        {
            // Measure rate of coding unit preamble
            auto m = h.template change<EstimateRate<void>>();
            if (m[transquant_bypass_enabled_flag()])
                m(cu_transquant_bypass_flag(), ae(v));
            if (m[slice_type()] != I)
                m(cu_skip_flag(cqt->x0, cqt->y0), ae(v));
            if (m[slice_type()] != I)
                m(pred_mode_flag(), ae(v));
        }

        // review: move inside "if (tryNxN) { }"
        CandidateStash<Sample> candidateNxN(**activeCandidate, *cqt, *stateReconstructionCache);
        candidateNxN.resetPieces();
        candidateNxN.copy(**activeCandidate, before, *cqt, true, true);
        candidateNxN.copyIntraReferenceSamplesLuma(**activeCandidate, before, *cqt);
        candidateNxN.copyIntraReferenceSamplesChroma(**activeCandidate, before, *cqt);
        candidateNxN.ContextsAndCost::copy(**activeCandidate);
        candidateNxN.StatePieces<Sample>::copyBefore(*cqt, **activeCandidate);
        candidateNxN.StateCodedData::copyBefore(*cqt, **activeCandidate);
        candidateNxN.codedCu.word0().part_mode = 1;


        {
            auto m = h.template change<EstimateRate<void>>();

            if (/*m[current(CuPredMode(cu.x0, cu.y0))] != MODE_INTRA || */cqt->log2CbSize == m[MinCbLog2SizeY()])
            {
                m(part_mode(), ae(v));
            }
            if (m[PartMode()] == PART_2Nx2N && m[pcm_enabled_flag()] &&
                cqt->log2CbSize >= m[Log2MinIpcmCbSizeY()] &&
                cqt->log2CbSize <= m[Log2MaxIpcmCbSizeY()])
            {
                m(pcm_flag(cqt->x0, cqt->y0), ae(v));
            }
        }

        static_cast<Candidate<Sample>*>(h)->codedCu.word1().raw = 0;
        static_cast<Candidate<Sample>*>(h)->codedCu.chromaOffset() = 0;
        static_cast<Candidate<Sample>*>(h)->transformTree = static_cast<Candidate<Sample>*>(h)->codedCu.firstTransformTree();

        h[MaxTrafoDepth()] = h[max_transform_hierarchy_depth_intra()];

        if (FORCE_PCM(h))
        {
            // Force this CU to be PCM
            // review: make this a function, place inside Reconstruct.h

            auto &pictureInput = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
            StateReconstructionCache<Sample> *stateReconstructionCache = h;

            // write CodedData for the IPCM CU - this just comprises an empty CodedCu with pcm_flag set - no actual IPCM data is stored in CodedData
            StateCodedData *stateCodedData = h;
            stateCodedData->codedCu.word0().pcm_flag = 1;
            stateCodedData->codedCu.chromaOffset() = 0;
            stateCodedData->codedDataAfter = stateCodedData->codedCu.firstTransformTreeChroma().p;

            // Create a "piece" for each colour component, copy source data direct to it and append it to the candidate's list
            blockData.intra.pcm = 1;
            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                auto const rc = residual_coding(cqt->x0, cqt->y0, cqt->log2CbSize - (cIdx ? 1 : 0), cIdx);
                auto const nCbS = 1 << rc.log2TrafoSize;

                auto const z = static_cast<Candidate<Sample>*>(h)->StatePieces<Sample>::zz[rc.cIdx];
                assert(z == (zPositionOf(rc) >> (rc.cIdx ? 2 : 0)));
                typename ReconstructionCache<Sample>::Piece piece = stateReconstructionCache->components[rc.cIdx].allocateBlock(rc.log2TrafoSize, z);
#ifdef DEBUG_PIECES
                piece.rc = rc;
#endif
                auto recSamples = stateReconstructionCache->components[rc.cIdx].get(piece);
                assert((int)recSamples.stride == 1 << rc.log2TrafoSize);

                auto in = pictureInput(cqt->x0, cqt->y0, cIdx);

                for (int y = 0; y < nCbS; ++y)
                {
                    memcpy(&recSamples(0, y), &in(0, y), nCbS * sizeof(Sample));
                }

                static_cast<Candidate<Sample>*>(h)->appendPiece(rc.cIdx, rc, rc.log2TrafoSize, piece.i);

                static_cast<NeighbourhoodEnc<Sample>*>(h)->snakeIntraReferenceSamples[cIdx].copyFrom2D(*cqt, pictureInput[cIdx], (cIdx ? 1 : 0));
            }

            static_cast<Snake<BlockData>::Cursor*>(h)->relocate(static_cast<Neighbourhood *>(h)->snake, *cqt, neighbourhood->MinCbLog2SizeYMinus1);
            static_cast<Snake<BlockData>::Cursor*>(h)->commit(*cqt, neighbourhood->MinCbLog2SizeYMinus1);
        }
        else
        {
            searchIntraPartition(h, IntraPartition(cqt->x0, cqt->y0, cqt->log2CbSize, 0, 0));
        }

        static_cast<Candidate<Sample>*>(h)->checkPieces(0, after, *cqt);

        const bool tryNxN = cqt->log2CbSize == h[MinCbLog2SizeY()] && !FORCE_PCM(h);
        if (tryNxN)
        {
            *activeCandidate = &candidateNxN;

            {
                auto m = h.template change<EstimateRate<void>>();
                m(part_mode(), ae(v));
            }

            candidateNxN.codedCu.chromaOffset() = 0;

            candidateNxN.codedCu.IntraPredModeY(0) = 0;
            candidateNxN.codedCu.IntraPredModeY(1) = 0;
            candidateNxN.codedCu.IntraPredModeY(2) = 0;
            candidateNxN.codedCu.IntraPredModeY(3) = 0;

            candidateNxN.codedDataAfter = candidateNxN.codedCu.firstTransformTree().p;

            h[MaxTrafoDepth()] = h[max_transform_hierarchy_depth_intra()] + 1;
            searchIntraPartition(h, IntraPartition(cqt->x0, cqt->y0, cqt->log2CbSize, 1, 0));
            searchIntraPartition(h, IntraPartition(cqt->x0, cqt->y0, cqt->log2CbSize, 1, 1));
            searchIntraPartition(h, IntraPartition(cqt->x0, cqt->y0, cqt->log2CbSize, 1, 2));
            searchIntraPartition(h, IntraPartition(cqt->x0, cqt->y0, cqt->log2CbSize, 1, 3));

            static_cast<Candidate<Sample>*>(h)->checkPieces(0, after, *cqt);

            *activeCandidate = challenger;

            if (candidateNxN.cost2() < challenger->cost2())
            {
                challenger->copy(candidateNxN, after, *cqt, false, true);
                challenger->copyIntraReferenceSamplesLuma(candidateNxN, after, *cqt);
                challenger->ContextsAndCost::copy(candidateNxN);
                challenger->StatePieces<Sample>::copyAfter(*cqt, candidateNxN, 0x1 /*= luma only */);
                challenger->StateCodedData::copyAfter(*cqt, candidateNxN);
            }
        }

        Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->searchIntraChroma);

        (*activeCandidate)->checkPieces(0, after, *cqt);

        if (!FORCE_PCM(h))
        {
            searchIntraChroma(h, *cqt);
        }

        (*activeCandidate)->checkPieces(1, after, *cqt);
        (*activeCandidate)->checkPieces(2, after, *cqt);

        *activeCandidate = originalCandidate;
    }
}

// Ignore PCM (for now at least) when 'Search'ing.
template <> struct Search<pcm_sample> : Null<pcm_sample> {};

// When an Element is encountered during 'Search'ing, ignore it...
template <class V> struct Search<Element<V, ae>> : Null<Element<V, ae>> {};
template <class V> struct Search<Element<V, f>> : Null<Element<V, f>> {};

template <class V>
struct SearchEstimateElementRate
{
    template <class H> static void go(const Element<V, ae> &e, H &h)
    {
        auto hEstimateRate = h.template change<EstimateRate<void>>();
        hEstimateRate(e);
    }
};

// ... with a few exceptions: the following are seen, estimate their bitrate.
template <> struct Search<Element<split_cu_flag, ae>> : SearchEstimateElementRate<split_cu_flag> {};


// Deleted CQTs occur at the picture edge when width or height is not a multiple of CTU size.
// During Seach, append null pieces so that when the reconstructed CTU is later assembled, such CQTs
// can be skipped with minimal logic.
template <class Direction>
template <class H>
void Search<Deleted<coding_quadtree, Direction>>::go(const Deleted<coding_quadtree, Direction> &cqt, H &h)
{
    using Sample = typename SampleType<H>::Type;

    Candidate<Sample> *candidate = h;

    // append "null" pieces - these will be skipped during assembly of the reconstructed planes
    candidate->appendPiece(0, residual_coding(cqt.x0, cqt.y0, cqt.log2CbSize, 0), uint16_t(cqt.log2CbSize), -1);
    candidate->appendPiece(1, residual_coding(cqt.x0, cqt.y0, cqt.log2CbSize - 1, 1), uint16_t(cqt.log2CbSize - 1), -1);
    candidate->appendPiece(2, residual_coding(cqt.x0, cqt.y0, cqt.log2CbSize - 1, 2), uint16_t(cqt.log2CbSize - 1), -1);

    candidate->check += widthOf(cqt)* heightOf(cqt);

    NeighbourhoodEnc<Sample> *neighbourhood = h;

    for (int cIdx = 0; cIdx < 3; ++cIdx)
    {
        auto &snake = neighbourhood->snakeIntraReferenceSamples[cIdx];
        auto const log2Resolution = cIdx ? 1 : 0;

        Sample sample;

        if (std::is_same<Direction, Down>::value)
            sample = snake.at(cqt.x0 + (1 << cqt.log2CbSize) - 1, cqt.y0 - 1, log2Resolution);
        else
            sample = snake.at(cqt.x0 - 1, cqt.y0 + (1 << cqt.log2CbSize) - 1, log2Resolution);

        // review: this commits two edges of the rectangle - really only need one	
        snake.commitRectangle(cqt, sample, log2Resolution);
    }
}

template <class H>
void Search<sao>::go(const sao &s, H &h)
{
    const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
    const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];
    auto h2 = h.template change<EstimateRate<void>>();
    if (s.rx > 0)
    {
        const bool leftCtbInSliceSeg = h[CtbAddrInRs()] > h[SliceAddrRs()];
        const bool leftCtbInTile = h[TileId(h[CtbAddrInTs()])] == h[TileId(h[CtbAddrRsToTs(h[CtbAddrInRs()] - 1)])];
        if (leftCtbInSliceSeg && leftCtbInTile)
            h2(sao_merge_left_flag(), ae(v));
    }
    if (s.ry > 0 && !h2[sao_merge_left_flag()])
    {
        const bool upCtbInSliceSeg = (h[CtbAddrInRs()] - h[PicWidthInCtbsY()]) >= h[SliceAddrRs()];
        const bool upCtbInTile = h[TileId(h[CtbAddrInTs()])] == h[TileId(h[CtbAddrRsToTs(h[CtbAddrInRs()] - h[PicWidthInCtbsY()])])];
        if (upCtbInSliceSeg && upCtbInTile)
            h2(sao_merge_up_flag(), ae(v));
    }
    if (!h2[sao_merge_up_flag()] && !h2[sao_merge_left_flag()])
    {
        if (h2[slice_sao_luma_flag()])
        {

            h2(sao_type_idx_luma(), ae(v));
            if (h2[SaoTypeIdx(0, s.rx, s.ry)] != 0)
            {
                for (int i = 0; i < 4; i++)
                    h2(sao_offset_abs(0, s.rx, s.ry, i), ae(v));
                if (h2[SaoTypeIdx(0, s.rx, s.ry)] == 1)
                {
                    for (int i = 0; i < 4; i++)
                        if (h2[sao_offset_abs(0, s.rx, s.ry, i)] != 0)
                            h2(sao_offset_sign(0, s.rx, s.ry, i), ae(v));
                    h2(sao_band_position(0, s.rx, s.ry), ae(v));
                }
                else
                {
                    h2(sao_eo_class_luma(), ae(v));
                }
            }
        }
        if (h2[slice_sao_chroma_flag()])
        {
            h2(sao_type_idx_chroma(), ae(v));
            if (h2[SaoTypeIdx(1, s.rx, s.ry)] != 0)
            {
                for (int i = 0; i < 4; i++)
                    h2(sao_offset_abs(1, s.rx, s.ry, i), ae(v));
                if (h2[SaoTypeIdx(1, s.rx, s.ry)] == 1)
                {
                    for (int i = 0; i < 4; i++)
                        if (h2[sao_offset_abs(1, s.rx, s.ry, i)] != 0)
                            h2(sao_offset_sign(1, s.rx, s.ry, i), ae(v));
                    h2(sao_band_position(1, s.rx, s.ry), ae(v));
                }
                else
                {
                    h2(sao_eo_class_chroma(), ae(v));
                }
            }
        }
    }
}

template <class H>
void Search<coding_quadtree>::go(const coding_quadtree &cqt, H &h)
{
    StateEncodePicture *stateEncodePicture = h;
    StateEncode *stateEncode = h;

    using Sample = typename SampleType<H>::Type;

    Candidate<Sample> **activeCandidate = h;
    Candidate<Sample> *originalCandidate = *activeCandidate;
    Neighbourhood *neighbourhood = h;
    //RCU depth:
    //RCUdepth can happen here
    if (cqt.cqtDepth == 0 && stateEncode->rcudepth)
    {
        originalCandidate->rcudepthstatus = 0;
        if (h[slice_type()] != I)
        {
            int depthsum = 0;
            if (cqt.x0 && cqt.y0)
            {
                int stepx = 32, stepy = 32;
                if (cqt.x0 + (1 << cqt.log2CbSize) > h[pic_width_in_luma_samples()])
                    stepx = 16;
                if (cqt.y0 + (1 << cqt.log2CbSize) > h[pic_height_in_luma_samples()])
                    stepy = 16;
                depthsum += neighbourhood->snakeMerge.get<Up>(cqt.x0, cqt.y0 - 1, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;
                depthsum += neighbourhood->snakeMerge.get<Up>(cqt.x0 + stepx, cqt.y0 - 1, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;

                depthsum += neighbourhood->snakeMerge.get<Left>(cqt.x0 - 1, cqt.y0, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;
                depthsum += neighbourhood->snakeMerge.get<Left>(cqt.x0 - 1, cqt.y0 + stepy, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;

                depthsum += neighbourhood->snakeMerge.get<Corner>(cqt.x0 - 1, cqt.y0 - 1, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;

                if (depthsum < 6)
                {
                    originalCandidate->rcudepthstatus = 1;
                }
                else if (depthsum < 14)
                {
                    originalCandidate->rcudepthstatus = 2;
                }
                else
                {
                    originalCandidate->rcudepthstatus = 3;
                }
            }
            else if (cqt.x0)
            {
                int stepx = 32;
                if (cqt.x0 + (1 << cqt.log2CbSize) > h[pic_width_in_luma_samples()])
                    stepx = 16;
                depthsum += neighbourhood->snakeMerge.get<Up>(cqt.x0, cqt.y0 - 1, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;
                depthsum += neighbourhood->snakeMerge.get<Up>(cqt.x0 + stepx, cqt.y0 - 1, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;
                if (depthsum < 4)
                {
                    originalCandidate->rcudepthstatus = 1;
                }
                else
                {
                    originalCandidate->rcudepthstatus = 2;
                }
            }
            else if (cqt.y0)
            {
                int stepy = 32;
                if (cqt.y0 + (1 << cqt.log2CbSize) > h[pic_height_in_luma_samples()])
                    stepy = 16;
                depthsum += neighbourhood->snakeMerge.get<Left>(cqt.x0 - 1, cqt.y0, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;
                depthsum += neighbourhood->snakeMerge.get<Left>(cqt.x0 - 1, cqt.y0 + stepy, neighbourhood->MinCbLog2SizeYMinus1).CtDepth;
                if (depthsum < 4)
                {
                    originalCandidate->rcudepthstatus = 1;
                }
                else
                {
                    originalCandidate->rcudepthstatus = 2;
                }
            }
        }
    }
    else if (cqt.cqtDepth == 0)
    {
        originalCandidate->rcudepthstatus = 0;
    }

    bool testFull = true;
    bool testSplit = true;


    if (stateEncode->rcudepth)
    {
        if (cqt.cqtDepth == 2 && originalCandidate->rcudepthstatus == 1)
            testSplit = false;
        if (cqt.cqtDepth == 0 && (originalCandidate->rcudepthstatus == 2 || originalCandidate->rcudepthstatus == 3))
            testFull = false;
        if (cqt.cqtDepth == 1 && originalCandidate->rcudepthstatus == 3)
            testFull = false;
    }

    originalCandidate->StateCodedData::startCu();

    bool mustSplit =
        cqt.x0 + (1 << cqt.log2CbSize) > h[pic_width_in_luma_samples()] ||
        cqt.y0 + (1 << cqt.log2CbSize) > h[pic_height_in_luma_samples()];

    if (FORCE_PCM(h) && cqt.log2CbSize > h[Log2MaxIpcmCbSizeY()])
    {
        mustSplit = true;
    }

    if (mustSplit)
    {
        h[split_cu_flag()] = 1;

        Syntax<coding_quadtree>::go(cqt, h);
    }
    else
    {
        // Copy initial state of contexts etc. before trying split=0
        StateReconstructionCache<Sample> *stateReconstructionCache = h;
        CandidateStash<Sample> splitCandidate(*originalCandidate, cqt, *stateReconstructionCache);
        splitCandidate.resetPieces();
        splitCandidate.copy(*originalCandidate, before, cqt, true, true);
        splitCandidate.copyIntraReferenceSamplesLuma(*originalCandidate, before, cqt);
        splitCandidate.copyIntraReferenceSamplesChroma(*originalCandidate, before, cqt);
        splitCandidate.ContextsAndCost::copy(*originalCandidate);
        splitCandidate.StatePieces<Sample>::copyBefore(cqt, *originalCandidate);
        splitCandidate.StateCodedData::copyBefore(cqt, *originalCandidate);
        splitCandidate.rcudepthstatus = originalCandidate->rcudepthstatus;

        if (testFull)
        {
            // single CU / split = 0
            h[split_cu_flag()] = 0;

            Syntax<coding_quadtree>::go(cqt, h);

            static_cast<Candidate<Sample> *>(h)->checkPieces(0, after, cqt);
            static_cast<Candidate<Sample> *>(h)->checkPieces(1, after, cqt);
            static_cast<Candidate<Sample> *>(h)->checkPieces(2, after, cqt);

            // resultant state of single CU / split=0 test is now in *originalCandidate
        }
        const bool trySplit = !testFull || (testSplit && (
            cqt.log2CbSize > h[MinCbLog2SizeY()] &&
            (!stateEncode->ecu || h[current(CuPredMode(cqt.x0, cqt.y0))] != MODE_SKIP) &&
            static_cast<Speed *>(h)->trySplit(cqt)));

        if (trySplit)
        {
            Candidate<Sample> &singleCandidate = *originalCandidate;

            // restore original states
            *activeCandidate = &splitCandidate;

            // split
            h[split_cu_flag()] = 1;

            Syntax<coding_quadtree>::go(cqt, h);

            static_cast<Candidate<Sample> *>(h)->checkPieces(0, after, cqt);
            static_cast<Candidate<Sample> *>(h)->checkPieces(1, after, cqt);
            static_cast<Candidate<Sample> *>(h)->checkPieces(2, after, cqt);

            *activeCandidate = originalCandidate;

            if (!testFull || splitCandidate.cost2() < singleCandidate.cost2())
            {
                // Copy
                originalCandidate->copy(splitCandidate, after, cqt, true, true);
                originalCandidate->copyIntraReferenceSamplesLuma(splitCandidate, after, cqt);
                originalCandidate->copyIntraReferenceSamplesChroma(splitCandidate, after, cqt);
                originalCandidate->ContextsAndCost::copy(splitCandidate);
                originalCandidate->StatePieces<Sample>::copyAfter(cqt, splitCandidate);
                originalCandidate->StateCodedData::copyAfter(cqt, splitCandidate);
            }
        }
    }
}


// returns true if a search of the specified partMode would be legal and is enabled
template <class H>
bool doSearchInterCuPartMode(H &h, coding_quadtree const &cqt, int partMode)
{
    if (partMode == PART_NxN)
    {
        if (cqt.log2CbSize > h[MinCbLog2SizeY()]) return false;
        if (cqt.log2CbSize == 3) return false;
        return false;
        // NxN is only used when  MinCbLog2SizeY > 3 and this is not a common use case.
    }
    else if (partMode)
    {
        // rectangular
        StateEncode *stateEncode = h;
        Speed *speed = h;
        if (partMode >= PART_2NxN && partMode <= PART_Nx2N)
        {
            if (stateEncode->smp) return true;
            if (stateEncode->nosmp) return false;
            if (!speed->useSmp(cqt.log2CbSize)) return false;
            return true;
        }
        if (partMode >= PART_2NxnU)
        {
            if (cqt.log2CbSize == h[MinCbLog2SizeY()]) return false;
            if (!h[amp_enabled_flag()]) return false;
        }
    }

    return true;
}


// Searches the current coding unit using a merged 2Nx2N partition
template <class H>
void searchMerge2Nx2N(H &h)
{
    coding_quadtree const *cqt = h;
    coding_unit const cu{ cqt->x0, cqt->y0, cqt->log2CbSize };

    auto h2 = h.template change<SearchMerge2Nx2N<typename H::Mode>>();
    Syntax<coding_unit>::go(cu, h2);
}


// Search the specified coding unit as MODE_INTER with the specified partMode.
template <typename Sample, class H>
void searchInterCuPartMode(H &h, coding_quadtree const &cqt, int partMode, CandidateStash<Sample> *&contender, CandidateStash<Sample> *&champion)
{
    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    StateEncode *stateEncode = h;
    StateReconstructionCache<Sample> *stateReconstructionCache = h;
    Neighbourhood *neighbourhood = h;
    Snake<BlockData>::Cursor *cursor = h;
    Candidate<Sample> **activeCandidate = h;

    // Review: set current blockData this only after searching a PU and then insert to snake
    BlockData &blockData = cursor->current(cqt.x0, cqt.y0, neighbourhood->MinCbLog2SizeYMinus1);
    blockData.refIdxPlus1[0] = 1;
    blockData.refIdxPlus1[1] = 1;
    blockData.skip = false;
    blockData.setup(&cqt, MODE_INTER);

    contender->~CandidateStash<Sample>();
    new (contender) CandidateStash<Sample>(**activeCandidate, cqt, *stateReconstructionCache);

    contender->resetPieces();

    contender->copy(**activeCandidate, before, cqt, true, true);
    contender->ContextsAndCost::copy(**activeCandidate);

    contender->StatePieces<Sample>::copyBefore(cqt, **activeCandidate);

    contender->StateCodedData::copyBefore(cqt, **activeCandidate);
    contender->StateCodedData::codedCu.init();
    contender->StateCodedData::codedCu.word0().CtDepth = cqt.cqtDepth;
    contender->StateCodedData::codedCu.word0().CuPredMode = MODE_INTER;
    contender->StateCodedData::codedCu.word0().part_mode = partMode < 0 ? 0 : partMode;
    contender->StateCodedData::firstPu();

    stateEncodeSubstream->partIdx = 0;

    stateEncodeSubstream->originalCandidate = *activeCandidate;
    *activeCandidate = contender;

    auto const zY = zPositionOf(cqt);
    auto const zC = zY >> 2;
    stateEncodeSubstream->interPieces[0][0] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[0].allocateBlock(cqt.log2CbSize, zY);
    stateEncodeSubstream->interPieces[1][0] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[1].allocateBlock(cqt.log2CbSize - 1, zC);
    stateEncodeSubstream->interPieces[2][0] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[2].allocateBlock(cqt.log2CbSize - 1, zC);
    stateEncodeSubstream->interPieces[0][1] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[0].allocateBlock(cqt.log2CbSize, zY);
    stateEncodeSubstream->interPieces[1][1] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[1].allocateBlock(cqt.log2CbSize - 1, zC);
    stateEncodeSubstream->interPieces[2][1] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[2].allocateBlock(cqt.log2CbSize - 1, zC);
    stateEncodeSubstream->interPieces[0][2] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[0].allocateBlock(cqt.log2CbSize, zY);
    stateEncodeSubstream->interPieces[1][2] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[1].allocateBlock(cqt.log2CbSize - 1, zC);
    stateEncodeSubstream->interPieces[2][2] = stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[2].allocateBlock(cqt.log2CbSize - 1, zC);
    if (partMode < 0)
    {
        contender->StateCodedData::codedCu.word1().merge[stateEncodeSubstream->partIdx] = partMode + h[MaxNumMergeCand()] + 1;
        if (stateEncode->fdm && champion->ContextsAndCost::rate != std::numeric_limits<Cost>::max() && champion->StateCodedData::codedCu.word0().CuPredMode == MODE_SKIP)
        {
            contender->noresidual = 1;
        }
        else
        {
            contender->noresidual = 0;
        }
        searchMerge2Nx2N(h);
    }
    else
    {
        // Descend call tree to perform PU search (ME, etc.) and TT search
        // CU is reconstructed during this call and distortion measured
        // Quantised coefficients and CBFs are decided.
        // Rate of PU is measured during PU decisions
        // Rate of TT is measured during second pass after quantization (and/or during RDOQ
        contender->noresidual = 0;
        if (stateEncode->fdam && !champion->StateCodedData::codedCu.word0().cbfWord)
        {
            contender->noresidual = 1;
        }
        else
        {
            contender->noresidual = 0;
        }
        Syntax<coding_unit>::go(coding_unit(cqt.x0, cqt.y0, cqt.log2CbSize), h);
    }

    *activeCandidate = stateEncodeSubstream->originalCandidate;

    if (contender->cost2() < champion->cost2())
    {
        std::swap(contender, champion);
    }
}


// Finds the best inter-coded configuration for the current CU.
template <typename Sample, class H>
bool searchInterCu(H &h, coding_quadtree const &cqt, CandidateStash<Sample> *&contender, CandidateStash<Sample> *&champion)
{
    StateReconstructionCache<Sample> *stateReconstructionCache = h;
    StateEncode *stateEncode = h;

    auto contextsBefore = *static_cast<Contexts *>(h);
    bool do2NxN = true, doNx2N = true;
    Aps<Sample> apsEngine;
    apsEngine.initApsModule(4);
    int maxNumMergeCandidates = h[MaxNumMergeCand()];
    contender->noresidual = 0;
    // Loop over all partition types
    for (int partMode = PART_2Nx2N - maxNumMergeCandidates; partMode <= PART_nRx2N; ++partMode)
    {
        if (doSearchInterCuPartMode(h, cqt, partMode))
        {
            if (stateEncode->aps && partMode == PART_2NxN && !do2NxN)
                continue;
            if (stateEncode->aps && partMode == PART_Nx2N && !doNx2N)
                continue;
            searchInterCuPartMode(h, cqt, partMode, contender, champion);

            if (stateEncode->aps && partMode == PART_2Nx2N)
                apsEngine.analyseResidueEnergy(champion, cqt, do2NxN, doNx2N);

            bool const cuHasResidual = !!champion->StateCodedData::codedCu.word0().cbfWord;
            if (partMode >= 0 && !cuHasResidual)
            {
                if (stateEncode->esd && (champion->StateCodedData::codedCu.word0().CuPredMode == MODE_SKIP)) break;
                if (stateEncode->cfm && (partMode > PART_2Nx2N) && (partMode == champion->StateCodedData::codedCu.word0().part_mode)) break;
            }
        }
    }

    if (champion->StateCodedData::codedCu.word0().part_mode == PART_2Nx2N)
    {
        // review: could do this later after inter/intra decision is made
        Neighbourhood *neighbourhood = champion;
        Snake<BlockData>::Cursor *cursor = champion;
        cursor->relocate(neighbourhood->snake, cqt, neighbourhood->MinCbLog2SizeYMinus1);
        auto const nCbS = 1 << cqt.log2CbSize;
        prediction_unit pu{ cqt.x0, cqt.y0, nCbS, nCbS };
        cursor->commit(pu, h[MinCbLog2SizeY()] - 1);
        neighbourhood->recordMerge(h, pu);
    }

    if (!champion->codedCu.word0().cbfWord)
    {
        // Cbfs are 0 - there is no residual (review, perhaps not best place to set this)
        champion->codedDataAfter = champion->codedPu.p;
    }

    // Return true if any component's CBF is set
    return !!champion->codedCu.word0().cbfWord;
}


template <>
struct Search<coding_unit>
{
    template <class H> static void go(const coding_unit &cu, H &h)
    {
        using Sample = typename SampleType<H>::Type;

        StateReconstructionCache<Sample> *stateReconstructionCache = h;
        StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
        coding_quadtree const *cqt = h;
        Speed *speed = h;

        Snake<BlockData>::Cursor *cursor = h;
        Neighbourhood *neighbourhood = h;
        cursor->relocate(neighbourhood->snake, cu, neighbourhood->MinCbLog2SizeYMinus1);
        cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1).intra.pcm = 0;

        Candidate<Sample> **activeCandidate = h;
        Candidate<Sample> *originalCandidate = *activeCandidate;

        // for checking bitrate
        CandidateStash<Sample> candidateCheckRate(*originalCandidate, *cqt, *stateReconstructionCache);
        candidateCheckRate.copy(*originalCandidate, before, *cqt, true, true);
        candidateCheckRate.ContextsAndCost::copy(*originalCandidate);
        candidateCheckRate.resetPieces();

        CandidateStash<Sample> candidateStash[2];
        CandidateStash<Sample> *challenger = &candidateStash[1];
        CandidateStash<Sample> *champion = &candidateStash[0];

        champion->ContextsAndCost::setMax();

        assert(h[SubWidthC()] == 2);
        assert(h[SubHeightC()] == 2);

        StateReconstructedPicture<Sample> *stateReconstructedPicture = h;
        Picture<Sample> &currPic = *stateReconstructedPicture->picture;

        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            stateEncodeSubstream->reconstructed0[cIdx] = currPic(0, 0, cIdx);
        }

        QpState *qpState = h;
        qpState->preCu(cu, h);

        if (h[cu_qp_delta_enabled_flag()])
        {
            const int qpOffset[4] = { 0, 3, -3, 1 };
            int row = (cu.y0 & (h[CtbSizeY()] - 1)) >> 3;
            int col = (cu.x0 & (h[CtbSizeY()] - 1)) >> 3;
            int qpValue = qpState->getQpInternal(row, col);
            // Add qp offset if adaptive qp is enabled
            StateEncode *stateEncode = h;
            if (stateEncode->useAq)
            {
                AdaptiveQuantisation &aqInfo = *static_cast<StateEncodePicture *>(h)->docket->aqInfo;
                int poc = static_cast<StateEncodePicture *>(h)->docket->poc;
                int qpOffset = aqInfo.getAqOffset(cu.y0, cu.x0, cqt->cqtDepth);
                qpValue += qpOffset;
                qpValue = Clip3(0, 51, qpValue);
            }

            qpState->setQpValue(qpValue);
        }

        bool doTryIntra = true;

        originalCandidate->startCu();
        originalCandidate->codedCu.init();
        originalCandidate->codedCu.word0().CtDepth = cqt->cqtDepth;

        if (h[slice_type()] != I && !FORCE_PCM(h))
        {
            static_cast<Candidate<Sample> *>(h)->checkPosition(before, cu);

            Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->searchInter);

            doTryIntra = searchInterCu(h, *cqt, challenger, champion);
            if (!(speed->doIntraInInter())) 
                doTryIntra = false;
        }

        static_cast<Candidate<Sample> *>(h)->checkPosition(before, cu);

        bool useIntra = false;

        if (doTryIntra)
        {
            static_cast<Candidate<Sample> *>(h)->checkPosition(before, cu);

            Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->searchIntra);

            searchIntraCu(h, challenger, champion);

            *activeCandidate = originalCandidate;

            if (challenger->cost2() < champion->cost2())
            {
                useIntra = true;

                originalCandidate->copy(*challenger, after, cu, false, true);
                originalCandidate->copyIntraReferenceSamplesLuma(*challenger, after, cu);
                originalCandidate->copyIntraReferenceSamplesChroma(*challenger, after, cu);
                originalCandidate->ContextsAndCost::copy(*challenger);
                originalCandidate->StatePieces<Sample>::copyAfter(cu, *challenger);
                originalCandidate->StateCodedData::copyAfter(cu, *challenger);
            }
        }

        if (!useIntra)
        {
            // Copy various information from champion Candidate<Sample> to h's Candidate...
            static_cast<StatePieces<Sample> *>(h)->copyAfter(cu, *champion);
            static_cast<StateCodedData *>(h)->copyAfter(*cqt, *champion);
            static_cast<Neighbourhood *>(h)->snake.copyBlockFrom(champion->snake, after, *cqt, static_cast<Neighbourhood *>(h)->MinCbLog2SizeYMinus1, 0, 0);
            static_cast<Neighbourhood *>(h)->snakeMerge.copyBlockFrom(champion->snakeMerge, after, *cqt, static_cast<Neighbourhood *>(h)->MinCbLog2SizeYMinus1, 0, 0);
            *static_cast<ContextsAndCost *>(h) = static_cast<ContextsAndCost &>(*champion);

            // update intra reference samples in snake buffers
            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                typename ReconstructionCache<Sample>::Piece piece = *(originalCandidate->pieces[cIdx][after] - 1);
                assert(piece.log2Size == cu.log2CbSize - (cIdx ? 1 : 0));
                auto recSamples = stateReconstructionCache->components[cIdx].get(piece).offset(-(cqt->x0 >> (cIdx ? 1 : 0)), -(cqt->y0 >> (cIdx ? 1 : 0)));
                originalCandidate->snakeIntraReferenceSamples[cIdx].copyFrom2D(*cqt, recSamples, cIdx ? 1 : 0);
            }
        }

        BlockData &blockData = cursor->current(cu.x0, cu.y0, neighbourhood->MinCbLog2SizeYMinus1);
        blockData.refIdxPlus1[0] = originalCandidate->codedCu.word0().CuPredMode == MODE_INTRA ? 0 : 88;
        blockData.refIdxPlus1[1] = originalCandidate->codedCu.word0().CuPredMode == MODE_INTRA ? 0 : 88;
        blockData.skip = originalCandidate->codedCu.word0().CuPredMode == MODE_SKIP;

        {
            coding_quadtree *cqt = h;
            Neighbourhood *neighbourhood = h;
            Snake<BlockData>::Cursor *cursor = h;
            BlockData &blockData = cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1);
            blockData.setup(cqt, originalCandidate->codedCu.word0().CuPredMode);
        }

        static_cast<QpState *>(h)->postCu(cu, h);

        neighbourhood->recordMerge(h, *cqt, useIntra);

        // estimate whole-CU bitrate again to make sure that computed during searches is correct
        if (!FORCE_PCM(h))
        {
            // copy coded data pointers
            static_cast<StateCodedData &>(candidateCheckRate) = static_cast<StateCodedData &>(*originalCandidate);

            *activeCandidate = &candidateCheckRate;

            // walk tree and estimate the rate of resultant bitstream
            auto w = h.template change<EstimateRate<void>>();
            w(cu);

            *activeCandidate = originalCandidate;

            candidateCheckRate.checkSameAs(*originalCandidate);
            ASSERT(candidateCheckRate.rate == originalCandidate->rate);
        }
    }
};


struct MvCandidate
{
    MvCandidate()
        :
        cost(std::numeric_limits<Cost>::max())
    {
    }

    template <class H>
    MvCandidate(H &h, int refList, MotionVector mv, MotionVector predictors[2])
    {
        auto const refIdx = 0;

        // Find the best predictor (0 or 1) for the given motion vector
        {
            MvCandidate &temp = *this;
            temp.mvpFlag = 0;

            // overloaded MotionVector operator-() does not optimise well in MSVC
            temp.mvd[0] = mv[0] - predictors[0][0];
            temp.mvd[1] = mv[1] - predictors[0][1];

            temp.cost = rateOf<H::Mode::value>(temp.mvd);

            EstimateRateBin<mvp_lX_flag> estimateRateMvp(h, 0);
            temp.cost += estimateRateMvp.rate(0);
        }

        {
            MvCandidate temp;
            temp.mvpFlag = 1;

            // overloaded MotionVector operator-() does not optimise well in MSVC
            temp.mvd[0] = mv[0] - predictors[1][0];
            temp.mvd[1] = mv[1] - predictors[1][1];

            temp.cost = rateOf<H::Mode::value>(temp.mvd);

            EstimateRateBin<mvp_lX_flag> estimateRateMvp(h, 0);
            temp.cost += estimateRateMvp.rate(1);

            this->consider(temp);
        }

        this->mv = mv;
    }

    bool consider(MvCandidate const &other)
    {
        bool const better = other.cost < this->cost;
        if (better) *this = other;
        return better;
    }

    MotionVector mv;
    MotionVector mvd;
    Cost cost;
    int mvpFlag;
};


template <class H> static void searchMotionUni(H &h, int refList)
{
    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    StateEncodePicture *stateEncodePicture = h;
    StateCodedData *stateCodedData = h;
    Speed *speed = h;
    coding_quadtree const *cqt = h;
    prediction_unit const *pu = h;
    mvd_coding const mvdc{ pu->x0, pu->y0, refList };

    //if (match(mvdc, h)) logEncode() << mvdc << " uni search (L" << mvdc.refList << ")\n";

    PuData puData;
    setPuDataMvpPredFlags(puData, h, !mvdc.refList, !!mvdc.refList);

    MvCandidate best;

    fullPelMotionEstimation(mvdc, h, best);

    MotionVector mvd = best.mvd;

    if (speed->doHalfPelRefinement())
    {
        using Sample = typename SampleType<H>::Type;

        auto &pictureInput = static_cast<PictureWrap<Sample>&>(*static_cast<StateEncodePicture2<Sample> *>(h)->docket->picture);

        auto const input = pictureInput(pu->x0, pu->y0, 0);
        auto &reconstructedPicture = static_cast<StateReconstructedPicture<Sample> &>(*h[RefPicList(mvdc.refList)][puData.refIdx(mvdc.refList)].dp->reconstructedPicture);
        Picture<Sample> &pictureReference = *reconstructedPicture.picture;
        auto const reference = pictureReference(pu->x0, pu->y0, 0);
        auto mv = best.mv;
        subPelRefinement(mv, mvd, mvdc, h, input, reference);
    }

    stateCodedData->codedPu.mvd(mvdc.refList) = mvd;
    stateCodedData->codedPu.word0().metadata[mvdc.refList].mvp_lX_flag = best.mvpFlag;
}


template <typename T, typename U>
T saturatedCast(U u)
{
    if (u > std::numeric_limits<T>::max()) return std::numeric_limits<T>::max();
    if (u < std::numeric_limits<T>::min()) return std::numeric_limits<T>::min();
    return static_cast<T>(u);
}


struct LimitFullPelMv
{
    template <class H>
    LimitFullPelMv(const prediction_unit &pu, H &h)
    {
        this->min[0] = -h[CtbSizeY()] - pu.x0;
        this->min[1] = -h[CtbSizeY()] - pu.y0;
        this->max[0] = h[pic_width_in_luma_samples()] + h[CtbSizeY()] - pu.x0 - pu.nPbW;
        this->max[1] = h[pic_height_in_luma_samples()] + h[CtbSizeY()] - pu.y0 - pu.nPbH;

        StateEncode *stateEncode = h;

        if (stateEncode->concurrentFrames > 1)
        {
            // Review: this constant limits how far to the right and downwards motion vectors can encroach on reference
            // picture data that may still be unavailable (reference picture encoding still in progress in another thread.
            // Needs review understanding effect of luma and chroma interpolation filters.  Perhaps it can be as low as 3.
            // Also, this limit needs to be applied to refinement (sub pixel and bi).
            const int howCloseDoYouDare = 15;

            int maxWavefront[2] =
            {
                h[xCtb()] + 3 * h[CtbSizeY()] - pu.x0 - pu.nPbW - howCloseDoYouDare,
                h[yCtb()] + 2 * h[CtbSizeY()] - pu.y0 - pu.nPbH - howCloseDoYouDare,
            };

            this->max[0] = std::min(this->max[0], MotionVector::ComponentType(maxWavefront[0]));
            this->max[1] = std::min(this->max[1], MotionVector::ComponentType(maxWavefront[1]));
        }
    }

    void operator()(MotionVector &mv) const
    {
        if (mv[0] < this->min[0]) mv[0] = this->min[0];
        if (mv[1] < this->min[1]) mv[1] = this->min[1];
        if (mv[0] > this->max[0]) mv[0] = this->max[0];
        if (mv[1] > this->max[1]) mv[1] = this->max[1];
    }
private:
    MotionVector min;
    MotionVector max;
};


template <typename Sample>
struct StateMeFullPel
{

    template <class H> StateMeFullPel(H &h, int refList, prediction_unit const &pu, MvCandidate &best) :
        best(best),
        tableSad(h),
        tableSadMultiref(h),
        limit(pu, h),
        refList(refList),
        rect(HAVOC_RECT(pu.nPbW, pu.nPbH))
    {
        this->functionSad = *havoc_get_sad(tableSad, pu.nPbW, pu.nPbH);
        this->functionSad4 = *havoc_get_sad_multiref(tableSadMultiref, 4, pu.nPbW, pu.nPbH);

        auto &picture = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
        this->src = Raster<Sample const>(picture[0], pu.x0, pu.y0);

        auto &referencePicture = static_cast<StateReconstructedPicture<Sample> &>(*h[RefPicList(refList)][0].dp->reconstructedPicture);
        this->ref = Raster<Sample const>((*referencePicture.picture)[0], pu.x0, pu.y0);

        StateEncode *stateEncode = h;
        if(stateEncode->useRateControl)
        {
            StateEncodePicture *stateEncodePicture = h;
            int segmentPoc = stateEncodePicture->docket->segmentPoc;
            stateEncode->rateControlParams->takeTokenOnRCLevel();
            double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalSqrtLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
            stateEncode->rateControlParams->releaseTokenOnRCLevel();
            this->lambda.set(value);
        }
        else
        {
            this->lambda.set(getReciprocalSqrtLambda(h));
        }
    }

    template <class H> bool considerPattern(H &h, MotionVector origin, MotionVector const *pattern, int n, int step, int dist, mvd_coding mvdc)
    {
        Mvp::Predictors *predictors = h;
        int const refIdx = 0;

        bool improved = false;
        for (int j = 0; j < n; j += 4 * step)
        {
            MotionVector mv[4];
            Sample const *refs[4];

            for (int i = 0; i < 4; ++i, pattern += step)
            {
                mv[i][0] = (origin[0] + dist * (*pattern)[0]) / 4;
                mv[i][1] = (origin[1] + dist * (*pattern)[1]) / 4;

                this->limit(mv[i]);


                refs[i] = &this->ref(mv[i][0], mv[i][1]);
            }

            int32_t sads[4];
            this->functionSad4(this->src.p, this->src.stride, refs, this->ref.stride, sads, this->rect);

            for (int i = 0; i < 4; ++i)
            {
                mv[i] *= 4;
                MvCandidate candidate(h, this->refList, mv[i], predictors->mvp[refIdx][mvdc.refList]);
                candidate.cost += lambda * sads[i];
                //if (match(mvdc, h)) logEncode() << mv[i] << " ";
                improved |= this->best.consider(candidate);
            }
        }
        return improved;
    }

    uint32_t rect;
    Lambda lambda;
    int refList;
    havoc_table_sad<Sample> *tableSad;
    havoc_table_sad_multiref<Sample> *tableSadMultiref;
    havoc_sad<Sample> *functionSad;
    havoc_sad_multiref<Sample> * functionSad4;
    MvCandidate &best;
    Raster<const Sample> src;
    Raster<const Sample> ref;
    LimitFullPelMv limit;
};


template <class H> static void searchMotionBi(H &h, int refList)
{
    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    StateCodedData *stateCodedData = h;
    coding_quadtree const *cqt = h;
    prediction_unit const *pu = h;
    mvd_coding const mvdc{ pu->x0, pu->y0, refList };
    Speed *speed = h;
    //if (match(mvdc, h)) logEncode() << mvdc << " BI search (L" << mvdc.refList << ")\n";

    PuData puData;
    setPuDataMvpPredFlags(puData, h, true, true);

    HAVOC_ALIGN(32, Sample, buffer[2][64 * 64]);
    Raster<Sample> predictionOther{ buffer[0], 64 };
    Raster<Sample> predictionIdeal{ buffer[1], 64 };
    MvCandidate best;
    StateMeFullPel<Sample> stateMeFullPel(h, mvdc.refList, *pu, best);
    {
        // predict from other reference

        auto refListOther = 1 - mvdc.refList;

        auto &referencePictureOther = static_cast<StateReconstructedPicture<Sample> &>(*h[RefPicList(1 - mvdc.refList)][puData.refIdx(refListOther)].dp->reconstructedPicture);
        auto reference = (*referencePictureOther.picture)(pu->x0, pu->y0, 0);

        auto mv = puData.mv(refListOther);
        auto const mvFrac = mv & 0x3;
        mv >>= 2;
        stateMeFullPel.limit(mv);
        reference += mv;
        HavocTablePredUni<Sample>  *tablePred = h;
        auto *f = *havocGetPredUni(tablePred, 8, pu->nPbW, pu->nPbH, mvFrac[0], mvFrac[1], h[BitDepthY()]);
        f(predictionOther.p, predictionOther.stride, reference.p, reference.stride, pu->nPbW, pu->nPbH, mvFrac[0], mvFrac[1], h[BitDepthY()]);
    }

    using Sample = typename SampleType<H>::Type;

    auto &pictureInput = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
    Raster<Sample const> input(pictureInput[0], pu->x0, pu->y0);

    // compute residual after considering other reference
    havoc::TableSubtractBi<Sample> *table = h;
    auto f = table->get();
    int constexpr bitDepth = 6 + 2 * sizeof(Sample);

    f(predictionIdeal.p, predictionIdeal.stride, predictionOther.p, predictionOther.stride, input.p, input.stride, pu->nPbW, pu->nPbH, bitDepth);

    stateMeFullPel.src = predictionIdeal;

    // start at best vector from unidirectional search
    auto startingMv = puData.mv(mvdc.refList);

    best.mv = (startingMv + 1) >> 2;
    stateMeFullPel.limit(best.mv);
    auto mv = best.mv;
    auto mv1 = best.mv;
    auto mv2 = best.mv;
    auto mv3 = best.mv;

    best.mv <<= 2;
    best.cost = std::numeric_limits<Cost>::max();

    MotionVector origin = best.mv;

    // exhaustive integer search over small area
    // review: check limits
    Lambda lambda;
    lambda.set(getReciprocalSqrtLambda(h) * 0.5);
    StateEncode *stateEncode = h;
    if(stateEncode->useRateControl)
    {
        StateEncodePicture *stateEncodePicture = h;
        int segmentPoc = stateEncodePicture->docket->segmentPoc;
        stateEncode->rateControlParams->takeTokenOnRCLevel();
        double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalSqrtLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
        stateEncode->rateControlParams->releaseTokenOnRCLevel();
        lambda.set(value * 0.5);
    }

    Mvp::Predictors *predictors = h;
    int const refIdx = 0;

    int const range = speed->useBiSmallSearchWindow() ? 1 : 5;

    for (int y = -range; y <= range; ++y)
    {
        Sample const *refs[4];
        int32_t sads[4];
        for (int x = -range; x <= range; ++x)
        {
            MotionVector mv = origin;
            mv[0] += 4 * x;
            mv[1] += 4 * y;
            mv >>= 2;
            stateMeFullPel.limit(mv);
            auto const i = (x + range) % 4;
            if (i == 0)
            {

                mv1 = mv;
                mv1[0] += 1;
                stateMeFullPel.limit(mv1);
                mv2 = mv;
                mv2[0] += 2;
                stateMeFullPel.limit(mv2);
                mv3 = mv;
                mv3[0] += 3;
                stateMeFullPel.limit(mv3);

                refs[0] = &stateMeFullPel.ref(mv[0], mv[1]);
                refs[1] = &stateMeFullPel.ref(mv1[0], mv1[1]);
                refs[2] = &stateMeFullPel.ref(mv2[0], mv2[1]);
                refs[3] = &stateMeFullPel.ref(mv3[0], mv3[1]);
                stateMeFullPel.functionSad4(stateMeFullPel.src.p, stateMeFullPel.src.stride, refs, stateMeFullPel.ref.stride, sads, stateMeFullPel.rect);

            }
            mv <<= 2;
            MotionVector mvNew = mv;
            MvCandidate candidate(h, mvdc.refList, mvNew, predictors->mvp[refIdx][refList]);
            candidate.cost += lambda * sads[i];
            best.consider(candidate);
        }
    }

    // fractional vector refinement
    auto &pictureReference = static_cast<StateReconstructedPicture<Sample> &>(*h[RefPicList(mvdc.refList)][puData.refIdx(int(mvdc.refList))].dp->reconstructedPicture);
    auto const reference = (*pictureReference.picture)(pu->x0, pu->y0, 0);

    if (speed->doHalfPelRefinement())
    {
        int refinement = (speed->doQuarterPelRefinement()) ? 1 : 2;
        for (int step = 2; step; step -= refinement)
        {
            MotionVector origin = best.mv;
            best.cost = std::numeric_limits<Cost>::max();
            for (int y = -step; y <= step; y += step)
            {
                for (int x = -step; x <= step; x += step)
                {
                    MotionVector mv = origin;
                    mv[0] += x;
                    mv[1] += y;

                    MvCandidate candidate(h, mvdc.refList, mv, predictors->mvp[refIdx][refList]);
                    candidate.cost += costDistortionMv(h, mv, predictionIdeal, reference, 0.5);
                    best.consider(candidate);
                }
            };
        }
    }
    stateCodedData->codedPu.mvd(mvdc.refList) = best.mvd;
    stateCodedData->codedPu.word0().metadata[mvdc.refList].mvp_lX_flag = best.mvpFlag;
}


template <class H>
Cost measurePuCost(H &h)
{
    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
    StatePicture *statePicture = h;
    Snake<BlockData>::Cursor *cursor = h;
    prediction_unit const pu = *static_cast<prediction_unit *>(h);

    BlockData &puData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);
    setPuData(h, puData);

    int32_t satd[3];
    {
        StateReconstructedPicture<Sample> *stateReconstructedPicture = h;

        // Reconstruct
        // review: duplication--we have already predicted luma (and measured its SATD) during sub-pixel refinement
        predictInter(*stateReconstructedPicture->picture, pu, h);

        // Measure distortion using SATD
        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            // review: call directly, don't use Compute<>?
            satd[cIdx] = Compute<H, Satd, prediction_unit>::go(h, pu, cIdx);
        }
    }

    {
        // measure PU's bitrate
        auto m = h.template change<Measure<void>>();
        Syntax<prediction_unit>::go(pu, m);
    }

    Lambda lambda;
    lambda.set(getReciprocalSqrtLambda(h));
    StateEncode *stateEncode = h;
    if(stateEncode->useRateControl)
    {
        StateEncodePicture *stateEncodePicture = h;
        int segmentPoc = stateEncodePicture->docket->segmentPoc;
        stateEncode->rateControlParams->takeTokenOnRCLevel();
        double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalSqrtLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
        stateEncode->rateControlParams->releaseTokenOnRCLevel();
        lambda.set(value);
    }
    StateEstimateRate *stateEstimateRate = h;

    return stateEstimateRate->rate + (satd[0] + satd[1] + satd[2]) * lambda;
}


// Finds best configuration for the given prediction unit
// - L0/L1/BI
// - reference pictures
// - consider merge / merge index
// - perform motion estimation
// Currently uses SATD as distortion metric
// Returns with
// - best configuration in CodedData sequence
// - contexts and rate updated accordingly
// - Snake<BlockData> updated accordingly
template <>
struct Search<prediction_unit>
{
    template <typename Sample>
    struct State
    {
        Neighbourhood *neighbourhood;
        Snake<BlockData>::Cursor *cursor;
        StateCodedData *stateCodedData;
        StateEncodeSubstream<Sample> *stateEncodeSubstream;
        coding_quadtree const *cqt;

        CodedData::Type bestCodedDataMerge;
        CodedData::Type bestCodedDataPuBuffer[10];

        Cost bestCost;
        Cost bestCostUni[2];
        CodedData::Type codedDataUniBest[2][3];

        ContextsAndCost originalContextsAndCost;
        ContextsAndCost bestContextsAndCost;

        template <class H> State(H &h)
            :
            neighbourhood(h),
            cursor(h),
            stateCodedData(h),
            stateEncodeSubstream(h),
            cqt(h)
        {
            bestCost = std::numeric_limits<Cost>::max();
            bestCostUni[0] = std::numeric_limits<Cost>::max();
            bestCostUni[1] = std::numeric_limits<Cost>::max();
        }

        template <class H> void searchMergeMode(int i, const prediction_unit &pu, H &h)
        {
            stateCodedData->codedCu.word1().merge[stateEncodeSubstream->partIdx] = i + 1;
            *static_cast<ContextsAndCost *>(h) = originalContextsAndCost;
            compare(h, measurePuCost(h));
        }

        template <class H> void searchMergeModes(const prediction_unit &pu, H &h)
        {
            populateMergeCandidates(h, pu);
            for (int i = 0; i < h[MaxNumMergeCand()]; ++i)
            {
                searchMergeMode(i, pu, h);
            }
        }

        template <class H> void searchUni(int refList, int refIdx, const prediction_unit &pu, H &h)
        {
            stateCodedData->codedCu.word1().merge[stateEncodeSubstream->partIdx] = 0;
            stateCodedData->codedPu.init();
            stateCodedData->codedPu.word0().metadata[refList].mvp_lX_flag = 0;
            stateCodedData->codedPu.word0().metadata[refList].ref_idx_lX = refIdx;
            stateCodedData->codedPu.word0().metadata[refList].predFlag = 1;
            stateCodedData->codedPu.mvd(refList) = MotionVector{ 0,0 };

            predictMvp(h, refList, refIdx);

            searchMotionUni(h, refList);

            *static_cast<ContextsAndCost *>(h) = originalContextsAndCost;
            Cost cost = measurePuCost(h);
            if (cost < bestCostUni[refList])
            {
                bestCostUni[refList] = cost;
                codedDataUniBest[refList][0] = stateCodedData->codedPu.p[0];
                codedDataUniBest[refList][1] = stateCodedData->codedPu.p[1];
                codedDataUniBest[refList][2] = stateCodedData->codedPu.p[2];
            }

            compare(h, cost);
        }

        template <class H> void searchBi(const prediction_unit &pu, H &h)
        {
            stateCodedData->codedCu.word1().merge[stateEncodeSubstream->partIdx] = 0;

            stateCodedData->codedPu.p[0] = codedDataUniBest[L0][0] | codedDataUniBest[L1][0];
            stateCodedData->codedPu.p[1] = codedDataUniBest[L0][1];
            stateCodedData->codedPu.p[2] = codedDataUniBest[L0][2];
            stateCodedData->codedPu.p[3] = codedDataUniBest[L1][1];
            stateCodedData->codedPu.p[4] = codedDataUniBest[L1][2];

            // review: prediction could be done more efficiently, one-off per prediction_unit?
            populatePuPredictors(h, pu, stateEncodeSubstream->partIdx);

            if (h[mvd_l1_zero_flag()])
            {
                // review: this decision based on the closest integer to each predictor - consider using accurate MVP costs
                auto const bestMvpL1Flag = stateEncodeSubstream->costMvdZero[L1][1] < stateEncodeSubstream->costMvdZero[L1][0];

                stateCodedData->codedPu.word0().metadata[L1].mvp_lX_flag = bestMvpL1Flag;
                stateCodedData->codedPu.mvd(L1) = MotionVector{ 0, 0 };

                searchMotionBi(h, L0);
            }
            else
            {
                searchMotionBi(h, L0);
                searchMotionBi(h, L1);
            }

            *static_cast<ContextsAndCost *>(h) = originalContextsAndCost;
            compare(h, measurePuCost(h));
        }

        template <class H> void compare(H &h, Cost cost)
        {
            if (cost < bestCost)
            {
                bestCost = cost;
                bestContextsAndCost = *static_cast<ContextsAndCost *>(h);
                bestCodedDataPuBuffer[0] = stateCodedData->codedPu.p[0];
                bestCodedDataPuBuffer[1] = stateCodedData->codedPu.p[1];
                bestCodedDataPuBuffer[2] = stateCodedData->codedPu.p[2];
                bestCodedDataPuBuffer[3] = stateCodedData->codedPu.p[3];
                bestCodedDataPuBuffer[4] = stateCodedData->codedPu.p[4];
                bestCodedDataMerge = stateCodedData->codedCu.word1().raw;
            }
        }

        template <class H> void go2(const prediction_unit &pu, H &h)
        {
            //if (match(pu, h)) logEncode() << pu << "\n";

            cursor->relocate(neighbourhood->snake, pu, neighbourhood->MinCbLog2SizeYMinus1);

            bool debug = false;//pu == prediction_unit(116, 0, 12, 16) && h[PicOrderCntVal()] == 8;

            originalContextsAndCost = *static_cast<ContextsAndCost *>(h);


            if (std::is_same<typename H::Tag, SearchMerge2Nx2N<void>>::value)
            {
                assert(h[PartMode()] == PART_2Nx2N);
                populateMergeCandidates(h, pu);
                StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
                Snake<BlockData>::Cursor *cursor = h;
                prediction_unit const pu = *static_cast<prediction_unit *>(h);
                BlockData &puData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);
                populatePuPredictors(h, pu, stateEncodeSubstream->partIdx);
                setPuData(h, puData);
                {
                    // measure PU's bitrate
                    auto m = h.template change<Measure<void>>();
                    Syntax<prediction_unit>::go(pu, m);
                }
            }
            else
            {
                if (h[PartMode()] != PART_2Nx2N)
                {
                    searchMergeModes(pu, h);
                }

                stateEncodeSubstream->costMvdZero[L0][0] = std::numeric_limits<Cost>::max();
                stateEncodeSubstream->costMvdZero[L0][1] = std::numeric_limits<Cost>::max();
                stateEncodeSubstream->costMvdZero[L1][0] = std::numeric_limits<Cost>::max();
                stateEncodeSubstream->costMvdZero[L1][1] = std::numeric_limits<Cost>::max();

                if (h[RefPicList(L0)][0].dp) searchUni(L0, 0, pu, h);
                if (h[RefPicList(L1)][0].dp) searchUni(L1, 0, pu, h);

                if (pu.nPbW + pu.nPbH != 12 &&
                    bestCostUni[L0] != std::numeric_limits<Cost>::max() &&
                    bestCostUni[L1] != std::numeric_limits<Cost>::max())
                {
                    searchBi(pu, h);
                }
                *static_cast<ContextsAndCost *>(h) = bestContextsAndCost;
                stateCodedData->codedPu.p[0] = bestCodedDataPuBuffer[0];
                stateCodedData->codedPu.p[1] = bestCodedDataPuBuffer[1];
                stateCodedData->codedPu.p[2] = bestCodedDataPuBuffer[2];
                stateCodedData->codedPu.p[3] = bestCodedDataPuBuffer[3];
                stateCodedData->codedPu.p[4] = bestCodedDataPuBuffer[4];
                stateCodedData->codedCu.word1().raw = bestCodedDataMerge;


                // review: need BI refinement here
            }



            cursor->relocate(neighbourhood->snake, pu, h[MinCbLog2SizeY()] - 1);
            BlockData &puData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);

            populatePuPredictors(h, pu, stateEncodeSubstream->partIdx);
            setPuData(h, puData);

            // Recreate the best prediction unit as just found in search

            struct PredAccessor
            {
                StateEncodeSubstream<Sample> *stateEncodeSubstream;
                Candidate<Sample> *candidate;
                coding_quadtree *cqt;
                Raster<Sample> operator()(int x0, int y0, int cIdx)
                {
                    auto predSamples = candidate->stateReconstructionCache->components[cIdx].get(stateEncodeSubstream->interPieces[cIdx][0]);
                    return predSamples.offset((x0 - cqt->x0) >> (cIdx ? 1 : 0), (y0 - cqt->y0) >> (cIdx ? 1 : 0));
                }
            };

            PredAccessor predAccessor;
            predAccessor.stateEncodeSubstream = h;
            predAccessor.candidate = h;
            predAccessor.cqt = h;

        StateReconstructedPicture<Sample> *stateReconstructedPicture = h;

            predictInter(*stateReconstructedPicture->picture, pu, h);

            copyBlock(predAccessor(cqt->x0, cqt->y0, 0), (*stateReconstructedPicture->picture)(cqt->x0, cqt->y0, 0), 1ull << cqt->log2CbSize, 1ull << cqt->log2CbSize);
            copyBlock(predAccessor(cqt->x0, cqt->y0, 1), (*stateReconstructedPicture->picture)(cqt->x0, cqt->y0, 1), 1ull << cqt->log2CbSize >> 1, 1ull << cqt->log2CbSize >> 1);
            copyBlock(predAccessor(cqt->x0, cqt->y0, 2), (*stateReconstructedPicture->picture)(cqt->x0, cqt->y0, 2), 1ull << cqt->log2CbSize >> 1, 1ull << cqt->log2CbSize >> 1);

            if (h[PartMode()] != PART_2Nx2N)
            {
                cursor->commit(pu, h[MinCbLog2SizeY()] - 1);
                neighbourhood->recordMerge(h, pu);
            }

            if (!h[merge_flag(pu.x0, pu.y0)])
            {
                stateCodedData->codedPu.p = stateCodedData->codedPu.next();
            }

            ++stateCodedData->partIdx;
            ++stateEncodeSubstream->partIdx;
        }
    };

    template <class H> static void go(const prediction_unit &pu, H &h)
    {
        using Sample = typename SampleType<H>::Type;
        State<Sample> s(h);
        s.go2(pu, h);
    }
};


// Review: remove k, instead adjust cost after function call
template <class H, typename Sample>
Cost costDistortionMv(H& h, MotionVector mv, Raster<Sample> input, Raster<Sample> reference, double k = 1.0)
{
    prediction_unit const *pu = h;

    HAVOC_ALIGN(32, Sample, buffer[64 * 64]);
    Raster<Sample> prediction{ buffer, 64 };

    auto const mvFrac = mv & 0x3;
    mv >>= 2;
    reference += mv;

    HavocTablePredUni<Sample>  *tablePred = h;
    auto *f = *havocGetPredUni(tablePred, 8, pu->nPbW, pu->nPbH, mvFrac[0], mvFrac[1], h[BitDepthY()]);
    f(prediction.p, prediction.stride, reference.p, reference.stride, pu->nPbW, pu->nPbH, mvFrac[0], mvFrac[1], h[BitDepthY()]);

    havoc_table_hadamard_satd<Sample>  *tableSatd = h;
    auto const distortion = measureSatd(tableSatd, input, prediction, pu->nPbW, pu->nPbH);

    Lambda lambda;
    lambda.set(getReciprocalSqrtLambda(h) * k);
    StateEncode *stateEncode = h;
    if(stateEncode->useRateControl)
    {
        StateEncodePicture *stateEncodePicture = h;
        int segmentPoc = stateEncodePicture->docket->segmentPoc;
        stateEncode->rateControlParams->takeTokenOnRCLevel();
        double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalSqrtLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
        stateEncode->rateControlParams->releaseTokenOnRCLevel();
        lambda.set(value*k);
    }

    return lambda * distortion;
}


// Performs luma motion compensation and measures SATD distortion, returns rate + distortion / lambda.
template <class H, typename Sample>
Cost costMv(H& h, MotionVector mv, MotionVector const &mvd, Raster<Sample> input, Raster<Sample> reference)
{
    return rateOf<H::Mode::value>(mvd) + costDistortionMv(h, mv, input, reference);
}


// Refinement search (half or quarter pixel)
template <typename Sample, class H, class P>
void patternSearch(H &h, mvd_coding mvdc, const P &pattern, bool tryOrigin, MotionVector &mv, MotionVector &mvd, Cost &bestCost, Raster<Sample> input, Raster<Sample> reference, int maxIterations = 0)
{
    //if (match(mvdc, h)) logEncode() << "pattern[" << boost::size(pattern) << "] ";

    int min = 0;
    int max = static_cast<int>(boost::size(pattern));
    bool improved;

    if (tryOrigin)
    {
        bestCost = costMv(h, mv, mvd, input, reference);
        //if (match(mvdc, h)) logEncode() << " [" << mvd << "]=" << bestCost;
    }

    do
    {
        Cost const prevCost = bestCost;

        int best;
        improved = false;

        for (int i = min; i < max; ++i)
        {
            MotionVector mvdTest = mvd;
            MotionVector mvTest = mv;
            mvdTest += pattern[i % boost::size(pattern)];
            mvTest += pattern[i % boost::size(pattern)];
            const Cost cost = costMv(h, mvTest, mvdTest, input, reference);
            //if (match(mvdc, h)) logEncode() << " " << mvdTest << "=" << cost;
            if (cost < bestCost)
            {
                best = i  % boost::size(pattern);
                bestCost = cost;
                improved = true;
            }
        }

        if (improved)
        {
            mvd += pattern[best];
            mv += pattern[best];
            min = best - 1;
            max = best + 2;
            //if (match(mvdc, h)) logEncode() << " (" << mvd << ")";
        }

    } while (improved && --maxIterations);

    //if (match(mvdc, h)) logEncode() << "\n";
}


// Integer motion estimation
template <class H> static void fullPelMotionEstimation(mvd_coding mvdc, H &h, MvCandidate &best)
{
    using Sample = typename SampleType<H>::Type;

    auto const refIdx = 0;

    Profiler::Timers *timers = h;
    StateCodedData *stateCodedData = h;
    Mvp::Predictors *predictors = h;
    auto *stateEncodeSubstream = &h[Concrete<StateSubstream>()];
    StateEncode *stateEncode = h;
    Speed *speed = h;
    Profiler::Scope scope(timers->searchMotionFullPel);

    Lambda lambda;
    lambda.set(getReciprocalSqrtLambda(h));

    if(stateEncode->useRateControl)
    {
        StateEncodePicture *stateEncodePicture = h;
        int segmentPoc = stateEncodePicture->docket->segmentPoc;
        stateEncode->rateControlParams->takeTokenOnRCLevel();
        double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalSqrtLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
        stateEncode->rateControlParams->releaseTokenOnRCLevel();
        lambda.set(value);
    }

    int searchWindow = (speed->useSmallSearchWindow()) ? 32 : 64;
    int maxCounter = (speed->useSmallSearchWindow()) ? 2 : 3;
    int rasterSearch = (speed->useSmallSearchWindow()) ? 120 : 240;
    const prediction_unit &pu = *static_cast<prediction_unit *>(h);
    const coding_quadtree *cqt = h;
    const coding_unit cu(cqt->x0, cqt->y0, cqt->log2CbSize);

    StateMeFullPel<Sample> stateMeFullPel(h, mvdc.refList, pu, best);

    EstimateRateBin<mvp_lX_flag> estimateRateMvp(h, 0);

    {
        // Consider zero MV as a starting point
        // Review: probably should do this *after* MVP's and/ integer2Nx2N candidates because they are more likely to trigger MET.
        MvCandidate candidate(h, mvdc.refList, { 0, 0 }, predictors->mvp[refIdx][mvdc.refList]);
        auto const sad = stateMeFullPel.functionSad(stateMeFullPel.src.p, stateMeFullPel.src.stride, &stateMeFullPel.ref(0, 0), stateMeFullPel.ref.stride, HAVOC_RECT(pu.nPbW, pu.nPbH));

        candidate.cost += lambda * sad;

        bool const better = best.consider(candidate);
        assert(better);
        if (better && stateEncode->met)
        {
            static const MotionVector diamond[4] = { { -4, 0 },{ 0, 4 },{ 4, 0 },{ 0, -4 } };
            bool triggerMet = !stateMeFullPel.considerPattern(h, best.mv, diamond, 4, 1, 1, mvdc);
            if (triggerMet && cu.log2CbSize >= 5)
            {
                // last two points of hexagon repeated because considerPattern likely only supports multiple-of-four patterns
                static const MotionVector hexagon[8] = { { 0, -8 },{ 8, -4 },{ 8, 4 },{ 0, 8 },{ -8, 4 },{ -8,-4 },{ -8, 4 },{ -8,-4 } };
                triggerMet = !stateMeFullPel.considerPattern(h, best.mv, hexagon, 8, 1, 1, mvdc);
            }
            if (triggerMet)
                return;
        }

        //if (match(mvdc, h)) logEncode() << candidate.mv << " ";
    }

    // Choose the best starting MVP from the two predictors
    MvCandidate candidate;
    for (candidate.mvpFlag = 0; candidate.mvpFlag < 2; ++candidate.mvpFlag)
    {
        MotionVector const &predictedVector = predictors->mvp[refIdx][mvdc.refList][candidate.mvpFlag];

        candidate.mv = (predictedVector + 1) >> 2;

        stateMeFullPel.limit(candidate.mv);

        candidate.mv <<= 2;

        // delta will be small, just the sub-pixel adjustment
        candidate.mvd = candidate.mv - predictedVector;

        candidate.cost = rateOf<H::Mode::value>(candidate.mvd);
        candidate.cost += estimateRateMvp.rate(candidate.mvpFlag);

        auto const sad = stateMeFullPel.functionSad(stateMeFullPel.src.p, stateMeFullPel.src.stride, &stateMeFullPel.ref(candidate.mv[0] >> 2, candidate.mv[1] >> 2), stateMeFullPel.ref.stride, HAVOC_RECT(pu.nPbW, pu.nPbH));

        candidate.cost += lambda * sad;

        stateEncodeSubstream->costMvdZero[mvdc.refList][candidate.mvpFlag] = candidate.cost;

        bool const better = best.consider(candidate);

        if (better && stateEncode->met)
        {
            static const MotionVector diamond[4] = { { -4, 0 },{ 0, 4 },{ 4, 0 },{ 0, -4 } };
            bool triggerMet = !stateMeFullPel.considerPattern(h, best.mv, diamond, 4, 1, 1, mvdc);
            if (triggerMet && cu.log2CbSize >= 5)
            {
                // last two points of hexagon repeated because considerPattern likely only supports multiple-of-four patterns
                static const MotionVector hexagon[8] = { { 0, -8 },{ 8, -4 },{ 8, 4 },{ 0, 8 },{ -8, 4 },{ -8,-4 },{ -8, 4 },{ -8,-4 } };
                triggerMet = !stateMeFullPel.considerPattern(h, best.mv, hexagon, 8, 1, 1, mvdc);
            }
            if (triggerMet)
                return;
        }
    }

    if (h[PartMode()] != PART_2Nx2N || cqt->cqtDepth != 0)
    {
        // Consider previous integer PART_2Nx2N best vector as starting vector
        auto mv = stateEncodeSubstream->mvPreviousInteger2Nx2N[mvdc.refList];
        assert(!(mv % 4));
        mv >>= 2;
        stateMeFullPel.limit(mv);
        mv <<= 2;
        MvCandidate candidate(h, mvdc.refList, mv, predictors->mvp[refIdx][mvdc.refList]);
        auto const sad = stateMeFullPel.functionSad(stateMeFullPel.src.p, stateMeFullPel.src.stride, &stateMeFullPel.ref(candidate.mv[0] >> 2, candidate.mv[1] >> 2), stateMeFullPel.ref.stride, HAVOC_RECT(pu.nPbW, pu.nPbH));
        candidate.cost += lambda * sad;
        bool const better = best.consider(candidate);

        if (better && stateEncode->met)
        {
            static const MotionVector diamond[4] = { { -4, 0 },{ 0, 4 },{ 4, 0 },{ 0, -4 } };
            bool triggerMet = !stateMeFullPel.considerPattern(h, best.mv, diamond, 4, 1, 1, mvdc);
            if (triggerMet && cu.log2CbSize >= 5)
            {
                // last two points of hexagon repeated because considerPattern likely only supports multiple-of-four patterns
                static const MotionVector hexagon[8] = { { 0, -8 },{ 8, -4 },{ 8, 4 },{ 0, 8 },{ -8, 4 },{ -8,-4 },{ -8, 4 },{ -8,-4 } };
                triggerMet = !stateMeFullPel.considerPattern(h, best.mv, hexagon, 8, 1, 1, mvdc);
            }
            if (triggerMet)
                return;
        }
    }

    //if (match(mvdc, h)) logEncode() << "starting point: " << best.mv << " mvp=" << best.mvpFlag << "\n";

    // HM style "star" search here

    MotionVector mvStart = stateMeFullPel.best.mv;

    int distBest = 0;
    int counter = 0;
    int step = 4;

    static const MotionVector diamond[16] = {
        { 0, -4 },{ 1, -3 },{ 2, -2 },{ 3, -1 },
        { 4, 0 },{ 3, 1 },{ 2, 2 },{ 1, 3 },
        { 0, 4 },{ -1, 3 },{ -2, 2 },{ -3, 1 },
        { -4, 0 },{ -3, -1 },{ -2, -2 },{ -1, -3 },
    };

    static const MotionVector square[16] = {
        { 4, -4 },{ 4, -2 },{ 4, 0 },{ 4, 2 },
        { 4, 4 },{ 2, 4 },{ 0, 4 },{ -2, 4 },
        { -4, 4 },{ -4, 2 },{ -4, 0 },{ -4, -2 },
        { -4, -4 },{ -2, -4 },{ 0, -4 },{ 2, -4 },
    };

    static const MotionVector square4[4] = { { -4, -4 },{ -4, 4 },{ 4, 4 },{ 4, -4 } };

    for (int dist = 1; dist <= searchWindow && counter < maxCounter; dist <<= 1)
    {
        if (dist == 2 || dist == 8)
        {
            step >>= 1;
        }

        if (stateMeFullPel.considerPattern(h, mvStart, diamond, 16, step, dist, mvdc))
        {
            distBest = dist;
            counter = 0;
        }
        else
        {
            ++counter;
        }

        //if (match(mvdc, h)) logEncode() << "initial dist " << dist << ": " << best.mv << " mvp=" << best.mvpFlag << "\n";
    }

    // fill some gaps
    if (distBest == 1)
    {
        // review: have checked two of these points already
        distBest = 0;
        stateMeFullPel.considerPattern(h, stateMeFullPel.best.mv, square4, 4, 1, 1, mvdc);
    }

    // raster refinement
    if (distBest > 5)
    {
        static const MotionVector line[4] = { { 0, 0 },{ 1, 0 },{ 2, 0 },{ 3, 0 } };
        MotionVector mv;
        for (mv[1] = -rasterSearch; mv[1] <= rasterSearch; mv[1] += 20)
        {
            for (mv[0] = -rasterSearch; mv[0] <= rasterSearch; mv[0] += 80)
            {
                stateMeFullPel.considerPattern(h, mv, line, 4, 1, 20, mvdc);
            }
        }
        distBest = 5;

        //if (match(mvdc, h)) logEncode() << "raster refinement: " << best.mv << " mvp=" << best.mvpFlag << "\n";
    }

    // star refinement
    while (distBest > 0)
    {
        MotionVector mvStart = stateMeFullPel.best.mv;

        distBest = 0;
        step = 4;
        for (int dist = 1; dist <= searchWindow; dist <<= 1)
        {
            if (dist == 2 || dist == 8)
            {
                step >>= 1;
            }

            if (stateMeFullPel.considerPattern(h, mvStart, diamond, 16, step, dist, mvdc))
            {
                distBest = dist;
            }

            //if (match(mvdc, h)) logEncode() << "star refinement " << dist << ": " << best.mv << " mvp=" << best.mvpFlag << "\n";
        }

        if (distBest == 1)
        {
            // review: two of these positions already tested
            stateMeFullPel.considerPattern(h, mvStart, square4, 4, 1, 1, mvdc);
            distBest = 0;
        }
    }
    if (!(speed->useSmallSearchWindow()))
    {
        int j;
        do
        {
            static const MotionVector diamond[4] = { { 0, -1 }, { -1, 0 }, { 0, 1 }, { 1, 0 } };

            MotionVector mv[4];
            Sample const *refs[4];
            for (int i = 0; i < 4; ++i)
            {
                mv[i] = best.mv / 4 + diamond[i];
                stateMeFullPel.limit(mv[i]);
                refs[i] = &stateMeFullPel.ref(mv[i][0], mv[i][1]);
            }

            int32_t sads[4];
            stateMeFullPel.functionSad4(stateMeFullPel.src.p, stateMeFullPel.src.stride, refs, stateMeFullPel.ref.stride, sads, HAVOC_RECT(pu.nPbW, pu.nPbH));

            j = -1;
            for (int i = 0; i < 4; ++i)
            {
                mv[i][0] *= 4;
                mv[i][1] *= 4;
                MvCandidate temp(h, mvdc.refList, mv[i], predictors->mvp[refIdx][mvdc.refList]);
                temp.cost += lambda * sads[i];
                if (best.consider(temp))
                {
                    j = i;
                }
            }
        } while (j >= 0);
    }

    if (h[PartMode()] == PART_2Nx2N)
    {
        stateEncodeSubstream->mvPreviousInteger2Nx2N[mvdc.refList] = best.mv;
    }
}


template <typename Sample, class H>
void subPelRefinement(MotionVector &mv, MotionVector &mvd, mvd_coding const &mvdc, H &h, Raster<Sample> input, Raster<Sample> reference)
{
    Profiler::Timers *timers = h;
    Cost bestCost;
    Speed *speed = h;

    {
        Profiler::Scope scope(timers->searchMotionHalfPel);
        const MotionVector diamond[8] = { { -2, -2 },{ 0, -2 },{ 2, -2 },{ -2, 0 },{ 2, 0 },{ -2, 2 },{ 0, 2 },{ 2, 2 } };
        patternSearch(h, mvdc, diamond, true, mv, mvd, bestCost, input, reference, 1);
    }
    if (speed->doQuarterPelRefinement())
    {
        Profiler::Scope scope(timers->searchMotionQuarterPel);
        const MotionVector diamond[8] = { { -1, -1 },{ 0, -1 },{ 1, -1 },{ -1, 0 },{ 1, 0 },{ -1, 1 },{ 0, 1 },{ 1, 1 } };
        patternSearch(h, mvdc, diamond, false, mv, mvd, bestCost, input, reference, 1);
    }
}


// when 'Searching', override the value of rqt_root_cbf and always call transform_tree
// review: should this code be shared with the intra search? (currently this is inter only)
template <>
struct Search<IfCbf<rqt_root_cbf, transform_tree>>
{
    template <class H> static void go(const IfCbf<rqt_root_cbf, transform_tree> &i, H &h)
    {
        using Sample = typename SampleType<H>::Type;

        transform_tree const &tt = i.f;
        StateEncodeSubstream<Sample> *stateEncodeSubstream = h;
        StateCodedData *stateCodedData = h;
        Candidate<Sample> *candidate = h;
        coding_quadtree const *cqt = h;

        assert(h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTER);

        // copy of coded data pointers as they are before encoding the transform tree.
        StateCodedData codedDataBefore = *stateCodedData;

        stateEncodeSubstream->ssdPrediction[0] = 0;
        stateEncodeSubstream->ssdPrediction[1] = 0;
        stateEncodeSubstream->ssdPrediction[2] = 0;
        // encode transform tree, measure distortion (SSD) and create corresponding CodedData
        auto const ssd = reconstructInter(tt, h);
        int rqtdepth = candidate->rqtdepth;
        Lambda reciprocalLambda = getReciprocalLambda(h);
        StateEncode *stateEncode = h;
        if(stateEncode->useRateControl)
        {
            StateEncodePicture *stateEncodePicture = h;
            int segmentPoc = stateEncodePicture->docket->segmentPoc;
            stateEncode->rateControlParams->takeTokenOnRCLevel();
            double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[ CtbAddrInRs()]);
            stateEncode->rateControlParams->releaseTokenOnRCLevel();
            reciprocalLambda.set(value);
        }
        candidate->lambdaDistortion += ssd * reciprocalLambda;
        StateCodedData codedDataAfter = *stateCodedData;

        // rewind CodedData pointers to start of CU in preparation for measuring rate
        *stateCodedData = codedDataBefore;
        stateCodedData->partIdx = 0;
        // review: sort out how snake is positioned / moves through this function
        // cu_skip_flag needs top/left CU neighbours
        // at entry to this function, snake is bottom/right
        // probably snake can remain top/left of CU throughout
        stateEncodeSubstream->originalCandidate->snake.checkPosition(before, *cqt);
        auto const snakePointerBefore = candidate->snake;
        candidate->snake = stateEncodeSubstream->originalCandidate->snake;
        candidate->relocate(candidate->snake, *cqt, candidate->MinCbLog2SizeYMinus1);

        bool haveResidual;
        auto m = h.template change<EstimateRate<void>>();
        // Measure rate of transform tree
        haveResidual = m[rqt_root_cbf()];

        bool const isMerge2Nx2N =
            h[PartMode()] == PART_2Nx2N &&
            h[merge_flag(cqt->x0, cqt->y0)];

        bool const encodeAsSkip = isMerge2Nx2N && !haveResidual;
        if (encodeAsSkip)
        {
            stateCodedData->codedCu.word0().CuPredMode = MODE_SKIP;

            Neighbourhood *neighbourhood = h;
            Snake<BlockData>::Cursor *cursor = h;
            BlockData &blockData = cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1);

            blockData.skip = true;

            // Reset contexts as before CU (review: overkill, only perhaps a couple of contexts were modified)
            *static_cast<Contexts *>(h) = *stateEncodeSubstream->originalCandidate;

            // Reset rate as before CU
            *static_cast<StateEstimateRate *>(h) = *stateEncodeSubstream->originalCandidate;

            // measure rate of skip syntax (review: possible duplication)
            if (m[MaxNumMergeCand()] > 1)
                m(merge_idx(cqt->x0, cqt->y0), ae(v));
        }

        // measure rate of CU preamble elements now that we know their values
        m(cu_skip_flag(cqt->x0, cqt->y0), ae(v));
        if (!h[current(cu_skip_flag(cqt->x0, cqt->y0))])
        {
            m(pred_mode_flag(), ae(v));
            m(part_mode(), ae(v));
        }

        ContextsAndCost snapshot;
        if (haveResidual && !isMerge2Nx2N) snapshot = *static_cast<ContextsAndCost *>(m);

        // mesaure rate of rqt_root_cbf() and transform_tree()
        if (!isMerge2Nx2N) m(rqt_root_cbf(), ae(v));
        if (haveResidual) m(tt);

        int const predicted = 0;
        int const reconstructed = 1;
        int keep = reconstructed;

        if (haveResidual)
        {
            Neighbourhood *neighbourhood = h;
            Snake<BlockData>::Cursor *cursor = h;
            BlockData &blockData = cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1);

            auto backupContextsAndCost = *static_cast<ContextsAndCost *>(candidate);
            auto backupCuWord0 = candidate->codedCu.word0().raw;
            auto backupCursor = *cursor;
            auto backupBlockData = blockData;

            assert(candidate->codedCu.word0().cbfWord);
            candidate->codedCu.word0().cbfWord = 0;
            assert(!h[rqt_root_cbf()]);
            assert(!h[cu_skip_flag(cqt->x0, cqt->y0)]);

            auto m = h.template change<EstimateRate<void>>();

            if (isMerge2Nx2N)
            {
                // measure rate of CU without residual (MODE_SKIP)
                *static_cast<ContextsAndCost *>(candidate) = *stateEncodeSubstream->originalCandidate;

                candidate->codedCu.word0().CuPredMode = MODE_SKIP;
                blockData.skip = true;
                assert(h[cu_skip_flag(cqt->x0, cqt->y0)]);

                // Skipped CU syntax is very simple:
                m(cu_skip_flag(cqt->x0, cqt->y0), ae(v));
                if (m[MaxNumMergeCand()] > 1)
                    m(merge_idx(cqt->x0, cqt->y0), ae(v));
            }
            else
            {
                // measure rate of CU without residual (rqt_root_cbf=0)
                *static_cast<ContextsAndCost *>(candidate) = snapshot;
                m(rqt_root_cbf(), ae(v));
            }

            // already have distortion measure
            Lambda lambda = getReciprocalLambda(h);
            if(stateEncode->useRateControl)
            {
                StateEncodePicture *stateEncodePicture = h;
                int segmentPoc = stateEncodePicture->docket->segmentPoc;
                stateEncode->rateControlParams->takeTokenOnRCLevel();
                double value = stateEncode->rateControlMap.find(segmentPoc)->second->getCtbReciprocalLambda(h[PicOrderCntVal()], h[CtbAddrInRs()]);
                stateEncode->rateControlParams->releaseTokenOnRCLevel();
                lambda.set(value);
            }
            candidate->lambdaDistortion = stateEncodeSubstream->originalCandidate->lambdaDistortion;
            candidate->lambdaDistortion += lambda * stateEncodeSubstream->ssdPrediction[0];
            candidate->lambdaDistortion += lambda * stateEncodeSubstream->ssdPrediction[1];
            candidate->lambdaDistortion += lambda * stateEncodeSubstream->ssdPrediction[2];

            if (candidate->cost2() < backupContextsAndCost.cost2())
            {
                // CU is cheaper with no residual: delete the coded residual
                candidate->codedDataAfter = candidate->codedCu.firstPredictionUnit().p;
                keep = predicted;
            }
            else
            {
                // restore everything as it was before trying with no residual
                blockData = backupBlockData;
                *cursor = backupCursor;
                *static_cast<ContextsAndCost *>(candidate) = backupContextsAndCost;
                candidate->codedCu.word0().raw = backupCuWord0;
                *static_cast<Snake<BlockData>::Cursor *>(candidate) = backupCursor;
            }
        }

        candidate->snake = snakePointerBefore;
        int keep3, flush1, flush2;
        if (keep == 0)
        {
            keep3 = 0;
            flush1 = 1;
            flush2 = 2;

        }
        else if(rqtdepth == 0)
        {
            flush1 = 0;
            keep3 = 1;
            flush2 = 2;
        }
        else if (rqtdepth == 1)
        {
            flush1 = 0;
            flush2 = 1;
            keep3 = 2;
        }
        stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[0].freeBlock(stateEncodeSubstream->interPieces[0][flush1]);
        stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[1].freeBlock(stateEncodeSubstream->interPieces[1][flush1]);
        stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[2].freeBlock(stateEncodeSubstream->interPieces[2][flush1]);

        stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[0].freeBlock(stateEncodeSubstream->interPieces[0][flush2]);
        stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[1].freeBlock(stateEncodeSubstream->interPieces[1][flush2]);
        stateEncodeSubstream->originalCandidate->stateReconstructionCache->components[2].freeBlock(stateEncodeSubstream->interPieces[2][flush2]);


        candidate->appendPiece(0, residual_coding((*cqt).x0, (*cqt).y0, (*cqt).log2CbSize, 0), uint16_t((*cqt).log2CbSize), stateEncodeSubstream->interPieces[0][keep3].i);
        candidate->appendPiece(1, residual_coding((*cqt).x0, (*cqt).y0, (*cqt).log2CbSize - 1, 1), uint16_t((*cqt).log2CbSize - 1), stateEncodeSubstream->interPieces[1][keep3].i);
        candidate->appendPiece(2, residual_coding((*cqt).x0, (*cqt).y0, (*cqt).log2CbSize - 1, 2), uint16_t((*cqt).log2CbSize - 1), stateEncodeSubstream->interPieces[2][keep3].i);
    }
};
