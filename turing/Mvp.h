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

#ifndef INCLUDED_Mvp_h
#define INCLUDED_Mvp_h

#pragma once

#include "Global.h"
#include "Dsp.h"
#include "BlockData.h"
#include "HevcMath.h"
#include "StatePicture.h"
#include "GlobalState.h"


static MotionVector distScale(MotionVector mv, int tD, int tB)
{
    int const td = Clip3(-128, 127, tD);
    int const tb = Clip3(-128, 127, tB);
    int const tx = (16384 + (Abs(td) >> 1)) / td;
    int const distScaleFactor = Clip3(-4096, 4095, (tb * tx + 32) >> 6);
    mv[0] = Clip3(-32768, 32767, Sign(distScaleFactor * mv[0])*((Abs(distScaleFactor * mv[0]) + 127) >> 8));
    mv[1] = Clip3(-32768, 32767, Sign(distScaleFactor * mv[1])*((Abs(distScaleFactor * mv[1]) + 127) >> 8));
    return mv;
}


// returns availableFlagLXCol
template <class H>
bool deriveCollocatedMotionVectors(H &h2, PuData &puData, StateCollocatedMotion *motion, int refList, int x, int y)
{
    MotionVector &mvLXCol = puData.mv(refList);
    int const refIdxLX = puData.refIdx(refList);

    if (!motion)
    {
        h2(Violation("7.4.7.1", "unavailable temporal motion predictor")); // CondCheck 7.4.7.1-AA, 7.4.7.1-AB
        throw Abort();
    }

    // Review: this is opaque and unnecessary?
    auto h = h2.extend(&*motion);

    x = (x >> 4) << 4;
    y = (y >> 4) << 4;

    const StatePicture &currPic = *static_cast<StatePicture *>(h);

    bool availableFlagLXCol;

    StateCollocatedMotion &stateCollocatedMotion = *static_cast<StateCollocatedMotion *>(h);
    const PuData &puDataCol = stateCollocatedMotion(x, y);

    const bool colPbIsIntra = !puDataCol.isAvailable() || (puDataCol.predFlag(L0) == 0 && puDataCol.predFlag(L1) == 0);
    if (colPbIsIntra)
    {
        mvLXCol = MotionVector{ 0, 0 };
        availableFlagLXCol = false;
    }
    else
    {
        int listCol;
        if (!puDataCol.predFlag(L0))
            listCol = L1;
        else if (!puDataCol.predFlag(L1))
            listCol = L0;
        else
        {
            assert(puDataCol.predFlag(L0) && puDataCol.predFlag(L1));
            if (static_cast<StatePicture *>(h)->allBackwards)
                listCol = refList;
            else
                listCol = h[collocated_from_l0_flag()] ? L1 : L0;
        }

        const MotionVector mvCol = puDataCol.mv(listCol);
        int const refIdxCol = listCol ? puDataCol.refIdx(L1) : puDataCol.refIdx(L0);

        const bool longTermRefPicCol = stateCollocatedMotion.getReference(puDataCol, listCol) == LONG_TERM;

        if (LongTermRefPic(h[RefPicList(refList)][refIdxLX]) != longTermRefPicCol)
        {
            availableFlagLXCol = false;
            mvLXCol = MotionVector{ 0, 0 };
        }
        else
        {
            availableFlagLXCol = true;

            int const colPocDiff = DiffPicOrderCnt(h, motion->poc, stateCollocatedMotion.getPoc(puDataCol, listCol));
            int const currPocDiff = DiffPicOrderCnt(h, h[PicOrderCntVal()], h[RefPicList(refList)][refIdxLX]);

            if (LongTermRefPic(h[RefPicList(refList)][refIdxLX]) || colPocDiff == currPocDiff)
                mvLXCol = mvCol;
            else
                mvLXCol = distScale(mvCol, colPocDiff, currPocDiff);
        }
    }
    if (!availableFlagLXCol)
        puData.setRefIdx(h, refList, -1);
    return availableFlagLXCol;
}


template <class H>
const StatePicture *getColPic(H &h, int *poc = 0)
{
    if (h[collocated_ref_idx()] >= 0 || h[collocated_ref_idx()] < 16)
    {
        auto &pic = h[RefPicList(h[collocated_from_l0_flag()] ? 0 : 1)][h[collocated_ref_idx()]];
        if (pic.dp)
        {
            if (poc) 
                *poc = (*pic.dp)[PicOrderCntVal()];
            return pic.dp.get();
        }
    }
    return 0;
}


// returns availableFlagLXCol
template <class H>
bool deriveTemporalLumaMotionVectorPredictors(H &h, PuData &puDataCol, prediction_unit const &pu, int refList, int refIdxLX)
{
    auto const &xPb = pu.x0;
    auto const &yPb = pu.y0;
    auto const &nPbW = pu.nPbW;
    auto const &nPbH = pu.nPbH;

    assert(h[slice_temporal_mvp_enabled_flag()]);

    int const poc = h[PicOrderCntVal()];

    // 1. The variable colPic, specifying the collocated picture, is derived as follows:
    const StatePicture *colPic = getColPic(h);

    //  2. The bottom right collocated motion vector is derived as follows:
    int const xColBr = xPb + nPbW;
    int const yColBr = yPb + nPbH;

    bool availableFlagLXCol = false; // omission in standard?

    if (((yPb >> h[CtbLog2SizeY()]) == (yColBr >> h[CtbLog2SizeY()]))
        && (yColBr < h[pic_height_in_luma_samples()])
        && (xColBr < h[pic_width_in_luma_samples()]))
    {
        puDataCol.setRefIdx(h, refList, refIdxLX);
        availableFlagLXCol = deriveCollocatedMotionVectors(h, puDataCol, colPic->motion.get(), refList, xColBr, yColBr);
    }

    // 3.
    if (!availableFlagLXCol)
    {
        // The central collocated motion vector is derived as follows
        int const xColCtr = xPb + (nPbW >> 1);
        int const yColCtr = yPb + (nPbH >> 1);

        puDataCol.setRefIdx(h, refList, refIdxLX);
        availableFlagLXCol = deriveCollocatedMotionVectors(h, puDataCol, colPic->motion.get(), refList, xColCtr, yColCtr);
    }

    return availableFlagLXCol;
}


struct MotionCand
{
    int refList;
    const char *name;
};

template <> struct Syntax<MotionCand> : Null<MotionCand> {};


// Computes the two MVP predictors
template <class H>
void predictMvp(H &h, int refList, int refIdxLX)
{
    int const X = refList;
    int const Y = X ? L0 : L1;

    prediction_unit const &pu = *static_cast<prediction_unit *>(h);
    auto const xPb = pu.x0;
    auto const yPb = pu.y0;

    // Derivation process for luma motion vector prediction
    // The motion vector predictor mvpLX is derived in the following ordered steps:
    // 1.
    MotionVector mvLXA;
    MotionVector mvLXB;
    int availableFlagLXA;
    int availableFlagLXB;
    int kA, kB;

    {
        // Derivation process for motion vector predictor candidates
        availableFlagLXA = 0;
        mvLXA = MotionVector{ 0, 0 };

        const PuData puDataA[2] = {
            neighbourPuData(h, xPb - 1, yPb + pu.nPbH), // A0
            neighbourPuData(h, xPb - 1, yPb + pu.nPbH - 1) // A1
        };

        const bool isScaledFlagLX = puDataA[0].isAvailable() || puDataA[1].isAvailable();

        for (int k = 0; k < 2; ++k)
        {
            if (puDataA[k].isAvailable() && availableFlagLXA == 0)
            {
                if (puDataA[k].predFlag(X) && DiffPicOrderCnt(h, h[RefPicList(X)][puDataA[k].refIdx(X)], h[RefPicList(X)][refIdxLX]) == 0)
                {
                    availableFlagLXA = 1;
                    mvLXA = puDataA[k].mv(X);
                    kA = k;
                }
                else if (puDataA[k].predFlag(Y) && DiffPicOrderCnt(h, h[RefPicList(Y)][puDataA[k].refIdx(Y)], h[RefPicList(X)][refIdxLX]) == 0)
                {
                    availableFlagLXA = 1;
                    mvLXA = puDataA[k].mv(Y);
                    kA = k;
                }
            }
        }

        // 7.
        for (int k = 0; availableFlagLXA == 0 && k < 2; ++k)
        {
            int refIdxA;
            RefPicList RefPicListA;
            if (puDataA[k].isAvailable())
            {
                if (puDataA[k].predFlag(X)
                    && LongTermRefPic(h[RefPicList(X)][refIdxLX]) == LongTermRefPic(h[RefPicList(X)][puDataA[k].refIdx(X)]))
                {
                    availableFlagLXA = 1;
                    mvLXA = puDataA[k].mv(X);
                    refIdxA = puDataA[k].refIdx(X);
                    RefPicListA = RefPicList(X);
                    kA = k;
                }
                else if (puDataA[k].predFlag(Y)
                    && LongTermRefPic(h[RefPicList(X)][refIdxLX]) == LongTermRefPic(h[RefPicList(Y)][puDataA[k].refIdx(Y)]))
                {
                    availableFlagLXA = 1;
                    mvLXA = puDataA[k].mv(Y);
                    refIdxA = puDataA[k].refIdx(Y);
                    RefPicListA = RefPicList(Y);
                    kA = k;
                }
            }
            if (availableFlagLXA == 1
                && DiffPicOrderCnt(h, h[RefPicListA][refIdxA], h[RefPicList(X)][refIdxLX]) != 0
                && h[RefPicListA][refIdxA].reference == SHORT_TERM
                && h[RefPicList(X)][refIdxLX].reference == SHORT_TERM)
            {
                int const tD = DiffPicOrderCnt(h, h[PicOrderCntVal()], h[RefPicListA][refIdxA]);
                int const tB = DiffPicOrderCnt(h, h[PicOrderCntVal()], h[RefPicList(X)][refIdxLX]);
                mvLXA = distScale(mvLXA, tD, tB);
            }
        }

        // The motion vector mvLXB and the availability flag availableFlagLXB are derived in the following ordered steps:
        availableFlagLXB = 0;
        mvLXB = MotionVector{ 0, 0 };

        int const xNbB[3] = { xPb + pu.nPbW, xPb + pu.nPbW - 1, xPb - 1 };
        int const yNbB[3] = { yPb - 1, yPb - 1, yPb - 1 };

        int refIdxB = 0;
        for (int k = 0; k < 3; ++k)
        {
            const PuData &puDataBk = neighbourPuData(h, xNbB[k], yNbB[k]);
            const bool availableBk = puDataBk.isAvailable();

            if (availableBk && availableFlagLXB == 0)
            {
                if (puDataBk.predFlag(X) && DiffPicOrderCnt(h, h[RefPicList(X)][puDataBk.refIdx(X)], h[RefPicList(X)][refIdxLX]) == 0)
                {
                    availableFlagLXB = 1;
                    mvLXB = puDataBk.mv(X);
                    refIdxB = puDataBk.refIdx(X);
                    kB = k;
                }
                else if (puDataBk.predFlag(Y) && DiffPicOrderCnt(h, h[RefPicList(Y)][puDataBk.refIdx(Y)], h[RefPicList(X)][refIdxLX]) == 0)
                {
                    availableFlagLXB = 1;
                    mvLXB = puDataBk.mv(Y);
                    refIdxB = puDataBk.refIdx(Y);
                    kB = k;
                }
            }
        }

        // 4.
        if (!isScaledFlagLX && availableFlagLXB == 1)
        {
            availableFlagLXA = 1;
            mvLXA = mvLXB;
            kA = kB + 2;
        }


        // 5.
        if (!isScaledFlagLX)
        {
            availableFlagLXB = 0;
            for (int k = 0; k < 3 && availableFlagLXB != 1; ++k)
            {
                const PuData &puDataBk = neighbourPuData(h, xNbB[k], yNbB[k]);
                const bool availableBk = puDataBk.isAvailable();

                RefPicList RefPicListB;

                if (availableBk && availableFlagLXB == 0)
                {
                    if (puDataBk.predFlag(X)
                        && LongTermRefPic(h[RefPicList(X)][refIdxLX]) == LongTermRefPic(h[RefPicList(X)][puDataBk.refIdx(X)]))
                    {
                        availableFlagLXB = 1;
                        mvLXB = puDataBk.mv(X);
                        refIdxB = puDataBk.refIdx(X);
                        RefPicListB = RefPicList(X);
                        kB = k;
                    }
                    else if (puDataBk.predFlag(Y)
                        && LongTermRefPic(h[RefPicList(X)][refIdxLX]) == LongTermRefPic(h[RefPicList(Y)][puDataBk.refIdx(Y)]))
                    {
                        availableFlagLXB = 1;
                        mvLXB = puDataBk.mv(Y);
                        refIdxB = puDataBk.refIdx(Y);
                        RefPicListB = RefPicList(Y);
                        kB = k;
                    }
                }

                if (availableFlagLXB == 1 && DiffPicOrderCnt(h, h[RefPicListB][refIdxB], h[RefPicList(X)][refIdxLX]) != 0
                    && h[RefPicListB][refIdxB].reference == SHORT_TERM
                    && h[RefPicList(X)][refIdxLX].reference == SHORT_TERM)
                {
                    int const tD = DiffPicOrderCnt(h, h[PicOrderCntVal()], h[RefPicListB][refIdxB]);
                    int const tB = DiffPicOrderCnt(h, h[PicOrderCntVal()], h[RefPicList(X)][refIdxLX]);
                    mvLXB = distScale(mvLXB, tD, tB);
                }
            }
        }
    }

    MotionVector mvLXCol;
    int availableFlagLXCol;

    // 2.
    if (availableFlagLXA == 1 && availableFlagLXB == 1 && mvLXA != mvLXB)
    {
        availableFlagLXCol = 0;
    }
    else
    {
        if (h[slice_temporal_mvp_enabled_flag()] == 0)
        {
            mvLXCol = MotionVector();
            availableFlagLXCol = 0;
        }
        else
        {
            PuData puDataCol; // temporary
            availableFlagLXCol = deriveTemporalLumaMotionVectorPredictors(h, puDataCol, pu, X, refIdxLX);
            mvLXCol = puDataCol.mv(X);
        }
    }

    MotionVector mvpListLX[3];

    // 3.
    int i = 0;
    if (availableFlagLXA)
    {
        mvpListLX[i++] = mvLXA;
        static const char *names[] = { "A0", "A1", "B0", "B1", "B2" };
        h(MotionCand{ X, names[kA] });
    }
    if (availableFlagLXB)
    {
        mvpListLX[i++] = mvLXB;

        // 4.
        if (i == 2 && mvpListLX[0] == mvpListLX[1])
        {
            i = 1;
        }
        else
        {
            static const char *names[] = { "B0", "B1", "B2" };
            h(MotionCand{ X, names[kB] });
        }
    }
    if (availableFlagLXCol)
    {
        mvpListLX[i++] = mvLXCol;
        h(MotionCand{ X, "Col" });
    }

    // 4.
    int numMvpCandLX = i;

    while (numMvpCandLX < 2)
    {
        mvpListLX[numMvpCandLX++] = MotionVector{ 0, 0 };
        h(MotionCand{ X, "Zero" });
    }

    Mvp::Predictors *predictors = h;
    predictors->mvp[refIdxLX][X][0] = mvpListLX[0];
    predictors->mvp[refIdxLX][X][1] = mvpListLX[1];
}


// populate Mvp::Predictors::mvp[][] for all active pictures in L0 and L1 
template <class H> 
void populateMvp(H &h)
{
    for (int refIdx = 0; refIdx <= h[num_ref_idx_l0_active_minus1()]; ++refIdx)
    {
        predictMvp(h, L0, refIdx);
    }
    for (int refIdx = 0; refIdx <= h[num_ref_idx_l1_active_minus1()]; ++refIdx)
    {
        predictMvp(h, L1, refIdx);
    }
}


template<class H>
MotionVector deriveLumaMotionVectorPrediction(H &h, int refList, const PuData &puData)
{
    predictMvp(h, refList, puData.refIdx(refList));

    Mvp::Predictors *predictors = h;

    prediction_unit const &pu = *static_cast<prediction_unit *>(h);
    auto const xPb = pu.x0;
    auto const yPb = pu.y0;

    // 5.
    auto const mvp_flag = refList ? h[mvp_l1_flag(xPb, yPb)] : h[mvp_l0_flag(xPb, yPb)];
    return predictors->mvp[puData.refIdx(refList)][refList][mvp_flag];
}


// Computes the two MVP predictors for refList
// Sets puData.mv[refList] based on mvd, mvp
template <class H>
void computeMv(H &h, PuData &puData, int refList)
{
    const prediction_unit *pu = h;
    int const xPb = pu->x0;
    int const yPb = pu->y0;

    const MotionVector mvd = h[Mvd(refList, xPb, yPb)];
    const MotionVector mvp = deriveLumaMotionVectorPrediction(h, refList, puData);

    static_assert(sizeof(MotionVector::ComponentType) == 2, "HEVC standard specifies an unsaturated (i.e. wrap-around) 16-bitaddition of MVP and MVD");

    puData.mv(refList) = mvd + mvp;
}


template <class H>
void populateMergeCandidatesInner(H &h, prediction_unit const &puOrig, PuData *const begin, PuData * const end)
{
    if (false)
    {
        // review: can use this or similar techniques to accelerate when actual predictors aren't needed
        Mvp::Predictors *predictors = h;
        PuData puData = {};
        puData.setRefIdx(h, L0, 0);
        for (int i = 0; i < 5; ++i)
            predictors->merge[i] = puData;
        return;
    }

    StateSubstream *stateSubstream = h;
    coding_quadtree const *cqt = h;

    // coding_unit information
    auto const &xCb = cqt->x0;
    auto const &yCb = cqt->y0;
    auto const nCbS = 1 << cqt->log2CbSize;

    // prediction_unit information
    auto pu = puOrig;
    auto &xPb = pu.x0;
    auto &yPb = pu.y0;
    auto &nPbW = pu.nPbW;
    auto &nPbH = pu.nPbH;
    auto partIdx = stateSubstream->partIdx;

    // Derivation process for luma motion vectors for merge mode
    auto const &xOrigP = puOrig.x0;
    auto const &yOrigP = puOrig.y0;;
    auto const &nOrigPbW = puOrig.nPbW;
    auto const &nOrigPbH = puOrig.nPbH;

    if (h[Log2ParMrgLevel()] > 2 && nCbS == 8)
    {
        xPb = xCb;
        yPb = yCb;
        nPbW = nCbS;
        nPbH = nCbS;
        partIdx = 0;
    }

    PuData puDataA1;
    PuData puDataB1;

    PuData *puData = begin;

    // The motion vectors mvL0 and mvL1, the reference indices refIdxL0 and refIdxL1, and the prediction utilization flags predFlagL0 and predFlagL1 are derived by the following ordered steps:
    // 1.
    {
        // Derivation process for spatial merging candidates

        int const xNbA1 = xPb - 1;
        int const yNbA1 = yPb + nPbH - 1;
        puDataA1.reset();
        bool const illegalA1 = partIdx && nPbW < nPbH;
        assert(illegalA1 == ((h[PartMode()] == PART_Nx2N || h[PartMode()] == PART_nLx2N || h[PartMode()] == PART_nRx2N) && partIdx == 1));
        if (!illegalA1)
        {
            PuMergeNeighbour<-1, 0>::get(h, puDataA1, xNbA1, yNbA1, xPb, yPb);
            if (puDataA1.isAvailable())
            {
                *puData = puDataA1;
                h(MotionCand{ -1, "A1" });
                if (++puData == end)
                    return;
            }
        }
    }

    {
        int const xNbB1 = xPb + nPbW - 1;
        int const yNbB1 = yPb - 1;
        puDataB1.reset();
        bool const illegalB1 = partIdx && nPbH < nPbW;
        assert(illegalB1 == ((h[PartMode()] == PART_2NxN || h[PartMode()] == PART_2NxnU || h[PartMode()] == PART_2NxnD) && partIdx == 1));
        if (!illegalB1)
        {
            PuMergeNeighbour<0, -1>::get(h, puDataB1, xNbB1, yNbB1, xPb, yPb);
            if (puDataB1.isAvailable() && puDataA1 != puDataB1)
            {
                *puData = puDataB1;
                h(MotionCand{ -1, "B1" });
                if (++puData == end)
                    return;
            }
        }
    }

    {
        auto const xNbB0 = xPb + nPbW;
        auto const yNbB0 = yPb - 1;
        puData->reset();
        PuMergeNeighbour<0, -1>::get(h, *puData, xNbB0, yNbB0, xPb, yPb);
        if (puData->isAvailable() && *puData != puDataB1)
        {
            h(MotionCand{ -1, "B0" });
            if (++puData == end)
                return;
        }
    }

    {
        auto const xNbA0 = xPb - 1;
        auto const yNbA0 = yPb + nPbH;
        puData->reset();
        PuMergeNeighbour<-1, 0>::get(h, *puData, xNbA0, yNbA0, xPb, yPb);

        if (puData->isAvailable() && *puData != puDataA1)
        {
            h(MotionCand{ -1, "A0" });
            if (++puData == end)
                return;
        }
    }

    if (puData != begin + 4)
    {
        auto const xNbB2 = xPb - 1;
        auto const yNbB2 = yPb - 1;
        puData->reset();
        PuMergeNeighbour<-1, -1>::get(h, *puData, xNbB2, yNbB2, xPb, yPb);
        if (puData->isAvailable() && *puData != puDataA1 && *puData != puDataB1)
        {
            h(MotionCand{ -1, "B2" });
            if (++puData == end)
                return;
        }
    }

    if (h[slice_temporal_mvp_enabled_flag()])
    {
        puData->reset();

        // 3.
        bool availableFlagCol = deriveTemporalLumaMotionVectorPredictors(h, *puData, pu, L0, 0);

        // 4.
        if (h[slice_type()] == B)
            availableFlagCol |= deriveTemporalLumaMotionVectorPredictors(h, *puData, pu, L1, 0);

        if (availableFlagCol)
        {
            h(MotionCand{ -1, "col" });
            if (++puData == end)
                return;
        }
    }

    // 7.
    // Derivation process for combined bi-predictive merging candidates
    if (h[slice_type()] == B)
    {
        auto const numOrigMergeCand = puData - begin;
        auto const combMax = "\x0\x0\x2\x6\xc"[numOrigMergeCand];
        assert(combMax == numOrigMergeCand * (numOrigMergeCand - 1));
        for (int combIdx = 0; combIdx < combMax; ++combIdx)
        {
            // 1.
            auto const l0CandIdx = "\0\1\0\2\1\2\0\3\1\3\2\3"[combIdx];
            auto const l1CandIdx = "\1\0\2\0\2\1\3\0\3\1\3\2"[combIdx];

            // 2.
            auto const &l0Cand = begin[l0CandIdx];
            auto const &l1Cand = begin[l1CandIdx];

            // 3.
            if (l0Cand.predFlag(L0) && l1Cand.predFlag(L1) &&
                ((DiffPicOrderCnt(h, h[RefPicList(L0)][l0Cand.refIdx(L0)], h[RefPicList(L1)][l1Cand.refIdx(L1)]) != 0) || (l0Cand.mv(L0) != l1Cand.mv(L1))))
            {
                *puData = l0Cand;
                puData->mv(L1) = l1Cand.mv(L1);
                puData->setRefIdx(h, L1, l1Cand.refIdx(L1));
                h(MotionCand{ -1, "Comb" });
                if (++puData == end)
                    return;
            }
        }
    }

    // 8.
    // Derivation process for zero motion vector merging candidates
    int numRefIdxMinus1 = h[num_ref_idx_l0_active_minus1()];
    if (h[slice_type()] == B && h[num_ref_idx_l1_active_minus1()] < numRefIdxMinus1)
        numRefIdxMinus1 = h[num_ref_idx_l1_active_minus1()];

    for (int zeroIdx = 0; zeroIdx <= numRefIdxMinus1; ++zeroIdx)
    {
        puData->reset();
        puData->setRefIdx(h, L0, zeroIdx);
        if (h[slice_type()] == B)
            puData->setRefIdx(h, L1, zeroIdx);
        h(MotionCand{ -1, "zero" });
        if (++puData == end)
            return;
    }

    do
    {
        puData->reset();
        puData->setRefIdx(h, L0, 0);
        if (h[slice_type()] == B) 
            puData->setRefIdx(h, L1, 0);
        h(MotionCand{ -1, "zero" });
    } while (++puData != end);
}


template <class H>
void populateMergeCandidates(H &h, prediction_unit const &pu, int mergeIdx=-1)
{
    Mvp::Predictors *predictors = h;

    ++mergeIdx;
    auto const max = mergeIdx ? mergeIdx : h[MaxNumMergeCand()];
    auto const end = predictors->merge + max;

    populateMergeCandidatesInner(h, pu, predictors->merge, end);

    auto const &nOrigPbW = pu.nPbW;
    auto const &nOrigPbH = pu.nPbH;

    // 10.
    if (nOrigPbW + nOrigPbH == 12)
        for (PuData *puData = predictors->merge; puData != end; ++puData)
            if (puData->predFlag(L0) && puData->predFlag(L1))
                puData->setRefIdx(h, L1, -1);
}
 

// Populates MVP predictors if merge_flag=0 or populate merge candidates if merge_flag=1
template <class H>
void populatePuPredictors(H &h, prediction_unit pu, int partIdx, int mergeIdx=-1)
{
    Mvp::Predictors *predictors = h;
    coding_quadtree const *cqt = h;

    PuData puData;

    int xPb = pu.x0;
    int yPb = pu.y0;

    int const xCb = cqt->x0;
    int const yCb = cqt->y0;
    int const nCbS = 1 << cqt->log2CbSize;
    int nPbW = pu.nPbW;
    int nPbH = pu.nPbH;

    int const merge = h[merge_flag(xPb, yPb)];
    if (merge)
    {
        auto const xOrigP = xPb;
        auto const yOrigP = yPb;
        populateMergeCandidates(h, pu, mergeIdx);
    }
    else
    {
        puData.reset();
        if (h[inter_pred_idc(xPb, yPb)] == PRED_L0 || h[inter_pred_idc(xPb, yPb)] == PRED_BI)
            puData.setRefIdx(h, L0, h[ref_idx_l0(xPb, yPb)]);
        if (h[inter_pred_idc(xPb, yPb)] == PRED_L1 || h[inter_pred_idc(xPb, yPb)] == PRED_BI)
            puData.setRefIdx(h, L1, h[ref_idx_l1(xPb, yPb)]);
        if (puData.predFlag(L0)) 
            computeMv(h, puData, L0);
        if (puData.predFlag(L1)) 
            computeMv(h, puData, L1);
    }
}


// Computes puData based on merge candidates and merge_idx
template <class H>
void setPuDataMerge(PuData &puData, H &h)
{
    prediction_unit const *pu = h;
    Mvp::Predictors *predictors = h;
    puData = predictors->merge[h[merge_idx(pu->x0, pu->y0)]];
}


// Computes puData based on MVP predictors and Mvd
template <class H>
void setPuDataMvpPredFlags(PuData &puData, H &h, bool predFlagL0, bool predFlagL1)
{
    Mvp::Predictors *predictors = h;
    prediction_unit const *pu = h;
    auto const &xPb = pu->x0;
    auto const &yPb = pu->y0;

    puData.reset();

    if (predFlagL0)
    {
        auto const refIdx = h[ref_idx_l0(xPb, yPb)];
        puData.setRefIdx(h, L0, refIdx);
        auto const mvpFlag = h[mvp_l0_flag(xPb, yPb)];
        auto const &mvp = predictors->mvp[refIdx][L0][mvpFlag];
        auto const &mvd = h[Mvd(L0, xPb, yPb)];
        puData.mv(L0) = mvp + mvd;
    }

    if (predFlagL1)
    {
        auto const refIdx = h[ref_idx_l1(xPb, yPb)];
        puData.setRefIdx(h, L1, refIdx);
        auto const mvpFlag = h[mvp_l1_flag(xPb, yPb)];
        auto const &mvp = predictors->mvp[refIdx][L1][mvpFlag];
        auto const &mvd = h[Mvd(L1, xPb, yPb)];
        puData.mv(L1) = mvp + mvd;
    }
}


// Computes puData based on MVP predictors and Mvd
template <class H>
void setPuDataMvp(PuData &puData, H &h)
{
    prediction_unit const *pu = h;
    auto const &xPb = pu->x0;
    auto const &yPb = pu->y0;

    setPuDataMvpPredFlags(
        puData, h,
        h[inter_pred_idc(xPb, yPb)] == PRED_L0 || h[inter_pred_idc(xPb, yPb)] == PRED_BI,
        h[inter_pred_idc(xPb, yPb)] == PRED_L1 || h[inter_pred_idc(xPb, yPb)] == PRED_BI);
}


// Computes puData based on merge_flag, merge candidates/merge_idx or MVP predictors/Mvd
template <class H>
void setPuData(H &h, PuData &puData)
{
    prediction_unit const *pu = h;
    auto const &xPb = pu->x0;
    auto const &yPb = pu->y0;
    int const merge = h[merge_flag(xPb, yPb)];

    if (merge)
        setPuDataMerge(puData, h);
    else
        setPuDataMvp(puData, h);
}


// Populates MVP predictors if merge_flag=0 or populate merge candidates if merge_flag=1
// then sets puData based on merge_flag
template <class H>
void processPredictionUnit(H &h, prediction_unit pu, PuData &puData, int partIdx, int mergeIdx=-1)
{
    populatePuPredictors(h, pu, partIdx, mergeIdx);
    setPuData(h, puData);
}

#endif
