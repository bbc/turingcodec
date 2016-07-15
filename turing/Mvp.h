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


namespace Mvp {

    struct Predictors
    {
        PuData merge[5];
        MotionVector mvp[15 /* refIdx */][2 /* refList */][2 /* mvp_flag */];
    };

}


static MotionVector distScale(MotionVector mv, int tD, int tB)
{
    int const td = Clip3(-128, 127, tD);
    int const tb = Clip3(-128, 127, tB);
    int const tx = (16384 + (Abs(td)>>1)) / td;
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
    auto h =  h2.extend(&*motion);

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
        if (puDataCol.predFlag(L0) == 0)
        {
            listCol = L1;
        }
        else if (puDataCol.predFlag(L1) == 0)
        {
            listCol = L0;
        }
        else
        {
            assert(puDataCol.predFlag(L0) == 1 && puDataCol.predFlag(L1) == 1);
            if (static_cast<StatePicture *>(h)->allBackwards)
            {
                listCol = refList;
            }
            else
            {
                listCol = h[collocated_from_l0_flag()] ? L1 : L0;
            }
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
            {
                mvLXCol = mvCol;
            }
            else
            {
                mvLXCol = distScale(mvCol, colPocDiff, currPocDiff);
            }
        }
    }
    if (!availableFlagLXCol)
    {
        puData.setRefIdx(h, refList, -1);
    }
    return availableFlagLXCol;
}


template <class H>
const StatePicture *getColPic(H &h, int *poc=0)
{
    if (h[collocated_ref_idx()] >= 0 || h[collocated_ref_idx()] < 16)
    {
        auto &pic = h[RefPicList(h[collocated_from_l0_flag()] ? 0 : 1)][h[collocated_ref_idx()]];
        if (pic.dp)
        {
            if (poc) *poc = (*pic.dp)[PicOrderCntVal()];
            return pic.dp.get();
        }
    }
    return 0;
}


// returns availableFlagLXCol
template <class H>
bool deriveTemporalLumaMotionVectorPredictors(H &h, PuData &puDataCol, MotionVector &mvLXCol, int refList, int refIdxLX, int xPb, int yPb, int nPbW, int nPbH)
{
    if (h[slice_temporal_mvp_enabled_flag()] == 0)
    {
        mvLXCol = MotionVector{ 0, 0 };
        return false;
    }

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
        mvLXCol = puDataCol.mv(refList);
    }

    // 3.
    if (!availableFlagLXCol)
    {
        // The central collocated motion vector is derived as follows
        int const xColCtr = xPb + (nPbW >> 1);
        int const yColCtr = yPb + (nPbH >> 1);

        puDataCol.setRefIdx(h, refList, refIdxLX);
        availableFlagLXCol = deriveCollocatedMotionVectors(h, puDataCol, colPic->motion.get(), refList, xColCtr, yColCtr);
        mvLXCol = puDataCol.mv(refList);
    }

    return availableFlagLXCol;
}


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
            availableFlagLXCol = deriveTemporalLumaMotionVectorPredictors(h, puDataCol, mvLXCol, X, refIdxLX, xPb, yPb, pu.nPbW, pu.nPbH);
        }
    }

    MotionVector mvpListLX[3];

    // 3.
    int i = 0;
    if (availableFlagLXA)
    {
        mvpListLX[i++] = mvLXA;
        static const char *names[] = { "A0", "A1", "B0", "B1", "B2" };
        h(MotionCand(X, names[kA]));
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
            h(MotionCand(X, names[kB]));
        }
    }
    if (availableFlagLXCol)
    {
        mvpListLX[i++] = mvLXCol;
        h(MotionCand(X, "Col"));
    }

    // 4.
    int numMvpCandLX = i;

    while (numMvpCandLX < 2)
    {
        mvpListLX[numMvpCandLX++] = MotionVector{ 0, 0 };
        h(MotionCand(X, "Zero"));
    }

    Mvp::Predictors *predictors = h;
    auto const refIdx = X ? h[ref_idx_l1(xPb, yPb)] : h[ref_idx_l0(xPb, yPb)];
    predictors->mvp[refIdx][X][0] = mvpListLX[0];
    predictors->mvp[refIdx][X][1] = mvpListLX[1];
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
void populateMergeCandidates(H &h, prediction_unit pu, int partIdx)
{
    coding_quadtree const *cqt = h;

    int xPb = pu.x0;
    int yPb = pu.y0;
    int const poc = h[PicOrderCntVal()];

    const coding_unit cu(cqt->x0, cqt->y0, cqt->log2CbSize);
    int const xCb = cu.x0;
    int const yCb = cu.y0;
    int const nCbS = 1 << cu.log2CbSize;
    int nPbW = pu.nPbW;
    int nPbH = pu.nPbH;

    // Derivation process for luma motion vectors for merge mode
    int const xOrigP = xPb;
    int const yOrigP = yPb;
    int const nOrigPbW = pu.nPbW;
    int const nOrigPbH = pu.nPbH;

    if (h[Log2ParMrgLevel()] > 2 && nCbS == 8)
    {
        xPb = xCb;
        yPb = yCb;
        nPbW = nCbS;
        nPbH = nCbS;
        partIdx = 0;
    }

    static int const A0 = 0;
    static int const A1 = 1;
    static int const B0 = 2;
    static int const B1 = 3;
    static int const B2 = 4;
    static int const Col = 5;
    static int const combCand = 6;
    static int const zeroCand = 20;

    bool available[6];

    // Review: work out actual needed array size
    bool availableFlag[60] = { false };
    PuData puDataN[60];

    // The motion vectors mvL0 and mvL1, the reference indices refIdxL0 and refIdxL1, and the prediction utilization flags predFlagL0 and predFlagL1 are derived by the following ordered steps:
    // 1.
    {
        // Derivation process for spatial merging candidates

        // A1
        int const xNbA1 = xPb - 1;
        int const yNbA1 = yPb + nPbH - 1;
        {
            puDataN[A1].reset();

            const bool illegalA1 = (h[PartMode()] == PART_Nx2N || h[PartMode()] == PART_nLx2N || h[PartMode()] == PART_nRx2N) && partIdx == 1;

            if (!illegalA1) PuMergeNeighbour<-1, 0>::get(h, puDataN[A1], xNbA1, yNbA1, xPb, yPb);

            available[A1] = puDataN[A1].isAvailable();
            availableFlag[A1] = available[A1] ? 1 : 0;
        }

        // B1
        int const xNbB1 = xPb + nPbW - 1;
        int const yNbB1 = yPb - 1;
        {
            puDataN[B1].reset();
            available[B1] = false;

            const bool illegalB1 = (h[PartMode()] == PART_2NxN || h[PartMode()] == PART_2NxnU || h[PartMode()] == PART_2NxnD) && partIdx == 1;

            if (!illegalB1)
            {
                PuMergeNeighbour<0, -1>::get(h, puDataN[B1], xNbB1, yNbB1, xPb, yPb);
                available[B1] = puDataN[B1].isAvailable();
                availableFlag[B1] = available[B1] ? 1 : 0;
            }

            if (puDataN[A1] == puDataN[B1]) availableFlag[B1] = 0;
        }

        // B0
        int const xNbB0 = xPb + nPbW;
        int const yNbB0 = yPb - 1;
        {
            puDataN[B0].reset();
            available[B0] = false;

            PuMergeNeighbour<0, -1>::get(h, puDataN[B0], xNbB0, yNbB0, xPb, yPb);
            available[B0] = puDataN[B0].isAvailable();
            availableFlag[B0] = available[B0] ? 1 : 0;

            if (puDataN[B0] == puDataN[B1]) availableFlag[B0] = 0;
        }

        // A0
        int const xNbA0 = xPb - 1;
        int const yNbA0 = yPb + nPbH;
        {
            puDataN[A0].reset();
            available[A0] = false;

            PuMergeNeighbour<-1, 0>::get(h, puDataN[A0], xNbA0, yNbA0, xPb, yPb);
            available[A0] = puDataN[A0].isAvailable();
            availableFlag[A0] = available[A0] ? 1 : 0;

            if (puDataN[A0] == puDataN[A1]) availableFlag[A0] = 0;
        }

        // B2
        int const xNbB2 = xPb - 1;
        int const yNbB2 = yPb - 1;
        {
            puDataN[B2].reset();
            available[B2] = false;
            availableFlag[B2] = 0;

            if (availableFlag[A0] + availableFlag[A1] + availableFlag[B0] + availableFlag[B1] != 4)
            {
                PuMergeNeighbour<-1, -1>::get(h, puDataN[B2], xNbB2, yNbB2, xPb, yPb);
                available[B2] = puDataN[B2].isAvailable();
                availableFlag[B2] = available[B2] ? 1 : 0;

                if (puDataN[B2] == puDataN[A1]) availableFlag[B2] = 0;
                if (puDataN[B2] == puDataN[B1]) availableFlag[B2] = 0;
            }
        }
    }

    // 2.
    puDataN[Col].reset();

    // 3.
    {
        MotionVector dummy;
        int dummy2 = 0;
        const bool availableFlagL0Col = deriveTemporalLumaMotionVectorPredictors(h, puDataN[Col], dummy, L0, dummy2, xPb, yPb, nPbW, nPbH);
        availableFlag[Col] = availableFlagL0Col;
    }

    // 4.
    if (h[slice_type()] == B)
    {
        MotionVector dummy;
        int dummy2 = 0;
        const bool availableFlagL1Col = deriveTemporalLumaMotionVectorPredictors(h, puDataN[Col], dummy, L1, dummy2, xPb, yPb, nPbW, nPbH);
        availableFlag[Col] = availableFlag[Col] || availableFlagL1Col;
    }


    // 5.
    int mergeCandList[60];

    int i = 0;
    if (availableFlag[A1])
    {
        mergeCandList[i++] = A1;
        h(MotionCand(-1, "A1"));
    }
    if (availableFlag[B1])
    {
        mergeCandList[i++] = B1;
        h(MotionCand(-1, "B1"));
    }
    if (availableFlag[B0])
    {
        mergeCandList[i++] = B0;
        h(MotionCand(-1, "B0"));
    }
    if (availableFlag[A0])
    {
        mergeCandList[i++] = A0;
        h(MotionCand(-1, "A0"));
    }
    if (availableFlag[B2])
    {
        mergeCandList[i++] = B2;
        h(MotionCand(-1, "B2"));
    }
    if (availableFlag[Col])
    {
        mergeCandList[i++] = Col;
        h(MotionCand(-1, "col"));
    }

    // 6.
    int const numOrigMergeCand = i;
    int numCurrMergeCand = i;

    // 7.
    if (h[slice_type()] == B)
    {
        // Derivation process for combined bi-predictive merging candidates

        if (numOrigMergeCand > 1 && numOrigMergeCand < h[MaxNumMergeCand()])
        {
            int numInputMergeCand = numCurrMergeCand;
            int combIdx = 0;
            bool combStop = false;
            do
            {
                // 1.
                int l0CandIdx = "\0\1\0\2\1\2\0\3\1\3\2\3"[combIdx];
                int l1CandIdx = "\1\0\2\0\2\1\3\0\3\1\3\2"[combIdx];

                // 2.
                int l0Cand = mergeCandList[l0CandIdx];
                int l1Cand = mergeCandList[l1CandIdx];

                // 3.
                if (puDataN[l0Cand].predFlag(L0) && puDataN[l1Cand].predFlag(L1) &&
                        ((DiffPicOrderCnt(h, h[RefPicList(L0)][puDataN[l0Cand].refIdx(L0)], h[RefPicList(L1)][puDataN[l1Cand].refIdx(L1)]) != 0) || (puDataN[l0Cand].mv(L0) != puDataN[l1Cand].mv(L1))))
                {
                    int k = numCurrMergeCand - numInputMergeCand;

                    mergeCandList[numCurrMergeCand] = combCand + k;
                    h(MotionCand(-1, "Comb"));

                    puDataN[combCand + k] = puDataN[l0Cand];
                    puDataN[combCand + k].mv(L1) = puDataN[l1Cand].mv(L1);
                    puDataN[combCand + k].setRefIdx(h, L1, puDataN[l1Cand].refIdx(L1));

                    numCurrMergeCand = numCurrMergeCand + 1;
                }

                // 4.
                ++combIdx;

                // 5.
                if (combIdx == (numOrigMergeCand * (numOrigMergeCand - 1)) || numCurrMergeCand == h[MaxNumMergeCand()])
                {
                    combStop = true;
                }
            } while (!combStop);
        }
    }

    // 8.
    {
        // Derivation process for zero motion vector merging candidates
        int numRefIdx;
        if (h[slice_type()] == P)
        {
            numRefIdx = h[num_ref_idx_l0_active_minus1()] + 1;
        }
        else
        {
            assert(h[slice_type()] == B);
            numRefIdx = std::min(h[num_ref_idx_l0_active_minus1()] + 1, h[num_ref_idx_l1_active_minus1()] + 1);
        }

        if (numCurrMergeCand < h[MaxNumMergeCand()])
        {
            int const numInputMergeCand = numCurrMergeCand;

            int zeroIdx = 0;
            do
            {
                // 1.
                if (h[slice_type()] == P)
                {
                    int m = numCurrMergeCand - numInputMergeCand;

                    mergeCandList[numCurrMergeCand] = zeroCand + m;

                    puDataN[zeroCand + m].reset();
                    puDataN[zeroCand + m].setRefIdx(h, L0, (zeroIdx < numRefIdx) ? zeroIdx : 0);

                    numCurrMergeCand = numCurrMergeCand + 1;
                    h(MotionCand(-1, "zero"));
                }
                else
                {
                    int m = numCurrMergeCand - numInputMergeCand;

                    mergeCandList[numCurrMergeCand] = zeroCand + m;

                    puDataN[zeroCand + m].reset();
                    puDataN[zeroCand + m].setRefIdx(h, L0, (zeroIdx < numRefIdx) ? zeroIdx : 0);
                    puDataN[zeroCand + m].setRefIdx(h, L1, (zeroIdx < numRefIdx) ? zeroIdx : 0);

                    numCurrMergeCand = numCurrMergeCand + 1;
                    h(MotionCand(-1, "zero"));
                }

                // 2.
                ++zeroIdx;
            } while (!(numCurrMergeCand == h[MaxNumMergeCand()]));
        }
    }

    for (int i = 0; i < h[MaxNumMergeCand()]; ++i)
    {
        const auto candType = mergeCandList[i];

        // 10. (appears before step 10 so that all values of Mvp::Predictors::merge[] are set correctly)
        if (puDataN[candType].predFlag(L0) && puDataN[candType].predFlag(L1) && (nOrigPbW + nOrigPbH) == 12)
        {
            puDataN[candType].setRefIdx(h, L1, -1);
        }

        Mvp::Predictors *predictors = h;
        predictors->merge[i] = puDataN[candType];
    }
}


// Populates MVP predictors if merge_flag=0 or populate merge candidates if merge_flag=1
template <class H>
void populatePuPredictors(H &h, prediction_unit pu, int partIdx)
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
        populateMergeCandidates(h, pu, partIdx);
        puData = predictors->merge[h[merge_idx(xOrigP, yOrigP)]];
    }
    else
    {
        puData.reset();
        if (h[inter_pred_idc(xPb, yPb)] == PRED_L0 || h[inter_pred_idc(xPb, yPb)] == PRED_BI)
        {
            puData.setRefIdx(h, L0, h[ref_idx_l0(xPb, yPb)]);
        }
        if (h[inter_pred_idc(xPb, yPb)] == PRED_L1 || h[inter_pred_idc(xPb, yPb)] == PRED_BI)
        {
            puData.setRefIdx(h, L1, h[ref_idx_l1(xPb, yPb)]);
        }
        if (puData.predFlag(L0)) computeMv(h, puData, L0);
        if (puData.predFlag(L1)) computeMv(h, puData, L1);
    }

    assert(puData.isAvailable());
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
    {
        setPuDataMerge(puData, h);
    }
    else
    {
        setPuDataMvp(puData, h);
    }
}


// Populates MVP predictors if merge_flag=0 or populate merge candidates if merge_flag=1
// then sets puData based on merge_flag
template <class H>
void processPredictionUnit(H &h, prediction_unit pu, PuData &puData, int partIdx)
{
    populatePuPredictors(h, pu, partIdx);
    setPuData(h, puData);
}

#endif
