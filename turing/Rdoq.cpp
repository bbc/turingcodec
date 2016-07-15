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

#include "Rdoq.h"
#include "Write.h"
#include "Binarization.h"


static inline int32_t estimateBits(int binValue, ContextModel &context)
{
    int state = context.getState();
    static const EntropyEstimate<Cost> bitEstimate;
    return (int32_t)(bitEstimate.table[state ^ binValue].value >> 1);
}



int Rdoq::runQuantisation(short *dst, const short *src, int quantiserScale, int quantiserShift,
                          int n, residual_coding rc,
                          int scanIdx, bool isIntra, bool isSdhEnabled)
{

    int log2Size = rc.log2TrafoSize;
    int blockSize = 1 << log2Size;
    int log2InCgUnits = log2Size - 2;
    int lastScanPos = -1;
    int lastCoeffGroupPos = -1;
    int contextSet = 0;
    int greater1CtxIdx = 1;
    int greater1Cnt = 0;
    int greater2Cnt = 0;
    int golombRiceParameter = 0;
    int cbf = 0;

    m_totalDistCoeff0.set(0.0);
    m_rdCostTu.set(0.0);

    int rateIncUp[32 * 32];
    int rateIncDown[32 * 32];
    int sigRateDelta[32 * 32];
    int deltaU[32 * 32];

    ::memset(rateIncUp, 0, sizeof(int) *  n);
    ::memset(rateIncDown, 0, sizeof(int) *  n);
    ::memset(sigRateDelta, 0, sizeof(int) *  n);
    ::memset(deltaU, 0, sizeof(int) *  n);

    int lastScanPosCgBased = -1;

    const int greater1CtxOff = (rc.cIdx > 0) ? 16 : 0;
    const int greater2CtxOff = (rc.cIdx > 0) ? 4 : 0;

    int     lastScanPosCoeffBased = -1;

    ::memset(m_codedSubBlockFlagArray, 0, sizeof(int) * 64);

    int totalCg = n >> 4;
    int scanPosTuBased;
    cgRdData rdData;

    //-----------------------------------------------------------------------------
    //| Step 1: Perform quantisation and adjust each level to minimise the RD cost
    //-----------------------------------------------------------------------------
    for (int scanPosCgBased = totalCg - 1; scanPosCgBased >= 0; scanPosCgBased--)
    {
        int cgCol = ScanOrder(log2InCgUnits, scanIdx, scanPosCgBased, 0);
        int cgRow = ScanOrder(log2InCgUnits, scanIdx, scanPosCgBased, 1);
        int cgPos = cgRow * (1 << log2InCgUnits) + cgCol;

        // Reset RD data stats for this CG
        rdData.nonZeroCoeffsBeforePos0 = 0;
        rdData.distCoeff0.set(0.0);
        rdData.rateCostSig.set(0.0);
        rdData.rateCostSigPos0.set(0.0);
        rdData.rdCostCoeff.set(0.0);

        const int prevCsbf = getPrevCsbf(m_codedSubBlockFlagArray, cgCol, cgRow, log2Size);
        for (int scanPosCoeffBased = m_cgSize - 1; scanPosCoeffBased >= 0; scanPosCoeffBased--)
        {
            scanPosTuBased = scanPosCgBased*m_cgSize + scanPosCoeffBased;

            int coeffCol = (cgCol << 2) + ScanOrder(2, scanIdx, scanPosCoeffBased, 0);
            int coeffRow = (cgRow << 2) + ScanOrder(2, scanIdx, scanPosCoeffBased, 1);
            int coeffPos = (coeffRow << log2Size) + coeffCol;

            int sourceAbsCoeff = abs(src[coeffPos]);
            int sourceAbsCoeffScaled = abs(src[coeffPos]) * quantiserScale;

            int quantLevel = (sourceAbsCoeffScaled + (1 << (quantiserShift - 1))) >> quantiserShift;

            int32_t zeroCoeffQuantErr = sourceAbsCoeff;
            m_distCoeff0[scanPosTuBased] = zeroCoeffQuantErr * zeroCoeffQuantErr * m_distortionScale;
            m_totalDistCoeff0 += m_distCoeff0[scanPosTuBased];
            dst[coeffPos] = quantLevel;

            if (quantLevel > 0 && lastScanPosCoeffBased < 0)
            {
                // Found one coefficient which is not rounded to zero.
                // It needs to be code so store its scanning position (for last sign. coeff pos)
                // Init the content set coeff_abs_level_greater1_flag and coeff_abs_level_greater2_flag
                lastScanPosCoeffBased = scanPosTuBased;
                contextSet = (scanPosTuBased < 16 || rc.cIdx != 0) ? 0 : 2;
                lastScanPosCgBased = scanPosCgBased;
            }

            if (lastScanPosCoeffBased >= 0)
            {
                // Given that a non zero level is obtained, RDOQ can start adjusting the coefficient level
                int  quantLevelAdjusted;
                int  greater1CtxInc = 4 * contextSet + greater1CtxIdx + greater1CtxOff;
                int  greater2CtxInc = contextSet + greater2CtxOff;

                int   coeffSigCtxInc = getCoeffSigCtxInc(prevCsbf, scanIdx, coeffCol, coeffRow, log2Size, rc.cIdx);
                quantLevelAdjusted = getAdjustedQuantLevel(scanPosTuBased,
                                                           sourceAbsCoeff, quantLevel, coeffSigCtxInc, greater1CtxInc, greater2CtxInc, golombRiceParameter,
                                                           greater1Cnt, greater2Cnt, quantiserShift, scanPosTuBased == lastScanPosCoeffBased);

                deltaU[coeffPos] = (sourceAbsCoeffScaled - (quantLevelAdjusted << quantiserShift)) >> (quantiserShift - 8);

                if (!(scanPosTuBased == lastScanPosCoeffBased))
                {
                    sigRateDelta[coeffPos] = estimateBits(1, m_contexts->get<sig_coeff_flag>(coeffSigCtxInc)) -
                            estimateBits(0, m_contexts->get<sig_coeff_flag>(coeffSigCtxInc));
                }

                if (quantLevelAdjusted > 0)
                {
                    int rateNow = getLevelRate(quantLevelAdjusted, greater1CtxInc, greater2CtxInc, golombRiceParameter, greater1Cnt, greater2Cnt);
                    rateIncUp[coeffPos] = getLevelRate(quantLevelAdjusted + 1, greater1CtxInc, greater2CtxInc, golombRiceParameter, greater1Cnt, greater2Cnt) - rateNow;
                    rateIncDown[coeffPos] = getLevelRate(quantLevelAdjusted - 1, greater1CtxInc, greater2CtxInc, golombRiceParameter, greater1Cnt, greater2Cnt) - rateNow;
                }
                else
                {
                    rateIncUp[coeffPos] = estimateBits(0, m_contexts->get<coeff_abs_level_greater1_flag>(greater1CtxInc));
                }

                dst[coeffPos] = quantLevelAdjusted;
                m_rdCostTu += m_rdCostCoeff[scanPosTuBased];

                updateEntropyCodingEngine(quantLevelAdjusted, scanPosTuBased, rc.cIdx,
                                          greater1Cnt, greater2Cnt, golombRiceParameter,
                                          greater1CtxIdx, contextSet);

            }
            else
            {
                // This coefficient is rounded to zero, add to the TU-based RD cost the associated distortion
                m_rdCostTu += m_distCoeff0[scanPosTuBased];
            }

            // Update the RD and rate costs for this CG:
            // 1) Add the rate cost for coefficient significance flag
            // 2) If we are processing the DC coefficient add its significance cost
            rdData.rateCostSig += m_rateCostCoeffSig[scanPosTuBased];
            if (scanPosCoeffBased == 0)
            {
                rdData.rateCostSigPos0 = m_rateCostCoeffSig[scanPosTuBased];
            }

            // If the level is different from zero then set the significance flag for this CG to 1
            // Update the RD cost for coded coefficients
            // Update the RD cost for coefficient rounded to zero
            // If it's not the last coefficient then increase the number of non-zero coefficients before last position
            if (dst[coeffPos])
            {
                m_codedSubBlockFlagArray[cgPos] = 1;
                rdData.rdCostCoeff += m_rdCostCoeff[scanPosTuBased] - m_rateCostCoeffSig[scanPosTuBased];
                rdData.distCoeff0 += m_distCoeff0[scanPosTuBased];
                if (scanPosCoeffBased != 0)
                {
                    rdData.nonZeroCoeffsBeforePos0++;
                }
            }
        }

        //----------------------------------------------------------------------------
        //| Step 2: Check whether setting all levels to zero provides a better RD cost
        //----------------------------------------------------------------------------
        if (lastScanPosCgBased >= 0)
        {
            if (scanPosCgBased)
            {
                if (m_codedSubBlockFlagArray[cgPos] == 0)
                {
                    // CG with all levels equal to zero. Adjust the TU-based RD cost:
                    // 1) Add the cost of sending the sign. flag for this CG
                    // 2) Subtract the cost associated with the significance of each coefficient
                    // 3) Set the rate cost for this CG equal to the cost of sending the sign. flag
                    int  csbfCtxInc = getCgSigCtxInc(m_codedSubBlockFlagArray, cgCol, cgRow, log2Size, rc.cIdx);
                    Cost costCsbfSigZero = m_lambda * estimateBits(0, m_contexts->get<coded_sub_block_flag>(csbfCtxInc));
                    m_rdCostTu += costCsbfSigZero - rdData.rateCostSig;
                    m_rateCostCgSig[scanPosCgBased] = costCsbfSigZero;
                }
                else
                {
                    // CG with some level different from zero. Check whether the all zero option is better
                    // Avoid last CG in scanning order. Its RD analysis will be dealt in Step 3 below
                    if (scanPosCgBased < lastScanPosCgBased)
                    {
                        if (rdData.nonZeroCoeffsBeforePos0 == 0)
                        {
                            // This CG has all coefficients rounded to zero but the last in reverse scan order.
                            // In this case the significance flag for this coefficient is inferred to be 1.
                            // Therefore adjust the RD cost accordingly:
                            // 1) Subtract the rate cost for this coefficient from the TU-based RD cost
                            // 2) Subtract the rate cost for this coefficient from the rate cost of this CG
                            m_rdCostTu -= rdData.rateCostSigPos0;
                            rdData.rateCostSig -= rdData.rateCostSigPos0;
                        }

                        // Compute the rate costs to send the significance flag equal to 0 or 1 for this CG.
                        // These will be needed later
                        int  csbfCtxInc = getCgSigCtxInc(m_codedSubBlockFlagArray, cgCol, cgRow, log2Size, rc.cIdx);
                        Cost rateCostCsbfSigZero = m_lambda * estimateBits(0, m_contexts->get<coded_sub_block_flag>(csbfCtxInc));
                        Cost rateCostCsbfSigOne = m_lambda * estimateBits(1, m_contexts->get<coded_sub_block_flag>(csbfCtxInc));

                        // Set the cost of all-zero to the TU-based RD cost and then adjust the costs:
                        // 1) Add the rate cost of sign. flag (equal to 1) to the TU-based RD cost
                        // 2) Add the rate cost of sign. flag (equal to 0) to the all-zero RD cost
                        // 3) Set the sign. flag rate cost of this CG to the  rate cost of sign. flag (equal to 1)
                        Cost rdCostAllZero = m_rdCostTu;
                        m_rdCostTu += rateCostCsbfSigOne;
                        rdCostAllZero += rateCostCsbfSigZero;
                        m_rateCostCgSig[scanPosCgBased] = rateCostCsbfSigOne;

                        // Further adjust the RD cost for the all-zero solution
                        // 1) Add the distortion associated with setting all levels to zero for this CG
                        // 2) Subtract the RD cost associated with coding all levels in this CG
                        // 3) Subtract the RD cost associated with the significance for all coeffs in this CG
                        rdCostAllZero += rdData.distCoeff0;
                        rdCostAllZero -= rdData.rdCostCoeff;
                        rdCostAllZero -= rdData.rateCostSig;

                        // Check whether the RD cost for the all-zero solution is lower
                        // If yes, do the following:
                        // 1) Set the significance flag for this CG to zero
                        // 2) Set the TU-based RD cost to the all zero solution
                        //    This cost accounts also for what done in previous CGs
                        // 3) Set the rate cost for sign. flag of this CG to the rate cost of sign. flag (equal to 0)
                        // 4) Set all levels of this CG to zero
                        if (rdCostAllZero < m_rdCostTu)
                        {
                            m_codedSubBlockFlagArray[cgPos] = 0;
                            m_rdCostTu = rdCostAllZero;
                            m_rateCostCgSig[scanPosCgBased] = rateCostCsbfSigZero;

                            for (int iScanPosinCG = m_cgSize - 1; iScanPosinCG >= 0; iScanPosinCG--)
                            {
                                int coeffCol = (cgCol << 2) + ScanOrder(2, scanIdx, iScanPosinCG, 0);
                                int coeffRow = (cgRow << 2) + ScanOrder(2, scanIdx, iScanPosinCG, 1);
                                int coeffPos = (coeffRow << log2Size) + coeffCol;
                                scanPosTuBased = scanPosCgBased * m_cgSize + iScanPosinCG;
                                int uiBlkPos = coeffPos;

                                if (dst[uiBlkPos])
                                {
                                    dst[uiBlkPos] = 0;
                                    m_rdCostCoeff[scanPosTuBased] = m_distCoeff0[scanPosTuBased];
                                    m_rateCostCoeffSig[scanPosTuBased].set(0.0);
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                // CG containing the DC coefficient has always significance equal to 1
                m_codedSubBlockFlagArray[cgPos] = 1;
            }
        }
    }


    if (lastScanPosCoeffBased < 0)
    {
        // All coefficients rounded to zero, don't bother doing RD search for last sign. position
        return 0;
    }

    Cost rdCostTuBest;     // Best RD cost for this TU
    int  cbfCtxInc = 0; // Context idx for coded block flag (cbf) associated with this TU
    int  lastPosIdx = 0; // Index for last significant coefficient position
    rdCostTuBest.set(0);

    //---------------------------------------------------------------------------------------
    //| Step 3: Estimate the last significant coefficient position by minimising the RD Cost
    //---------------------------------------------------------------------------------------
    // Initialise the best RD cost for this TU as if all coefficients where to be rounded to 0
    // Adjust the TU-based RD cost to include the cost of signalling the cbf for this TU
    if (!isIntra && rc.cIdx == 0 && true /*pcCU->getTransformIdx( uiAbsPartIdx ) == 0 //TODO extend for RQT*/)
    {
        cbfCtxInc = 0;
        rdCostTuBest = m_totalDistCoeff0 + m_lambda * estimateBits(0, m_contexts->get<rqt_root_cbf>(cbfCtxInc));
        m_rdCostTu += m_lambda * estimateBits(1, m_contexts->get<rqt_root_cbf>(cbfCtxInc));
    }
    else
    {
        cbfCtxInc = getCbfCtxIdx(rc.cIdx == 0, 0); // TODO: extend for RQT
        int32_t rateZero, rateOne;
        switch (rc.cIdx)
        {
            case 0:
                // Luma
                rateZero = estimateBits(0, m_contexts->get<cbf_luma>(cbfCtxInc));
                rateOne = estimateBits(1, m_contexts->get<cbf_luma>(cbfCtxInc));
                break;
            case 1:
                // Cb
                rateZero = estimateBits(0, m_contexts->get<typename Context<cbf_cb>::Type>(cbfCtxInc));
                rateOne = estimateBits(1, m_contexts->get<typename Context<cbf_cb>::Type>(cbfCtxInc));
                break;
            case 2:
                // Cr
                rateZero = estimateBits(0, m_contexts->get<typename Context<cbf_cr>::Type>(cbfCtxInc));
                rateOne = estimateBits(1, m_contexts->get<typename Context<cbf_cr>::Type>(cbfCtxInc));
                break;
            default:
                assert(0);
        }
        rdCostTuBest = m_totalDistCoeff0 + m_lambda * rateZero;
        m_rdCostTu += m_lambda * rateOne;
    }

    // Search for the last significance coefficient position which minimises the RD cost
    bool foundLast = false;
    for (int scanPosCgBased = lastScanPosCgBased; scanPosCgBased >= 0; scanPosCgBased--)
    {
        int cgCol = ScanOrder(log2InCgUnits, scanIdx, scanPosCgBased, 0);
        int cgRow = ScanOrder(log2InCgUnits, scanIdx, scanPosCgBased, 1);
        int cgPos = cgRow * (1 << log2InCgUnits) + cgCol;

        // Adjust TU-based RD cost by subtracting the rate cost for the significance flag for this CG
        // In fact
        // a) If this CG contains the last sign. coefficient then its significance is inferred to be 1
        // b) Otherwise if this CG is not significant, then then last sign. coeff will be above the current scanning position
        m_rdCostTu -= m_rateCostCgSig[scanPosCgBased];
        if (m_codedSubBlockFlagArray[cgPos])
        {
            for (int scanPosCoeffBased = m_cgSize - 1; scanPosCoeffBased >= 0; scanPosCoeffBased--)
            {
                scanPosTuBased = scanPosCgBased*m_cgSize + scanPosCoeffBased;
                if (scanPosTuBased > lastScanPosCoeffBased) continue;

                int coeffCol = (cgCol << 2) + ScanOrder(2, scanIdx, scanPosCoeffBased, 0);
                int coeffRow = (cgRow << 2) + ScanOrder(2, scanIdx, scanPosCoeffBased, 1);
                int coeffPos = (coeffRow << log2Size) + coeffCol;

                if (dst[coeffPos])
                {

                    // Compute the RD cost if this coefficient were the last significant one
                    // 1) Compute the rate cost associated with coding this position
                    // 2) Total RD cost will be TU-based cost plus position cost and minus the significance flag cost since its flag is inferred to be 1
                    Cost rateCostLastSig = scanIdx == 2 ? getLastSigCoeffPosRateCost(coeffRow, coeffCol, rc.cIdx, log2Size) :
                            getLastSigCoeffPosRateCost(coeffCol, coeffRow, rc.cIdx, log2Size);
                    Cost totalCost = m_rdCostTu + rateCostLastSig - m_rateCostCoeffSig[scanPosTuBased];

                    // Check whether setting the last significance coeff. pos here minimises the RD cost
                    // If yes, do the following:
                    // 1) Set last sign. coefficient position to the current one
                    // 2) Set the best RD for this TU to the total RD cost
                    if (totalCost < rdCostTuBest)
                    {
                        lastPosIdx = scanPosTuBased + 1;
                        rdCostTuBest = totalCost;
                    }

                    // Check whether this coefficient has a high level value. If yes, there's no point
                    // checking other positions. In fact, the more we move towards low frequencies the
                    // higher the level values and the less likely the RD cost will improve
                    if (dst[coeffPos] > 1)
                    {
                        foundLast = true;
                        break;
                    }
                    m_rdCostTu -= m_rdCostCoeff[scanPosTuBased];
                    m_rdCostTu += m_distCoeff0[scanPosTuBased];
                }
                else
                {
                    m_rdCostTu -= m_rateCostCoeffSig[scanPosTuBased];
                }
            }
            if (foundLast)
            {
                break;
            }
        }
    }

    int *scanArray = new int[n];

    for (int cgScan = 0, idx = 0; cgScan < totalCg; cgScan++)
    {
        int cgCol = ScanOrder(log2InCgUnits, scanIdx, cgScan, 0);
        int cgRow = ScanOrder(log2InCgUnits, scanIdx, cgScan, 1);
        for (int coeffScan = 0; coeffScan < m_cgSize; coeffScan++)
        {
            int coeffCol = (cgCol << 2) + ScanOrder(2, scanIdx, coeffScan, 0);
            int coeffRow = (cgRow << 2) + ScanOrder(2, scanIdx, coeffScan, 1);
            int coeffPos = (coeffRow << log2Size) + coeffCol;
            scanArray[idx++] = coeffPos;
        }
    }

    // Recover the sign for coded levels
    int absSumLevel = 0;
    for (int scanPos = 0; scanPos < lastPosIdx; scanPos++)
    {
        int blkPos = scanArray[scanPos];
        int level = dst[blkPos];
        absSumLevel += level;
        dst[blkPos] = (src[blkPos] < 0) ? -level : level;
        cbf |= level;
    }

    // Set to zero uncoded levels
    for (int scanPos = lastPosIdx; scanPos <= lastScanPosCoeffBased; scanPos++)
    {
        dst[scanArray[scanPos]] = 0;
    }

    //	Sign data hiding
    if (isSdhEnabled && absSumLevel >= 2)
    {
        signDataHiding(totalCg, dst, src, scanArray, rateIncUp, rateIncDown, sigRateDelta, deltaU);
    }

    delete[] scanArray;

    return cbf;
}

int Rdoq::getAdjustedQuantLevel(int        scanPosTuBased,
                                int        srcCoeff,
                                int        quantLevel,
                                int        coeffSigCtxInc,
                                int        greater1CtxInc,
                                int        greater2CtxInc,
                                int        golombRiceParameter,
                                int        greater1Cnt,
                                int        greater2Cnt,
                                int        quantizerShift,
                                bool       lastScannedCoeff)
{
    Cost currSigCoeffCost; currSigCoeffCost.set(0);
    int  optimalQuantLevel = 0;

    if (!lastScannedCoeff && quantLevel < 3)
    {
        // Compute RD costs in case the function has to return
        m_rateCostCoeffSig[scanPosTuBased] = m_lambda * estimateBits(0, m_contexts->get<sig_coeff_flag>(coeffSigCtxInc));
        m_rdCostCoeff[scanPosTuBased] = m_distCoeff0[scanPosTuBased] + m_rateCostCoeffSig[scanPosTuBased];
        if (quantLevel == 0)
        {
            return optimalQuantLevel;
        }
    }
    else
    {
        m_rdCostCoeff[scanPosTuBased].set(numeric_limits<double>::max());
    }

    // Significance for last coefficient is always inferred to be 1.
    // Associated cost will be zero in this case
    if (!lastScannedCoeff)
    {
        currSigCoeffCost = m_lambda * estimateBits(1, m_contexts->get<sig_coeff_flag>(coeffSigCtxInc));
    }

    // Perform RD optimisation for level adjustment
    int minQuantLevel = (quantLevel > 1 ? quantLevel - 1 : 1);
    for (int currQuantLevel = quantLevel; currQuantLevel >= minQuantLevel; currQuantLevel--)
    {
        //int64_t quantError = srcCoeffScaled  - ( currQuantLevel << quantizerShift );
        int reconLevel = getInvQuantisedLevel(currQuantLevel);
        int32_t quantError = srcCoeff - reconLevel;
        Cost currRdCost = quantError * quantError * m_distortionScale +
                getLevelRateCost(currQuantLevel, greater1CtxInc, greater2CtxInc,
                                 golombRiceParameter, greater1Cnt, greater2Cnt);
        currRdCost += currSigCoeffCost;

        if (currRdCost < m_rdCostCoeff[scanPosTuBased])
        {
            optimalQuantLevel = currQuantLevel;
            m_rdCostCoeff[scanPosTuBased] = currRdCost;
            m_rateCostCoeffSig[scanPosTuBased] = currSigCoeffCost;
        }
    }

    return optimalQuantLevel;
}

int Rdoq::getCoeffSigCtxInc(int   prevCsbf,
                            int   scanIdx,
                            int   xC,
                            int   yC,
                            int   log2TrafoSize,
                            int   cIdx)
{
    int coeffSigCtxInc;
    if (log2TrafoSize == 2)
    {
        const int ctxIdxMap[16] =
        {
                0, 1, 4, 5,
                2, 3, 4, 5,
                6, 6, 8, 8,
                7, 7, 8, 8
        };
        coeffSigCtxInc = ctxIdxMap[(yC << 2) + xC];
    }
    else if (xC + yC == 0)
    {
        coeffSigCtxInc = 0;
    }
    else
    {
        const int xS = xC >> 2;
        const int yS = yC >> 2;

        const int xP = xC & 3;
        const int yP = yC & 3;

        if (prevCsbf == 0)
        {
            coeffSigCtxInc = (xP + yP == 0) ? 2 : (xP + yP < 3) ? 1 : 0;
        }
        else if (prevCsbf == 1)
        {
            coeffSigCtxInc = (yP == 0) ? 2 : (yP == 1) ? 1 : 0;
        }
        else if (prevCsbf == 2)
        {
            coeffSigCtxInc = (xP == 0) ? 2 : (xP == 1) ? 1 : 0;
        }
        else
        {
            assert(prevCsbf == 3);
            coeffSigCtxInc = 2;
        }

        if (cIdx == 0)
        {
            if (xS + yS > 0)
            {
                coeffSigCtxInc += 3;
            }

            if (log2TrafoSize == 3)
            {
                coeffSigCtxInc += (scanIdx == 0) ? 9 : 15;
            }
            else
            {
                coeffSigCtxInc += 21;
            }
        }
        else
        {
            if (log2TrafoSize == 3)
            {
                coeffSigCtxInc += 9;
            }
            else
            {
                coeffSigCtxInc += 12;
            }
        }
    }

    assert(coeffSigCtxInc >= 0);
    assert(coeffSigCtxInc < 54);

    if (cIdx == 0)
    {
        return coeffSigCtxInc;
    }
    else
    {
        return 27 + coeffSigCtxInc;
    }
}


int  Rdoq::getPrevCsbf(const int* codedsubBlockFlagArray, int xS, int yS, int log2TrafoSize)
{
    int prevCsbf = 0;
    int widthInCgUnits = 1 << (log2TrafoSize - 2);

    if (xS < (1 << (log2TrafoSize - 2)) - 1)
    {
        prevCsbf += codedsubBlockFlagArray[yS * widthInCgUnits + xS + 1];
    }
    if (yS < (1 << (log2TrafoSize - 2)) - 1)
    {
        prevCsbf += (codedsubBlockFlagArray[(yS + 1) * widthInCgUnits + xS] << 1);
    }

    return prevCsbf;
}

Cost  Rdoq::getLevelRateCost(int  quantLevel,
                             int  greater1CtxInc,
                             int  greater2CtxInc,
                             int  golombRiceParameter,
                             int  greater1Cnt,
                             int  greater2Cnt)
{
    int32_t currentRate = 32768;
    int baseLevel = (greater1Cnt < MAX_REG_BINS_GREATER1) ? (2 + (greater2Cnt < MAX_REG_BINS_GREATER2)) : 1;

    if (quantLevel >= baseLevel)
    {
        int symbol = quantLevel - baseLevel;
        int length;
        if (symbol < (3 << golombRiceParameter))
        {
            length = symbol >> golombRiceParameter;
            currentRate += (length + 1 + golombRiceParameter) << 15;
        }
        else
        {
            length = golombRiceParameter;
            symbol = symbol - (3 << golombRiceParameter);
            while (symbol >= (1 << length))
            {
                symbol -= (1 << (length++));
            }
            currentRate += (3 + length + 1 - golombRiceParameter + length) << 15;
        }
        if (greater1Cnt < MAX_REG_BINS_GREATER1)
        {
            currentRate += estimateBits(1, m_contexts->get<coeff_abs_level_greater1_flag>(greater1CtxInc));

            if (greater2Cnt < MAX_REG_BINS_GREATER2)
            {
                currentRate += estimateBits(1, m_contexts->get<coeff_abs_level_greater2_flag>(greater2CtxInc));
            }
        }
    }
    else
        if (quantLevel == 1)
        {
            currentRate += estimateBits(0, m_contexts->get<coeff_abs_level_greater1_flag>(greater1CtxInc));
        }
        else if (quantLevel == 2)
        {
            currentRate += estimateBits(1, m_contexts->get<coeff_abs_level_greater1_flag>(greater1CtxInc));
            currentRate += estimateBits(0, m_contexts->get<coeff_abs_level_greater2_flag>(greater2CtxInc));
        }
        else
        {
            assert(0);
        }
    return m_lambda * currentRate;
}

int Rdoq::getCgSigCtxInc(const int  *codedsubBlockFlagArray,
                         const int   xS,
                         const int   yS,
                         int         log2TrafoSize,
                         int         cIdx)
{
    const int widthInCgUnits = 1 << (log2TrafoSize - 2);
    int csbfCtx = 0;
    if (xS < (1 << (log2TrafoSize - 2)) - 1)
    {
        csbfCtx += codedsubBlockFlagArray[yS * widthInCgUnits + xS + 1];
    }
    if (yS < (1 << (log2TrafoSize - 2)) - 1)
    {
        csbfCtx += codedsubBlockFlagArray[(yS + 1) * widthInCgUnits + xS];
    }

    if (cIdx == 0)
        return min(csbfCtx, 1);
    else
        return 2 + min(csbfCtx, 1);
}

int Rdoq::getCbfCtxIdx(bool isLuma, int rqtDepth)
{
    if (!isLuma)
    {
        return rqtDepth;
    }
    else
    {
        const int uiCtx = (rqtDepth == 0 ? 1 : 0);
        return uiCtx;
    }
}
Cost Rdoq::getLastSigCoeffPosRateCost(const int    xC,
                                      const int    yC,
                                      const int    cIdx,
                                      const int    log2TrasfoSize)
{
    const int32_t binarisationLengthForPosition[32] = { 0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9 };

    // Get the length of the binary string associated with each coordinate (xC, yC)
    int32_t stringLengthX = binarisationLengthForPosition[xC];
    int32_t stringLengthY = binarisationLengthForPosition[yC];

    int32_t rateX = 0, rateY = 0;
    int ctxX, ctxY;

    // Compute the rate of each associated binary string
    for (int i = 0; i < stringLengthX; i++)
    {
        ctxX = getLastSignCoeffPosPrefixCtxInc(i, cIdx, log2TrasfoSize);
        rateX += estimateBits(1, m_contexts->get<last_sig_coeff_x_prefix>(ctxX));
    }
    if (stringLengthX < 9)
    {
        ctxX = getLastSignCoeffPosPrefixCtxInc(stringLengthX, cIdx, log2TrasfoSize);
        rateX += estimateBits(0, m_contexts->get<last_sig_coeff_x_prefix>(ctxX));
    }

    for (int i = 0; i < stringLengthY; i++)
    {
        ctxY = getLastSignCoeffPosPrefixCtxInc(i, cIdx, log2TrasfoSize);
        rateY += estimateBits(1, m_contexts->get<last_sig_coeff_y_prefix>(ctxY));
    }
    if (stringLengthY < 9)
    {
        ctxY = getLastSignCoeffPosPrefixCtxInc(stringLengthY, cIdx, log2TrasfoSize);
        rateY += estimateBits(0, m_contexts->get<last_sig_coeff_y_prefix>(ctxY));
    }

    // Add the additional bits (if any) required by coding the string suffix
    int32_t rateTotal = rateX + rateY;
    if (stringLengthX > 3)
    {
        rateTotal += 32768 * ((stringLengthX - 2) >> 1);
    }
    if (stringLengthY > 3)
    {
        rateTotal += 32768 * ((stringLengthY - 2) >> 1);
    }
    return m_lambda * rateTotal;
}

int Rdoq::getLastSignCoeffPosPrefixCtxInc(int binIdx, int cIdx, int log2TrafoSize)
{
    const int ctxOffset = cIdx ? 15 : (3 * (log2TrafoSize - 2) + ((log2TrafoSize - 1) >> 2));
    const int ctxShift = cIdx ? (log2TrafoSize - 2) : ((log2TrafoSize + 1) >> 2);
    const int ctxInc = (binIdx >> ctxShift) + ctxOffset;
    return Clip3(0, 17, ctxInc);
}

void Rdoq::updateEntropyCodingEngine(int quantisedLevel, int scanPosTuBased, int cIdx,
                                     int &greater1Cnt, int &greater2Cnt,
                                     int &golombRiceParameter, int &greater1CtxIdx,
                                     int &contextSet)
{
    // Update Golomb-Rice parameter and associated counters
    const int baseLevel = (greater1Cnt < MAX_REG_BINS_GREATER1) ? (2 + (greater2Cnt < MAX_REG_BINS_GREATER2)) : 1;
    if (quantisedLevel >= baseLevel)
    {
        if (quantisedLevel > 3 * (1 << golombRiceParameter))
        {
            golombRiceParameter = min<int>(golombRiceParameter + 1, 4);
        }
    }
    if (quantisedLevel >= 1)
    {
        greater1Cnt++;
    }
    if (quantisedLevel > 1)
    {
        // Honours the condition: lastGreater1Flag is equal to 1 in clause 9.3.4.2.6
        greater1CtxIdx = 0;
        greater2Cnt++;
    }
    else if ((greater1CtxIdx < 3) && (greater1CtxIdx > 0) && quantisedLevel)
    {
        // Honours the condition: lastGreater1Flag is equal to 0 in clause 9.3.4.2.6
        greater1CtxIdx++;
    }

    // New coefficient group starts. Update as follows:
    // Update context set
    // Reset Golomb-Rice parameter
    // Reset greater1 and greater2 counters
    // Reset greater1CtxInc to 1
    if ((scanPosTuBased % 16 == 0) && (scanPosTuBased > 0))
    {
        golombRiceParameter = 0;
        greater1Cnt = 0;
        greater2Cnt = 0;
        contextSet = (scanPosTuBased == 16 || cIdx != 0) ? 0 : 2;
        if (greater1CtxIdx == 0)
        {
            contextSet++;
        }
        greater1CtxIdx = 1;
    }
}

int  Rdoq::getLevelRate(int             quantLevel,
                        int             greater1CtxInc,
                        int             greater2CtxInc,
                        int             golombRiceParameter,
                        int             greater1Cnt,
                        int             greater2Cnt)
{
    int rate = 0;
    int baseLevel = (greater1Cnt < MAX_REG_BINS_GREATER1) ? (2 + (greater2Cnt < MAX_REG_BINS_GREATER2)) : 1;

    const int golombRiceRange[5] =
    {
            7, 14, 26, 46, 78
    };

    const int golombRicePrefixLen[5] =
    {
            8, 7, 6, 5, 4
    };

    if (quantLevel >= baseLevel)
    {
        int  symbol = quantLevel - baseLevel;
        int  maxVlc = golombRiceRange[golombRiceParameter];
        bool expGolomb = (symbol > maxVlc);

        if (expGolomb)
        {
            quantLevel = symbol - maxVlc;
            int iEGS = 1;  for (int uiMax = 2; quantLevel >= uiMax; uiMax <<= 1, iEGS += 2);
            rate += iEGS << 15;
            symbol = min<int>(symbol, (maxVlc + 1));
        }

        int prefixLen = symbol >> (golombRiceParameter + 1);
        int ui16NumBins = min<int>(prefixLen, golombRicePrefixLen[golombRiceParameter]) + golombRiceParameter;

        rate += ui16NumBins << 15;

        if (greater1Cnt < MAX_REG_BINS_GREATER1)
        {
            rate += estimateBits(1, m_contexts->get<coeff_abs_level_greater1_flag>(greater1CtxInc));

            if (greater2Cnt < MAX_REG_BINS_GREATER2)
            {
                rate += estimateBits(1, m_contexts->get<coeff_abs_level_greater2_flag>(greater2CtxInc));
            }
        }
    }
    else
        if (quantLevel == 0)
        {
            return 0;
        }
        else if (quantLevel == 1)
        {
            rate += estimateBits(0, m_contexts->get<coeff_abs_level_greater1_flag>(greater1CtxInc));
        }
        else if (quantLevel == 2)
        {
            rate += estimateBits(1, m_contexts->get<coeff_abs_level_greater1_flag>(greater1CtxInc));
            rate += estimateBits(0, m_contexts->get<coeff_abs_level_greater2_flag>(greater2CtxInc));
        }
        else
        {
            assert(0);
        }
    return rate;
}

void Rdoq::signDataHiding(
        int totalCg,
        short *dst,
        const short *src,
        int *scanArray,
        int *rateIncUp,
        int *rateIncDown,
        int *sigRateDelta,
        int *deltaU
)
{
    int lastCG = -1;
    int currAbsSum = 0;
    int coeffIdx;

    for (int subSet = totalCg - 1; subSet >= 0; subSet--)
    {
        int  subPos = subSet << 4;
        int  firstNZPosInCG = m_cgSize, lastNZPosInCG = -1;
        currAbsSum = 0;

        for (coeffIdx = m_cgSize - 1; coeffIdx >= 0; --coeffIdx)
        {
            if (dst[scanArray[coeffIdx + subPos]])
            {
                lastNZPosInCG = coeffIdx;
                break;
            }
        }

        for (coeffIdx = 0; coeffIdx < m_cgSize; coeffIdx++)
        {
            if (dst[scanArray[coeffIdx + subPos]])
            {
                firstNZPosInCG = coeffIdx;
                break;
            }
        }

        for (coeffIdx = firstNZPosInCG; coeffIdx <= lastNZPosInCG; coeffIdx++)
        {
            currAbsSum += dst[scanArray[coeffIdx + subPos]];
        }

        if (lastNZPosInCG >= 0 && lastCG == -1)
        {
            lastCG = 1;
        }

        if (lastNZPosInCG - firstNZPosInCG >= 4)
        {
            int signbit = (dst[scanArray[subPos + firstNZPosInCG]] > 0 ? 0 : 1);
            if (signbit != (currAbsSum & 0x1))  // sign will be hidden but levels need to be changed accordingly
            {
                // calculate the cost associated with the change of the level
                int minCostInc = std::numeric_limits<int>::max(), currCost = std::numeric_limits<int>::max();
                int minPos = -1, finalChange = 0, currChange = 0;

                for (coeffIdx = (lastCG == 1 ? lastNZPosInCG : m_cgSize - 1); coeffIdx >= 0; --coeffIdx)
                {
                    int uiBlkPos = scanArray[coeffIdx + subPos];
                    if (dst[uiBlkPos] != 0)
                    {
                        int costForLevelUp = m_shdRdFactor * (-deltaU[uiBlkPos]) + rateIncUp[uiBlkPos];
                        int costForLevelDown = m_shdRdFactor * (deltaU[uiBlkPos]) + rateIncDown[uiBlkPos]
                                                                                                - (abs(dst[uiBlkPos]) == 1 ? ((1 << 15) + sigRateDelta[uiBlkPos]) : 0);

                        if (lastCG == 1 && lastNZPosInCG == coeffIdx && abs(dst[uiBlkPos]) == 1)
                        {
                            costForLevelDown -= (4 << 15);
                        }

                        if (costForLevelUp < costForLevelDown)
                        {
                            currCost = costForLevelUp;
                            currChange = 1;
                        }
                        else
                        {
                            currChange = -1;
                            if (coeffIdx == firstNZPosInCG && abs(dst[uiBlkPos]) == 1)
                            {
                                currCost = std::numeric_limits<int>::max();
                            }
                            else
                            {
                                currCost = costForLevelDown;
                            }
                        }
                    }
                    else
                    {
                        currCost = m_shdRdFactor * (-(abs(deltaU[uiBlkPos]))) + (1 << 15) + rateIncUp[uiBlkPos] + sigRateDelta[uiBlkPos];
                        currChange = 1;

                        if (coeffIdx < firstNZPosInCG)
                        {
                            int thissignbit = (src[uiBlkPos] >= 0 ? 0 : 1);
                            if (thissignbit != signbit)
                            {
                                currCost = std::numeric_limits<int>::max();
                            }
                        }
                    }

                    if (currCost < minCostInc)
                    {
                        minCostInc = currCost;
                        finalChange = currChange;
                        minPos = uiBlkPos;
                    }
                }

                if (dst[minPos] == 32767 || dst[minPos] == -32768)
                {
                    finalChange = -1;
                }

                if (src[minPos] >= 0)
                {
                    dst[minPos] += finalChange;
                }
                else
                {
                    dst[minPos] -= finalChange;
                }
            }
        }

        if (lastCG == 1)
        {
            lastCG = 0;
        }
    }
}
