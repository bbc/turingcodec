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

#ifndef INCLUDED_Rdoq_h
#define INCLUDED_Rdoq_h
#include "ScanOrder.h"
#include "Global.h"
#include "Cabac.h"
#include "Syntax.h"
#include "Cost.h"

using namespace std;

//TODO: Revise all of these macros
#define DIST_SCALE_PRECISION       15
#define MAX_TR_DYNAMIC_RANGE       15
#define MAX_REG_BINS_GREATER1       8
#define MAX_REG_BINS_GREATER2       1

typedef struct
{
    int    nonZeroCoeffsBeforePos0; // Number of non zero coefficients before last position in
    // reverse scan order
    Cost rdCostCoeff;             // RD cost for all coded coefficients in a CG
    Cost distCoeff0;              // Distortion when all coefficients are rounded to zero in a CG
    Cost rateCostSig;             // Lambda-rate cost for sign. flag of all coeffs. in this CG
    Cost rateCostSigPos0;         // Lambda-rate cost for sign. flag of last coefficient in reverse scan order
} cgRdData;

class Rdoq
{
private:

    int m_bitDepth;
    const int m_cgSize   = 16;
    int m_invQuantiserScale;
    int m_invQuantiserShift;
    int m_invQuantiserOffset;
    int m_codedSubBlockFlagArray[64];

    Lambda m_lambda;
    int    m_shdRdFactor;

    FixedPoint<int32_t, 16> m_distortionScale;

    Cost m_totalDistCoeff0;          // Distortion when all coefficients are rounded to zero
    Cost m_rdCostCoeff[32 * 32] = { { 0 } };   // RD cost for each coefficient which is coded
    Cost m_rateCostCoeffSig[32 * 32] = { { 0 } }; // Lambda-rate cost associated with the sig flag
    // of each coefficient (coded or rounded to zero)
    Cost m_distCoeff0[32 * 32] = { { 0 } };       // Distortion associated with each coefficient
    // rounded to zero
    Cost m_rdCostTu;                    // RD cost associated with current TU
    Cost m_rateCostCgSig[64] = { { 0 } };       // Lambda-rate cost associated with the sig flag for
    // each CG

    Contexts *m_contexts;

    int getAdjustedQuantLevel(int            scanPosTuBased,
                              int            srcCoeff,
                              int            quantLevel,
                              int            coeffSigCtxInc,
                              int            greater1CtxInc,
                              int            greater2CtxInc,
                              int            golombRiceParameter,
                              int            greater1Cnt,
                              int            greater2Cnt,
                              int            quantizerShift,
                              bool           lastScannedCoeff);

    void updateEntropyCodingEngine(int             quantisedLevel,
                                   int             scanPosTuBased,
                                   int             cIdx,
                                   int            &greater1Cnt,
                                   int            &greater2Cnt,
                                   int            &golombRiceParameter,
                                   int            &greater1CtxIdx,
                                   int            &contextSet);

    Cost getLevelRateCost(int             quantLevel,
                          int             greater1CtxInc,
                          int             greater2CtxInc,
                          int             golombRiceParameter,
                          int             greater1Cnt,
                          int             greater2Cnt);

    int  getLevelRate(int             quantLevel,
                      int             greater1CtxInc,
                      int             greater2CtxInc,
                      int             golombRiceParameter,
                      int             greater1Cnt,
                      int             greater2Cnt);

    Cost getLastSigCoeffPosRateCost(const int       xC,
                                    const int       yC,
                                    const int       cIdx,
                                    const int       log2TrasfoSize);

    int getCgSigCtxInc(const int      *codedsubBlockFlagArray,
                       const int       xS,
                       const int       yS,
                       int             log2TrafoSize,
                       int             cIdx);

    int getCoeffSigCtxInc(int             previousCsbf,
                          int             scanIdx,
                          int             posX,
                          int             posY,
                          int             log2BlockSize,
                          int             textureType);

    int  getPrevCsbf(const int*      codedsubBlockFlagArray,
                     int             xS,
                     int             yS,
                     int             log2TrafoSize);

    int getCbfCtxIdx(bool            isLuma,
                     int             rqtDepth);

    int getLastSignCoeffPosPrefixCtxInc(int             binIdx,
                                        int             cIdx,
                                        int             log2TrafoSize);
    int getInvQuantisedLevel(int             quantLevel)
    {
        int clippedLevel = Clip3(-32768, 32767, quantLevel);
        int reconLevel = (clippedLevel * m_invQuantiserScale + m_invQuantiserOffset) >> m_invQuantiserShift;
        return Clip3(-32768, 32767, reconLevel);
    }

    void signDataHiding(
            int   totalCg,
            short *dst,
            const short *src,
            int *scanArray,
            int *rateIncUp,
            int *rateIncDown,
            int *sigRateDelta,
            int *deltaU
    );

public:
    Rdoq() : 	m_invQuantiserScale(0),
    m_invQuantiserShift(0),
    m_invQuantiserOffset(0),
    m_contexts(0)
    {
        m_totalDistCoeff0.set(0.0);
        m_rdCostTu.set(0.0);
        m_lambda.set(1.0);
        m_distortionScale.set(1.0);
        m_shdRdFactor = 1;
        m_bitDepth = 8;
    }
    Rdoq(double lambda, Contexts *contexts, int quantiserScale, int invQuantScale, int log2Size, int bitDepth)
    {
        m_lambda.set(lambda);
        m_shdRdFactor = (int)(invQuantScale * invQuantScale / lambda / 16 + 0.5);
        m_contexts = contexts;

        m_bitDepth = bitDepth;

        m_totalDistCoeff0.set(0.0);
        m_rdCostTu.set(0.0);

        //TODO: Consider pre-computing this value
        int transformShift       = MAX_TR_DYNAMIC_RANGE - m_bitDepth - log2Size;
        int distortionScaleShift = DIST_SCALE_PRECISION - 2 * transformShift - 2 * (m_bitDepth - 8);
        m_distortionScale.set(1 << distortionScaleShift);
        m_invQuantiserScale      = invQuantScale;
        m_invQuantiserShift      = 20 - 14 - transformShift;
        m_invQuantiserOffset     = 1 << (m_invQuantiserShift - 1);
    }
    int runQuantisation(short           *dst,
                        const short     *src,
                        int              quantiserScale,
                        int              quantiserShift,
                        int              n,
                        residual_coding  rc,
                        int              scanIdx,
                        bool             isIntra,
                        bool             isSdhEnabled);
};

#endif /* TURING_QUANTIZE_RDOQ_H_ */
