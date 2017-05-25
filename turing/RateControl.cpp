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
/*
 * RateControl.cpp
 *
 *  Created on: 4 May 2016
 *      Author: Matteo Naccari
 */

#include "RateControl.h"
#include "EstimateIntraComplexity.h"
                   
PictureController::PictureController(std::shared_ptr<DataStorage> dataStorage, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctbSize, int targetBits, double averageBpp, std::shared_ptr<InputQueue::Docket> docket) :
    m_targetBits(targetBits),
    m_pixelsPerPicture(pictureHeight*pictureWidth),
    m_averageBpp(averageBpp),
    m_dataStorageAccess(dataStorage),
    m_pictureHeightInCtbs(pictureHeightInCtbs),
    m_pictureWidthInCtbs(pictureWidthInCtbs),
    m_pictureSizeInCtbs(pictureHeightInCtbs*pictureWidthInCtbs),
    m_pictureLambdaIni(NON_VALID_LAMBDA),
    m_pictureLambdaFin(NON_VALID_LAMBDA),
    m_pictureQp(NON_VALID_QP),
    m_pictureAlpha((docket->sliceType == I) ? ALPHA_INTRA : INITIAL_ALPHA),
    m_pictureBeta((docket->sliceType == I) ? BETA_INTRA : INITIAL_BETA),
    m_pictureBpp(0.0),
    m_isIntra(docket->sliceType == I),
    m_sopLevel(docket->sliceType == I ? 0 : docket->sopLevel),
    m_totalCostIntra(0),
    m_paramsUpdated(false),
    m_hierarchyLevel(docket->hierarchyLevel),
    m_numLeftSameHierarchyLevel(docket->numSameHierarchyLevel),
    m_sopSize(docket->currentGopSize),
    m_poc(docket->poc),
    m_codingBits(0),
    m_intraPoc(docket->intraFramePoc),
    m_segmentPoc(docket->segmentPoc)
{
    m_ctbControllerEngine = new CtbController[m_pictureSizeInCtbs];

    for(int r = 0, ctbIdx = 0; r < m_pictureHeightInCtbs; r++)
    {
        const int pixelsInCurrentRow = r == (pictureHeightInCtbs - 1) ? pictureHeight - ctbSize*r : ctbSize;

        for(int c = 0; c < pictureWidthInCtbs; c++, ctbIdx++)
        {
            const int pixelsInCurrentCol = c == (pictureWidthInCtbs - 1) ? pictureWidth - ctbSize*c : ctbSize;
            const int ctbPixels = pixelsInCurrentRow * pixelsInCurrentCol;
            m_ctbControllerEngine[ctbIdx].setNumberOfPixels(ctbPixels);
        }
    }

    {
        if (averageBpp < 0.005)
        {
            m_alphaUpdateStep = 0.005;
            m_betaUpdateStep = 0.0025;
        }
        else if (averageBpp < 0.01)
        {
            m_alphaUpdateStep = 0.01;
            m_betaUpdateStep = 0.005;
        }
        else if (averageBpp < 0.05)
        {
            m_alphaUpdateStep = 0.05;
            m_betaUpdateStep = 0.025;
        }
        else if (averageBpp < 0.1)
        {
            m_alphaUpdateStep = 0.1;
            m_betaUpdateStep = 0.05;
        }
        else if (averageBpp < 0.2)
        {
            m_alphaUpdateStep = 0.15;
            m_betaUpdateStep = 0.075;
        }
        else if (averageBpp < 0.5)
        {
            m_alphaUpdateStep = 0.2;
            m_betaUpdateStep = 0.1;
        }
        else
        {
            m_alphaUpdateStep = 0.4;
            m_betaUpdateStep = 0.2;
        }
    }

    if (m_isIntra)
    {
        m_totalCostIntra = docket->icInfo->getSatdSum();
        // Initialise the array of CTB controllers for intra
        for (int r = 0, ctbIdx = 0; r < m_pictureHeightInCtbs; r++)
        {
            for (int c = 0; c < m_pictureWidthInCtbs; c++, ctbIdx++)
            {
                m_ctbControllerEngine[ctbIdx].setCtbQp(NON_VALID_QP);
                int costIntraCtb = docket->icInfo->getSatdCtu(r, c);
                assert(costIntraCtb != -1);
                m_ctbControllerEngine[ctbIdx].setCostIntra((double)costIntraCtb);
                const int estimatedBits = (int)((double)m_targetBits * (double)costIntraCtb / (double)m_totalCostIntra);
                m_ctbControllerEngine[ctbIdx].setTargetBitsEst(estimatedBits);
            }
        }
    }
}

double PictureController::estimateLambda(double previousLambdaIni)
{
    CodedPicture &currentPicture = m_dataStorageAccess->getPictureAtLevel(m_sopLevel);
    double bpp = (double)m_targetBits / (double)m_pixelsPerPicture;
    double alpha = currentPicture.getAlpha();
    double beta  = currentPicture.getBeta();
    double lambdaForClipping = currentPicture.getLambda();
    //	Compute lambda according to the model for intra or inter
    double lambda;

    if (m_isIntra)
    {
        
        auto costIntra = static_cast<int>(m_totalCostIntra);
        double madPixelBased = pow((double)costIntra / (double)m_pixelsPerPicture, BETA_INTRA_MAD);
        lambda = (alpha / 256.0) * pow(madPixelBased / (double)bpp, beta);
    }
    else
    {
        lambda = alpha * pow(bpp, beta);
    }

#if !DISABLE_CLIP_PL
    if (lambdaForClipping > 0.0 )
    {
        lambdaForClipping = Clip3(0.1, 10000.0, lambdaForClipping);
        lambda = Clip3(lambdaForClipping * pow(2.0, -2.0 / 3.0), lambdaForClipping * pow(2.0, 2.0 / 3.0), lambda);
    }

    if (previousLambdaIni > 0.0)
    {
        previousLambdaIni = Clip3(0.1, 2000.0, previousLambdaIni);
        lambda = Clip3(previousLambdaIni * pow(2.0, -10.0 / 3.0), previousLambdaIni * pow(2.0, 10.0 / 3.0), lambda);
    }
    else
    {
        lambda = Clip3( 0.1, 10000.0, lambda );
    }
#else
    lambda = Clip3(0.1, 10000.0, lambda);
#endif

    m_pictureLambdaIni = lambda;
    CodedCtb *codedCtbAtLevel = m_dataStorageAccess->getCodedCtbAtLevel(m_sopLevel);
    for (int ctuIdx = 0; ctuIdx < m_pictureSizeInCtbs; ctuIdx++)
    {
        const double alpha = codedCtbAtLevel[ctuIdx].getAlpha();
        const double beta  = codedCtbAtLevel[ctuIdx].getBeta();
        m_ctbControllerEngine[ctuIdx].setCtbAlpha(alpha);
        m_ctbControllerEngine[ctuIdx].setCtbBeta(beta);
    }
    
    if (!m_isIntra)
    {
        double totalWeight = 0.0;
        std::vector<double> ctbWeights(m_pictureSizeInCtbs, 0.0);

        // Compute the weight of each CTU based on the picture level lambda
        for (int ctbIdx = 0; ctbIdx < m_pictureSizeInCtbs; ctbIdx++)
        {
            double ctbWeight = pow(lambda / alpha, 1.0 / beta) * m_ctbControllerEngine[ctbIdx].getNumberOfPixels();
            if (ctbWeight < 0.01)
            {
                ctbWeight = 0.01;
            }
            ctbWeights[ctbIdx] = ctbWeight;
            totalWeight += ctbWeight;
        }

        // Compute the estimated target bits for each CTU based on its weight
        for (int ctbIdx = 0; ctbIdx < m_pictureSizeInCtbs; ctbIdx++)
        {
            const int ctbBitsEst = static_cast<int>(ctbWeights[ctbIdx] / totalWeight * static_cast<double>(m_targetBits));
            m_ctbControllerEngine[ctbIdx].setTargetBitsEst(ctbBitsEst);
        }
    }
    return lambda;
}

int PictureController::estimateQp(int sopLevel, int previousQp)
{
    int qp = int(4.2005 * log(m_pictureLambdaIni) + 13.7122 + 0.5);

    CodedPicture &pictureAtCurrentLevel = m_dataStorageAccess->getPictureAtLevel(sopLevel);

    int qpForClipping = pictureAtCurrentLevel.getQp();

    if (qpForClipping > NON_VALID_QP)
    {
        qp = Clip3(qpForClipping - 3, qpForClipping + 3, qp);
    }

    if (previousQp > NON_VALID_QP)
    {
        qp = Clip3(previousQp - 5, previousQp + 5, qp);
    }

    qp = Clip3(2, 46, qp);
    m_pictureQp = qp;
    return qp;
}

void PictureController::getCodedInfoFromWavefront(bool wpp, int ctbAddrInRs, int &lastValidQp, double& lastValidLambda, int128bit &cumulativeTargetBits, int128bit &cumulativeBitsSpent, int64_t &cumulativeCtbsCoded, bool currFrame, int callerPoc)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    int rowStart = ctbAddrInRs / m_pictureWidthInCtbs;
    int colStart = (ctbAddrInRs % m_pictureWidthInCtbs);

    for (int r = rowStart; r >= 0; r--)
    {
        for (int c = colStart; c >= 0; c--)
        {
            if(!(currFrame && c==colStart && r == rowStart))
            {
                int currentIdx = r*m_pictureWidthInCtbs + c;
                assert(0 <= currentIdx && currentIdx < m_pictureSizeInCtbs);
                assert(m_ctbControllerEngine[currentIdx].getFinishedFlag());
                if(lastValidQp == NON_VALID_QP)
                    lastValidQp = m_ctbControllerEngine[currentIdx].getCtbQp();
                if (lastValidLambda == NON_VALID_LAMBDA)
                    lastValidLambda = m_ctbControllerEngine[currentIdx].getCtbLambda();
                cumulativeBitsSpent += static_cast<int64_t>(m_ctbControllerEngine[currentIdx].getCodedBits());
                cumulativeTargetBits += static_cast<int64_t>(m_ctbControllerEngine[currentIdx].getTargetBitsEst());
                cumulativeCtbsCoded++;
            }
        }
        colStart = wpp ? (min<int>(m_pictureWidthInCtbs - 1, colStart +  1)) : (m_pictureWidthInCtbs - 1);
    }
}

int PictureController::estimateCtbLambdaAndQp(int ctbAddrInRs)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    if(m_isIntra)
    {
        return estimateCtbLambdaAndQpIntra(ctbAddrInRs);
    }
    else
    {
        return estimateCtbLambdaAndQpInter(ctbAddrInRs);
    }
}

int PictureController::estimateCtbLambdaAndQpInter(int ctbAddrInRs)
{
    // Step 1: Estimate the lambda from lambda = alpha * (bpp)^beta;
    const double alpha = m_ctbControllerEngine[ctbAddrInRs].getCtbAlpha();
    const double beta  = m_ctbControllerEngine[ctbAddrInRs].getCtbBeta();

    const double bpp   = m_ctbControllerEngine[ctbAddrInRs].getCtbBpp();

    double estLambdaCtb = alpha * pow(bpp, beta);

    const int pictureQp = m_pictureQp;

#if !DISABLE_CLIP_CTL
    if (m_pictureLambdaIni > 0.0 )
    {
       estLambdaCtb = Clip3(m_pictureLambdaIni * pow( 2.0, -4.0/2.0 ), m_pictureLambdaIni * pow( 2.0, 4.0/2.0 ), estLambdaCtb );
    }
    else
    {
        estLambdaCtb = Clip3(0.1, 10000.0, estLambdaCtb);
    }


#endif
    if (estLambdaCtb < 0.1)
    {
        estLambdaCtb = 0.1;
    }

    m_ctbControllerEngine[ctbAddrInRs].setCtbLambda(estLambdaCtb);
    m_ctbControllerEngine[ctbAddrInRs].setCtbReciprocalLambda(1.0 / estLambdaCtb);
    m_ctbControllerEngine[ctbAddrInRs].setCtbReciprocalSqrtLambda(1.0 / sqrt(estLambdaCtb));

    // Step 2: Estimate the QP from lambda using the nonlinear relationship
    int ctbEstQp = (int)( 4.2005 * log( estLambdaCtb ) + 13.7122 + 0.5 );

    const int lastValidQp = m_ctbControllerEngine[ctbAddrInRs].getLastValidQp();

    int minQP = m_pictureQp - 3;
    int maxQP = m_pictureQp + 3;
    if (lastValidQp != NON_VALID_QP)
    {
        minQP = max(lastValidQp - 2, minQP);
        maxQP = min(lastValidQp + 2, maxQP);
    }
    ctbEstQp = Clip3(minQP, maxQP, ctbEstQp );
    ctbEstQp = Clip3(2, 46, ctbEstQp);

    m_ctbControllerEngine[ctbAddrInRs].setCtbQp(ctbEstQp);
    return ctbEstQp;
}

int PictureController::estimateCtbLambdaAndQpIntra(int ctbAddrInRs)
{
    // Step 1: Estimate the lambda from the model lambda = (alpha/256) * (cost/bpp)^beta
    const double alpha = m_ctbControllerEngine[ctbAddrInRs].getCtbAlpha();
    const double beta = m_ctbControllerEngine[ctbAddrInRs].getCtbBeta();

    double costPixelBased = m_ctbControllerEngine[ctbAddrInRs].getCostIntra() /
            static_cast<double>(m_ctbControllerEngine[ctbAddrInRs].getNumberOfPixels());
    const double bpp = m_ctbControllerEngine[ctbAddrInRs].getCtbBpp();
    costPixelBased = pow(costPixelBased, BETA_INTRA_MAD);
    double lambda = (alpha/256.0) * pow( costPixelBased/bpp, beta );
    const int lastValidQp = m_ctbControllerEngine[ctbAddrInRs].getLastValidQp();

    int minQP = m_pictureQp - 4;
    int maxQP = m_pictureQp + 4;
    if (lastValidQp != NON_VALID_QP)
    {
        minQP = max(lastValidQp - 2, minQP);
        maxQP = min(lastValidQp + 2, maxQP);
    }
#if !DISABLE_CLIP_CTL


    double maxLambda = exp(((double)(maxQP+0.49)-13.7122)/4.2005);
    double minLambda = exp(((double)(minQP-0.49)-13.7122)/4.2005);

    lambda = Clip3(minLambda, maxLambda, lambda);
#endif
    lambda = Clip3(0.1, 10000.0, lambda);

    // Step 2: Estimate the QP using the nonlinear relationship
    int qp = int( 4.2005 * log(lambda) + 13.7122 + 0.5 );
    qp = Clip3(minQP, maxQP, qp);
    qp = Clip3(2, 46, qp);

    m_ctbControllerEngine[ctbAddrInRs].setCtbQp(qp);
    m_ctbControllerEngine[ctbAddrInRs].setCtbLambda(lambda);
    m_ctbControllerEngine[ctbAddrInRs].setCtbReciprocalLambda(1.0 / lambda);
    m_ctbControllerEngine[ctbAddrInRs].setCtbReciprocalSqrtLambda(1.0 / sqrt(lambda));

    return qp;
}

void PictureController::storeCtbParameters(int codingBits, bool isIntra, int ctbAddrInRs, int pictureLevel, std::shared_ptr<CpbInfo> cpbInfo)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    m_ctbControllerEngine[ctbAddrInRs].setCodedBits(codingBits);
    m_ctbControllerEngine[ctbAddrInRs].setFinishedFlag(true);
    bool isValidCtb = m_ctbControllerEngine[ctbAddrInRs].getValidityFlag(); // i.e. is not a skipped CTB
    m_ctbControllerEngine[ctbAddrInRs].setIsIntra(isIntra);
}

void PictureController::updateModelParameters()
{
    CodedPicture &codedPictureAtLevel = m_dataStorageAccess->getPictureAtLevel(m_sopLevel);
    double currentAlpha = codedPictureAtLevel.getAlpha();
    double currentBeta = codedPictureAtLevel.getBeta();
    double targetBpp = static_cast<double>(m_targetBits) / static_cast<double>(m_pixelsPerPicture);
    double actualBpp = m_pictureBpp;
    
    codedPictureAtLevel.updateCounter();
    codedPictureAtLevel.updateBits(targetBpp, actualBpp);
    double averageTargetBpp = codedPictureAtLevel.getAverageTargetBpp();
    double averageActualBpp = codedPictureAtLevel.getAverageActualBpp();

    double initialLambda = currentAlpha * pow(averageTargetBpp, currentBeta);
    double finalLambda = currentAlpha * pow(averageActualBpp, currentBeta);

    if (m_isIntra)
    {
        double lnbpp = log(pow(static_cast<double>(m_totalCostIntra) / static_cast<double>(m_pixelsPerPicture), BETA_INTRA_MAD));
        double diffLambda = currentBeta*(log(static_cast<double>(m_codingBits)) - log(static_cast<double>(m_targetBits)));

        diffLambda = Clip3(-0.125, 0.125, 0.25*diffLambda);
        currentAlpha = currentAlpha * exp(diffLambda);
        currentBeta = currentBeta + diffLambda / lnbpp;
    }
    else
    {
        if (initialLambda < 0.01 || finalLambda < 0.01 || averageActualBpp < 0.0001)
        {
            currentAlpha *= (1.0 - m_alphaUpdateStep / 2.0);
            currentBeta *= (1.0 - m_betaUpdateStep / 2.0);
        }
        else
        {
            currentAlpha += m_alphaUpdateStep * (log(initialLambda) - log(finalLambda)) * currentAlpha;
            double lnbpp = log(averageActualBpp);
            lnbpp = Clip3(-5.0, -0.1, lnbpp);
            currentBeta += m_betaUpdateStep * (log(initialLambda) - log(finalLambda)) * lnbpp;

        }
        currentAlpha = Clip3(ALPHA_MIN, ALPHA_MAX, currentAlpha);
        currentBeta = Clip3(BETA_MIN, BETA_MAX, currentBeta);
    }

    codedPictureAtLevel.setAlpha(currentAlpha);
    codedPictureAtLevel.setBeta(currentBeta);
    codedPictureAtLevel.setLambda(m_pictureLambdaIni);
    codedPictureAtLevel.setQp(m_pictureQp);
    const int headerBits = m_codingBits - getTotalCodingBitsCtbs();
    codedPictureAtLevel.setHeaderBits(headerBits);
    
    if (!m_isIntra)
    {
        for (int ctbIdx = 0; ctbIdx < m_pictureSizeInCtbs; ctbIdx++)
        {
            CodedCtb *codedCtbAtLevel = m_dataStorageAccess->getCodedCtbAtLevel(m_sopLevel);
            double alpha = codedCtbAtLevel[ctbIdx].getAlpha();
            double beta  = codedCtbAtLevel[ctbIdx].getBeta();

            int ctbNumberOfPixels = m_ctbControllerEngine[ctbIdx].getNumberOfPixels();
            const int codingBits  = m_ctbControllerEngine[ctbIdx].getCodedBits();

            double ctbTargetBpp = m_ctbControllerEngine[ctbIdx].getCtbBpp();
            double ctbLambda = m_ctbControllerEngine[ctbIdx].getCtbLambda();//alpha * pow(ctbTargetBpp, beta);

            double bpp = static_cast<double>(codingBits) / static_cast<double>(ctbNumberOfPixels);
            double computedLambda = alpha * pow(bpp, beta);

            if (ctbLambda < 0.01 || computedLambda < 0.01 || bpp < 0.0001)
            {
                alpha *= (1.0 - m_alphaUpdateStep / 2.0);
                beta *= (1.0 - m_betaUpdateStep / 2.0);
            }
            else
            {
                computedLambda = Clip3(ctbLambda / 10.0, ctbLambda * 10.0, computedLambda);
                alpha += m_alphaUpdateStep * (log(ctbLambda) - log(computedLambda)) * alpha;
                double lnbpp = log(bpp);
                lnbpp = Clip3(-5.0, -0.1, lnbpp);
                beta += m_betaUpdateStep * (log(ctbLambda) - log(computedLambda)) * lnbpp;
            }
            alpha = Clip3(ALPHA_MIN, ALPHA_MAX, alpha);
            beta = Clip3(BETA_MIN, BETA_MAX, beta);

            codedCtbAtLevel[ctbIdx].setAlpha(alpha);
            codedCtbAtLevel[ctbIdx].setBeta(beta);
        }
    }
    
}

void PictureController::computeCtbTargetBits(int ctbAddrInRs, int lastValidQp, double lastValidLambda, int128bit cumulativeTargetBits, int128bit cumulativeBitsSpent, int128bit cumulativeDeltaBits, int64_t ctbsLeftInPicture, int64_t cumulativeCtbsCoded, int bitdepth, std::shared_ptr<CpbInfo> cpbInfo)
{
    int avgBits = 0;
    int influence = CTB_SMOOTH_WINDOW;
    int ctbEstBits = m_ctbControllerEngine[ctbAddrInRs].getTargetBitsEst();
    cumulativeCtbsCoded = std::max<int64_t>(1, cumulativeCtbsCoded);
    int128bit tempDiff = (cumulativeTargetBits - cumulativeBitsSpent);
    int maxBitsPerCtb = (bitdepth * m_ctbControllerEngine[ctbAddrInRs].getNumberOfPixels()* 3) >>1;
    int delta = static_cast<int>(tempDiff.to64());

    avgBits = Clip3(1, maxBitsPerCtb, ctbEstBits + delta);

    m_ctbControllerEngine[ctbAddrInRs].setTargetBits(avgBits);
    const double bpp = static_cast<double>(avgBits) / static_cast<double>(m_ctbControllerEngine[ctbAddrInRs].getNumberOfPixels());
    m_ctbControllerEngine[ctbAddrInRs].setCtbBpp(bpp);
    m_ctbControllerEngine[ctbAddrInRs].setLastValidQp(lastValidQp);
    m_ctbControllerEngine[ctbAddrInRs].setLastValidLambda(lastValidLambda);
}

int PictureController::getFinishedCtbs()
{
    int counter = 0;
    for(int ctbIdx = 0; ctbIdx < m_pictureSizeInCtbs; ctbIdx++)
    {
        counter += m_ctbControllerEngine[ctbIdx].getFinishedFlag();
    }
    return counter;
}

int PictureController::getTotalCodingBitsCtbs()
{
    int total = 0;
    for(int ctbIdx = 0; ctbIdx < m_pictureSizeInCtbs; ctbIdx++)
    {
        assert(m_ctbControllerEngine[ctbIdx].getFinishedFlag());
        total += m_ctbControllerEngine[ctbIdx].getCodedBits();
    }
    return total;
}

void PictureController::updateAfterEncoding(const int codingBits)
{
    m_codingBits = codingBits;
    m_pictureBpp = static_cast<double>(codingBits) / static_cast<double>(m_pixelsPerPicture);
    m_pictureLambdaFin = m_pictureAlpha * pow(m_pictureBpp, m_pictureBeta);
}


SOPController::SOPController(std::shared_ptr<DataStorage> dataStorageAccess, int sopId, int size, int *bitrateWeight, int targetBits, bool intraInSop) :
    m_size(size),
    m_targetBits(targetBits),
    m_sopId(sopId),
    m_framesCoded(0),
    m_bitsCoded(0),
    m_dataStorageAccess(dataStorageAccess)
{
    m_weight = new int[m_size];
    int totalSum = 0;
    int offset = (int)intraInSop;
    for(int i = 0; i < (m_size - offset); i++)
    {
        totalSum += bitrateWeight[i];
    }

    for(int i = 0; i < (m_size - offset); i++)
    {
        m_weight[i] = (bitrateWeight[i] * targetBits) / totalSum;
    }
}

SOPController::SOPController(std::shared_ptr<DataStorage> dataStorageAccess, int sopId, int size, int *bitrateWeight, int *maxbitrateWeight, bool intraInSop, int totalIntraPeriodBits, int numFramesInIP, int maxSopSize) :
    m_size(size),
    m_sopId(sopId),
    m_framesCoded(0),
    m_bitsCoded(0),
    m_dataStorageAccess(dataStorageAccess)
{
    if(size==1)
    {    
        m_targetBits = (int)(static_cast<double>(totalIntraPeriodBits) / static_cast<double>(numFramesInIP - 1));
        return;
    }
    int totalSum = 0;
    const int numSOPsInIP = static_cast<int>(static_cast<double>(numFramesInIP) / static_cast<double>(maxSopSize));
    for (int i = 0; i < numFramesInIP; i++)
    {
        const int j = (i % maxSopSize);
        const int k = static_cast<int>((static_cast<double>(i) / maxSopSize));
        if (k != (numSOPsInIP - 1))
        {
            totalSum += maxbitrateWeight[j];
        }
        else if (!intraInSop || j != (maxSopSize - 1))
        {
            totalSum += bitrateWeight[j];
        }
    }
    int offset = (int)intraInSop;
    int sopTargetBits = 0;
    m_weight = new int[m_size];
    for (int j = 0; j <  (m_size - offset); j++)
    {
        m_weight[j] = (bitrateWeight[j] * totalIntraPeriodBits) / totalSum;
        sopTargetBits += m_weight[j];
    }
    m_targetBits = sopTargetBits;
}

int SOPController::getRateCurrentPicture(std::shared_ptr<InputQueue::Docket> docket, std::shared_ptr<DataStorage> dataStorage)
{
    const int pocInSop = docket->pocInSop;
    assert(pocInSop <= m_size);

    const int pictureLevel = docket->sopLevel;

    CodedPicture &pictureAtLevel = m_dataStorageAccess->getPictureAtLevel(pictureLevel);
    const int headerBits = pictureAtLevel.getHeaderBitsEstimate();

    int pictureRate = m_weight[pocInSop - 1] - headerBits;
    if(pictureRate < 100)
        pictureRate = 100;
    return pictureRate;
}

void SOPController::updateSopController(int bitsCoded)
{
    m_framesCoded++;
    m_bitsCoded += bitsCoded;
    assert(m_framesCoded <= m_size);
}

void DataStorage::initCtbStorage(int numberOfUnits)
{
    for(int level = 0; level < m_codedCtbsPerLevel.size(); level++)
    {
        m_codedCtbsPerLevel[level] = new CodedCtb[numberOfUnits];
    }
    for (int ctbIdx = 0; ctbIdx < numberOfUnits; ctbIdx++)
    {
        m_codedCtbsPerLevel[0][ctbIdx].setAlpha(ALPHA_INTRA);
        m_codedCtbsPerLevel[0][ctbIdx].setBeta(BETA_INTRA);
    }

}

void DataStorage::initPictureStorage()
{
    m_codedPicturesPerLevel[0].setAlpha(ALPHA_INTRA);
    m_codedPicturesPerLevel[0].setBeta(BETA_INTRA);
    for(int level = 1; level < m_codedPicturesPerLevel.size(); level++)
    {
        m_codedPicturesPerLevel[level].setAlpha(INITIAL_ALPHA);
        m_codedPicturesPerLevel[level].setBeta(INITIAL_BETA);
    }
}

SequenceController::SequenceController(std::shared_ptr<RateControlParameters> rateControlParams)
{
    m_targetRate = rateControlParams->m_targetRate * 1000; // kbps to bps conversion
    m_smoothingWindow = static_cast<int>(rateControlParams->m_frameRate);
    int totalFrames = rateControlParams->m_intraPeriod != 1 ? rateControlParams->m_intraPeriod : ((static_cast<int>(rateControlParams->m_frameRate) + 4) / 8) * 8;
    m_totalFrames = totalFrames;
    m_frameRate = rateControlParams->m_frameRate;
    m_sopSize = rateControlParams->m_sopSize;

    m_averageRate = m_targetRate / m_frameRate;
    m_targetBits = static_cast<int64_t>(m_averageRate * totalFrames);
    m_bitsPerPicture = (int)(m_targetBits / static_cast<double>(totalFrames));
    m_pixelsPerPicture = rateControlParams->m_picHeight * rateControlParams->m_picWidth;
    m_baseQp = rateControlParams->m_baseQp;

    //	Compute the bit per pixel
    m_averageBpp = m_averageRate / (double)m_pixelsPerPicture;
    m_picHeight = rateControlParams->m_picHeight;
    m_picWidth = rateControlParams->m_picWidth;
    m_ctbSize = rateControlParams->m_ctbSize;

    // Initialise the array of CTB controllers
    int picHeightInCtbsY = m_picHeight % m_ctbSize == 0 ? m_picHeight / m_ctbSize : m_picHeight / m_ctbSize + 1;
    int picWidthInCtbsY = m_picWidth  % m_ctbSize == 0 ? m_picWidth / m_ctbSize : m_picWidth / m_ctbSize + 1;
    m_picSizeInCtbsY = picHeightInCtbsY * picWidthInCtbsY;
    m_picHeightInCtbs = picHeightInCtbsY;
    m_picWidthInCtbs = picWidthInCtbsY;

    m_totalCtbs = m_totalFrames * m_picSizeInCtbsY;
    m_averageBitsPerCtb = m_bitsPerPicture / m_picSizeInCtbsY;

    m_useWpp = rateControlParams->m_useWpp;
    m_concurrentFrames = rateControlParams->m_concurrentFrames;
    m_dataStorage.reset(new DataStorage());
    m_dataStorage->initCtbStorage(m_picSizeInCtbsY);
    m_dataStorage->initPictureStorage();
}


SequenceController::SequenceController(double targetRate,
                                       double frameRate,
                                       int    intraPeriod,
                                       int    sopSize,
                                       int    picHeight,
                                       int    picWidth,
                                       int    ctbSize,
                                       int    baseQp,
                                       bool   useWpp,
                                       int    concurrentFrames
                                       )
{
    m_targetRate        = targetRate * 1000; // kbps to bps conversion
    m_smoothingWindow   = static_cast<int>(frameRate);
    int totalFrames     = intraPeriod != 1 ? intraPeriod : ((static_cast<int>(frameRate) + 4)/8)*8;
    m_totalFrames       = totalFrames;
    m_frameRate         = frameRate;
    m_sopSize           = sopSize;

    m_averageRate       = m_targetRate / m_frameRate;
    m_targetBits        = static_cast<int64_t>(m_averageRate * totalFrames);
    m_bitsPerPicture    = (int)(m_targetBits / static_cast<double>(totalFrames));
    m_pixelsPerPicture  = picHeight * picWidth;
    m_baseQp            = baseQp;

    //	Compute the bit per pixel
    m_averageBpp = m_averageRate / (double)m_pixelsPerPicture;
    m_picHeight = picHeight;
    m_picWidth  = picWidth;
    m_ctbSize   = ctbSize;

    // Initialise the array of CTB controllers
    int picHeightInCtbsY = picHeight % ctbSize == 0 ? picHeight / ctbSize : picHeight / ctbSize + 1;
    int picWidthInCtbsY  = picWidth  % ctbSize == 0 ? picWidth  / ctbSize : picWidth  / ctbSize + 1;
    m_picSizeInCtbsY     = picHeightInCtbsY * picWidthInCtbsY;
    m_picHeightInCtbs    = picHeightInCtbsY;
    m_picWidthInCtbs     = picWidthInCtbsY;

    m_totalCtbs          = m_totalFrames * m_picSizeInCtbsY;
    m_averageBitsPerCtb = m_bitsPerPicture / m_picSizeInCtbsY;

    m_useWpp = useWpp;
    m_concurrentFrames = concurrentFrames;
}

void SequenceController::initNewSop(std::shared_ptr<InputQueue::Docket> docket)

{
    // Initialise the frame-based weights according to the value of bpp
    int *bitrateProfile = 0;
    int *maxbitrateProfile = 0;
    int sopOneSizeProfile = 1;
    int sopTargetBits = 0;
    bool intraInSop = false;

    int64_t totalBitsIP = 0;

    if(docket->currentGopSize == 1)
    {
        bitrateProfile = &sopOneSizeProfile;
    }
    else
    {
        if(m_averageBpp <= 0.05)
        {
            bitrateProfile = m_sopWeight[docket->currentGopSize-2][0];
            maxbitrateProfile = m_sopWeight[m_sopSize - 2][0];
        }
        else if(0.05 < m_averageBpp && m_averageBpp <= 0.1)
        {
            bitrateProfile = m_sopWeight[docket->currentGopSize-2][1];
            maxbitrateProfile = m_sopWeight[m_sopSize - 2][1];
        }
        else if(0.1 < m_averageBpp && m_averageBpp <= 0.2)
        {
            bitrateProfile = m_sopWeight[docket->currentGopSize-2][2];
            maxbitrateProfile = m_sopWeight[m_sopSize - 2][2];
        }
        else
        {
            bitrateProfile = m_sopWeight[docket->currentGopSize-2][3];
            maxbitrateProfile = m_sopWeight[m_sopSize - 2][3];
        }
    }

    {
        auto currentIpData = m_ipDataEngine.find(docket->intraFramePoc);
        assert(currentIpData != m_ipDataEngine.end());
        int framesEncoded = 0;
        int64_t bitsSpent = 0;
        int64_t targetBits = 0;
        int bitsPerPictureIp = 0;

        if(docket->poc != docket->intraFramePoc)
        {
            {
                //get the CURRENT intra period:
                auto currentIpData = m_ipDataEngine.find(docket->intraFramePoc);
                assert(currentIpData != m_ipDataEngine.end());
                bitsPerPictureIp = (int)(static_cast<double>(currentIpData->second->m_targetBitsEst) / static_cast<double>(m_totalFrames - 1));
                totalBitsIP = currentIpData->second->m_targetBitsEst;
            }
        }
        else
        {
            intraInSop = true;
            int previousIntraPoc = -1;
            {
                // Get the PREVIOUS intra period:
                int currIntraPoc = docket->intraFramePoc;
                for (auto& ipController : m_ipDataEngine)
                {
                    if (ipController.first > previousIntraPoc && ipController.first < currIntraPoc)
                    {
                        previousIntraPoc = ipController.first;
                    }
                }
                if (previousIntraPoc != -1)
                {
                    auto previousIpData = m_ipDataEngine.find(previousIntraPoc);
                    bitsPerPictureIp = (int)(static_cast<double>(previousIpData->second->m_targetBitsEst) / static_cast<double>(m_totalFrames - 1));
                    totalBitsIP = previousIpData->second->m_targetBitsEst;
                }
            }
        }
        int sopSize = docket->currentGopSize - (int)intraInSop;
        sopTargetBits = max<int>(200, bitsPerPictureIp*sopSize);
    }

    {
        assert(m_sopControllerEngine.find(docket->sopId) == m_sopControllerEngine.end());
        m_sopControllerEngine[docket->sopId] = new SOPController(m_dataStorage, docket->sopId, docket->currentGopSize, bitrateProfile, sopTargetBits, intraInSop);
    }
}
void SequenceController::decreaseNumLeftSameHierarchyLevel(std::shared_ptr<InputQueue::Docket> docket)
{
    int currentPictureLevel = docket->sopLevel;
    int pocInSop = docket->pocInSop;
    int poc = docket->poc;

    {
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        for (auto &pictureController : m_pictureControllerEngine)
        {
            if (pictureController.second->getHierarchyLevel() == currentPictureController->second->getHierarchyLevel() && pictureController.first != poc)
            {
                pictureController.second->decreaseNumLeftSameHierarchyLevel();
                currentPictureController->second->decreaseNumLeftSameHierarchyLevel();
            }
        }

    }
}

void SequenceController::pictureRateAllocation(std::shared_ptr<InputQueue::Docket> docket, std::shared_ptr<CpbInfo> cpbInfo)
{
    int currentPictureLevel = docket->sopLevel;
    int pocInSop            = docket->pocInSop;
    int poc                 = docket->poc;
    int picTargetBits;
    CodedPicture &codedPictureAtLevel = m_dataStorage->getPictureAtLevel(currentPictureLevel);
    codedPictureAtLevel.setLevel(currentPictureLevel);
    
    {
        int sopId = docket->sopId;
        auto currentSopController = m_sopControllerEngine.find(sopId);
        assert(currentSopController != m_sopControllerEngine.end());
        picTargetBits = currentSopController->second->getRateCurrentPicture(docket, m_dataStorage);
    }

    // Cpb correction
    cpbInfo->adjustAllocatedBits(picTargetBits);

    {
        assert(m_pictureControllerEngine.find(poc) == m_pictureControllerEngine.end());

        double bpp = static_cast<double>(picTargetBits) / static_cast<double>(m_pixelsPerPicture);
        double actualbpp = codedPictureAtLevel.getAverageActualBpp();
        double targetbpp = codedPictureAtLevel.getAverageTargetBpp();
        double error = codedPictureAtLevel.isVectorFull() ? abs(actualbpp/targetbpp - 1) : 0.0;
       
        PictureController *currentPictureController = new PictureController(m_dataStorage, m_picHeight, m_picWidth, m_picHeightInCtbs, m_picWidthInCtbs, m_ctbSize, picTargetBits, bpp, docket);
        currentPictureController->setPictureAlpha(codedPictureAtLevel.getAlpha());
        currentPictureController->setPictureBeta(codedPictureAtLevel.getBeta());

        for (auto &pictureController : m_pictureControllerEngine)
        {
            const bool notUpdated = !pictureController.second->getParamsUpdates();
            const bool correctHierarchyLevel = pictureController.second->getHierarchyLevel() == (docket->hierarchyLevel - m_concurrentFrames);
            
            if(notUpdated && correctHierarchyLevel)
            {
                pictureController.second->updateModelParameters();
                const int codingBits = pictureController.second->getCodingBits();
                cpbInfo->updateCpbStatus(codingBits);
                pictureController.second->setParamsUpdates(true);
            }
        }
        m_pictureControllerEngine[poc] = currentPictureController;
    }
}

double SequenceController::estimatePictureLambda(int poc)
{
    double estLambda;
    {
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        double previousLambdaIni = NON_VALID_LAMBDA;
        for (auto& pictureController : m_pictureControllerEngine)
        {
            if (pictureController.second->getHierarchyLevel() == (currentPictureController->second->getHierarchyLevel() - 1))
            {
                previousLambdaIni = pictureController.second->getPictureLambdaIni();
                break;
            }

        }
        estLambda = currentPictureController->second->estimateLambda(previousLambdaIni);
    }
    return estLambda;
}

int SequenceController::deriveQpFromLambda(double lambda, bool isIntra, int sopLevel, int poc)
{
    int estQp;
    int level = isIntra ? 0 : sopLevel;
    {
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        int  previousQp = NON_VALID_QP;
        for (auto& pictureController : m_pictureControllerEngine)
        {
            if (pictureController.second->getHierarchyLevel() == (currentPictureController->second->getHierarchyLevel() - 1))
            {
                previousQp = pictureController.second->getPictureQp();
            }

        }
        estQp = currentPictureController->second->estimateQp(level, previousQp);
    }
    return estQp;
}

void SequenceController::pictureRateAllocationIntra(std::shared_ptr<InputQueue::Docket> docket, std::shared_ptr<CpbInfo> cpbInfo)
{
    double cost = static_cast<double>(docket->icInfo->getSatdSum());
    CodedPicture &pictureAtLevel = m_dataStorage->getPictureAtLevel(0);
    pictureAtLevel.setLevel(0);

    int64_t availableBits = 0;
    {
        auto currentIntraPeriod = m_ipDataEngine.find(docket->intraFramePoc);
        assert(currentIntraPeriod != m_ipDataEngine.end());
        availableBits = currentIntraPeriod->second->m_targetBitsEst;
    }
    int averageBitsLeft = static_cast<int>(static_cast<double>(availableBits) / static_cast<double>(m_totalFrames));
    if (averageBitsLeft < 200)
    {
        averageBitsLeft = 200; 
    }
    double alpha = 0.25, beta = 0.5582;

    if (averageBitsLeft * 40 < m_pixelsPerPicture)
    {
        alpha = 0.25;
    }
    else
    {
        alpha = 0.30;
    }

    int currentBitsPerPicture = (int)(alpha* pow(cost*4.0/(double)averageBitsLeft, beta)*(double)averageBitsLeft + 0.5);

    if(currentBitsPerPicture > (int)(static_cast<double>(availableBits)*0.5))
        currentBitsPerPicture = (int)(static_cast<double>(availableBits)*0.5);

    if(currentBitsPerPicture < 200)
    {
        currentBitsPerPicture = 200;
    }
    
    if(docket->poc != 0)
        cpbInfo->adjustAllocatedBits(currentBitsPerPicture);


    {
        auto currentIntraPeriod = m_ipDataEngine.find(docket->intraFramePoc);
        currentIntraPeriod->second->m_targetBitsEst -= currentBitsPerPicture;
    }
    
    {
        assert(m_pictureControllerEngine.find(docket->poc) == m_pictureControllerEngine.end());
        double bpp = static_cast<double>(currentBitsPerPicture) / static_cast<double>(m_pixelsPerPicture);
        CodedPicture &codedPictureAtLevel = m_dataStorage->getPictureAtLevel(0);
        double actualbpp = codedPictureAtLevel.getAverageActualBpp();
        double targetbpp = codedPictureAtLevel.getAverageTargetBpp();
        double error = codedPictureAtLevel.isVectorFull() ? abs(actualbpp / targetbpp - 1) : 0.0;
        PictureController *currentPictureController = new PictureController(m_dataStorage, m_picHeight, m_picWidth, m_picHeightInCtbs, m_picWidthInCtbs, m_ctbSize, currentBitsPerPicture, bpp, docket);
        currentPictureController->setPictureAlpha(pictureAtLevel.getAlpha());
        currentPictureController->setPictureBeta(pictureAtLevel.getBeta());
        m_pictureControllerEngine[docket->poc] = currentPictureController;
        for (auto &pictureController : m_pictureControllerEngine)
        {
            const bool notUpdated = !pictureController.second->getParamsUpdates();
            const bool correctHierarchyLevel = pictureController.second->getHierarchyLevel() == (docket->hierarchyLevel - m_concurrentFrames);
            if (notUpdated && correctHierarchyLevel)
            {
                pictureController.second->updateModelParameters();
                const int codingBits = pictureController.second->getCodingBits();
                cpbInfo->updateCpbStatus(codingBits);
                pictureController.second->setParamsUpdates(true);
            }
        }
    }
}

int SequenceController::estimateCtbLambdaAndQp(int ctbAddrInRs, int poc)
{
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());
    return currentPictureController->second->estimateCtbLambdaAndQp(ctbAddrInRs);
}

void SequenceController::storeCtbParameters(int codingBitsCtb, bool isIntra, int ctbAddrInRs, int currentPictureLevel, int poc, std::shared_ptr<CpbInfo> cpbInfo)
{
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    currentPictureController->second->storeCtbParameters(codingBitsCtb, isIntra, ctbAddrInRs, currentPictureLevel, cpbInfo);

}

void SequenceController::initNewIntraPeriod(std::shared_ptr<InputQueue::Docket> docket)
{
    int bitsSpent = 0;
    int64_t deltaBits = 0;
    int128bit totalBitsSpent = 0;
    int128bit totalTargetBits = 0;
    int totalFramesCoded = 0;
    const bool carryForward = (docket->absolutePoc != docket->segmentPoc);
    
    if (carryForward)
    { 
        int previousIntraPoc = -1;
        {
            // Get the PREVIOUS intra period:
            int currIntraPoc = docket->intraFramePoc;
            for (auto& ipController : m_ipDataEngine)
            {
                if (ipController.first > previousIntraPoc && ipController.first < currIntraPoc)
                {
                    previousIntraPoc = ipController.first;
                }
            }
        }
        if (previousIntraPoc != -1)
        {
            auto previousIpData = m_ipDataEngine.find(previousIntraPoc);
            totalBitsSpent = previousIpData->second->m_bitsSpent;
            totalTargetBits = previousIpData->second->m_targetBits;
            totalFramesCoded = previousIpData->second->m_framesEncoded;

            for (auto& pictureController : m_pictureControllerEngine)
            {
                const bool correctHierarchyLevel = pictureController.second->getHierarchyLevel() <= (docket->hierarchyLevel - m_concurrentFrames);
                const bool previousIntraPeriod = pictureController.second->getIntraPoc() == previousIntraPoc;
                if (correctHierarchyLevel && previousIntraPeriod)
                {
                    totalBitsSpent += static_cast<int64_t>(pictureController.second->getCodingBits());
                    totalTargetBits += static_cast<int64_t>(pictureController.second->getPictureTargetBits());
                    totalFramesCoded++;
                }
            }
            int128bit tempDelta = totalTargetBits - totalBitsSpent;
            deltaBits = tempDelta.to64();
        }
    }

    {
        assert(m_ipDataEngine.find(docket->intraFramePoc) == m_ipDataEngine.end());
        IntraPeriodData *currentIpData = new IntraPeriodData();
        currentIpData->m_deltaBits = (totalTargetBits - totalBitsSpent);
        currentIpData->m_framesLeft = m_totalFrames;
        currentIpData->m_targetBitsEst = m_targetBits + deltaBits;

        m_ipDataEngine[docket->intraFramePoc] = currentIpData;
    }
}
