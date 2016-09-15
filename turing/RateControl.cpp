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

PictureController::PictureController(DataStorage *storageAccess, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctuSize, int targetBits, double averageBpp, bool isIntra, int sopLevel) :
    m_targetBits(targetBits),
    m_pixelsPerPicture(pictureHeight*pictureWidth),
    m_averageBpp(averageBpp),
    m_dataStorageAccess(storageAccess),
    m_pictureHeightInCtbs(pictureHeightInCtbs),
    m_pictureWidthInCtbs(pictureWidthInCtbs),
    m_pictureSizeInCtbs(pictureHeightInCtbs*pictureWidthInCtbs),
    m_ctusCoded(0),
    m_isIntra(isIntra),
    m_sopLevel(sopLevel),
    m_ctusLeft(m_pictureSizeInCtbs),
    m_totalCostIntra(0)
{
    m_ctuControllerEngine = new CtuController[m_pictureSizeInCtbs];

    for(int r = 0, ctuIdx = 0; r < m_pictureHeightInCtbs; r++)
    {
        const int pixelsInCurrentRow = r == (pictureHeightInCtbs - 1) ? pictureHeight - ctuSize*r : ctuSize;

        for(int c = 0; c < pictureWidthInCtbs; c++, ctuIdx++)
        {
            const int pixelsInCurrentCol = c == (pictureWidthInCtbs - 1) ? pictureWidth - ctuSize*c : ctuSize;
            const int ctuPixels = pixelsInCurrentRow * pixelsInCurrentCol;
            m_ctuControllerEngine[ctuIdx].setNumberOfPixels(ctuPixels);
        }
    }

    if ( averageBpp < 0.03 )
    {
        m_alphaUpdateStep = 0.01;
        m_betaUpdateStep  = 0.005;
    }
    else if ( averageBpp < 0.08 )
    {
        m_alphaUpdateStep = 0.05;
        m_betaUpdateStep  = 0.025;
    }
    else if ( averageBpp < 0.2 )
    {
        m_alphaUpdateStep = 0.1;
        m_betaUpdateStep  = 0.05;
    }
    else if ( averageBpp < 0.5 )
    {
        m_alphaUpdateStep = 0.2;
        m_betaUpdateStep  = 0.1;
    }
    else
    {
        m_alphaUpdateStep = 0.4;
        m_betaUpdateStep  = 0.2;
    }
}

PictureController::PictureController(DataStorage *storageAccess, EstimateIntraComplexity &icInfo, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctuSize, int targetBits, double averageBpp, bool isIntra, int sopLevel) :
    m_targetBits(targetBits),
    m_pixelsPerPicture(pictureHeight*pictureWidth),
    m_averageBpp(averageBpp),
    m_dataStorageAccess(storageAccess),
    m_pictureHeightInCtbs(pictureHeightInCtbs),
    m_pictureWidthInCtbs(pictureWidthInCtbs),
    m_pictureSizeInCtbs(pictureHeightInCtbs*pictureWidthInCtbs),
    m_ctusCoded(0),
    m_isIntra(isIntra),
    m_sopLevel(sopLevel),
    m_ctusLeft(m_pictureSizeInCtbs)
{
    int cost = icInfo.getSatdSum();
    m_totalCostIntra = cost;

    m_ctuControllerEngine = new CtuController[m_pictureSizeInCtbs];

    // Initialise the array of CTU controllers
    for(int r = 0, ctuIdx = 0; r  < m_pictureHeightInCtbs; r++)
    {
        for(int c = 0; c < m_pictureWidthInCtbs; c++, ctuIdx++)
        {
            m_ctuControllerEngine[ctuIdx].setCtuQp(NON_VALID_QP);

            int costIntraCtu = icInfo.getSatdCtu(r, c);
            assert(costIntraCtu != -1);
            m_ctuControllerEngine[ctuIdx].setCostIntra((double)costIntraCtu);

            const int estimatedBits = (int)((double)m_targetBits * (double)costIntraCtu / (double)cost);
            m_ctuControllerEngine[ctuIdx].setTargetBitsEst(estimatedBits);
        }
    }

    m_ctusLeft = m_pictureSizeInCtbs;
    m_ctusCoded = 0;

    if ( averageBpp < 0.03 )
    {
        m_alphaUpdateStep = 0.01;
        m_betaUpdateStep  = 0.005;
    }
    else if ( averageBpp < 0.08 )
    {
        m_alphaUpdateStep = 0.05;
        m_betaUpdateStep  = 0.025;
    }
    else if ( averageBpp < 0.2 )
    {
        m_alphaUpdateStep = 0.1;
        m_betaUpdateStep  = 0.05;
    }
    else if ( averageBpp < 0.5 )
    {
        m_alphaUpdateStep = 0.2;
        m_betaUpdateStep  = 0.1;
    }
    else
    {
        m_alphaUpdateStep = 0.4;
        m_betaUpdateStep  = 0.2;
    }
}

void PictureController::updateLambdaModelIntra(double &alpha,
                                               double &beta,
                                               CodedPicture &picture,
                                               int actualBits)
{
    double lnbpp = log(pow(picture.getCostIntra() / (double)m_pixelsPerPicture, BETA_INTRA));
    double diffLambda = beta*(log((double)actualBits)-log((double)m_targetBits));

    diffLambda = Clip3(-0.125, 0.125, 0.25*diffLambda);
    alpha    =  alpha * exp(diffLambda);
    beta     =  beta + diffLambda / lnbpp;
}

double PictureController::estimateLambda(int level, bool isIntra)
{
    double bpp = (double)m_targetBits / (double)m_pixelsPerPicture;

    CodedPicture &currentPicture = m_dataStorageAccess->getPictureAtLevel(level);
    list<CodedPicture> &previousPictures = m_dataStorageAccess->getPreviousCodedPictures();

    double alpha = currentPicture.getAlpha();
    double beta  = currentPicture.getBeta();

    //	Compute lambda according to the model for intra or inter
    double lambda;
    if(isIntra)
    {
        auto costIntra = static_cast<int>(currentPicture.getCostIntra());
        double madPixelBased = pow((double)costIntra/(double)m_pixelsPerPicture, BETA_INTRA_MAD);
        lambda = (alpha/256.0) * pow( madPixelBased/ (double)bpp, beta );
    }
    else
    {
        lambda = alpha * pow(bpp, beta);
    }

    list<CodedPicture>::iterator it;
    double lastLevelLambda = -1.0;
    double lastPicLambda   = -1.0;
    double lastValidLambda = -1.0;

    for ( it = previousPictures.begin(); it != previousPictures.end(); it++ )
    {
        if ( (*it).getLevel() == level )
        {
            lastLevelLambda = (*it).getLambda();
        }
        lastPicLambda     = (*it).getLambda();

        if ( lastPicLambda > 0.0 )
        {
            lastValidLambda = lastPicLambda;
        }
    }

    if ( lastLevelLambda > 0.0 )
    {
        lastLevelLambda = Clip3( 0.1, 10000.0, lastLevelLambda );
        lambda = Clip3( lastLevelLambda * pow( 2.0, -3.0/3.0 ), lastLevelLambda * pow( 2.0, 3.0/3.0 ), lambda );
    }

    if ( lastPicLambda > 0.0 )
    {
        lastPicLambda = Clip3( 0.1, 2000.0, lastPicLambda );
        lambda = Clip3( lastPicLambda * pow( 2.0, -10.0/3.0 ), lastPicLambda * pow( 2.0, 10.0/3.0 ), lambda );
    }
    else if ( lastValidLambda > 0.0 )
    {
        lastValidLambda = Clip3( 0.1, 2000.0, lastValidLambda );
        lambda = Clip3( lastValidLambda * pow(2.0, -10.0/3.0), lastValidLambda * pow(2.0, 10.0/3.0), lambda );
    }
    else
    {
        lambda = Clip3( 0.1, 10000.0, lambda );
    }

    if ( lambda < 0.1 )
    {
        lambda = 0.1;
    }

    currentPicture.setLambda(lambda);

    // Initialise the CTU weights with the estimated lambda
    CodedCtu *codedCtuAtLevel = m_dataStorageAccess->getCodedCtuAtLevel(level);
    double totalWeight = 0.0;
    for(int ctuIdx = 0; ctuIdx < m_pictureSizeInCtbs; ctuIdx++)
    {
        const double alpha = codedCtuAtLevel[ctuIdx].getAlpha();
        const double beta  = codedCtuAtLevel[ctuIdx].getBeta();

        double ctuWeight = pow(lambda/alpha, 1.0/beta) * m_ctuControllerEngine[ctuIdx].getNumberOfPixels();
        if(ctuWeight < 0.01)
        {
            ctuWeight = 0.01;
        }
        m_ctuControllerEngine[ctuIdx].setCtuWeight(ctuWeight);
        totalWeight += ctuWeight;
    }

    // Weight normalisation
    for(int ctuIdx = 0; ctuIdx < m_pictureSizeInCtbs; ctuIdx++)
    {
        const double ctuWeightNormalised = m_ctuControllerEngine[ctuIdx].getCtuWeight() / totalWeight * m_targetBits;
        m_ctuControllerEngine[ctuIdx].setCtuWeight(ctuWeightNormalised);
    }

    return lambda;
}

void PictureController::updatePictureController(int currentLevel,
                                                bool isIntra,
                                                double &lastCodedLambda)
{
    CodedPicture &currentPicture = m_dataStorageAccess->getPictureAtLevel(currentLevel);

    const int bitsSpent = getTotalCodingBits();

    double currentAlpha = currentPicture.getAlpha();
    double currentBeta  = currentPicture.getBeta();

    if(isIntra)
    {
        updateLambdaModelIntra(currentAlpha, currentBeta, currentPicture, bitsSpent);
    }
    else
    {
        double bpp          = (double)bitsSpent / (double)m_pixelsPerPicture;

        // Compute lambda using the model and the actual bits spent
        double lambdaComputed  = currentAlpha * pow(bpp, currentBeta);
        double lambdaEstimated = currentPicture.getLambda();

        if ( lambdaEstimated < 0.01 || lambdaComputed < 0.01 || bpp < 0.0001 )
        {
            currentAlpha *= ( 1.0 - m_alphaUpdateStep / 2.0 );
            currentBeta  *= ( 1.0 - m_betaUpdateStep / 2.0 );

            currentAlpha = Clip3( ALPHA_MIN, ALPHA_MAX, currentAlpha );
            currentBeta  = Clip3( BETA_MIN,  BETA_MAX,  currentBeta  );

            currentPicture.setAlpha(currentAlpha);
            currentPicture.setBeta(currentBeta);

            return;
        }

        lambdaComputed = Clip3( lambdaEstimated / 10.0, lambdaEstimated * 10.0, lambdaComputed );
        currentAlpha += m_alphaUpdateStep * ( log( lambdaEstimated ) - log( lambdaComputed ) ) * currentAlpha;
        double lnbpp = log( bpp );
        lnbpp = Clip3( -5.0, -0.1, lnbpp );

        currentBeta  += m_betaUpdateStep * ( log( lambdaEstimated ) - log( lambdaComputed ) ) * lnbpp;

        currentAlpha = Clip3( ALPHA_MIN, ALPHA_MAX, currentAlpha );
        currentBeta  = Clip3( BETA_MIN,  BETA_MAX,  currentBeta  );
    }

    currentPicture.setAlpha(currentAlpha);
    currentPicture.setBeta(currentBeta);
    currentPicture.setCodingBits(bitsSpent - currentPicture.getHeaderBits());

    if ( currentPicture.getLevel() == 1 )
    {
        double currLambda = Clip3( 0.1, 10000.0, currentPicture.getLambda() );
        lastCodedLambda   = LAMBDA_WEIGHT * lastCodedLambda + (1 - LAMBDA_WEIGHT) * currLambda;
    }
    int counter = getFinishedCtus();
    assert(m_pictureSizeInCtbs == counter);
}

void PictureController::getCodedInfoFromWavefront(bool wpp, int ctbAddrInRs, int &cumulativeMad, int &cumulativeBitsEstimated, int &cumulativeBitsSpent, int &cumulativeCtusCoded, double &cumulativeWeight, bool currFrame)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    int rowStart = ctbAddrInRs / m_pictureWidthInCtbs;
    int colStart = (ctbAddrInRs % m_pictureWidthInCtbs) - (int)currFrame;
    cumulativeMad = cumulativeBitsEstimated = cumulativeBitsSpent = cumulativeCtusCoded = 0;
    cumulativeWeight = 0.0;



    for (int r = rowStart; r >= 0; r--)
    {
        for (int c = colStart; c >= 0; c--)
        {
            int currentIdx = r*m_pictureWidthInCtbs + c;
            assert(0 <= currentIdx && currentIdx < m_pictureSizeInCtbs);
            assert(m_ctuControllerEngine[currentIdx].getFinishedFlag());
            cumulativeMad += static_cast<int>(m_ctuControllerEngine[currentIdx].getCostIntra());
            cumulativeBitsEstimated += m_ctuControllerEngine[currentIdx].getTargetBitsEst();
            cumulativeBitsSpent += m_ctuControllerEngine[currentIdx].getCodedBits();
            cumulativeWeight += m_ctuControllerEngine[currentIdx].getCtuWeight();
            cumulativeCtusCoded++;
        }
        colStart = wpp ? (min<int>(m_pictureWidthInCtbs - 1, colStart + 2)) : (m_pictureWidthInCtbs - 1);
    }

}

int PictureController::estimateCtuLambdaAndQp(bool isIntra, int ctbAddrInRs, int pictureLevel, int sliceQp)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    if(isIntra)
    {
        return estimateCtuLambdaAndQpIntra(sliceQp, ctbAddrInRs);
    }
    else
    {
        return estimateCtuLambdaAndQpInter(ctbAddrInRs, pictureLevel);
    }
}

int PictureController::estimateCtuLambdaAndQpInter(int ctbAddrInRs, int pictureLevel)
{
    // Step 1: Estimate the lambda from lambda = alpha * (bpp)^beta;
    CodedCtu *codedCtuAtLevel = m_dataStorageAccess->getCodedCtuAtLevel(pictureLevel);
    const double alpha = codedCtuAtLevel[ctbAddrInRs].getAlpha();
    const double beta  = codedCtuAtLevel[ctbAddrInRs].getBeta();
    const double bpp   = m_ctuControllerEngine[ctbAddrInRs].getCtuBpp();

    double estLambdaCtu = alpha * pow(bpp, beta);

    CodedPicture &currentPicture = m_dataStorageAccess->getPictureAtLevel(pictureLevel);
    const double pictureLambda = currentPicture.getLambda();

    if ( pictureLambda > 0.0 )
    {
        estLambdaCtu = Clip3( pictureLambda * pow( 2.0, -2.0/3.0 ), pictureLambda * pow( 2.0, 2.0/3.0 ), estLambdaCtu );
    }
    else
    {
        estLambdaCtu = Clip3( 10.0, 1000.0, estLambdaCtu );
    }

    if ( estLambdaCtu < 0.1 )
    {
        estLambdaCtu = 0.1;
    }

    // Step 2: Estimate the QP from lambda using the nonlinear relationship
    m_ctuControllerEngine[ctbAddrInRs].setCtuLambda(estLambdaCtu);
    m_ctuControllerEngine[ctbAddrInRs].setCtuReciprocalLambda(1.0 / estLambdaCtu);
    m_ctuControllerEngine[ctbAddrInRs].setCtuReciprocalSqrtLambda(1.0 / sqrt(estLambdaCtu));
    int ctuEstQp = (int)( 4.2005 * log( estLambdaCtu ) + 13.7122 + 0.5 );

    int pictureQp = currentPicture.getQp();

    ctuEstQp = Clip3( pictureQp - 2, pictureQp + 2, ctuEstQp );

    ctuEstQp = Clip3(0, 51, ctuEstQp);

    m_ctuControllerEngine[ctbAddrInRs].setCtuQp(ctuEstQp);

    return ctuEstQp;
}

int PictureController::estimateCtuLambdaAndQpIntra(int sliceQp, int ctbAddrInRs)
{
    // Step 1: Estimate the lambda from the model lambda = (alpha/256) * (cost/bpp)^beta
    CodedPicture &pictureAtLevel = m_dataStorageAccess->getPictureAtLevel(0); // Level 0 for intra pictures only

    double   alpha = pictureAtLevel.getAlpha();
    double   beta  = pictureAtLevel.getBeta();

    double costPixelBased = m_ctuControllerEngine[ctbAddrInRs].getCostIntra() /
            static_cast<double>(m_ctuControllerEngine[ctbAddrInRs].getNumberOfPixels());
    const double bpp = m_ctuControllerEngine[ctbAddrInRs].getCtuBpp();
    costPixelBased = pow(costPixelBased, BETA_INTRA_MAD);
    double lambda = (alpha/256.0) * pow( costPixelBased/bpp, beta );

    const int minQP = sliceQp - 2;
    const int maxQP = sliceQp + 2;

    double maxLambda = exp(((double)(maxQP+0.49)-13.7122)/4.2005);
    double minLambda = exp(((double)(minQP-0.49)-13.7122)/4.2005);

    lambda = Clip3(minLambda, maxLambda, lambda);

    // Step 2: Estimate the QP using the nonlinear relationship
    int qp = int( 4.2005 * log(lambda) + 13.7122 + 0.5 );
    qp = Clip3(minQP, maxQP, qp);
    qp = Clip3(0, 51, qp);
    m_ctuControllerEngine[ctbAddrInRs].setCtuQp(qp);
    m_ctuControllerEngine[ctbAddrInRs].setCtuLambda(lambda);
    m_ctuControllerEngine[ctbAddrInRs].setCtuReciprocalLambda(1.0 / lambda);
    m_ctuControllerEngine[ctbAddrInRs].setCtuReciprocalLambda(1.0 / sqrt(lambda));
    return qp;
}

void PictureController::updateCtuController(int codingBits, bool isIntra, int ctbAddrInRs, int pictureLevel)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    m_ctuControllerEngine[ctbAddrInRs].setCodedBits(codingBits);
    m_ctuControllerEngine[ctbAddrInRs].setFinishedFlag(true);

    bool isValidCtu = m_ctuControllerEngine[ctbAddrInRs].getValidityFlag(); // i.e. is not a skipped CTU
    if(!isValidCtu)
    {
        m_ctuControllerEngine[ctbAddrInRs].setCtuQp(NON_VALID_QP);
    }

    if ( isIntra )
    {
        m_ctusLeft--;
        m_ctusCoded++;
        return;
    }

    CodedCtu *codedCtuAtLevel = m_dataStorageAccess->getCodedCtuAtLevel(pictureLevel);
    double alpha = codedCtuAtLevel[ctbAddrInRs].getAlpha();
    double beta  = codedCtuAtLevel[ctbAddrInRs].getBeta();

    int ctuNumberOfPixels = m_ctuControllerEngine[ctbAddrInRs].getNumberOfPixels();
    const double bpp      = ( double )codingBits/( double )ctuNumberOfPixels;
    double computedLambda = alpha * pow( bpp, beta );
    double ctuLambda      = m_ctuControllerEngine[ctbAddrInRs].getCtuLambda();

    if( ctuLambda < 0.01 || computedLambda < 0.01 || bpp < 0.0001 )
    {
        alpha *= ( 1.0 - m_alphaUpdateStep / 2.0 );
        beta  *= ( 1.0 - m_betaUpdateStep  / 2.0 );

        alpha = Clip3( ALPHA_MIN, ALPHA_MAX, alpha );
        beta  = Clip3( BETA_MIN,  BETA_MAX,  beta  );

        codedCtuAtLevel[ctbAddrInRs].setAlpha(alpha);
        codedCtuAtLevel[ctbAddrInRs].setBeta(beta);

        m_ctusLeft--;
        m_ctusCoded++;
        return;
    }

    computedLambda = Clip3( ctuLambda / 10.0, ctuLambda * 10.0, computedLambda );
    alpha += m_alphaUpdateStep * ( log( ctuLambda ) - log( computedLambda ) ) * alpha;
    double lnbpp = log( bpp );
    lnbpp = Clip3( -5.0, -0.1, lnbpp );
    beta  += m_betaUpdateStep * ( log( ctuLambda ) - log( computedLambda ) ) * lnbpp;

    alpha = Clip3( ALPHA_MIN, ALPHA_MAX, alpha );
    beta  = Clip3( BETA_MIN,  BETA_MAX,  beta  );

    codedCtuAtLevel[ctbAddrInRs].setAlpha(alpha);
    codedCtuAtLevel[ctbAddrInRs].setBeta(beta);

    m_ctusLeft--;
    m_ctusCoded++;
}

void PictureController::getAveragePictureQpAndLambda(int &averageQp, double &averageLambda)
{
    int counter = getFinishedCtus();
    assert(counter == m_pictureSizeInCtbs);
    int ctuCnt = 0;
    int sumQp  = 0;
    double sumLambda = 0.0;
    for(int ctuIdx = 0; ctuIdx < m_pictureSizeInCtbs; ctuIdx++)
    {
        int currentCtuQp = m_ctuControllerEngine[ctuIdx].getCtuQp();
        double currentCtuLambda = m_ctuControllerEngine[ctuIdx].getCtuLambda();
        if(currentCtuQp != NON_VALID_QP)
        {
            ctuCnt++;
            sumQp += currentCtuQp;
            sumLambda += currentCtuLambda;
        }
    }
    if(ctuCnt)
    {
        averageQp     = sumQp / ctuCnt;
        averageQp     = Clip3(0, 51, averageQp);
        averageLambda = sumLambda / ctuCnt;
    }
    else
    {
        averageQp     = NON_VALID_QP;
        averageLambda = NON_VALID_LAMBDA;
    }
}


void   PictureController::computeCtuTargetBits(bool isIntraSlice, int ctbAddrInRs, int cumulativeMad, int cumulativeBitsEstimated, int cumulativeBitsSpent, int cumulativeCtusCoded, double cumulativeWeight)
{

    int avgBits = 0;
    int currentCtusLeft = m_pictureSizeInCtbs - cumulativeCtusCoded;
    int influence = std::min<int>(CTU_SMOOTH_WINDOW, currentCtusLeft);

    if (isIntraSlice)
    {
        double MAD = m_ctuControllerEngine[ctbAddrInRs].getCostIntra();
        int remainingCostIntra = m_totalCostIntra - cumulativeMad;

        if (remainingCostIntra > 0.1)
        {
            int estimateBitsCtu = static_cast<int>((MAD / static_cast<double>(cumulativeMad + MAD)) * static_cast<double>(m_targetBits));
            avgBits = estimateBitsCtu + (cumulativeBitsEstimated - cumulativeBitsSpent) / influence;
        }
        else
        {
            avgBits = static_cast<int>((m_targetBits - cumulativeBitsSpent) / currentCtusLeft);
        }
    }
    else
    {
        const double currentCtuWeight = m_ctuControllerEngine[ctbAddrInRs].getCtuWeight();
        const double totalWeight = cumulativeWeight + currentCtuWeight;
        avgBits = static_cast<int>(currentCtuWeight + (totalWeight - static_cast<double>(cumulativeBitsSpent)) / static_cast<double>(influence));
    }

    if (avgBits < 1)
    {
        avgBits = 1;
    }
    m_ctuControllerEngine[ctbAddrInRs].setTargetBits(avgBits);
    const double bpp = static_cast<double>(avgBits) / static_cast<double>(m_ctuControllerEngine[ctbAddrInRs].getNumberOfPixels());
    m_ctuControllerEngine[ctbAddrInRs].setCtuBpp(bpp);
}

int PictureController::getFinishedCtus()
{
    int counter = 0;
    for(int ctuIdx = 0; ctuIdx < m_pictureSizeInCtbs; ctuIdx++)
    {
        counter += m_ctuControllerEngine[ctuIdx].getFinishedFlag();
    }
    return counter;
}

int PictureController::getTotalCodingBits()
{
    int total = 0;
    for(int ctuIdx = 0; ctuIdx < m_pictureSizeInCtbs; ctuIdx++)
    {
        assert(m_ctuControllerEngine[ctuIdx].getFinishedFlag());
        total += m_ctuControllerEngine[ctuIdx].getCodedBits();
    }
    return total;
}

SOPController::SOPController(DataStorage *dataAccess, int sopId, int size, int *bitrateWeight, int targetBits) :
    m_size(size),
    m_targetBits(targetBits),
    m_dataStorageAccess(dataAccess),
    m_sopId(sopId),
    m_framesCoded(0),
    m_bitsCoded(0)
{
    m_weight = new int[m_size];
    int totalSum = 0;
    for(int i = 0; i < m_size; i++)
    {
        totalSum += bitrateWeight[i];
    }

    for(int i = 0; i < m_size; i++)
    {
        m_weight[i] = (bitrateWeight[i] * targetBits) / totalSum;
    }
}

int SOPController::getEstimatedHeaderBits(int level)
{
    int totalPreviousPics = 0;
    int totalPreviousBits = 0;

    list<CodedPicture> &previousPictures = m_dataStorageAccess->getPreviousCodedPictures();

    list<CodedPicture>::iterator it;
    for ( it = previousPictures.begin(); it != previousPictures.end(); it++ )
    {
        if ( (*it).getLevel() == level )
        {
            totalPreviousBits += (*it).getHeaderBits();
            totalPreviousPics++;
        }
    }

    int estHeaderBits = 0;
    if ( totalPreviousPics > 0 )
    {
        estHeaderBits = totalPreviousBits / totalPreviousPics;
    }

    return estHeaderBits;
}

int SOPController::getRateCurrentPicture(int pocInSop)
{
    assert(pocInSop <= m_size);
    const int headerBits = 0; // To be updated with the right function call
    int pictureRate = m_weight[pocInSop - 1] - headerBits;
    if(pictureRate < 100)
        pictureRate = 100; // As in HM
    return pictureRate;
}

void SOPController::updateSopController(int bitsCoded)
{
    m_framesCoded++;
    m_bitsCoded += bitsCoded;
    assert(m_framesCoded <= m_size);
}

void DataStorage::addCodedPicture(int level)
{
    assert(level < MAX_NUM_LEVELS);
    if(m_previousCodedPictures.size() > MAX_LIST_SIZE)
    {
        m_previousCodedPictures.pop_front();
    }
    m_previousCodedPictures.push_back(m_codedPicturesPerLevel[level]);
}

void DataStorage::initCtuStorage(int numberOfUnits)
{
    for(int level = 0; level < m_codedCtusPerLevel.size(); level++)
    {
        m_codedCtusPerLevel[level] = new CodedCtu[numberOfUnits];
    }
}

void DataStorage::resetPreviousCodedPicturesBuffer()
{
    while(m_previousCodedPictures.size() > 0)
    {
        m_previousCodedPictures.pop_front();
    }
}

void DataStorage::resetCodedPicturesAtLevelMemory()
{
    m_codedPicturesPerLevel[0].setAlpha(ALPHA_INTRA);
    m_codedPicturesPerLevel[0].setBeta(BETA_INTRA);
    for(int level = 1; level < m_codedPicturesPerLevel.size(); level++)
    {
        m_codedPicturesPerLevel[level].setAlpha(INITIAL_ALPHA);
        m_codedPicturesPerLevel[level].setBeta(INITIAL_BETA);
    }
}

void DataStorage::resetCodedCtuAtLevelMemory(int totalCtusPerFrame)
{
    for(int level = 0; level < m_codedCtusPerLevel.size(); level++)
    {
        CodedCtu *codedCtuAtLevel = m_codedCtusPerLevel[level];
        for(int ctuIdx = 0; ctuIdx < totalCtusPerFrame; ctuIdx++)
        {
            codedCtuAtLevel[ctuIdx].setAlpha(INITIAL_ALPHA);
            codedCtuAtLevel[ctuIdx].setBeta(INITIAL_BETA);
        }
    }
}

SequenceController::SequenceController(double targetRate,
                                       double frameRate,
                                       int intraPeriod,
                                       int sopSize,
                                       int picHeight,
                                       int picWidth,
                                       int ctuSize,
                                       int baseQp)
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
    m_bitsLeft          = m_targetBits;
    m_pixelsPerPicture  = picHeight * picWidth;
    m_baseQp            = baseQp;

    //	Compute the bit per pixel
    m_averageBpp = m_averageRate / (double)m_pixelsPerPicture;

    m_lastCodedPictureLambda = 0.0;
    m_picHeight = picHeight;
    m_picWidth  = picWidth;
    m_ctuSize   = ctuSize;

    // Initialise the array of CTU controllers
    int picHeightInCtbsY = picHeight % ctuSize == 0 ? picHeight / ctuSize : picHeight / ctuSize + 1;
    int picWidthInCtbsY  = picWidth  % ctuSize == 0 ? picWidth  / ctuSize : picWidth  / ctuSize + 1;
    m_picSizeInCtbsY     = picHeightInCtbsY * picWidthInCtbsY;
    m_picHeightInCtbs    = picHeightInCtbsY;
    m_picWidthInCtbs     = picWidthInCtbsY;

    m_dataStorageEngine       = new DataStorage();
    m_dataStorageEngine->initCtuStorage(m_picSizeInCtbsY);
    m_averageBitsPerCtb = m_bitsPerPicture / m_picSizeInCtbsY;

    CodedPicture &pictureIntra = m_dataStorageEngine->getPictureAtLevel(0);
    pictureIntra.setAlpha(ALPHA_INTRA);
    pictureIntra.setBeta(BETA_INTRA);

    m_codedFrames = 0;
    m_bitsSpent = 0;

#if WRITE_RC_LOG
    m_logFile.open("rcLogFileCpb.txt", ios::out);
    m_logFile<<"----------------------------------------------------------------------------------\n";
    m_logFile<<"|    POC |     Target |    Lambda |   QP | Coded bits |       Total | CPB Status |\n";
    m_logFile<<"----------------------------------------------------------------------------------\n";
    m_logFile.flush();
#endif
}

void SequenceController::initNewSop(int sopId, int sopSize)
{
    // Initialise the frame-based weights according to the value of bpp
    int *bitrateProfile = 0;
    int sopOneSizeProfile = 1;
    if(sopSize == 1)
    {
        bitrateProfile = &sopOneSizeProfile;
    }
    else
    {
        if(m_averageBpp <= 0.05)
        {
            bitrateProfile = m_sopWeight[sopSize-2][0];
        }
        else if(0.05 < m_averageBpp && m_averageBpp <= 0.1)
        {
            bitrateProfile = m_sopWeight[sopSize-2][1];
        }
        else if(0.1 < m_averageBpp && m_averageBpp <= 0.2)
        {
            bitrateProfile = m_sopWeight[sopSize-2][2];
        }
        else
            bitrateProfile = m_sopWeight[sopSize-2][3];
    }

    unique_lock<mutex> lock(m_sopControllerMutex);
    assert(m_sopControllerEngine.find(sopId) == m_sopControllerEngine.end());
    m_sopControllerEngine[sopId] = new SOPController(m_dataStorageEngine, sopId, sopSize, bitrateProfile, m_bitsPerPicture*sopSize);
}
void SequenceController::pictureRateAllocation(std::shared_ptr<InputQueue::Docket> docket)
{
    int currentPictureLevel = docket->sopLevel;
    int pocInSop            = docket->pocInSop;
    int poc                 = docket->poc;
    int picTargetBits;
    CodedPicture &codedPictureAtLevel = m_dataStorageEngine->getPictureAtLevel(currentPictureLevel);
    codedPictureAtLevel.setLevel(currentPictureLevel);
    {
        int sopId = docket->sopId;
        unique_lock<mutex> lock(m_sopControllerMutex);
        auto currentSopController = m_sopControllerEngine.find(sopId);
        assert(currentSopController != m_sopControllerEngine.end());
        picTargetBits = currentSopController->second->getRateCurrentPicture(pocInSop);
    }

    // Cpb correction
    m_cpbControllerEngine.adjustAllocatedBits(picTargetBits);

    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        assert(m_pictureControllerEngine.find(poc) == m_pictureControllerEngine.end());
        m_pictureControllerEngine[poc] = new PictureController(m_dataStorageEngine, m_picHeight, m_picWidth, m_picHeightInCtbs, m_picWidthInCtbs, m_ctuSize, picTargetBits, m_averageBpp,false,currentPictureLevel);
    }

}


void SequenceController::updateSequenceController(bool isIntra, int sopLevel, int poc, int sopId)
{
    int codingBits, pictureTargetBits;
    // Take the lock on the picture controller and update it
    {
        unique_lock<mutex> lockPicture(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        codingBits = currentPictureController->second->getTotalCodingBits();
        pictureTargetBits = currentPictureController->second->getPictureTargetBits();
        m_bitsLeft -= static_cast<int64_t>(codingBits);
        m_codedFrames++;
        m_bitsSpent += codingBits;
        double lambda;
        int qp;

        currentPictureController->second->getAveragePictureQpAndLambda(qp, lambda);
        int currentPictureLevel = isIntra ? 0 : sopLevel;
        CodedPicture &pictureAtLevel = m_dataStorageEngine->getPictureAtLevel(currentPictureLevel);
        if(!isIntra)
        {
            pictureAtLevel.setQp(qp);
            pictureAtLevel.setLambda(lambda);
        }
        pictureAtLevel.setPoc(poc);
        m_dataStorageEngine->addCodedPicture(currentPictureLevel);
        currentPictureController->second->updatePictureController(currentPictureLevel,
                                                                  isIntra,
                                                                  m_lastCodedPictureLambda);
        assert(m_pictureControllerEngine.erase(poc) == 1);
        // Lock released
    }

    {
        unique_lock<mutex> lockSop(m_sopControllerMutex);
        auto currentSopController = m_sopControllerEngine.find(sopId);
        assert(currentSopController != m_sopControllerEngine.end());

        if(!isIntra)
        {
            currentSopController->second->updateSopController(codingBits);
        }
        else
        {
            currentSopController->second->updateSopController(pictureTargetBits);
        }
        // Remove this SOP controller if all frames have been encoded
        if(currentSopController->second->finished())
        {
            assert(m_sopControllerEngine.erase(sopId) == 1);
        }
    }

    // Update Cpb status
    m_cpbControllerEngine.updateCpbStatus(codingBits);
}

double SequenceController::estimatePictureLambda(bool isIntra, int sopLevel, int poc)
{
    int currentPictureLevel = isIntra ? 0 : sopLevel;

    double estLambda;
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());

        estLambda = currentPictureController->second->estimateLambda(currentPictureLevel, isIntra);
    }

    return estLambda;
}

int SequenceController::deriveQpFromLambda(double lambda, bool isIntra, int sopLevel)
{
    int qp = int( 4.2005 * log( lambda ) + 13.7122 + 0.5 );

    int level = isIntra ? 0 : sopLevel;

    int lastLevelQp = NON_VALID_QP;
    int lastPicQp   = NON_VALID_QP;
    int lastValidQp = NON_VALID_QP;
    list<CodedPicture>::iterator it;
    list<CodedPicture> &previousCodedPictures = m_dataStorageEngine->getPreviousCodedPictures();
    for ( it = previousCodedPictures.begin(); it != previousCodedPictures.end(); it++ )
    {
        if ( (*it).getLevel() == level )
        {
            lastLevelQp = (*it).getQp();
        }
        lastPicQp = (*it).getQp();
        if ( lastPicQp > NON_VALID_QP )
        {
            lastValidQp = lastPicQp;
        }
    }

    if ( lastLevelQp > NON_VALID_QP )
    {
        qp = Clip3( lastLevelQp - 3, lastLevelQp + 3, qp );
    }

    if( lastPicQp > NON_VALID_QP )
    {
        qp = Clip3( lastPicQp - 10, lastPicQp + 10, qp );
    }
    else if( lastValidQp > NON_VALID_QP )
    {
        qp = Clip3( lastValidQp - 10, lastValidQp + 10, qp );
    }

    qp = Clip3(0, 51, qp);

    CodedPicture &pictureAtLevel = m_dataStorageEngine->getPictureAtLevel(level);
    pictureAtLevel.setQp(qp);

    return qp;
}

void SequenceController::pictureRateAllocationIntra(EstimateIntraComplexity &icInfo, int poc)
{
    int cost = icInfo.getSatdSum();
    CodedPicture &pictureAtLevel = m_dataStorageEngine->getPictureAtLevel(0);
    pictureAtLevel.setCostIntra((double)cost);
    pictureAtLevel.setLevel(0);

    // Get average bits left
    int averageBitsLeft = static_cast<int>(m_bitsLeft / static_cast<int64_t>(m_totalFrames - m_codedFrames));

    // Refine the bits estimate for intra
    double alpha=0.25, beta=0.5582;

    if (averageBitsLeft*40 < m_pixelsPerPicture)
    {
        alpha=0.25;
    }
    else
    {
        alpha=0.30;
    }

    int currentBitsPerPicture = (int)(alpha* pow(cost*4.0/(double)averageBitsLeft, beta)*(double)averageBitsLeft + 0.5);

    // Cpb correction
    m_cpbControllerEngine.adjustAllocatedBits(currentBitsPerPicture);

    if(currentBitsPerPicture < 200)
    {
        currentBitsPerPicture = 200; // As in HM RC
    }

    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        assert(m_pictureControllerEngine.find(poc) == m_pictureControllerEngine.end());
        m_pictureControllerEngine[poc] = new PictureController(m_dataStorageEngine, icInfo, m_picHeight, m_picWidth, m_picHeightInCtbs, m_picWidthInCtbs, m_ctuSize, currentBitsPerPicture, m_averageBpp,false,0);
    }


}

void SequenceController::setHeaderBits(int bits, bool isIntra, int sopLevel, int poc)
{
    // First check if the current POC has already been stored in the previous coded pictures memory (it is done in a different thread)
    const bool found = insertHeaderBitsData(bits, poc);
    if(!found)
    {
        // Not found in the previous coded pictures, which means it is still being processed by the substream thread, so it's in the picture at level memory
        int currentPictureLevel = isIntra ? 0 : sopLevel;
        CodedPicture &pictureAtLevel = m_dataStorageEngine->getPictureAtLevel(currentPictureLevel);
        pictureAtLevel.setHeaderBits(bits);
    }
}

int SequenceController::estimateCtuLambdaAndQp(bool isIntra, int ctbAddrInRs, int currentPictureLevel, int poc, int sliceQp)
{
    unique_lock<mutex> lock(m_pictureControllerMutex);
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());
    return currentPictureController->second->estimateCtuLambdaAndQp(isIntra, ctbAddrInRs, currentPictureLevel, sliceQp);
}

void SequenceController::updateCtuController(int codingBitsCtu, bool isIntra, int ctbAddrInRs, int currentPictureLevel, int poc)
{
    unique_lock<mutex> lock(m_pictureControllerMutex);
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    currentPictureController->second->updateCtuController(codingBitsCtu, isIntra, ctbAddrInRs, currentPictureLevel);

}

void SequenceController::getAveragePictureQpAndLambda(int &averageQp, double &averageLambda, int poc)
{
    unique_lock<mutex> lock(m_pictureControllerMutex);
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    currentPictureController->second->getAveragePictureQpAndLambda(averageQp, averageLambda);
}

void SequenceController::resetSequenceControllerMemory()
{
    // Reset the buffers used for prediction and parameter smoothing
    m_dataStorageEngine->resetPreviousCodedPicturesBuffer();
    m_dataStorageEngine->resetCodedPicturesAtLevelMemory();
    m_dataStorageEngine->resetCodedCtuAtLevelMemory(m_picSizeInCtbsY);
}

bool SequenceController::insertHeaderBitsData(const int headerBits, const int poc)
{
    list<CodedPicture> &previousCodedPictures = m_dataStorageEngine->getPreviousCodedPictures();
    list<CodedPicture>::iterator it;
    for ( it = previousCodedPictures.begin(); it != previousCodedPictures.end(); it++ )
    {
       if((*it).getPoc() == poc)
       {
           (*it).setHeaderBits(headerBits);
           return true;
       }
    }
    return false;
}
