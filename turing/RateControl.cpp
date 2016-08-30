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

PictureController::PictureController(DataStorage *storageAccess, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctuSize, int targetBits, double averageBpp) :
    m_targetBits(targetBits),
    m_pixelsPerPicture(pictureHeight*pictureWidth),
    m_averageBpp(averageBpp),
    m_dataStorageAccess(storageAccess),
    m_pictureHeightInCtbs(pictureHeightInCtbs),
    m_pictureWidthInCtbs(pictureWidthInCtbs),
    m_pictureSizeInCtbs(pictureHeightInCtbs*pictureWidthInCtbs),
    m_ctusCoded(0),
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

PictureController::PictureController(DataStorage *storageAccess, EstimateIntraComplexity &icInfo, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctuSize, int targetBits, double averageBpp) :
    m_targetBits(targetBits),
    m_pixelsPerPicture(pictureHeight*pictureWidth),
    m_averageBpp(averageBpp),
    m_dataStorageAccess(storageAccess),
    m_pictureHeightInCtbs(pictureHeightInCtbs),
    m_pictureWidthInCtbs(pictureWidthInCtbs),
    m_pictureSizeInCtbs(pictureHeightInCtbs*pictureWidthInCtbs),
    m_ctusCoded(0),
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

double PictureController::getCtuTargetBits(bool isIntraSlice, int ctbAddrInRs)
{
    int avgBits = 0;
    assert(ctbAddrInRs < m_pictureSizeInCtbs);

    int cumulativeMad, cumulativeBitsEstimated, cumulativeBitsSpent, cumulativeCtusCoded;
    double cumulativeWeight;
    getCodedInfoFromWavefront(ctbAddrInRs, cumulativeMad, cumulativeBitsEstimated, cumulativeBitsSpent, cumulativeCtusCoded, cumulativeWeight);
    int currentCtusLeft = m_pictureSizeInCtbs - cumulativeCtusCoded;
    int influence = std::min<int>(CTU_SMOOTH_WINDOW, currentCtusLeft);

    if(isIntraSlice)
    {
        double MAD        = m_ctuControllerEngine[ctbAddrInRs].getCostIntra();
        int remainingCostIntra = m_totalCostIntra - cumulativeMad;

        if (remainingCostIntra > 0.1 )
        {
            int estimateBitsCtu = static_cast<int>((MAD / static_cast<double>(cumulativeMad + MAD)) * static_cast<double>(m_targetBits));
            avgBits = estimateBitsCtu + (cumulativeBitsEstimated - cumulativeBitsSpent) / influence;
        }
        else
        {
            avgBits = static_cast<int>( (m_targetBits - cumulativeBitsSpent) / currentCtusLeft );
        }
    }
    else
    {
        const double currentCtuWeight = m_ctuControllerEngine[ctbAddrInRs].getCtuWeight();
        const double totalWeight = cumulativeWeight + currentCtuWeight;
        avgBits = static_cast<int>(currentCtuWeight + (totalWeight - static_cast<double>(cumulativeBitsSpent))/static_cast<double>(influence));
    }

    if(avgBits < 1)
    {
        avgBits = 1;
    }

    const double bpp = (double)avgBits / m_ctuControllerEngine[ctbAddrInRs].getNumberOfPixels();
    m_ctuControllerEngine[ctbAddrInRs].setTargetBits(avgBits);
    return bpp;
}

void PictureController::getCodedInfoFromWavefront(int ctbAddrInRs, int &cumulativeMad, int &cumulativeBitsEstimated, int &cumulativeBitsSpent, int &cumulativeCtusCoded, double &cumulativeWeight)
{

    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    int rowStart = ctbAddrInRs / m_pictureWidthInCtbs;
    int colStart = (ctbAddrInRs % m_pictureWidthInCtbs) - 1;
    cumulativeMad = cumulativeBitsEstimated = cumulativeBitsSpent = cumulativeCtusCoded = 0;
    cumulativeWeight = 0.0;

    for(int r = rowStart; r >= 0; r--)
    {
        for(int c = colStart; c >= 0; c--)
        {
            int currentIdx = r*m_pictureWidthInCtbs + c;
            assert(0 <= currentIdx && currentIdx < m_pictureSizeInCtbs);
            assert(m_ctuControllerEngine[currentIdx].getFinishedFlag());
            cumulativeMad           += static_cast<int>(m_ctuControllerEngine[currentIdx].getCostIntra());
            cumulativeBitsEstimated += m_ctuControllerEngine[currentIdx].getTargetBitsEst();
            cumulativeBitsSpent     += m_ctuControllerEngine[currentIdx].getCodedBits();
            cumulativeWeight        += m_ctuControllerEngine[currentIdx].getCtuWeight();
            cumulativeCtusCoded++;
        }
        colStart = min<int>(m_pictureWidthInCtbs - 1, colStart + 1);
    }

}

double PictureController::getCtuEstLambda(double bpp, int ctbAddrInRs, int pictureLevel)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    CodedCtu *codedCtuAtLevel = m_dataStorageAccess->getCodedCtuAtLevel(pictureLevel);
    const double alpha = codedCtuAtLevel[ctbAddrInRs].getAlpha();
    const double beta  = codedCtuAtLevel[ctbAddrInRs].getBeta();

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

    m_ctuControllerEngine[ctbAddrInRs].setCtuLambda(estLambdaCtu);
    m_ctuControllerEngine[ctbAddrInRs].setCtuReciprocalLambda(1.0 / estLambdaCtu);
    m_ctuControllerEngine[ctbAddrInRs].setCtuReciprocalSqrtLambda(1.0 / sqrt(estLambdaCtu));
    return estLambdaCtu;

}

int PictureController::getCtuEstQp(int ctbAddrInRs, int pictureLevel)
{
    assert(ctbAddrInRs < m_pictureSizeInCtbs);

    const double ctuLambda = m_ctuControllerEngine[ctbAddrInRs].getCtuLambda();
    int ctuEstQp = (int)( 4.2005 * log( ctuLambda ) + 13.7122 + 0.5 );

    CodedPicture &currentPicture = m_dataStorageAccess->getPictureAtLevel(pictureLevel);
    int pictureQp = currentPicture.getQp();

    ctuEstQp = Clip3( pictureQp - 2, pictureQp + 2, ctuEstQp );

    ctuEstQp = Clip3(0, 51, ctuEstQp);

    m_ctuControllerEngine[ctbAddrInRs].setCtuQp(ctuEstQp);

    return ctuEstQp;
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

void PictureController::getCtuLambdaAndQp(double bpp, int sliceQp, int ctbAddrInRs, double &lambda, int &qp)
{

    assert(ctbAddrInRs < m_pictureSizeInCtbs);
    int currentLevel = 0; // only called for intra slices!
    CodedPicture &pictureAtLevel = m_dataStorageAccess->getPictureAtLevel(currentLevel);

    double   alpha = pictureAtLevel.getAlpha();
    double   beta  = pictureAtLevel.getBeta();

    double costPixelBased = m_ctuControllerEngine[ctbAddrInRs].getCostIntra() /
            (double)m_ctuControllerEngine[ctbAddrInRs].getNumberOfPixels();
    costPixelBased = pow(costPixelBased, BETA_INTRA_MAD);
    lambda = (alpha/256.0) * pow( costPixelBased/bpp, beta );

    const int minQP = sliceQp - 2;
    const int maxQP = sliceQp + 2;

    double maxLambda = exp(((double)(maxQP+0.49)-13.7122)/4.2005);
    double minLambda = exp(((double)(minQP-0.49)-13.7122)/4.2005);

    lambda = Clip3(minLambda, maxLambda, lambda);

    qp = int( 4.2005 * log(lambda) + 13.7122 + 0.5 );
    qp = Clip3(minQP, maxQP, qp);
    qp = Clip3(0, 51, qp);
    m_ctuControllerEngine[ctbAddrInRs].setCtuQp(qp);
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

void SOPController::initSOPController(int size,
                                      int targetBits,
                                      int *weight)
{
    if(m_weight == 0)
    {
        m_weight = new int[size];
    }
    else
    {
        // Reallocate the memory in case SOP changes its size
        delete[] m_weight;
        m_weight = new int[size];
    }
    m_size = size;

    m_targetBits = targetBits;

    int ratioSum = 0;
    for(int i = 0; i < m_size; i++)
    {
        ratioSum += weight[i];
    }

    for(int i = 0; i < m_size; i++)
    {
        m_weight[i] = (int)(((double)targetBits * weight[i])/ratioSum);
    }

    m_bitsLeft = targetBits;
    m_framesLeft = m_size;
}

void SOPController::updateSopController(int bitsSpent)
{
    m_bitsLeft -= bitsSpent;
    m_framesLeft--;
}

int SOPController::allocateRateCurrentPicture(int level)
{
    int currentPicPosition = m_size - m_framesLeft;
    int currentPicRatio    = m_weight[currentPicPosition];
    int sumPicRatio        = 0;
    int headerBits         = getEstimatedHeaderBits(level);

    for(int picIdx = currentPicPosition; picIdx < m_size; picIdx++)
    {
        sumPicRatio += m_weight[picIdx];
    }

    int picTargetBits = (int)((double)m_bitsLeft * (double)currentPicRatio / (double)sumPicRatio);

    if(picTargetBits < 100)
    {
        picTargetBits = 100; // Same as in HM
    }

    picTargetBits = static_cast<int>(0.1 * (double)picTargetBits + 0.9 * (double)m_weight[currentPicPosition]);

    if( picTargetBits < headerBits + 100 )
    {
        picTargetBits = headerBits + 100;
    }

    picTargetBits -= headerBits;

    return picTargetBits;
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
                                       int totalLevels,
                                       int baseQp)
{
    m_targetRate       = targetRate * 1000; // kbps to bps conversion
    m_smoothingWindow  = static_cast<int>(frameRate);
    int totalFrames    = intraPeriod != 1 ? intraPeriod : ((static_cast<int>(frameRate) + 4)/8)*8;
    m_totalFrames      = totalFrames;
    m_frameRate        = frameRate;
    m_sopSize          = sopSize;

    m_averageRate      = m_targetRate / m_frameRate;
    m_targetBits       = static_cast<int64_t>(m_averageRate * totalFrames);
    m_bitsPerPicture   = (int)(m_targetBits / (double)totalFrames);
    m_bitsLeft         = m_targetBits;
    m_pixelsPerPicture = picHeight * picWidth;
    m_totalLevels      = totalLevels;
    m_baseQp           = baseQp;

    //	Compute the bit per pixel
    m_averageBpp = m_averageRate / (double)m_pixelsPerPicture;

    // Initialise the frame-based weights according to the value of bpp
    if(m_averageBpp <= 0.05)
    {
        m_frameWeight = m_sopWeight[m_sopSize-2][0];
    }
    else if(0.05 < m_averageBpp && m_averageBpp <= 0.1)
    {
        m_frameWeight = m_sopWeight[m_sopSize-2][1];
    }
    else if(0.1 < m_averageBpp && m_averageBpp <= 0.2)
    {
        m_frameWeight = m_sopWeight[m_sopSize-2][2];
    }
    else
        m_frameWeight = m_sopWeight[m_sopSize-2][3];

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

    m_sopControllerEngine     = new SOPController(m_dataStorageEngine);
    m_dataStorageEngine->initCtuStorage(m_picSizeInCtbsY);

    CodedPicture &pictureIntra = m_dataStorageEngine->getPictureAtLevel(0);
    pictureIntra.setAlpha(ALPHA_INTRA);
    pictureIntra.setBeta(BETA_INTRA);

    m_codedFrames = 0;
    m_bitsSpent = 0;

#if WRITE_RC_LOG
    m_logFile.open("rcLogFileCpb.txt", ios::out);
    m_logFile<<"---------------------------------------------------------------------------------\n";
    m_logFile<<"|    POC |     Target |    Lambda |   QP | Coded bits |      Total | CPB Status |\n";
    m_logFile<<"---------------------------------------------------------------------------------\n";
    m_logFile.flush();
#endif
}

void SequenceController::initNewSop()
{
    int realInfluencePicture  = m_smoothingWindow;
    int currentBitsPerPicture = (m_bitsPerPicture * (m_codedFrames + m_smoothingWindow) - m_bitsSpent) / realInfluencePicture;
    int currentBitsSop = currentBitsPerPicture * m_sopSize;

    if(currentBitsSop < 200)
    {
        currentBitsSop = 200; // As in HM RC
    }

#if SOP_ADAPTIVE
    if(m_lastCodedPictureLambda < 0.1)
    {
        m_sopControllerEngine->initSOPController(m_sopSize, currentBitsSop, m_frameWeight);
    }
    else
    {
        // Adaptive computation of frame weights
        double *lambdaRatio = new double[m_sopSize];
        double *coefficient = new double[m_sopSize];
        double *exponent    = new double[m_sopSize];
        int    *frameWeight = new int[m_sopSize];
        double targetBpp    = (double)currentBitsSop / (double)m_pixelsPerPicture;

        if ( m_lastCodedPictureLambda < 90.0 )
        {
            double lambdaP4 = 0.725 * log( m_lastCodedPictureLambda ) + 0.7963;
            double sopLambdaRatio[7][8] = {{1.0, lambdaP4, 1.0, 1.0,  1.0     , 1.0, 1.0,   1.0},
                                           {1.0, lambdaP4, 1.3, 1.0,  1.0     , 1.0, 1.0,   1.0},
                                           {1.0, lambdaP4, 1.3, 1.3,  1.0     , 1.0, 1.0,   1.0},
                                           {1.0, lambdaP4, 1.3, 1.3,  lambdaP4, 1.0, 1.0,   1.0},
                                           {1.0, lambdaP4, 1.3, 1.3,  lambdaP4, 1.3, 1.0,   1.0},
                                           {1.0, lambdaP4, 1.3, 3.25, 3.25    , 1.3, 3.25,  1.0},
                                           {1.0, lambdaP4, 1.3, 3.25, 3.25    , 1.3, 3.25, 3.25}};
            for(int idx = 0; idx < m_sopSize; idx++)
            {
                lambdaRatio[idx] = sopLambdaRatio[m_sopSize-2][idx];
                if(idx > 1)
                {
                    lambdaRatio[idx] *= lambdaRatio[1];
                }
            }
        }
        else
        {
            double sopLambdaRatio[7][8] = {{1.0, 4.0, 1.0,  1.0,  1.0, 1.0,  1.0,  1.0},
                                           {1.0, 4.0, 5.0,  1.0,  1.0, 1.0,  1.0,  1.0},
                                           {1.0, 4.0, 5.0,  5.0,  1.0, 1.0,  1.0,  1.0},
                                           {1.0, 4.0, 5.0,  5.0,  4.0, 1.0,  1.0,  1.0},
                                           {1.0, 4.0, 5.0,  5.0,  4.0, 5.0,  1.0,  1.0},
                                           {1.0, 4.0, 5.0, 12.3, 12.3, 5.0, 12.3,  1.0},
                                           {1.0, 4.0, 5.0, 12.3, 12.3, 5.0, 12.3, 12.3}};
            for(int idx = 0; idx < m_sopSize; idx++)
            {
                lambdaRatio[idx] = sopLambdaRatio[m_sopSize-2][idx];
            }
        }
        computeEquationCoefficients(coefficient, exponent, lambdaRatio);

        double adaptiveLambda = solveEquationWithBisection(coefficient, exponent, targetBpp);

        for(int sopIdx = 0; sopIdx < m_sopSize; sopIdx++)
        {
            frameWeight[sopIdx] = (int)(coefficient[sopIdx] * pow(adaptiveLambda, exponent[sopIdx]) * m_pixelsPerPicture);
        }
        m_sopControllerEngine->initSOPController(m_sopSize, currentBitsSop, frameWeight);

        delete[] lambdaRatio;
        delete[] coefficient;
        delete[] exponent;
        delete[] frameWeight;
    }
#else
    m_sopControllerEngine->initSOPController(m_sopSize, currentBitsSop, m_frameWeight);
#endif

}
void SequenceController::pictureRateAllocation(int currentPictureLevel, int poc)
{
    CodedPicture &codedPictureAtLevel = m_dataStorageEngine->getPictureAtLevel(currentPictureLevel);
    codedPictureAtLevel.setLevel(currentPictureLevel);
    int picTargetBits = m_sopControllerEngine->allocateRateCurrentPicture(currentPictureLevel);

    // Cpb correction
    m_cpbControllerEngine.adjustAllocatedBits(picTargetBits);

    assert(m_pictureControllerEngine.find(poc) == m_pictureControllerEngine.end());
    m_pictureControllerEngine[poc] = new PictureController(m_dataStorageEngine, m_picHeight, m_picWidth, m_picHeightInCtbs, m_picWidthInCtbs, m_ctuSize, picTargetBits, m_averageBpp);

}

void SequenceController::updateSequenceController(bool isIntra, int sopLevel, int poc)
{
    // Update the bit budget and frames to be encoded at sequence and SOP levels
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());
    const int codingBits = currentPictureController->second->getTotalCodingBits();
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

    if(!isIntra)
    {
        m_sopControllerEngine->updateSopController(codingBits);
    }
    else
    {
        m_sopControllerEngine->updateSopController(currentPictureController->second->getPictureTargetBits());
    }

    m_dataStorageEngine->addCodedPicture(currentPictureLevel);

    currentPictureController->second->updatePictureController(currentPictureLevel,
                                                              isIntra,
                                                              m_lastCodedPictureLambda);

    // Update Cpb status
    m_cpbControllerEngine.updateCpbStatus(codingBits);

    // Delete the picture controller associated with this POC
    assert(m_pictureControllerEngine.erase(poc) == 1);
}

double SequenceController::estimatePictureLambda(bool isIntra, int sopLevel, int poc)
{
    int currentPictureLevel = isIntra ? 0 : sopLevel;

    double estLambda;
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    estLambda = currentPictureController->second->estimateLambda(currentPictureLevel, isIntra);

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

    assert(m_pictureControllerEngine.find(poc) == m_pictureControllerEngine.end());
    m_pictureControllerEngine[poc] = new PictureController(m_dataStorageEngine, icInfo, m_picHeight, m_picWidth, m_picHeightInCtbs, m_picWidthInCtbs, m_ctuSize, currentBitsPerPicture, m_averageBpp);


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

double SequenceController::getCtuTargetBits(bool isIntraSlice, int ctbAddrInRs, int poc)
{
    int avgBits = 0;
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    const double bpp = currentPictureController->second->getCtuTargetBits(isIntraSlice, ctbAddrInRs);
    return bpp;
}

double SequenceController::getCtuEstLambda(double bpp, int ctbAddrInRs, int currentPictureLevel, int poc)
{
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());
    double estLambdaCtu = currentPictureController->second->getCtuEstLambda(bpp, ctbAddrInRs, currentPictureLevel);

    return estLambdaCtu;
}

int SequenceController::getCtuEstQp(int ctbAddrInRs, int currentPictureLevel, int poc)
{
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    int ctuEstQp = currentPictureController->second->getCtuEstQp(ctbAddrInRs, currentPictureLevel);
    return ctuEstQp;
}

void SequenceController::updateCtuController(int codingBitsCtu, bool isIntra, int ctbAddrInRs, int currentPictureLevel, int poc)
{
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    currentPictureController->second->updateCtuController(codingBitsCtu, isIntra, ctbAddrInRs, currentPictureLevel);

}

void SequenceController::getAveragePictureQpAndLambda(int &averageQp, double &averageLambda, int poc)
{
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());

    currentPictureController->second->getAveragePictureQpAndLambda(averageQp, averageLambda);
}

#if SOP_ADAPTIVE
void SequenceController::computeEquationCoefficients(double *coefficient, double *exponent, double *lambdaRatio)
{
    int    sopPosition2Level[7][8] = {{1, 2, 0, 0, 0, 0, 0, 0},  //[SOP size][SOP position]
                                        {1, 2, 3, 0, 0, 0, 0, 0},
                                        {1, 2, 3, 3, 0, 0, 0, 0},
                                        {1, 2, 3, 3, 2, 0, 0, 0},
                                        {1, 2, 3, 3, 2, 3, 0, 0},
                                        {1, 2, 3, 4, 4, 3, 4, 0},
                                        {1, 2, 3, 4, 4, 3, 4, 4}};
    for ( int sopIdx = 0; sopIdx < m_sopSize; sopIdx++ )
    {
        int level    = sopPosition2Level[m_sopSize-2][sopIdx];
        CodedPicture &pictureAtLevel = m_dataStorageEngine->getPictureAtLevel(level);
        double alpha = pictureAtLevel.getAlpha();
        double beta  = pictureAtLevel.getBeta();
        coefficient[sopIdx] = pow( 1.0/alpha, 1.0/beta ) * pow( lambdaRatio[sopIdx], 1.0/beta );
        exponent[sopIdx]    = 1.0/beta;
    }
}

double SequenceController::solveEquationWithBisection(double *coefficient, double *exponent, double targetBpp)
{
    double lambdaSolution = 100.0;
    double m = 0.1;
    double M = 10000.0;
    for ( int iteration = 0; iteration < BISECTION_MAX_IT; iteration++ )
    {
        double currentSopBpp = 0.0;
        for ( int sopIdx = 0; sopIdx < m_sopSize; sopIdx++ )
        {
            currentSopBpp += coefficient[sopIdx] * pow( lambdaSolution, exponent[sopIdx] );
        }

        if ( fabs( currentSopBpp - targetBpp ) < 0.000001 )
        {
            break;
        }

        if ( currentSopBpp > targetBpp )
        {
            m = lambdaSolution;
            lambdaSolution = ( lambdaSolution + M ) / 2.0;
        }
        else
        {
            M = lambdaSolution;
            lambdaSolution = ( lambdaSolution + m ) / 2.0;
        }
    }

    lambdaSolution = Clip3( 0.1, 10000.0, lambdaSolution );
    return lambdaSolution;
}
#endif

void SequenceController::getCtuEstLambdaAndQp(double bpp, int sliceQp, int ctbAddrInRs, double &lambda, int &qp, int poc)
{
    auto currentPictureController = m_pictureControllerEngine.find(poc);
    assert(currentPictureController != m_pictureControllerEngine.end());
    currentPictureController->second->getCtuLambdaAndQp(bpp, sliceQp, ctbAddrInRs, lambda, qp);
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
