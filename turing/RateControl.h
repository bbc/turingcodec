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
 * RateControl.h
 *
 *  Created on: 11 Feb 2016
 *      Author: Matteo Naccari
 *      Copyright: BBC R&D 2016
 */

#ifndef Included_RateControl_h
#define Included_RateControl_h

#include <iostream>
#include <list>
#include <array>
#include "HevcMath.h"
#include "EstimateIntraComplexity.h"

#define NON_VALID_LAMBDA  -1.0
#define NON_VALID_QP      -1000
#define INITIAL_ALPHA      3.2001
#define INITIAL_BETA      -1.367
#define ALPHA_UPDATE_STEP  0.1
#define BETA_UPDATE_STEP   0.05
#define ALPHA_MIN          0.05
#define ALPHA_MAX          500.0
#define BETA_MIN          -3.0
#define BETA_MAX          -0.1
#define LAMBDA_WEIGHT      0.5
#define MAX_LIST_SIZE      32
#define ALPHA_INTRA        6.7542
#define BETA_INTRA         1.7860
#define MAX_NUM_LEVELS     6
#define CTU_SMOOT_WINDOW   4
#define BISECTION_MAX_IT   40
#define SOP_ADAPTIVE       0
#define BETA_INTRA_MAD     1.2517

using namespace std;

struct EstimateIntraComplexity;

class CpbInfo
{
private:
    int m_cpbStatus;
    int m_cpbSize;
    int m_cpbBufferingRate;

public:
    CpbInfo() : m_cpbStatus(0),
    m_cpbSize(0),
    m_cpbBufferingRate(0) {}

    CpbInfo(int cpbSize, int targetRate, double frameRate, double initialFullness)
    {
        m_cpbSize          = cpbSize;
        m_cpbStatus        = static_cast<int>(static_cast<double>(m_cpbSize * initialFullness));
        m_cpbBufferingRate = static_cast<int>(static_cast<double>(targetRate) / frameRate);
    }
    void updateCpbStatus(int codingBits)
    {
        m_cpbStatus += m_cpbBufferingRate - codingBits;
    }

    void setCpbInfo(int cpbSize, int targetRate, double frameRate, double initialFullness)
    {
        m_cpbSize          = cpbSize;
        m_cpbStatus        = static_cast<int>(static_cast<double>(m_cpbSize * initialFullness));
        m_cpbBufferingRate = static_cast<int>(static_cast<double>(targetRate) / frameRate);
    }

    int getCpbStatus()
    {
        return m_cpbStatus;
    }

    int getCpbBufferRate()
    {
        return m_cpbBufferingRate;
    }

    int getCpbSize()
    {
        return m_cpbSize;
    }

    void adjustAllocatedBits(int &currentBitsPerPicture)
    {
        // Cpb correction
        int estimatedCpbFullness = m_cpbStatus + m_cpbBufferingRate;

        // Check whether the allocated bits will lead to cpb overflow
        int overflowLevel = m_cpbSize * 0.9;
        if(estimatedCpbFullness - currentBitsPerPicture > overflowLevel)
        {
            currentBitsPerPicture = estimatedCpbFullness - overflowLevel;
        }

        // Check whether the allocated bits will lead to cpb underflow
        estimatedCpbFullness -= m_cpbBufferingRate;
        int underflowLevel = m_cpbSize * 0.1;
        if(estimatedCpbFullness - currentBitsPerPicture < underflowLevel)
        {
            currentBitsPerPicture = estimatedCpbFullness - underflowLevel;
        }
    }
};

class CodedCtu
{
private:
    double m_alpha;
    double m_beta;

public:
    CodedCtu() : m_alpha(INITIAL_ALPHA), m_beta(INITIAL_BETA) {}
    void   setAlpha(double value) { m_alpha = value; }
    void   setBeta (double value) { m_beta  = value; }

    double getAlpha() { return m_alpha; }
    double getBeta()  { return m_beta;  }
};

class CodedPicture
{
private:
    int    m_codingBits;
    int    m_headerBits;
    double m_lambda;
    int    m_qp;
    double m_alpha;
    double m_beta;
    int    m_level;
    double m_costIntra;

public:
    CodedPicture() : m_codingBits(0),
    m_headerBits(0),
    m_lambda(NON_VALID_LAMBDA),
    m_qp(NON_VALID_QP),
    m_alpha(INITIAL_ALPHA),
    m_beta(INITIAL_BETA),
    m_level(-1),
    m_costIntra(0.0) {}

    void setCodingBits(int bits)   { m_codingBits = bits; }
    void setHeaderBits(int bits)   { m_headerBits = bits; }
    void setLambda(double lambda)  { m_lambda = lambda;   }
    void setQp(int qp)             { m_qp = qp;           }
    void setAlpha(double alpha)    { m_alpha = alpha;     }
    void setBeta(double  beta)     { m_beta  = beta;      }
    void setLevel(int level)       { m_level = level;     }
    void setCostIntra(double cost) { m_costIntra = cost;  }

    double getAlpha()              { return m_alpha;      }
    double getBeta()               { return m_beta;       }
    double getLambda()             { return m_lambda;     }
    int    getQp()                 { return m_qp;         }
    int    getCodingBits()         { return m_codingBits; }
    int    getHeaderBits()         { return m_headerBits; }
    int    getLevel()              { return m_level;      }
    double getCostIntra()          { return m_costIntra;  }
};

class DataStorage
{
private:
    list<CodedPicture> m_previousCodedPictures;
    array<CodedPicture, MAX_NUM_LEVELS> m_codedPicturesPerLevel;
    array<CodedCtu*, MAX_NUM_LEVELS> m_codedCtusPerLevel;

public:

    ~DataStorage()
    {
        while(m_previousCodedPictures.size() > 0)
        {
            m_previousCodedPictures.pop_front();
        }
        for(int level = 0; level < m_codedCtusPerLevel.size(); level++)
        {
            delete[] m_codedCtusPerLevel[level];
        }
    }

    list<CodedPicture>& getPreviousCodedPictures()
	{
        return m_previousCodedPictures;
	}

    CodedPicture& getPictureAtLevel(int level)
    {
        assert(level < m_codedPicturesPerLevel.size());
        return m_codedPicturesPerLevel[level];
    }

    CodedCtu *getCodedCtuAtLevel(int level)
    {
        assert(level < m_codedCtusPerLevel.size());
        return m_codedCtusPerLevel[level];
    }

    void addCodedPicture(int level);

    void initCtuStorage(int numberOfUnits);

    void resetPreviousCodedPicturesBuffer();

    void resetCodedPicturesAtLevelMemory();

    void resetCodedCtuAtLevelMemory(int totalCtusPerFrame);

};

class CtuController
{
private:
    int m_codedBits;
    int m_ctuQp;
    int m_targetBits;
    double m_ctuLambda;
    double m_ctuWeight;
    int m_numberOfPixels;
    double m_costIntra;
    int m_targetBitsLeft;
    bool m_isValid;       // non fully skipped CTU

public:

    CtuController() : m_codedBits(0),
    m_ctuQp(NON_VALID_QP),
    m_costIntra(0),
    m_ctuLambda(NON_VALID_LAMBDA),
    m_ctuWeight(1.0),
    m_targetBits(0),
    m_numberOfPixels(0),
    m_targetBitsLeft(0),
    m_isValid(false) {}
    void setCodedBits(int bits)        { m_codedBits      = bits;   }
    void setCtuQp(int qp)              { m_ctuQp          = qp;     }
    void setTargetBits(int bits)       { m_targetBits     = bits;   }
    void setTargetBitsLeft(int bits)   { m_targetBitsLeft = bits;   }
    void setCtuLambda(double lambda)   { m_ctuLambda      = lambda; }
    void setCtuWeight(double weight)   { m_ctuWeight      = weight; }
    void setNumberOfPixels(int pixels) { m_numberOfPixels = pixels; }
    void setCostIntra(double cost)     { m_costIntra      = cost;   }
    void setValidityFlag(bool flag)    { m_isValid        = flag;   }
    void updateValidityFlag(bool flag)
    {
        m_isValid       |= flag;
    }

    int     getCodedBits()             { return m_codedBits;        }
    int     getCtuQp()                 { return m_ctuQp;            }
    int     getTargetBits()            { return m_targetBits;       }
    int     getTargetBitsLeft()        { return m_targetBitsLeft;   }
    double  getCtuLambda()             { return m_ctuLambda;        }
    double  getCtuWeight()             { return m_ctuWeight;        }
    int     getNumberOfPixels()        { return m_numberOfPixels;   }
    double  getCostIntra()             { return m_costIntra;        }
    bool    getValidityFlag()          { return m_isValid;          }
};

class PictureController
{
private:
    int m_targetBits;
    int m_bitsLeft;
    int m_pixelsPerPicture;
    double m_averageBpp;
    double m_alphaUpdateStep;
    double m_betaUpdateStep;
    DataStorage * m_dataStorageAccess;

public:
    PictureController() :
    m_targetBits(0),
    m_bitsLeft(0),
    m_pixelsPerPicture(0),
    m_averageBpp(0.0),
    m_alphaUpdateStep(ALPHA_UPDATE_STEP),
    m_betaUpdateStep(BETA_UPDATE_STEP),
    m_dataStorageAccess(0) {}

    PictureController(DataStorage *storageAccess) :
    m_targetBits(0),
    m_bitsLeft(0),
    m_pixelsPerPicture(0),
    m_averageBpp(0.0),
    m_alphaUpdateStep(ALPHA_UPDATE_STEP),
    m_betaUpdateStep(BETA_UPDATE_STEP),
    m_dataStorageAccess(storageAccess) {}

    void updateLambdaModelIntra(double &alpha,
                                double &beta,
                                CodedPicture &picture,
                                int actualBits);

    void initPictureController(int targetBits,
                               int pixelsPerPicture,
                               double averageBpp);

    double estimateLambda(int level, bool isIntra);

    void updatePictureController(int bitsSpent,
                                 int currentLevel,
                                 bool isIntra,
                                 double &lastCodedLambda);

    void setBitsIntra(int bits,
                      int pixelsCount);
    int    getBitsLeft()            { return m_bitsLeft;        }
    double getAlphaUpdateStep()     { return m_alphaUpdateStep; }
    double getBetaUpdateStep()      { return m_betaUpdateStep;  }
    void   updateBitsLeft(int bits) { m_bitsLeft -= bits;       }
    int    getPictureTargetBits()   { return m_targetBits;      }
};

class SOPController
{
private:
    int  m_size;
    int  m_targetBits;
    int  m_bitsLeft;
    int  m_framesLeft;
    int *m_weight;
    double m_averageBpp;
    DataStorage *m_dataStorageAccess;

    int getEstimatedHeaderBits(int level);

public:
    SOPController() : m_size(0),
    m_targetBits(0),
    m_bitsLeft(0),
    m_framesLeft(0),
    m_weight(0),
    m_averageBpp(0.0),
    m_dataStorageAccess(0) {}

    SOPController(DataStorage *dataAccess) :
        m_size(0),
        m_targetBits(0),
        m_bitsLeft(0),
        m_framesLeft(0),
        m_weight(0),
        m_averageBpp(0.0),
        m_dataStorageAccess(dataAccess) {}
    ~SOPController()
    {
        if(m_weight)
            delete[] m_weight;
    }

    void initSOPController(int size,
                           int targetBits,
                           int *weight);

    void updateSopController(int bitsSpent);

    int allocateRateCurrentPicture(int sequenceLevelFramesLeft, int level);

    int getFramesLeft() { return m_framesLeft; }
};

class SequenceController
{
private:
    double  m_targetRate;
    int     m_smoothingWindow;
    int     m_totalFrames;
    int     m_framesLeft;
    double  m_frameRate;
    int     m_sopSize;
    double  m_averageRate;
    double  m_averageBpp;
    int     m_bitsPerPicture;
    int64_t m_bitsLeft;
    int64_t m_targetBits;
    int     m_pixelsPerPicture;
    int     m_sopWeight[7][4][8] = {{{30, 8, 0, 0, 0, 0, 0, 0}, {25, 7, 0, 0, 0, 0, 0, 0}, {20, 6, 0, 0, 0, 0, 0, 0}, {15, 5, 0, 0, 0, 0, 0, 0}},
                                    {{30, 8, 4, 0, 0, 0, 0, 0}, {25, 7, 4, 0, 0, 0, 0, 0}, {20, 6, 4, 0, 0, 0, 0, 0}, {15, 5, 4, 0, 0, 0, 0, 0}},
                                    {{30, 8, 4, 4, 0, 0, 0, 0}, {25, 7, 4, 4, 0, 0, 0, 0}, {20, 6, 4, 4, 0, 0, 0, 0}, {15, 5, 4, 4, 0, 0, 0, 0}},
                                    {{30, 8, 4, 4, 8, 0, 0, 0}, {25, 7, 4, 4, 7, 0, 0, 0}, {20, 6, 4, 4, 6, 0, 0, 0}, {15, 5, 4, 4, 5, 0, 0, 0}},
                                    {{30, 8, 4, 4, 8, 4, 0, 0}, {25, 7, 4, 4, 7, 4, 0, 0}, {20, 6, 4, 4, 6, 4, 0, 0}, {15, 5, 4, 4, 5, 4, 0, 0}},
                                    {{30, 8, 4, 1, 1, 4, 1, 0}, {25, 7, 4, 1, 1, 4, 1, 0}, {20, 6, 4, 1, 1, 4, 1, 0}, {15, 5, 4, 1, 1, 4, 1, 0}},
                                    {{30, 8, 4, 1, 1, 4, 1, 1}, {25, 7, 4, 1, 1, 4, 1, 1}, {20, 6, 4, 1, 1, 4, 1, 1}, {15, 5, 4, 1, 1, 4, 1, 1}}}; // As in HM

    int   *m_frameWeight;
    int    m_totalLevels;
    int    m_baseQp;
    SOPController     *m_sopControllerEngine;
    DataStorage       *m_dataStorageEngine;
    PictureController *m_pictureControllerEngine;
    CtuController     *m_ctuControllerEngine;
    CpbInfo            m_cpbControllerEngine;
    double m_lastCodedPictureLambda;
    int    m_currentCodingBits;
    int    m_picSizeInCtbsY;
    int    m_picHeightInCtbs;
    int    m_picWidhtInCtbs;
    int    m_ctusLeft;
    int    m_ctusCoded;
    double m_remainingCostIntra;

    void   resetCtuController();
#if SOP_ADAPTIVE
    void   computeEquationCoefficients(double *coefficient, double *exponent, double *lambdaRatio);
    double solveEquationWithBisection (double *coefficient, double *exponent, double targetBpp);
#endif

public:
    SequenceController() : m_targetRate(0),
    m_smoothingWindow(40),
    m_framesLeft(0),
    m_frameRate(0),
    m_sopSize(8),
    m_averageRate(0),
    m_targetBits(0),
    m_bitsPerPicture(0),
    m_bitsLeft(0),
    m_pixelsPerPicture(0),
    m_frameWeight(0),
    m_totalFrames(0),
    m_totalLevels(0),
    m_baseQp(0),
    m_lastCodedPictureLambda(0.0),
    m_averageBpp(0.0),
    m_currentCodingBits(0),
    m_dataStorageEngine(0),
    m_sopControllerEngine(0),
    m_pictureControllerEngine(0),
    m_ctuControllerEngine(0),
    m_picSizeInCtbsY(0),
    m_ctusCoded(0),
    m_ctusLeft(0),
    m_picHeightInCtbs(0),
    m_picWidhtInCtbs(0),
    m_remainingCostIntra(0.0)
    {
    }

    SequenceController(double targetRate,
                       int totalFrames,
                       double frameRate,
                       int sopSize,
                       int picHeight,
                       int picWidth,
                       int ctuSize,
                       int totalLevels,
                       int baseQp);

    ~SequenceController()
    {
        if(m_dataStorageEngine)
            delete m_dataStorageEngine;
        if(m_sopControllerEngine)
            delete m_sopControllerEngine;
        if(m_pictureControllerEngine)
            delete m_pictureControllerEngine;
        if(m_ctuControllerEngine)
            delete[] m_ctuControllerEngine;
    }

    void initNewSop();

    void pictureRateAllocation(int currentPictureLevel);

    void updateSequenceController(int bitsSpent, int qp, double lambda, bool isIntra, int sopLevel);

    int getBaseQp() { return m_baseQp; }

    double estimatePictureLambda(bool isIntra, int sopLevel);

    int deriveQpFromLambda(double lambda, bool isIntra, int sopLevel);

    void setCodingBits(int codingBits)
    {
        m_currentCodingBits = codingBits;
    }

    void pictureRateAllocationIntra(EstimateIntraComplexity &icInfo);

    void setHeaderBits(int bits, bool isIntra, int sopLevel);

    double getCtuTargetBits(bool isIntraSlice, int ctbAddrInRs);

    double getCtuEstLambda(double bpp, int ctbAddrInRs, int currentPictureLevel);
    void   getCtuEstLambdaAndQp(double bpp, int sliceQp, int ctbAddrInRs, double &lambda, int &qp);

    int    getCtuEstQp(int ctbAddrInRs, int currentPictureLevel);

    void   updateCtuController(int codingBitsCtu, bool isIntra, int ctbAddrInRs, int currentPictureLevel);

    void   updateValidityFlag(bool flag, int ctuIdx)
    {
        m_ctuControllerEngine[ctuIdx].updateValidityFlag(flag);
    }

    void setValidityFlag(bool flag, int ctuIdx)
    {
        m_ctuControllerEngine[ctuIdx].setValidityFlag(flag);
    }

    int getCtuStoredQp(int ctuIdx)
    {
        return m_ctuControllerEngine[ctuIdx].getCtuQp();
    }

    void getAveragePictureQpAndLambda(int &averageQp, double &averageLambda);

    void setSopSize(int size)
    {
        m_sopSize = size;

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
    }

    void reset();

    void initCpbInfo(int cpbMaxSize)
    {
        m_cpbControllerEngine.setCpbInfo(cpbMaxSize, (int)m_targetRate, m_frameRate, 0.8);
    }
};


#endif /* RATECONTROL_H_ */
