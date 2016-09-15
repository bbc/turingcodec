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
#include <fstream>
#include <list>
#include <array>
#include <map>
#include <mutex>
#include "HevcMath.h"
#include "EstimateIntraComplexity.h"
#include "InputQueue.h"

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
#define CTU_SMOOTH_WINDOW  4
#define BISECTION_MAX_IT   40
#define BETA_INTRA_MAD     1.2517
#define WRITE_RC_LOG       1

using namespace std;

struct EstimateIntraComplexity;

class CpbInfo
{
private:
    int m_cpbStatus;
    int m_cpbSize;
    int m_cpbBufferingRate;
    mutex m_cpbStatusMutex;

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
        unique_lock<mutex> lock(m_cpbStatusMutex);
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
        unique_lock<mutex> lock(m_cpbStatusMutex);
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
        unique_lock<mutex> lock(m_cpbStatusMutex);

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
            currentBitsPerPicture = max<int>(200, estimatedCpbFullness - underflowLevel);
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
    int    m_poc;

public:
    CodedPicture() : m_codingBits(0),
    m_headerBits(0),
    m_lambda(NON_VALID_LAMBDA),
    m_qp(NON_VALID_QP),
    m_alpha(INITIAL_ALPHA),
    m_beta(INITIAL_BETA),
    m_level(-1),
    m_costIntra(0.0),
    m_poc(-1) {}

    void setCodingBits(int bits)   { m_codingBits = bits; }
    void setHeaderBits(int bits)   { m_headerBits = bits; }
    void setLambda(double lambda)  { m_lambda = lambda;   }
    void setQp(int qp)             { m_qp = qp;           }
    void setAlpha(double alpha)    { m_alpha = alpha;     }
    void setBeta(double  beta)     { m_beta  = beta;      }
    void setLevel(int level)       { m_level = level;     }
    void setCostIntra(double cost) { m_costIntra = cost;  }
    void setPoc(int poc)           { m_poc = poc;         }

    double getAlpha()              { return m_alpha;      }
    double getBeta()               { return m_beta;       }
    double getLambda()             { return m_lambda;     }
    int    getQp()                 { return m_qp;         }
    int    getCodingBits()         { return m_codingBits; }
    int    getHeaderBits()         { return m_headerBits; }
    int    getLevel()              { return m_level;      }
    double getCostIntra()          { return m_costIntra;  }
    int    getPoc()                { return m_poc;        }
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
    int m_targetBits;
    int m_targetBitsEst;
    int m_ctuQp;
    double m_ctuLambda;
    double m_ctuReciprocalLambda;
    double m_ctuReciprocalSqrtLambda;
    double m_ctuWeight;
    double m_costIntra;
    double m_bpp;
    int m_numberOfPixels;
    bool m_isValid;       // non fully skipped CTU
    bool m_finished;

public:

    CtuController() : m_codedBits(0),
    m_ctuQp(NON_VALID_QP),
    m_costIntra(0),
    m_ctuLambda(NON_VALID_LAMBDA),
    m_ctuWeight(1.0),
    m_targetBits(0),
    m_numberOfPixels(0),
    m_isValid(false),
    m_targetBitsEst(0),
    m_finished(false),
    m_ctuReciprocalLambda(NON_VALID_LAMBDA),
    m_ctuReciprocalSqrtLambda(NON_VALID_LAMBDA),
    m_bpp(0.0) {}

    void setCodedBits(int bits)                   { m_codedBits      = bits;           }
    void setCtuQp(int qp)                         { m_ctuQp          = qp;             }
    void setTargetBits(int bits)                  { m_targetBits     = bits;           }
    void setTargetBitsEst(int bits)               { m_targetBitsEst  = bits;           }
    void setCtuLambda(double lambda)              { m_ctuLambda      = lambda;         }
    void setCtuReciprocalLambda(double value)     { m_ctuReciprocalLambda = value;     }
    void setCtuReciprocalSqrtLambda(double value) { m_ctuReciprocalSqrtLambda = value; }
    void setCtuWeight(double weight)              { m_ctuWeight      = weight;         }
    void setNumberOfPixels(int pixels)            { m_numberOfPixels = pixels;         }
    void setCostIntra(double cost)                { m_costIntra      = cost;           }
    void setValidityFlag(bool flag)               { m_isValid        = flag;           }
    void updateValidityFlag(bool flag)
    {
        m_isValid       |= flag;
    }
    void setFinishedFlag(bool flag)               { m_finished       = flag;           }
    void setCtuBpp(double bpp)                    { m_bpp            = bpp;            }

    int     getCodedBits()               { return m_codedBits;               }
    int     getCtuQp()                   { return m_ctuQp;                   }
    int     getTargetBits()              { return m_targetBits;              }
    int     getTargetBitsEst()           { return m_targetBitsEst;           }
    double  getCtuLambda()               { return m_ctuLambda;               }
    double  getCtuReciprocalLambda()     { return m_ctuReciprocalLambda;     }
    double  getCtuReciprocalSqrtLambda() { return m_ctuReciprocalSqrtLambda; }
    double  getCtuWeight()               { return m_ctuWeight;               }
    int     getNumberOfPixels()          { return m_numberOfPixels;          }
    double  getCostIntra()               { return m_costIntra;               }
    bool    getValidityFlag()            { return m_isValid;                 }
    bool    getFinishedFlag()            { return m_finished;                }
    double  getCtuBpp()                  { return m_bpp;                     }
};

class PictureController
{
private:
    int m_targetBits;
    int m_pixelsPerPicture;
    int m_pictureHeightInCtbs;
    int m_pictureWidthInCtbs;
    int m_pictureSizeInCtbs;
    int m_ctusLeft;
    int m_ctusCoded;
    int m_totalCostIntra;
    double m_averageBpp;
    double m_alphaUpdateStep;
    double m_betaUpdateStep;
    DataStorage   *m_dataStorageAccess;
    CtuController *m_ctuControllerEngine;

    void getCodedInfoFromWavefront(int ctbAddrInRs, int &cumulativeMad, int &cumulativeBitsEstimated, int &cumulativeBitsSpent, int &cumulativeCtusCoded, double &cumulativeWeigth);
    int  estimateCtuLambdaAndQpIntra(int sliceQp, int ctbAddrInRs);
    int  estimateCtuLambdaAndQpInter(int ctbAddrInRs, int pictureLevel);
    void updateLambdaModelIntra(double &alpha,
                                double &beta,
                                CodedPicture &picture,
                                int actualBits);

public:
    PictureController() :
    m_targetBits(0),
    m_pixelsPerPicture(0),
    m_averageBpp(0.0),
    m_alphaUpdateStep(ALPHA_UPDATE_STEP),
    m_betaUpdateStep(BETA_UPDATE_STEP),
    m_dataStorageAccess(0),
    m_ctuControllerEngine(0),
    m_pictureHeightInCtbs(0),
    m_pictureWidthInCtbs(0),
    m_pictureSizeInCtbs(0),
    m_ctusCoded(0),
    m_ctusLeft(0),
    m_totalCostIntra(0) {}

    PictureController(DataStorage *storageAccess, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctuSize, int targetBits, double averageBpp);

    PictureController(DataStorage *storageAccess, EstimateIntraComplexity &icInfo, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctuSize, int targetBits, double averageBpp);

    ~PictureController()
    {
        if(m_ctuControllerEngine)
            delete[] m_ctuControllerEngine;
    }

    double estimateLambda(int level, bool isIntra);

    void updatePictureController(int currentLevel,
                                 bool isIntra,
                                 double &lastCodedLambda);

    int    getPictureTargetBits()   { return m_targetBits;      }

    void   updateCtuValidityFlag(bool flag, int ctuIdx)
    {
        assert(ctuIdx < m_pictureSizeInCtbs);
        m_ctuControllerEngine[ctuIdx].updateValidityFlag(flag);
    }
    void   setCtuValidityFlag(bool flag, int ctuIdx)
    {
        assert(ctuIdx < m_pictureSizeInCtbs);
        m_ctuControllerEngine[ctuIdx].setValidityFlag(flag);
    }
    int    getCtuStoredQp(int ctuIdx)
    {
        assert(ctuIdx < m_pictureSizeInCtbs);
        return m_ctuControllerEngine[ctuIdx].getCtuQp();
    }
    void   computeCtuTargetBits(bool isIntraSlice, int ctbAddrInRs);
    int    estimateCtuLambdaAndQp(bool isIntra, int ctbAddrInRs, int pictureLevel, int sliceQp);
    void   updateCtuController(int codingBits, bool isIntra, int ctbAddrInRs, int pictureLevel);
    void   getAveragePictureQpAndLambda(int &averageQp, double &averageLambda);
    int    getFinishedCtus();
    double getCtuLambda(int ctbAddrInRs)
    {
        assert(ctbAddrInRs < m_pictureSizeInCtbs);
        return m_ctuControllerEngine[ctbAddrInRs].getCtuLambda();
    }
    double getCtuReciprocalLambda(int ctbAddrInRs)
    {
        assert(ctbAddrInRs < m_pictureSizeInCtbs);
        return m_ctuControllerEngine[ctbAddrInRs].getCtuReciprocalLambda();
    }
    double getCtuReciprocalSqrtLambda(int ctbAddrInRs)
    {
        assert(ctbAddrInRs < m_pictureSizeInCtbs);
        return m_ctuControllerEngine[ctbAddrInRs].getCtuReciprocalSqrtLambda();
    }
    int getTotalCodingBits();
};

class SOPController
{
private:
    int  m_size;
    int  m_targetBits;
    int  m_sopId;
    int  m_framesCoded;
    int  m_bitsCoded;
    int *m_weight;
    DataStorage *m_dataStorageAccess;

    int getEstimatedHeaderBits(int level);

public:
    SOPController() : m_size(0),
    m_targetBits(0),
    m_weight(0),
    m_dataStorageAccess(0),
    m_sopId(-1),
    m_framesCoded(0),
    m_bitsCoded(0) {}

    SOPController(DataStorage *dataAccess, int sopId, int size, int *bitrateWeight, int targetBits);

    ~SOPController()
    {
        if(m_weight)
        {
            delete[] m_weight;
        }
    }

    int getRateCurrentPicture(int pocInSop);

    void updateSopController(int bitsCoded);

    bool finished()
    {
        return m_framesCoded == m_size;
    }

};

class SequenceController
{
private:
    double  m_targetRate;
    int     m_smoothingWindow;
    int     m_totalFrames;
    int     m_codedFrames;
    int     m_bitsSpent;
    double  m_frameRate;
    int     m_sopSize;
    double  m_averageRate;
    double  m_averageBpp;
    int     m_bitsPerPicture;
    int64_t m_bitsLeft;
    int64_t m_targetBits;
    int     m_pixelsPerPicture;
//    int     m_sopWeight[7][4][8] = {{{30, 8, 0, 0, 0, 0, 0, 0}, {25, 7, 0, 0, 0, 0, 0, 0}, {20, 6, 0, 0, 0, 0, 0, 0}, {15, 5, 0, 0, 0, 0, 0, 0}},
//                                    {{30, 8, 4, 0, 0, 0, 0, 0}, {25, 7, 4, 0, 0, 0, 0, 0}, {20, 6, 4, 0, 0, 0, 0, 0}, {15, 5, 4, 0, 0, 0, 0, 0}},
//                                    {{30, 8, 4, 4, 0, 0, 0, 0}, {25, 7, 4, 4, 0, 0, 0, 0}, {20, 6, 4, 4, 0, 0, 0, 0}, {15, 5, 4, 4, 0, 0, 0, 0}},
//                                    {{30, 8, 4, 4, 8, 0, 0, 0}, {25, 7, 4, 4, 7, 0, 0, 0}, {20, 6, 4, 4, 6, 0, 0, 0}, {15, 5, 4, 4, 5, 0, 0, 0}},
//                                    {{30, 8, 4, 4, 8, 4, 0, 0}, {25, 7, 4, 4, 7, 4, 0, 0}, {20, 6, 4, 4, 6, 4, 0, 0}, {15, 5, 4, 4, 5, 4, 0, 0}},
//                                    {{30, 8, 4, 1, 1, 4, 1, 0}, {25, 7, 4, 1, 1, 4, 1, 0}, {20, 6, 4, 1, 1, 4, 1, 0}, {15, 5, 4, 1, 1, 4, 1, 0}},
//                                    {{30, 8, 4, 1, 1, 4, 1, 1}, {25, 7, 4, 1, 1, 4, 1, 1}, {20, 6, 4, 1, 1, 4, 1, 1}, {15, 5, 4, 1, 1, 4, 1, 1}}}; // As in HM
    int     m_sopWeight[7][4][8] = {{{8, 30, 0, 0, 0, 0, 0, 0}, {7, 25, 0, 0, 0, 0, 0, 0}, {6, 20, 0, 0, 0, 0, 0, 0}, {5, 15, 0, 0, 0, 0, 0, 0}},
                                    {{4, 8, 30, 0, 0, 0, 0, 0}, {4, 7, 25, 0, 0, 0, 0, 0}, {4, 6, 20, 0, 0, 0, 0, 0}, {4, 5, 15, 0, 0, 0, 0, 0}},
                                    {{4, 4, 8, 30, 0, 0, 0, 0}, {4, 4, 7, 25, 0, 0, 0, 0}, {4, 4, 6, 20, 0, 0, 0, 0}, {4, 4, 5, 15, 0, 0, 0, 0}},
                                    {{4, 4, 8, 8, 30, 0, 0, 0}, {4, 4, 7, 7, 25, 0, 0, 0}, {4, 4, 6, 6, 20, 0, 0, 0}, {4, 4, 5, 5, 15, 0, 0, 0}},
                                    {{4, 4, 8, 4, 8, 30, 0, 0}, {4, 4, 7, 4, 7, 25, 0, 0}, {4, 4, 6, 4, 6, 20, 0, 0}, {4, 4, 5, 6, 5, 15, 0, 0}},
                                    {{1, 4, 1, 8, 1, 4, 30, 0}, {1, 4, 1, 7, 1, 4, 25, 0}, {1, 4, 1, 6, 1, 4, 20, 0}, {1, 4, 1, 6, 1, 4, 15, 0}},
                                    {{1, 4, 1, 8, 1, 4, 1, 30}, {1, 4, 1, 7, 1, 4, 1, 25}, {1, 4, 1, 6, 1, 4, 1, 20}, {1, 4, 1, 5, 1, 4, 1, 15}}}; // As in HM

    int    m_baseQp;

    DataStorage                 *m_dataStorageEngine;
    map<int, PictureController*> m_pictureControllerEngine;
    map<int, SOPController*>     m_sopControllerEngine;
    CpbInfo                      m_cpbControllerEngine;

    mutex  m_pictureControllerMutex;
    mutex  m_sopControllerMutex;

    double m_lastCodedPictureLambda;
    int    m_picSizeInCtbsY;
    int    m_picHeightInCtbs;
    int    m_picWidthInCtbs;
    int    m_picHeight;
    int    m_picWidth;
    int    m_ctuSize;
    int    m_averageBitsPerCtb;
#if WRITE_RC_LOG
    ofstream m_logFile;
#endif

    bool insertHeaderBitsData(const int headerBits, const int poc);

public:
    SequenceController() : m_targetRate(0),
    m_smoothingWindow(40),
    m_frameRate(0),
    m_sopSize(8),
    m_averageRate(0),
    m_targetBits(0),
    m_bitsPerPicture(0),
    m_bitsLeft(0),
    m_pixelsPerPicture(0),
    m_totalFrames(0),
    m_baseQp(0),
    m_lastCodedPictureLambda(0.0),
    m_averageBpp(0.0),
    m_dataStorageEngine(0),
    m_picSizeInCtbsY(0),
    m_picHeightInCtbs(0),
    m_picWidthInCtbs(0),
    m_codedFrames(0),
    m_bitsSpent(0),
    m_picHeight(0),
    m_picWidth(0),
    m_ctuSize(0),
    m_averageBitsPerCtb(0)
    {
    }

    SequenceController(double targetRate,
                       double frameRate,
                       int intraPeriod,
                       int sopSize,
                       int picHeight,
                       int picWidth,
                       int ctuSize,
                       int baseQp);

    ~SequenceController()
    {
        if(m_dataStorageEngine)
            delete m_dataStorageEngine;
        m_pictureControllerEngine.clear();
        m_sopControllerEngine.clear();
#if WRITE_RC_LOG
        if(m_logFile)
            m_logFile.close();
#endif
    }

    void initNewSop(int sopId, int sopSize);

    void pictureRateAllocation(std::shared_ptr<InputQueue::Docket> docket);

    void updateSequenceController(bool isIntra, int sopLevel, int poc, int sopId);

    int getBaseQp() { return m_baseQp; }

    double estimatePictureLambda(bool isIntra, int sopLevel, int poc);

    int deriveQpFromLambda(double lambda, bool isIntra, int sopLevel);

    void pictureRateAllocationIntra(EstimateIntraComplexity &icInfo, int poc);

    void setHeaderBits(int bits, bool isIntra, int sopLevel, int poc);

    void   computeCtuTargetBits(bool isIntraSlice, int ctbAddrInRs, int poc);

    int    estimateCtuLambdaAndQp(bool isIntra, int ctbAddrInRs, int currentPictureLevel, int poc, int sliceQp);

    void   updateCtuController(int codingBitsCtu, bool isIntra, int ctbAddrInRs, int currentPictureLevel, int poc);

    void   updateValidityFlag(bool flag, int ctuIdx, int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        currentPictureController->second->updateCtuValidityFlag(flag, ctuIdx);
    }

    void setValidityFlag(bool flag, int ctuIdx, int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        currentPictureController->second->setCtuValidityFlag(flag, ctuIdx);
    }

    int getCtuStoredQp(int ctuIdx, int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtuStoredQp(ctuIdx);
    }

    void getAveragePictureQpAndLambda(int &averageQp, double &averageLambda, int poc);

    void resetSequenceControllerMemory();

    void initCpbInfo(int cpbMaxSize)
    {
        m_cpbControllerEngine.setCpbInfo(cpbMaxSize, (int)m_targetRate, m_frameRate, 0.8);
    }

#if WRITE_RC_LOG
    void writetoLogFile(string s)
    {
        m_logFile<<s;
        m_logFile.flush();
    }
    int getPictureTargetBits(int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getPictureTargetBits();
    }
#endif

    int getCpbFullness()
    {
        return m_cpbControllerEngine.getCpbStatus();
    }

    void resetSequenceControllerState(bool carryForward)
    {
        m_codedFrames = 0;
        m_bitsSpent = 0;
        if(carryForward)
        {
            m_bitsLeft = m_targetBits + (int)((double)m_bitsLeft * 0.1);
        }
        else
        {
            m_bitsLeft = m_targetBits;
        }
    }

    int getTotalFrames()
    {
        return m_totalFrames;
    }

    double getCtuLambda(int poc, int ctbAddrInRs)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtuLambda(ctbAddrInRs);
    }
    double getCtuReciprocalLambda(int poc, int ctbAddrInRs)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtuReciprocalLambda(ctbAddrInRs);
    }
    double getCtuReciprocalSqrtLambda(int poc, int ctbAddrInRs)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtuReciprocalSqrtLambda(ctbAddrInRs);
    }
};


#endif /* RATECONTROL_H_ */
