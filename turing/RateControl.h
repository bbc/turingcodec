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
#define CTB_SMOOTH_WINDOW  4
#define BETA_INTRA_MAD     1.2517
#define WRITE_RC_LOG       0

using namespace std;

struct EstimateIntraComplexity;
struct StateEncodePicture;
struct StatePicture;

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
        unique_lock<mutex> lockCpbInfo(m_cpbStatusMutex);
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
        unique_lock<mutex> lockCpbInfo(m_cpbStatusMutex);
        return m_cpbStatus;
    }

    void adjustAllocatedBits(int &currentBitsPerPicture)
    {
        unique_lock<mutex> lockCpbInfo(m_cpbStatusMutex);

        // Cpb correction
        int estimatedCpbFullness = m_cpbStatus + m_cpbBufferingRate;

        // Check whether the allocated bits will lead to cpb overflow
        int overflowLevel = static_cast<int>(static_cast<double>(m_cpbSize) * 0.9);
        if(estimatedCpbFullness - currentBitsPerPicture > overflowLevel)
        {
            currentBitsPerPicture = estimatedCpbFullness - overflowLevel;
        }

        // Check whether the allocated bits will lead to cpb underflow
        estimatedCpbFullness -= m_cpbBufferingRate;
        int underflowLevel = static_cast<int>(static_cast<double>(m_cpbSize) * 0.1);
        if(estimatedCpbFullness - currentBitsPerPicture < underflowLevel)
        {
            currentBitsPerPicture = max<int>(200, estimatedCpbFullness - underflowLevel);
        }
    }
};

class CodedCtb
{
private:
    double m_alpha;
    double m_beta;

public:
    CodedCtb() : m_alpha(INITIAL_ALPHA), m_beta(INITIAL_BETA) {}
    void   setAlpha(double value) { m_alpha = value; }
    void   setBeta (double value) { m_beta  = value; }

    double getAlpha()             { return m_alpha;  }
    double getBeta()              { return m_beta;   }
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
    int    m_poc;
    int    m_updateCounter;

public:
    CodedPicture() : m_codingBits(0),
    m_headerBits(0),
    m_lambda(NON_VALID_LAMBDA),
    m_qp(NON_VALID_QP),
    m_alpha(INITIAL_ALPHA),
    m_beta(INITIAL_BETA),
    m_level(-1),
    m_poc(-1),
    m_updateCounter(0) {}

    void setCodingBits(int bits)   { m_codingBits = bits; }
    void setHeaderBits(int bits)
    {
        m_headerBits = bits;
        m_updateCounter++;
    }
    void setLambda(double lambda)  { m_lambda = lambda;   }
    void setQp(int qp)             { m_qp = qp;           }
    void setAlpha(double alpha)    { m_alpha = alpha;     }
    void setBeta(double  beta)     { m_beta  = beta;      }
    void setLevel(int level)       { m_level = level;     }
    void setPoc(int poc)           { m_poc = poc;         }

    double getAlpha()              { return m_alpha;      }
    double getBeta()               { return m_beta;       }
    double getLambda()             { return m_lambda;     }
    int    getQp()                 { return m_qp;         }
    int    getCodingBits()         { return m_codingBits; }
    int    getHeaderBits()         { return m_headerBits; }
    int    getLevel()              { return m_level;      }
    int    getPoc()                { return m_poc;        }
    int    getHeaderBitsEstimate()
    {
        int estimate = 0;
        if(m_updateCounter)
            estimate = m_headerBits / m_updateCounter;
        return estimate;
    }
};

class DataStorage
{
private:
    array<CodedPicture, MAX_NUM_LEVELS> m_codedPicturesPerLevel;
    array<CodedCtb*, MAX_NUM_LEVELS>    m_codedCtbsPerLevel;
    mutex m_codedPicsLevelToken;
    mutex m_codedCtbsLevelToken;

public:

    ~DataStorage()
    {
        for(int level = 0; level < MAX_NUM_LEVELS; level++)
        {
            delete[] m_codedCtbsPerLevel[level];
        }
    }

    CodedPicture& getPictureAtLevel(int level)
    {
        assert(level < m_codedPicturesPerLevel.size());
        return m_codedPicturesPerLevel[level];
    }

    CodedCtb *getCodedCtbAtLevel(int level)
    {
        assert(level < m_codedCtbsPerLevel.size());
        return m_codedCtbsPerLevel[level];
    }

    void initCtbStorage(int numberOfUnits);

    void initPictureStorage();

    void takeTokenOnPictureLevel()    { m_codedPicsLevelToken.lock();   }

    void releaseTokenOnPictureLevel() { m_codedPicsLevelToken.unlock(); }

    void takeTokenOnCodedCtbs()       { m_codedCtbsLevelToken.lock();   }

    void releaseTokenOnCodedCtbs()    { m_codedCtbsLevelToken.unlock(); }

};

class CtbController
{
private:
    int    m_codedBits;
    int    m_targetBits;
    int    m_targetBitsEst;
    int    m_ctbQp;
    int    m_lastValidQp;
    double m_ctbLambda;
    double m_ctbAlpha;
    double m_ctbBeta;
    double m_ctbReciprocalLambda;
    double m_ctbReciprocalSqrtLambda;
    double m_costIntra;
    double m_bpp;
    int    m_numberOfPixels;
    bool   m_isValid;       // non fully skipped CTB
    bool   m_finished;
    bool   m_isIntra;

public:

    CtbController() : m_codedBits(0),
    m_ctbQp(NON_VALID_QP),
    m_lastValidQp(NON_VALID_QP),
    m_costIntra(0),
    m_ctbLambda(NON_VALID_LAMBDA),
    m_ctbAlpha(INITIAL_ALPHA),
    m_ctbBeta(INITIAL_BETA),
    m_targetBits(0),
    m_numberOfPixels(0),
    m_isValid(false),
    m_targetBitsEst(0),
    m_finished(false),
    m_ctbReciprocalLambda(NON_VALID_LAMBDA),
    m_ctbReciprocalSqrtLambda(NON_VALID_LAMBDA),
    m_isIntra(false),
    m_bpp(0.0) {}

    void setCodedBits(int bits)                   { m_codedBits      = bits;           }
    void setCtbQp(int qp)                         { m_ctbQp          = qp;             }
    void setLastValidQp(int qp)                   { m_lastValidQp = qp;                }
    void setTargetBits(int bits)                  { m_targetBits     = bits;           }
    void setTargetBitsEst(int bits)               { m_targetBitsEst  = bits;           }
    void setCtbLambda(double lambda)              { m_ctbLambda      = lambda;         }
    void setCtbAlpha(double alpha)                { m_ctbAlpha       = alpha;          }
    void setCtbBeta(double beta)                  { m_ctbBeta        = beta;           }
    void setCtbReciprocalLambda(double value)     { m_ctbReciprocalLambda = value;     }
    void setCtbReciprocalSqrtLambda(double value) { m_ctbReciprocalSqrtLambda = value; }
    void setNumberOfPixels(int pixels)            { m_numberOfPixels = pixels;         }
    void setCostIntra(double cost)                { m_costIntra      = cost;           }
    void setValidityFlag(bool flag)               { m_isValid        = flag;           }
    void updateValidityFlag(bool flag)
    {
        m_isValid       |= flag;
    }
    void setFinishedFlag(bool flag)               { m_finished       = flag;           }
    void setCtbBpp(double bpp)                    { m_bpp            = bpp;            }
    void setIsIntra(bool isIntra)                 { m_isIntra        = isIntra;        }  

    int     getCodedBits()                        { return m_codedBits;                }
    int     getCtbQp()                            { return m_ctbQp;                    }
    int     getLastValidQp()                      { return m_lastValidQp;              }
    int     getTargetBits()                       { return m_targetBits;               }
    int     getTargetBitsEst()                    { return m_targetBitsEst;            }
    double  getCtbLambda()                        { return m_ctbLambda;                }
    double  getCtbAlpha()                         { return m_ctbAlpha;                 }
    double  getCtbBeta()                          { return m_ctbBeta;                  }
    double  getCtbReciprocalLambda()              { return m_ctbReciprocalLambda;      }
    double  getCtbReciprocalSqrtLambda()          { return m_ctbReciprocalSqrtLambda;  }
    int     getNumberOfPixels()                   { return m_numberOfPixels;           }
    double  getCostIntra()                        { return m_costIntra;                }
    bool    getValidityFlag()                     { return m_isValid;                  }
    bool    getFinishedFlag()                     { return m_finished;                 }
    double  getCtbBpp()                           { return m_bpp;                      }
    bool    getIsIntra()                          { return m_isIntra;                  }
};

class PictureController
{
private:
    int m_targetBits;
    int m_pixelsPerPicture;
    int m_pictureHeightInCtbs;
    int m_pictureWidthInCtbs;
    int m_pictureSizeInCtbs;
    int m_totalCostIntra;
    bool m_isIntra;
    int m_sopLevel;
    int m_hierarchyLevel;
    int m_numLeftSameHierarchyLevel;
    bool m_paramsUpdated;
    int m_sopSize;
    int m_poc;
    double m_averageBpp;
    double m_alphaUpdateStep;
    double m_betaUpdateStep;
    double m_pictureLambdaIni;
    double m_pictureLambdaFin;
    int    m_pictureQp;
    double m_pictureAlpha;
    double m_pictureBeta;
    double m_pictureBpp;
    DataStorage   *m_dataStorageAccess;
    CtbController *m_ctbControllerEngine;
    int    m_codingBits;
    int    m_intraPoc;
    int    m_segmentPoc;

    int  estimateCtbLambdaAndQpIntra(int ctbAddrInRs);
    int  estimateCtbLambdaAndQpInter(int ctbAddrInRs);

public:
    PictureController() :
    m_targetBits(0),
    m_pixelsPerPicture(0),
    m_averageBpp(0.0),
    m_alphaUpdateStep(ALPHA_UPDATE_STEP),
    m_betaUpdateStep(BETA_UPDATE_STEP),
    m_dataStorageAccess(0),
    m_pictureLambdaIni(NON_VALID_LAMBDA),
    m_pictureLambdaFin(NON_VALID_LAMBDA),
    m_pictureQp(NON_VALID_QP),
    m_pictureAlpha(0),
    m_pictureBeta(0),
    m_pictureBpp(0.0),
    m_ctbControllerEngine(0),
    m_pictureHeightInCtbs(0),
    m_pictureWidthInCtbs(0),
    m_pictureSizeInCtbs(0),
    m_totalCostIntra(0),
    m_sopLevel(0),
    m_isIntra(false),
    m_hierarchyLevel(-1),
    m_numLeftSameHierarchyLevel(-1),
    m_paramsUpdated(false),
    m_sopSize(-1),
    m_poc(-1),
    m_codingBits(0),
    m_segmentPoc(-1),
    m_intraPoc(-1) {}

    PictureController(DataStorage *storageAccess, int pictureHeight, int pictureWidth, int pictureHeightInCtbs, int pictureWidthInCtbs, int ctbSize, int targetBits, double averageBpp, std::shared_ptr<InputQueue::Docket> docket);

    ~PictureController()
    {
        if(m_ctbControllerEngine)
            delete[] m_ctbControllerEngine;
    }

    int getPictureQp()                 { return m_pictureQp;                 }
    double getPictureLambdaIni()       { return m_pictureLambdaIni;          }
    double getPictureAlpha()           { return m_pictureAlpha;              }
    double getPictureBeta()            { return m_pictureBeta;               }
    bool getIsIntra()                  { return m_isIntra;                   }
    int getSopLevel()                  { return m_sopLevel;                  }
    int getPictureSizeInCtbs()         { return m_pictureSizeInCtbs;         }
    int getHierarchyLevel()            { return m_hierarchyLevel;            }
    int getNumLeftSameHierarchyLevel() { return m_numLeftSameHierarchyLevel; }
    bool getParamsUpdates()            { return m_paramsUpdated;             }
    int getPoc()                       { return m_poc;                       }
    const int getCodingBits()          { return m_codingBits;                }
    const int getIntraPoc()            { return m_intraPoc;                  }
    const int getSegmentPoc()          { return m_segmentPoc;                }
    DataStorage* getStorageAccess()    { return m_dataStorageAccess;         }
    int getCtbStoredQp(int ctbIdx)
    {
        assert(ctbIdx < m_pictureSizeInCtbs);
        return m_ctbControllerEngine[ctbIdx].getCtbQp();
    }
    double getCtbLambda(int ctbAddrInRs)
    {
        assert(ctbAddrInRs < m_pictureSizeInCtbs);
        return m_ctbControllerEngine[ctbAddrInRs].getCtbLambda();
    }
    double getCtbReciprocalLambda(int ctbAddrInRs)
    {
        assert(ctbAddrInRs < m_pictureSizeInCtbs);
        return m_ctbControllerEngine[ctbAddrInRs].getCtbReciprocalLambda();
    }
    double getCtbReciprocalSqrtLambda(int ctbAddrInRs)
    {
        assert(ctbAddrInRs < m_pictureSizeInCtbs);
        return m_ctbControllerEngine[ctbAddrInRs].getCtbReciprocalSqrtLambda();
    }

    bool getCtbFinishedFlag(int ctbAddrInRs)
    {
        assert(ctbAddrInRs < m_pictureSizeInCtbs);
        return m_ctbControllerEngine[ctbAddrInRs].getFinishedFlag();
    }

    void setPictureAlpha(double alpha)        { m_pictureAlpha = alpha;          }
    void setPictureBeta(double beta)          { m_pictureBeta = beta;            }
    void decreaseNumLeftSameHierarchyLevel()  { m_numLeftSameHierarchyLevel--;   }
    void setParamsUpdates(bool paramsUpdates) { m_paramsUpdated = paramsUpdates; }
    void updateCtbValidityFlag(bool flag, int ctbIdx)
    {
        assert(ctbIdx < m_pictureSizeInCtbs);
        m_ctbControllerEngine[ctbIdx].updateValidityFlag(flag);
    }
    void setCtbValidityFlag(bool flag, int ctbIdx)
    {
        assert(ctbIdx < m_pictureSizeInCtbs);
        m_ctbControllerEngine[ctbIdx].setValidityFlag(flag);
    }

    void updateAfterEncoding(const int codingBits);

    void computeCtbTargetBits(int ctbAddrInRs, int lastValidQp, int64_t cumulativeTargetBits, int64_t cumulativeBitsSpent, int64_t ctbsLeftInPicture);

    void updateModelParameters();

    void getCodedInfoFromWavefront(bool wpp, int ctbAddrInRs, int& lastValidQp, int64_t &cumulativeTargetBits, int64_t &cumulativeBitsSpent, int64_t &cumulativeCtbsCoded, bool currFrame);

    double estimateLambda(double previousLambdaIni);

    int estimateQp(int sopLevel, int previousQp);

    int getPictureTargetBits()   { return m_targetBits;  }

    int estimateCtbLambdaAndQp(int ctbAddrInRs);

    void storeCtbParameters(int codingBits, bool isIntra, int ctbAddrInRs, int pictureLevel);

    void getAveragePictureQpAndLambda(int &averageQp, double &averageLambda);

    int getFinishedCtbs();

    int getTotalCodingBitsCtbs();

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

public:
    SOPController() : m_size(0),
    m_targetBits(0),
    m_weight(0),
    m_sopId(-1),
    m_framesCoded(0),
    m_bitsCoded(0),
    m_dataStorageAccess(nullptr) {}

    SOPController(DataStorage *dataStorageAccess, int sopId, int size, int *bitrateWeight, int targetBits);

    ~SOPController()
    {
        if(m_weight)
        {
            delete[] m_weight;
        }
    }

    int getRateCurrentPicture(std::shared_ptr<InputQueue::Docket> docket);

    void updateSopController(int bitsCoded);

    bool finished()
    {
        return m_framesCoded == m_size;
    }

};

struct IntraPeriodData
{
    int m_framesLeft;
    int64_t m_bitsSpent;
    int64_t m_targetBits;
    int64_t m_targetBitsEst;
    int m_totalCtbsCoded;
    int m_framesEncoded;
    IntraPeriodData() :
        m_framesLeft(0),
        m_bitsSpent(0),
        m_targetBits(0),
        m_targetBitsEst(0),
        m_totalCtbsCoded(0),
        m_framesEncoded(0) {};
};

struct SegmentData
{
    DataStorage    *m_dataStorageEngine;
    SegmentData(int picSizeInCtbsY)
    {
        m_dataStorageEngine = new DataStorage();
        m_dataStorageEngine->initCtbStorage(picSizeInCtbsY);
        m_dataStorageEngine->initPictureStorage();
    };

    ~SegmentData()
    {
        if (m_dataStorageEngine)
            delete m_dataStorageEngine;
    }
};


class SequenceController
{
private:
    double  m_targetRate;
    int     m_smoothingWindow;
    int     m_totalFrames;
    int     m_totalCtbs;
    double  m_frameRate;
    int     m_sopSize;
    double  m_averageRate;
    double  m_averageBpp;
    int     m_bitsPerPicture;
    int64_t m_targetBits;
    int     m_pixelsPerPicture;
    int     m_sopWeight[7][4][8] = {{{8, 30, 0, 0, 0, 0, 0, 0}, {7, 25, 0, 0, 0, 0, 0, 0}, {6, 20, 0, 0, 0, 0, 0, 0}, {5, 15, 0, 0, 0, 0, 0, 0}},
                                    {{4, 8, 30, 0, 0, 0, 0, 0}, {4, 7, 25, 0, 0, 0, 0, 0}, {4, 6, 20, 0, 0, 0, 0, 0}, {4, 5, 15, 0, 0, 0, 0, 0}},
                                    {{4, 4, 8, 30, 0, 0, 0, 0}, {4, 4, 7, 25, 0, 0, 0, 0}, {4, 4, 6, 20, 0, 0, 0, 0}, {4, 4, 5, 15, 0, 0, 0, 0}},
                                    {{4, 4, 8, 8, 30, 0, 0, 0}, {4, 4, 7, 7, 25, 0, 0, 0}, {4, 4, 6, 6, 20, 0, 0, 0}, {4, 4, 5, 5, 15, 0, 0, 0}},
                                    {{4, 4, 8, 4, 8, 30, 0, 0}, {4, 4, 7, 4, 7, 25, 0, 0}, {4, 4, 6, 4, 6, 20, 0, 0}, {4, 4, 5, 6, 5, 15, 0, 0}},
                                    {{1, 4, 1, 8, 1, 4, 30, 0}, {1, 4, 1, 7, 1, 4, 25, 0}, {1, 4, 1, 6, 1, 4, 20, 0}, {1, 4, 1, 6, 1, 4, 15, 0}},
                                    {{1, 4, 1, 8, 1, 4, 1, 30}, {1, 4, 1, 7, 1, 4, 1, 25}, {1, 4, 1, 6, 1, 4, 1, 20}, {1, 4, 1, 5, 1, 4, 1, 15}}}; // As in HM
    int     m_baseQp;

    map<int, PictureController*> m_pictureControllerEngine;
    map<int, SOPController*>     m_sopControllerEngine;
    CpbInfo                      m_cpbControllerEngine;
    map<int, IntraPeriodData*>   m_ipDataEngine;
    map<int, SegmentData*>       m_segmentData;

    mutex  m_pictureControllerMutex;
    mutex  m_sopControllerMutex;
    mutex  m_ipDataMutex;
    mutex  m_segmentDataMutex;

    int    m_picSizeInCtbsY;
    int    m_picHeightInCtbs;
    int    m_picWidthInCtbs;
    int    m_picHeight;
    int    m_picWidth;
    int    m_ctbSize;
    int    m_averageBitsPerCtb;
    bool   m_useWpp;
    int    m_concurrentFrames;
#if WRITE_RC_LOG
    ofstream m_logFile;
#endif

public:
    SequenceController() : m_targetRate(0),
    m_smoothingWindow(40),
    m_frameRate(0),
    m_sopSize(8),
    m_averageRate(0),
    m_targetBits(0),
    m_bitsPerPicture(0),
    m_pixelsPerPicture(0),
    m_totalFrames(0),
    m_totalCtbs(0),
    m_baseQp(0),
    m_averageBpp(0.0),
    m_picSizeInCtbsY(0),
    m_picHeightInCtbs(0),
    m_picWidthInCtbs(0),
    m_picHeight(0),
    m_picWidth(0),
    m_ctbSize(0),
    m_averageBitsPerCtb(0),
    m_useWpp(false),
    m_concurrentFrames(0)
    {
    }

    SequenceController(double targetRate,
                       double frameRate,
                       int    intraPeriod,
                       int    sopSize,
                       int    picHeight,
                       int    picWidth,
                       int    ctbSize,
                       int    baseQp,
                       bool   useWpp,
                       int    concurrentFrames);

    ~SequenceController()
    {
        m_pictureControllerEngine.clear();
        m_sopControllerEngine.clear();
        m_ipDataEngine.clear();
        m_segmentData.clear();
#if WRITE_RC_LOG
        if(m_logFile)
            m_logFile.close();
#endif
    }

    int  numLeftSameHierarchyLevel(int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        PictureController *currentPictureController = m_pictureControllerEngine.find(poc)->second;
        return currentPictureController->getNumLeftSameHierarchyLevel();
    }

    template <class H>
    void   computeCtbTargetBits(H &h, bool isIntraSlice, int ctbAddrInRs, int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());

        int avgBits = 0;
        int lastValidQp = NON_VALID_QP;
        int64_t  cumulativeTargetBits = 0,  cumulativeBitsSpent = 0, ctbsLeftInPicture = 0;

        getCodedInfoCtbs(h, ctbAddrInRs, lastValidQp, cumulativeTargetBits, cumulativeBitsSpent, ctbsLeftInPicture);

        if (!currentPictureController->second->getIsIntra())
        {
            unique_lock<mutex> lock(m_ipDataMutex);
            auto currentIpData = m_ipDataEngine.find(currentPictureController->second->getIntraPoc());
            assert(currentIpData != m_ipDataEngine.end());
            cumulativeBitsSpent += currentIpData->second->m_bitsSpent;
            cumulativeTargetBits += currentIpData->second->m_targetBits;
        }
        currentPictureController->second->computeCtbTargetBits(ctbAddrInRs, lastValidQp, cumulativeTargetBits, cumulativeBitsSpent, ctbsLeftInPicture);
    }

    template <class H>
    void getCodedInfoCtbs(H &h, int ctbAddrInRs, int &lastValidQp, int64_t &cumulativeTargetBits, int64_t &cumulativeBitsSpent, int64_t &ctbsLeftInPicture)
    {
        auto currentPictureController = m_pictureControllerEngine.find(h[PicOrderCntVal()]);
        assert(currentPictureController != m_pictureControllerEngine.end());
        int64_t cumulativeCtbsCoded = 0;
        currentPictureController->second->getCodedInfoFromWavefront(m_useWpp, ctbAddrInRs, lastValidQp, cumulativeTargetBits, cumulativeBitsSpent, cumulativeCtbsCoded, true);
        ctbsLeftInPicture = currentPictureController->second->getPictureSizeInCtbs() - cumulativeCtbsCoded;
        if (currentPictureController->second->getIsIntra())
            return;

        const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
        const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];

        for(auto &pictureController : m_pictureControllerEngine)
        {
            if(pictureController.second->getIntraPoc() == currentPictureController->second->getIntraPoc())
            {
                if(pictureController.second->getHierarchyLevel() < currentPictureController->second->getHierarchyLevel() && (pictureController.second->getHierarchyLevel() >= currentPictureController->second->getHierarchyLevel() - (m_concurrentFrames - 1)))
                {
                    int offset = currentPictureController->second->getSopLevel() - pictureController.second->getSopLevel();
                    offset = offset < 0 ? 0 : offset;
                    int depX = std::min(rx + (2 * offset), h[PicWidthInCtbsY()] - 1);
                    int depY = std::min(ry + (1 * offset), h[PicHeightInCtbsY()] - 1);
                    int ctbAddrInRsRef = ctbAddrInRs;// h[PicWidthInCtbsY()]* depY + depX;
                    int dummyLastQp;
                    pictureController.second->getCodedInfoFromWavefront(m_useWpp, ctbAddrInRsRef, dummyLastQp, cumulativeTargetBits, cumulativeBitsSpent, cumulativeCtbsCoded, false);
                }
                else if (pictureController.second->getHierarchyLevel() < currentPictureController->second->getHierarchyLevel())
                {
                    cumulativeBitsSpent += pictureController.second->getCodingBits();
                    cumulativeTargetBits += pictureController.second->getPictureTargetBits();
                }
            }
        }
    }

    void updateAfterEncoding(std::shared_ptr<InputQueue::Docket> docket, const int codingBits)
    {
        {
            unique_lock<mutex> lockPicture(m_pictureControllerMutex);
            auto currentPictureController = m_pictureControllerEngine.find(docket->poc);
            assert(currentPictureController != m_pictureControllerEngine.end());
            currentPictureController->second->updateAfterEncoding(codingBits);
        }

        {
            unique_lock<mutex> lockSopController(m_sopControllerMutex);
            auto currentSopController = m_sopControllerEngine.find(docket->sopId);
            assert(currentSopController != m_sopControllerEngine.end());
            if(!currentSopController->second->finished())
            {
                currentSopController->second->updateSopController(codingBits);
            }
        }
    }

    void updateSequenceControllerFinishedFrames(std::shared_ptr<InputQueue::Docket> docket)
    {
        // Take the token on both picture controller and ip data
        unique_lock<mutex> lockPicture(m_pictureControllerMutex);
        unique_lock<mutex> lockIpData(m_ipDataMutex);

        auto currentIpData = m_ipDataEngine.find(docket->intraFramePoc);
        assert(currentIpData != m_ipDataEngine.end());
        int poc = docket->poc;
        if(docket->numSameHierarchyLevel > 1)
        {
            for (auto &currentPictureController : m_pictureControllerEngine)
            {
                if (currentPictureController.second->getNumLeftSameHierarchyLevel()>1 ||( currentPictureController.second->getHierarchyLevel() == docket->hierarchyLevel  && currentPictureController.first != docket->poc   && currentPictureController.second->getFinishedCtbs() != m_picSizeInCtbsY))
                {
                    return;
                }
            }
        }
        bool removeIntraPeriod = true;
        bool removeSegment = true;
        for (auto currentPictureController = m_pictureControllerEngine.cbegin(); currentPictureController != m_pictureControllerEngine.cend() /* not hoisted */; /* no increment */)
        {
            const bool updated = currentPictureController->second->getParamsUpdates();
            const bool correctHierarchyLevel = currentPictureController->second->getHierarchyLevel() <= (docket->hierarchyLevel - m_concurrentFrames + 1);
            const bool correctIntraPeriod = currentPictureController->second->getIntraPoc() == docket->intraFramePoc;
            const bool previousIntraPeriod = currentPictureController->second->getIntraPoc() < docket->intraFramePoc;
            const bool previousSegment     = currentPictureController->second->getSegmentPoc() < docket->segmentPoc;
            if (updated && correctHierarchyLevel && correctIntraPeriod)
            {
                currentIpData->second->m_bitsSpent += currentPictureController->second->getCodingBits();
                currentIpData->second->m_targetBits += currentPictureController->second->getPictureTargetBits();
                currentIpData->second->m_totalCtbsCoded += currentPictureController->second->getPictureSizeInCtbs();
                currentIpData->second->m_framesEncoded++;
                m_pictureControllerEngine.erase(currentPictureController++);
                
            }
            else
            {
                if (previousIntraPeriod && currentPictureController->second->getFinishedCtbs() != m_picSizeInCtbsY)
                    removeIntraPeriod = false;

                if (previousSegment     && currentPictureController->second->getFinishedCtbs() != m_picSizeInCtbsY)
                    removeSegment = false;

                ++currentPictureController;
            }
        }

        if(removeIntraPeriod)
        {
            for (auto currentPictureController = m_pictureControllerEngine.cbegin(); currentPictureController != m_pictureControllerEngine.cend() /* not hoisted */; /* no increment */)
            {
                const bool previousIntraPeriod = currentPictureController->second->getIntraPoc() < docket->intraFramePoc;
                if (previousIntraPeriod)
                {
                    m_pictureControllerEngine.erase(currentPictureController++);
                }
                    
                else
                    ++currentPictureController;
            }

            for (auto currentIpData = m_ipDataEngine.cbegin(); currentIpData != m_ipDataEngine.cend() /* not hoisted */; /* no increment */)
            {
                const bool previousIntraPeriod = currentIpData->first < docket->intraFramePoc;
                if (previousIntraPeriod)
                    m_ipDataEngine.erase(currentIpData++);
                else
                    ++currentIpData;
            }
        }

        if (removeSegment)
        {
            for (auto currentSegment = m_segmentData.cbegin(); currentSegment != m_segmentData.cend() /* not hoisted */; /* no increment */)
            {
                const bool previousSegment = currentSegment->first < docket->segmentPoc;
                if (previousSegment)
                    m_segmentData.erase(currentSegment++);
                else
                    ++currentSegment;
            }
        }

        {
            // Remove also finished SOP controllers
            unique_lock<mutex> lockSopController(m_sopControllerMutex);
            for(auto currentSopController = m_sopControllerEngine.cbegin(); currentSopController != m_sopControllerEngine.cend(); )
            {
                if(currentSopController->second->finished())
                    m_sopControllerEngine.erase(currentSopController++);
                else
                    ++currentSopController;
            }
        }
    }

    void initNewSop(std::shared_ptr<InputQueue::Docket> docket);

    void pictureRateAllocation(std::shared_ptr<InputQueue::Docket> docket);

    void decreaseNumLeftSameHierarchyLevel(std::shared_ptr<InputQueue::Docket> docket);

    int getBaseQp() { return m_baseQp; }

    double estimatePictureLambda(int poc);

    int deriveQpFromLambda(double lambda, bool isIntra, int sopLevel, int poc);

    void pictureRateAllocationIntra(std::shared_ptr<InputQueue::Docket> docket);

    int estimateCtbLambdaAndQp(int ctbAddrInRs, int poc);

    void storeCtbParameters(int codingBitsCtb, bool isIntra, int ctbAddrInRs, int currentPictureLevel, int poc);

    void updateValidityFlag(bool flag, int ctbIdx, int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        currentPictureController->second->updateCtbValidityFlag(flag, ctbIdx);
    }

    void setValidityFlag(bool flag, int ctbIdx, int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        currentPictureController->second->setCtbValidityFlag(flag, ctbIdx);
    }

    int getCtbStoredQp(int ctbIdx, int poc)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtbStoredQp(ctbIdx);
    }

    void getAveragePictureQpAndLambda(int &averageQp, double &averageLambda, int poc);

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

    double getCtbLambda(int poc, int ctbAddrInRs)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtbLambda(ctbAddrInRs);
    }
    double getCtbReciprocalLambda(int poc, int ctbAddrInRs)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtbReciprocalLambda(ctbAddrInRs);
    }
    double getCtbReciprocalSqrtLambda(int poc, int ctbAddrInRs)
    {
        unique_lock<mutex> lock(m_pictureControllerMutex);
        auto currentPictureController = m_pictureControllerEngine.find(poc);
        assert(currentPictureController != m_pictureControllerEngine.end());
        return currentPictureController->second->getCtbReciprocalSqrtLambda(ctbAddrInRs);
    }

    void initNewIntraPeriod(std::shared_ptr<InputQueue::Docket> docket);
};


#endif /* RATECONTROL_H_ */
