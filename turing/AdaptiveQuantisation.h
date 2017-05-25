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

#ifndef INCLUDED_AdaptiveQuantisation_h
#define INCLUDED_AdaptiveQuantisation_h

#include "Picture.h"
#include "StatePicture.h"
#include <limits>
#include <fstream>

namespace {

    int log2Size(int size)
    {
        int log2 = 0;
        while (size != 1)
        {
            size >>= 1;
            ++log2;
        }
        return log2;
    }

};

class AdaptiveQuantisationUnit
{
private:
    double m_activity;

public:
    AdaptiveQuantisationUnit() : m_activity(0.0) {}
    ~AdaptiveQuantisationUnit() {}
    void   setActivity(double value) { m_activity = value; }
    double getActivity()             { return m_activity;  }
};

class AdaptiveQuantisationLayer
{
private:
    int                       m_layerPartitionHeight;
    int                       m_layerPartitionWidth;
    int                       m_picHeightInPartitionUnits;
    int                       m_picWidthInPartitionUnits;
    double                    m_averageActivity;
    AdaptiveQuantisationUnit *m_aqUnitArray;

public:
    AdaptiveQuantisationLayer() : m_aqUnitArray(0),
    m_layerPartitionHeight(0),
    m_layerPartitionWidth (0),
    m_averageActivity(0.0),
    m_picHeightInPartitionUnits(0),
    m_picWidthInPartitionUnits(0) {}

    void create(int picHeight, int picWidth, int unitHeight, int unitWidth)
    {
        m_layerPartitionHeight      = unitHeight;
        m_layerPartitionWidth       = unitWidth;
        m_averageActivity           = 0.0;
        m_picHeightInPartitionUnits = (picHeight + unitHeight - 1) / unitHeight;
        m_picWidthInPartitionUnits  = (picWidth  + unitWidth  - 1) / unitWidth;
        int n                       = m_picHeightInPartitionUnits * m_picWidthInPartitionUnits;
        m_aqUnitArray               = new AdaptiveQuantisationUnit[n];
    }

    ~AdaptiveQuantisationLayer()
    {
        if(m_aqUnitArray)
            delete[] m_aqUnitArray;
    }

    AdaptiveQuantisationUnit *getUnitArray() { return m_aqUnitArray; }

    int getLayerPartitionHeight()    { return m_layerPartitionHeight;      }
    int getLayerPartitionWidth()     { return m_layerPartitionWidth;       }
    int getPicHeightInPartUnits()    { return m_picHeightInPartitionUnits; }
    int getPicWidthInPartUnits()     { return m_picWidthInPartitionUnits;  }
    double getAverageActivity()      { return m_averageActivity;           }
    void   setAverageActivity(int v) { m_averageActivity = v; }
};

class AdaptiveQuantisation
{
private:
    int                        m_aqDepth;
    int                        m_aqRange;
    int                        m_pictureHeight;
    int                        m_pictureWidth;
    AdaptiveQuantisationLayer *m_aqLayerArray;

public:
    AdaptiveQuantisation()  : m_aqLayerArray(0),
    m_aqDepth(3),
    m_aqRange(6),
    m_pictureHeight(0),
    m_pictureWidth(0) {}
    AdaptiveQuantisation(int depth, int range, int picHeight, int picWidth, int maxCuSize)
    {
        m_aqDepth = depth;
        m_aqRange = range;
        m_pictureHeight = picHeight;
        m_pictureWidth  = picWidth;

        m_aqLayerArray = new AdaptiveQuantisationLayer[m_aqDepth+1];

        for(int d = 0; d <= m_aqDepth; d++)
        {
            m_aqLayerArray[d].create(picHeight, picWidth, maxCuSize >> d, maxCuSize >> d);
        }
    }
    ~AdaptiveQuantisation()
    {
        if(m_aqLayerArray)
            delete[] m_aqLayerArray;
    }

    AdaptiveQuantisationLayer *getAqLayerArray() { return m_aqLayerArray; }

    int getAqOffset( int row, int col, int depth)
    {
        depth                                = std::min<int>(depth, m_aqDepth);
        int picWidthInPartUnits              = m_aqLayerArray[depth].getPicWidthInPartUnits();
        int layerPartHeight                  = m_aqLayerArray[depth].getLayerPartitionHeight();
        int layerPartWidth                   = m_aqLayerArray[depth].getLayerPartitionWidth();
        int rowPartLayer                     = row / layerPartHeight;
        int colPartLayer                     = col / layerPartWidth;
        int n                                = rowPartLayer * picWidthInPartUnits + colPartLayer;
        AdaptiveQuantisationLayer &currLayer = m_aqLayerArray[depth];
        AdaptiveQuantisationUnit *currUnit   = currLayer.getUnitArray();
        double unitActivity                  = currUnit[n].getActivity();
        double layerActivity                 = currLayer.getAverageActivity();
        double activityScale                 = pow(2.0, m_aqRange/6.0);
        double normalisedActivity            = (activityScale * unitActivity + layerActivity) / (unitActivity + activityScale * layerActivity);
        int qpOffset                         = (int)floor(log(normalisedActivity) / log(2.0) * 6.0 + 0.49999);
        return qpOffset;
    }

    template<typename Sample>
    void preAnalysis(std::shared_ptr<PictureWrapper> picture)
    {
        auto &pictureInput = static_cast<PictureWrap<Sample> &>(*picture);

        for(int d = 0; d <= m_aqDepth; d++)
        {
            int layerPartHeight = m_aqLayerArray[d].getLayerPartitionHeight();
            int layerPartWidth  = m_aqLayerArray[d].getLayerPartitionWidth();
            AdaptiveQuantisationUnit *currentAqUnit = m_aqLayerArray[d].getUnitArray();
            double sumActivity = 0.0;

            for(int rowPicture = 0; rowPicture < m_pictureHeight; rowPicture += layerPartHeight)
            {
                int currentUnitHeight = std::min<int>(layerPartHeight, m_pictureHeight - rowPicture);
                int log2UnitHeight = log2Size(currentUnitHeight);
                for(int colPicture = 0; colPicture < m_pictureWidth; colPicture += layerPartWidth, currentAqUnit++)
                {
                    int currentUnitWidth = std::min<int>(layerPartWidth, m_pictureWidth - colPicture);
                    int log2UnitWidth = log2Size(currentUnitWidth);
                    uint64_t sum[4]       = {0, 0, 0, 0};
                    uint64_t sumSquare[4] = {0, 0, 0, 0};

                    // Block 0
                    for(int r = 0; r < currentUnitHeight >> 1; r++)
                    {
                        for(int c = 0; c < currentUnitWidth >> 1; c++)
                        {
                            int lumaPixelValue = pictureInput[0](colPicture + c, rowPicture + r);
                            sum[0] += lumaPixelValue;
                            sumSquare[0] = lumaPixelValue*lumaPixelValue;
                        }
                    }

                    // Block 1
                    for(int r = 0; r < currentUnitHeight >> 1; r++)
                    {
                        for(int c = currentUnitWidth >> 1; c < currentUnitWidth; c++)
                        {
                            int lumaPixelValue = pictureInput[0](colPicture + c, rowPicture + r);
                            sum[1] += lumaPixelValue;
                            sumSquare[1] += lumaPixelValue;
                        }
                    }

                    // Block 2
                    for(int r = currentUnitHeight >> 1; r < currentUnitHeight; r++)
                    {
                        for(int c = 0; c < currentUnitWidth >> 1; c++)
                        {
                            int lumaPixelValue = pictureInput[0](colPicture + c, rowPicture + r);
                            sum[2] += lumaPixelValue;
                            sumSquare[2] += lumaPixelValue*lumaPixelValue;
                        }
                    }

                    // Block 3
                    for(int r = currentUnitHeight >> 1; r < currentUnitHeight; r++)
                    {
                        for(int c = currentUnitWidth >> 1; c < currentUnitWidth; c++)
                        {
                            int lumaPixelValue = pictureInput[0](colPicture + c, rowPicture + r);
                            sum[3] += lumaPixelValue;
                            sumSquare[3] += lumaPixelValue*lumaPixelValue;
                        }
                    }

                    double minVar = std::numeric_limits<double>::max();
                    int numPixelUnit = (currentUnitHeight * currentUnitWidth) >> 2;
                    if(numPixelUnit)
                        for(int blockIdx = 0; blockIdx < 4; blockIdx++)
                        {
                            auto average  = static_cast<double>(sum[blockIdx] / numPixelUnit);
                            double variance = sumSquare[blockIdx] / numPixelUnit - average*average;
                            if(variance < minVar)
                                minVar = variance;
                        }
                    else
                        minVar = 0.0;

                    double activity = 1.0 + minVar;
                    currentAqUnit->setActivity(activity);
                    sumActivity += activity;
                }
            }
            double averageActivity = sumActivity / (m_aqLayerArray[d].getPicWidthInPartUnits() * m_aqLayerArray[d].getPicHeightInPartUnits());
            m_aqLayerArray[d].setAverageActivity(static_cast<int>(averageActivity));
        }
    }
};



#endif /* ADAPTIVEQUANTISATION_H_ */
