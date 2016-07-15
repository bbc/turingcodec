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

#ifndef INCLUDED_Aps_h
#define INCLUDED_Aps_h

#pragma once

#include "StateEncode.h"
#include "Global.h"

template <typename Sample>
class Aps
{
private:
    int m_maxApsDepth;
    int m_sadRatioThreshold;
    int m_sadThreshold;

public:
    Aps(): m_maxApsDepth(0),
    m_sadRatioThreshold(2), // 25% tolerance
    m_sadThreshold(4) {}
    ~Aps() {}
    void initApsModule( int maxApsDepth )
    {
        m_maxApsDepth   = maxApsDepth;
    }
    void analyseResidueEnergy (const CandidateStash<Sample> *champion,
                               const coding_quadtree cqt,
                               bool& do2NxN,
                               bool& doNx2N)
    {

        if(cqt.cqtDepth >= m_maxApsDepth)
            return;

        int height = 1 << (cqt.log2CbSize - 1);
        int  width = 1 << (cqt.log2CbSize - 1);

        do2NxN = true;
        doNx2N = true;

        //  Check whether the SAD ratio for 2Nx2N is inside tolerance
        int numerator   = (champion->sadResidueQuad[0][0] + champion->sadResidueQuad[0][1]);
        int denominator = (champion->sadResidueQuad[1][0] + champion->sadResidueQuad[1][1]);
        int threshold = m_sadThreshold * height * width * 2;
        if (numerator < threshold && denominator < threshold)
        {
            do2NxN = false;
        }
        else
        {
            int delta = denominator>>m_sadRatioThreshold;
            do2NxN = !(denominator - delta < numerator && numerator < denominator + delta);
        }

        numerator   = (champion->sadResidueQuad[0][0] + champion->sadResidueQuad[1][0]);
        denominator = (champion->sadResidueQuad[0][1] + champion->sadResidueQuad[1][1]);
        if (numerator < threshold && denominator < threshold)
        {
            doNx2N = false;
        }
        else
        {
            int delta = denominator>>m_sadRatioThreshold;
            doNx2N = !(denominator - delta < numerator && numerator < denominator + delta);
        }
    }

};

#endif
