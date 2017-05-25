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

#ifndef INCLUDED_CandModeList_H
#define INCLUDED_CandModeList_H

#include <algorithm>
#include <array>


struct CandModeList
{
    typedef Left A;
    typedef Up B;

    template <class Direction, class H>
    static int getCandidate(H &h, int xPb, int yPb)
    {
        const int xNb = xPb - (std::is_same<Direction, A>::value ? 1 : 0);
        const int yNb = yPb - (std::is_same<Direction, B>::value ? 1 : 0);

    AvailabilityCtu *availabilityCtu = h;
    
    if (!availabilityCtu->available(xPb, yPb, xNb, yNb, h[CtbLog2SizeY()]))
            return INTRA_DC;
        else if (h[CuPredMode(xNb, yNb)] != MODE_INTRA || h[pcm_flag(xNb, yNb)] == 1)
            return INTRA_DC;
        else if (std::is_same<Direction, B>::value && yPb - 1 < ((yPb >> h[CtbLog2SizeY()]) << h[CtbLog2SizeY()]))
            return INTRA_DC;
        else if (std::is_same<Direction, Left>::value)
            return h[left(-1, 0, IntraPredModeY(xNb, yNb))];
        else
            return h[up(0, -1, IntraPredModeY(xNb, yNb))];
    }

    template <class H>
    void populate(int dummy, H &h, int xPb, int yPb)
    {
        const int candIntraPredModeA = getCandidate<Left>(h, xPb, yPb);
        const int candIntraPredModeB = getCandidate<Up>(h, xPb, yPb);

        if (candIntraPredModeA == candIntraPredModeB)
        {
            this->neighbourModes = 1;

            if (candIntraPredModeA < 2)
            {
                candModeList[0] = INTRA_PLANAR;
                candModeList[1] = INTRA_DC;
                candModeList[2] = INTRA_ANGULAR26;
            }
            else
            {
                candModeList[0] = candIntraPredModeA;
                candModeList[1] = ((candIntraPredModeA + 29) % 32) + 2;
                candModeList[2] = ((candIntraPredModeA - 1) % 32) + 2;
            }
        }
        else
        {
            this->neighbourModes = 2;

            candModeList[0] = candIntraPredModeA;
            candModeList[1] = candIntraPredModeB;
            if (candModeList[0] != INTRA_PLANAR && candModeList[1] != INTRA_PLANAR)
            {
                candModeList[2] = INTRA_PLANAR;
            }
            else if (candModeList[0] != INTRA_DC && candModeList[1] != INTRA_DC)
            {
                candModeList[2] = INTRA_DC;
            }
            else
            {
                candModeList[2] = INTRA_ANGULAR26;
            }
        }
    }

    void sort()
    {
        std::sort(std::begin(this->candModeList), std::end(this->candModeList));
    }

    int operator[](int n)
    {
        return candModeList[n];
    }
    std::array<int, 3> candModeList;
    int neighbourModes;
};


static std::ostream &operator<<(std::ostream &o, const CandModeList &candModeList)
{
    o << candModeList.candModeList[0] << ",";
    o << candModeList.candModeList[1] << ",";
    o << candModeList.candModeList[2];
    return o;
}

#endif

