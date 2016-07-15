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

#ifndef INCLUDED_Levels_h
#define INCLUDED_Levels_h

#pragma once

// Information and code concerning the HEVC Level definitions (i.e. Level 5.0, etc.)

#include <string>
#include <sstream>


#define LEVEL_PARAMETERS \
        X(1, MaxLumaPs, "Max luma picture size", "samples") \
        X(1, MaxDpbSize, "Max DPB size", "pictures") \
        X(1, MaxCPB, "Max CPB size", "1000 bits") \
        X(1, MaxSliceSegmentsPerPicture, "Max slice segments per picture", "") \
        X(1, MaxTileRows, "Max # of tile rows", "") \
        X(1, MaxTileCols, "Max # of tile columns", "") \
        X(0, MaxLumaSr, "Max luma sample rate", "samples/sec") \
        X(0, MaxBR, "Max bit rate", "1000 bits/s") \
        X(0, MinCr, "Min Compression Ratio", "")  \

static const bool isMaximum(const char *description)
{
    return description[1] == 'a';
}

static const bool isMinimum(const char *description)
{
    return description[1] == 'i';
}

static std::string friendly(std::string units)
{
    if (units == "1000 bits") return "kbit";
    if (units == "1000 bits/s") return "kbit/s";
    return units;
}

struct Level
{
    enum Parameters
    {
#define X(general, name, description, units) name,
        LEVEL_PARAMETERS
#undef X
    };

    static const int nParameters =
#define X(general,name, description, units) 1 +
            LEVEL_PARAMETERS 0;
#undef X

    int units, tenths, tierFlag;
    double parameters[nParameters];

    int level_idc() const
    {
        return 30 * units + 3 * tenths;
    }

    std::string str() const
    {
        std::ostringstream os;
        os << "Level " << units;
        if (tenths) os << "." << tenths;
        os << (tierFlag ? " High" : " Main") << " tier";
        return os.str();
    }
};


static const Level levels[] = {
        { 1, 0, 0, { 36864, 0, 350, 16, 1, 1, 552960, 128, 2 } },
        { 2, 0, 0, { 122880, 0, 1500, 16, 1, 1, 3686400, 1500, 2 } },
        { 2, 1, 0, { 245760, 0, 3000, 20, 1, 1, 7372800, 3000, 2 } },
        { 3, 0, 0, { 552960, 0, 6000, 30, 2, 2, 16588800, 6000, 2 } },
        { 3, 1, 0, { 983040, 0, 10000, 40, 3, 3, 33177600, 10000, 2 } },
        { 4, 0, 0, { 2228224, 0, 12000, 75, 5, 5, 66846720, 12000, 4 } },
        { 4, 0, 1, { 2228224, 0, 30000, 75, 5, 5, 66846720, 30000, 4 } },
        { 4, 1, 0, { 2228224, 0, 20000, 75, 5, 5, 133693440, 20000, 4 } },
        { 4, 1, 1, { 2228224, 0, 50000, 75, 5, 5, 133693440, 50000, 4 } },
        { 5, 0, 0, { 8912896, 0, 25000, 200, 11, 10, 267386880, 25000, 6 } },
        { 5, 0, 1, { 8912896, 0, 100000, 200, 11, 10, 267386880, 100000, 6 } },
        { 5, 1, 0, { 8912896, 0, 40000, 200, 11, 10, 534773760, 40000, 8 } },
        { 5, 1, 1, { 8912896, 0, 160000, 200, 11, 10, 534773760, 160000, 8 } },
        { 5, 2, 0, { 8912896, 0, 60000, 200, 11, 10, 1069547520, 60000, 8 } },
        { 5, 2, 1, { 8912896, 0, 240000, 200, 11, 10, 1069547520, 240000, 8 } },
        { 6, 0, 0, { 35651584, 0, 60000, 600, 22, 20, 1069547520, 60000, 8 } },
        { 6, 0, 1, { 35651584, 0, 240000, 600, 22, 20, 1069547520, 240000, 8 } },
        { 6, 1, 0, { 35651584, 0, 120000, 600, 22, 20, 2139095040, 120000, 8 } },
        { 6, 1, 1, { 35651584, 0, 480000, 600, 22, 20, 2139095040, 480000, 8 } },
        { 6, 2, 0, { 35651584, 0, 240000, 600, 22, 20, 4278190080, 240000, 6 } },
        { 6, 2, 1, { 35651584, 0, 800000, 600, 22, 20, 4278190080, 800000, 6 } }
};

static const Level *getLevel(int level_idc, int tierFlag)
{
    for (auto &level : levels)
    {
        if (level.level_idc() == level_idc && level.tierFlag == tierFlag) return &level;
    }
    return 0;
}

#endif
