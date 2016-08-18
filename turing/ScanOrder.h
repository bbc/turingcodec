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

#ifndef INCLUDED_ScanOrder_h
#define INCLUDED_ScanOrder_h

#pragma once

#include <array>
#include <vector>
#include <cassert>
#include <cstdint>


template <class Scan>
static void diagScan(Scan &scan, int blkSize)
{
    int i = 0;
    int x = 0;
    int y = 0;
    bool stopLoop = false;
    while (!stopLoop)
    {
        while (y >= 0)
        {
            if (x < blkSize  &&  y < blkSize)
            {
                scan[i][0] = x;
                scan[i][1] = y;
                i++;
            }
            y--;
            x++;
        }
        y = x;
        x = 0;
        if (i >= blkSize  *  blkSize)
        {
            stopLoop = true;
        }
    }
}

template <class Scan>
static void horizScan(Scan &horScan, int blkSize)
{
    int i = 0;
    int y = 0;
    while (y < blkSize)
    {
        int x = 0;
        while (x < blkSize)
        {
            horScan[i][0] = x;
            horScan[i][1] = y;
            x++;
            i++;
        }
        y++;
    }
}

template <class Scan>
static void vertScan(Scan &verScan, int blkSize)
{
    int i = 0;
    int x = 0;
    while (x < blkSize)
    {
        int y = 0;
        while (y < blkSize)
        {
            verScan[i][0] = x;
            verScan[i][1] = y;
            y++;
            i++;
        }
        x++;
    }
}


struct ScanPos
{
    union
    {
        struct
        {
            uint8_t x;
            uint8_t y;
        };

        uint16_t both;
    };

    inline uint8_t const &operator[](int sComp) const
    {
        return sComp ? this->y : this->x;
    }

    inline uint8_t &operator[](int sComp)
    {
        return sComp ? this->y : this->x;
    }

    inline ScanPos operator<<=(int n)
    {
        this->both <<= n;
        return *this;
    }

    inline ScanPos operator<<(int n) const
    {
        ScanPos scanPos = *this;
        scanPos <<= n;
        return scanPos;
    }

    inline ScanPos operator+=(ScanPos rhs)
    {
        this->both += rhs.both;
        return *this;
    }

    inline ScanPos operator+(ScanPos rhs) const
    {
        ScanPos scanPos = *this;
        scanPos += rhs;
        return scanPos;
    }
};


typedef std::vector<ScanPos > Scan;


template <int blkSize>
struct ScanArray
{
    ScanPos lookup[3][blkSize * blkSize];
    ScanArray()
    {
        diagScan(this->lookup[0], blkSize);
        horizScan(this->lookup[1], blkSize);
        vertScan(this->lookup[2], blkSize);
    }
};


template <int log2BlockSize>
ScanArray<1 << log2BlockSize> const &scanPos()
{
    static const ScanArray<1 << log2BlockSize> scanArray;
    return scanArray;
}


static inline int ScanOrder(int log2BlockSize, int scanIdx, int sPos, int sComp)
{
    switch (log2BlockSize)
    {
    case 5: return scanPos<5>().lookup[scanIdx][sPos][sComp];
    case 4: return scanPos<4>().lookup[scanIdx][sPos][sComp];
    case 3: return scanPos<3>().lookup[scanIdx][sPos][sComp];
    case 2: return scanPos<2>().lookup[scanIdx][sPos][sComp];
    case 1: return scanPos<1>().lookup[scanIdx][sPos][sComp];
    default: return 0;
    }
}


template <int blkSize>
struct ScanArrayInverse
{
    uint8_t lookup[3][blkSize * blkSize];
    ScanArrayInverse()
    {
        ScanArray<blkSize> scanArray;
        for (int scanIdx = 0; scanIdx < 3; ++scanIdx)
            for (int i = 0; i < blkSize * blkSize; ++i)
            {
                auto const scanPos = scanArray.lookup[scanIdx][i];
                auto const rasterPos = (scanPos.y * blkSize) | scanPos.x;
                this->lookup[scanIdx][rasterPos] = i;
            }
    }
};


template <int log2BlockSize>
static inline uint8_t const &scanPosInverse(int scanIdx, int rasterPos)
{
    static const ScanArrayInverse<1 << log2BlockSize> scanArrayInverse;
    return scanArrayInverse.lookup[scanIdx][rasterPos];
}

#endif
