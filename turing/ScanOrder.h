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
    int x  =  0;
    int y  =  0;
    bool stopLoop  =  false;
    while( !stopLoop )
    {
        while( y  >=  0 )
        {
            if( x  <  blkSize  &&  y  <  blkSize )
            {
                scan[ i ][ 0 ]  =  x;
                scan[ i ][ 1 ]  =  y;
                i++;
            }
            y--;
            x++;
        }
        y  =  x;
        x  =  0;
        if( i  >=  blkSize  *  blkSize )
        {
            stopLoop  =  true;
        }
    }
}

template <class Scan>
static void horizScan(Scan &horScan, int blkSize)
{
    int i  =  0;
    int y  =  0;
    while( y  < blkSize )
    {
        int x  =  0;
        while( x  <  blkSize )
        {
            horScan[ i ][ 0 ]  = x;
            horScan[ i ][ 1 ]  = y;
            x++;
            i++;
        }
        y++;
    }
}

template <class Scan>
static void vertScan(Scan &verScan, int blkSize)
{
    int i  =  0;
    int x  =  0;
    while( x  < blkSize )
    {
        int y  =  0;
        while( y  <  blkSize )
        {
            verScan[ i ][ 0 ]  = x;
            verScan[ i ][ 1 ]  = y;
            y++;
            i++;
        }
        x++;
    }
}


typedef std::vector<std::array<uint8_t, 2> > Scan;


template <int blkSize>
struct ScanArray
{
    typedef short Type[3][blkSize*blkSize][2];
    Type buffer;
    ScanArray()
    {
        diagScan(this->buffer[0], blkSize);
        horizScan(this->buffer[1], blkSize);
        vertScan(this->buffer[2], blkSize);
    }
};

static const ScanArray<32> scanArray32;
static const ScanArray<16> scanArray16;
static const ScanArray<8> scanArray8;
static const ScanArray<4> scanArray4;
static const ScanArray<2> scanArray2;

static inline int ScanOrder(int log2BlockSize, int scanIdx, int sPos, int sComp)
{
    switch (log2BlockSize)
    {
        case 5: return scanArray32.buffer[scanIdx][sPos][sComp];
        case 4: return scanArray16.buffer[scanIdx][sPos][sComp];
        case 3: return scanArray8.buffer[scanIdx][sPos][sComp];
        case 2: return scanArray4.buffer[scanIdx][sPos][sComp];
        case 1: return scanArray2.buffer[scanIdx][sPos][sComp];
        case 0: return 0;
    }
    assert(0); return 0;
}

#endif
