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

#ifndef INCLUDED_StateSpatial_h
#define INCLUDED_StateSpatial_h

#pragma once

#include "Global.h"
#include "Snake.h"
#include "StateCollocatedMotion.h"
#include "Picture.h"


struct StatePictures;

struct SaoCtuData :
    ValueHolder<SaoTypeIdx>,
    ValueHolder<sao_offset_abs>,
    ValueHolder<sao_offset_sign>,
    ValueHolder<sao_band_position>,
    ValueHolder<SaoEoClass>,
    AccessOperators<SaoCtuData>
    {
        bool operator==(const SaoCtuData &other) const
        {
            return true;
        }
    };


struct Neighbourhood :
    Snake<BlockData>::Cursor
{
    Snake<BlockData>::Pointer snake;
    Snake<BlockData>::Pointer snakeMerge;
    int MinCbLog2SizeYMinus1;

    template <class H>
    void recordMerge(H &h, const prediction_unit &pu)
    {
        const auto parMrgLevel = 1 << h[Log2ParMrgLevel()];
        coding_quadtree *cqt = h;
        const auto nCbS = 1 << cqt->log2CbSize;

        if (parMrgLevel > 4 && nCbS == 8)
        {
            // Record merge data at CU resolution instead
        }
        else
        {
            // dimensions rounded
            const int x = pu.x0 & ~(parMrgLevel - 1);
            const int y = pu.y0 & ~(parMrgLevel - 1);
            int xRight = (pu.x0 + pu.nPbW) & ~(parMrgLevel - 1);
            int yBottom = (pu.y0 + pu.nPbH) & ~(parMrgLevel - 1);

            if (pu.x0 + pu.nPbW == h[pic_width_in_luma_samples()]) xRight = h[pic_width_in_luma_samples()];
            if (pu.y0 + pu.nPbH == h[pic_height_in_luma_samples()]) yBottom = h[pic_height_in_luma_samples()];

            if (x != xRight && y != yBottom)
            {
                for (int xx = x; xx < xRight; xx += (1 << this->MinCbLog2SizeYMinus1))
                {
                    this->snakeMerge.commit(this->snake.at(xx, pu.y0 + pu.nPbH - 1, this->MinCbLog2SizeYMinus1), xx, yBottom - 1, this->MinCbLog2SizeYMinus1);
                }
                for (int yy = y; yy < yBottom; yy += (1 << this->MinCbLog2SizeYMinus1))
                {
                    this->snakeMerge.commit(this->snake.at(pu.x0 + pu.nPbW - 1, yy, this->MinCbLog2SizeYMinus1), xRight - 1, yy, this->MinCbLog2SizeYMinus1);
                }
            }
        }
    }

    // review: this can be optimised because a cqt is always square, power-of-two size and aligned.
    template <class H>
    void recordMerge(H &h, const coding_quadtree &cqt, bool isIntra)
    {
        const auto parMrgLevel = 1 << h[Log2ParMrgLevel()];
        const auto nCbS = 1 << cqt.log2CbSize;

        if (isIntra || (parMrgLevel > 4 && nCbS == 8))
        {
            // dimensions rounded
            const int x = cqt.x0 & ~(parMrgLevel - 1);
            const int y = cqt.y0 & ~(parMrgLevel - 1);
            int xRight = (cqt.x0 + nCbS) & ~(parMrgLevel - 1);
            int yBottom = (cqt.y0 + nCbS) & ~(parMrgLevel - 1);

            if (cqt.x0 + nCbS == h[pic_width_in_luma_samples()]) xRight = h[pic_width_in_luma_samples()];
            if (cqt.y0 + nCbS == h[pic_height_in_luma_samples()]) yBottom = h[pic_height_in_luma_samples()];

            if (x != xRight && y != yBottom)
            {
                for (int xx = x; xx < xRight; xx += (1 << this->MinCbLog2SizeYMinus1))
                {
                    this->snakeMerge.commit(this->snake.at(xx, cqt.y0 + nCbS - 1, this->MinCbLog2SizeYMinus1), xx, yBottom - 1, this->MinCbLog2SizeYMinus1);
                }
                for (int yy = y; yy < yBottom; yy += (1 << this->MinCbLog2SizeYMinus1))
                {
                    this->snakeMerge.commit(this->snake.at(cqt.x0 + nCbS - 1, yy, this->MinCbLog2SizeYMinus1), xRight - 1, yy, this->MinCbLog2SizeYMinus1);
                }
            }
        }
    }

    template <class Rectangular>
    void print(std::ostream &o, Rectangular const &rectangular)
    {
        o << "Neighbourhood at " << rectangular << "\n";
        std::intptr_t diff = this->Snake<BlockData>::Cursor::p - this->snake.origin;
        o << "cursor offset: " << diff << "\n";
        o << "snake:\n";
        this->snake.print(o, rectangular, this->MinCbLog2SizeYMinus1, 2, 2);
        o << "snakeMerge:\n";
        this->snakeMerge.print(o, rectangular, this->MinCbLog2SizeYMinus1, 2, 2);
        o << "\n";
    }
};


template <typename Sample>
struct NeighbourhoodEnc :
    Neighbourhood
{
    typename Snake<Sample>::Pointer snakeIntraReferenceSamples[3];
};

struct StateSpatial
{
    template <class H>
    void resize(H &h, Neighbourhood &neighbourhood)
    {
        Turing::Rectangle rectangle;
        rectangle.x0 = 0;
        rectangle.y0 = 0;
        rectangle.width = h[PicWidthInCtbsY()];
        rectangle.height = h[PicHeightInCtbsY()];

        this->snakeSaoCtuData.resize(rectangle, 0);

        rectangle.width <<= h[CtbLog2SizeY()];
        rectangle.height <<= h[CtbLog2SizeY()];

        neighbourhood.MinCbLog2SizeYMinus1 = h[MinCbLog2SizeY()] - 1;

        this->snakeVector.resize(rectangle, neighbourhood.MinCbLog2SizeYMinus1);
        this->snakeVectorMerge.resize(rectangle, neighbourhood.MinCbLog2SizeYMinus1);

        neighbourhood.snake = this->snakeVector;
        neighbourhood.snakeMerge = this->snakeVectorMerge;
        neighbourhood.p = this->snakeVector.origin;
    }

    Snake<SaoCtuData>::Vector<0, 0> snakeSaoCtuData;
    Snake<BlockData>::Vector<1, 1> snakeVector;
    Snake<BlockData>::Vector<1, 1> snakeVectorMerge;
    Snake<uint8_t>::Vector<32, 32> snakeVectorIntraReferenceSamples8[3];
    Snake<uint16_t>::Vector<32, 32> snakeVectorIntraReferenceSamples16[3];
};

template <class V>
struct AccessSao :
    Access<V, SaoCtuData>
    {
        typedef typename Access<V, SaoCtuData>::Type Type;

        static Type get(V v, Snake<SaoCtuData>::Pointer &snake)
        {
            return Access<V, SaoCtuData>::get(v, const_cast<SaoCtuData&>(snake.at(v.rx, v.ry, 0)));
        }
    };

template <class V, class S>
struct Access<V, S, typename std::enable_if<Accessible<V, SaoCtuData>::value && std::is_base_of<StateSpatial, S>::value>::type>
{
    typedef typename ValueType<V>::Type& Type;
    static Type get(V v, StateSpatial &s)
    {
        return AccessSao<V>::get(v, s.snakeSaoCtuData);
    }
};

struct Unavailable
{
    Unavailable() { this->puData.reset(); }
    PuData puData;
};

template <class H>
const PuData neighbourPuData(H &h, int xN, int yN)
{
    int poc = h[PicOrderCntVal()];

    static const Unavailable unavailable;

    const int maskLow = h[CtbSizeY()] - 1;
    const int maskHigh = ~maskLow;

    const prediction_unit *pu = h;
    const int yCtbCurr = pu->y0 & maskHigh;
    const int yCtbN = yN & maskHigh;

    if (yCtbN > yCtbCurr)
    {
        // neighbour is in next CTU row
        return unavailable.puData;
    }
    else
    {
        // bottom-left sample of prediction unit
        const int xCurr = pu->x0 + pu->nPbW - 1;
        const int yCurr = pu->y0 + pu->nPbH - 1;

        AvailabilityCtu *availabilityCtu = h;

        if (!availabilityCtu->available(xCurr, yCurr, xN, yN, h[CtbLog2SizeY()]))
        {
            return unavailable.puData;
        }

        if (!compareZ(xCurr, yCurr - yCtbCurr, xN, yN - yCtbCurr))
        {
            return unavailable.puData;
        }
    }

    Snake<BlockData>::Cursor *cursor = h;
    return cursor->offset<Anywhere>(xN - pu->x0, yN - pu->y0, h[MinCbLog2SizeY()] - 1);
}

template <int dx, int dy>
struct PuMergeNeighbour
{
    template <class H>
    static void get(H &h, PuData &dest, int xNb, int yNb, int xPb, int yPb)
    {
        assert(!dest.isAvailable());

        assert(dx || xNb >= xPb);
        assert(dy || yNb >= yPb);

        const bool sameX = (xPb >> h[Log2ParMrgLevel()]) == (xNb >> h[Log2ParMrgLevel()]);
        const bool sameY = (yPb >> h[Log2ParMrgLevel()]) == (yNb >> h[Log2ParMrgLevel()]);

        if (sameX && sameY)
        {
            return;
        }

        AvailabilityCtu &availabilityCtu = *static_cast<AvailabilityCtu *>(h);

        Neighbourhood *neighbourhood = h;

        PuData &neighbour = dest;
        if (dx == -1)
        {
            if (dy == -1)
            {
                if (availabilityCtu.available<>(xPb, yPb, xNb, yNb, h[CtbLog2SizeY()]))
                {
                    if (sameX)
                    {
                        neighbour = neighbourhood->snakeMerge.get<Up>(xNb, yNb, neighbourhood->MinCbLog2SizeYMinus1);
                    }
                    else if (sameY)
                    {
                        neighbour = neighbourhood->snakeMerge.get<Left>(xNb, yNb, neighbourhood->MinCbLog2SizeYMinus1);
                    }
                    else
                    {
                        neighbour = neighbourhood->snakeMerge.get<Corner>(xNb, yNb, neighbourhood->MinCbLog2SizeYMinus1);
                    }
                }
            }
            else
            {
                if (availabilityCtu.available<'L'>(xPb, yPb, xNb, yNb, h[CtbLog2SizeY()]))
                {
                    neighbour = neighbourhood->snakeMerge.get<Left>(xNb, yNb, neighbourhood->MinCbLog2SizeYMinus1);
                }
            }
        }
        else
        {
            assert(dy == -1);
            if (availabilityCtu.available<'A'>(xPb, yPb, xNb, yNb, h[CtbLog2SizeY()]))
            {
                neighbour = neighbourhood->snakeMerge.get<Up>(xNb, yNb, neighbourhood->MinCbLog2SizeYMinus1);
            }
        }
    }
};

#endif
