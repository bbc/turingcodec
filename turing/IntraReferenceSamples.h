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

#ifndef INCLUDED_IntraReference_h
#define INCLUDED_IntraReference_h

#pragma once

#include "Global.h"
#include "Picture.h"
#include "BlockData.h"
#include "Snake.h"
#include <array>


// availability for intra prediction
template <class Direction, class H>
bool available(H &h, int cIdx, int xCurr, int yCurr, int dx, int dy)
{
    if (cIdx)
    {
        if (dx > 0) dx *= h[SubWidthC()];
        if (dy > 0) dy *= h[SubHeightC()];
    }

    bool a = h[availableX(xCurr, yCurr, xCurr + dx, yCurr + dy)];

    if (a && h[constrained_intra_pred_flag()])
    {
        const coding_quadtree *cqt = h;
        if ((xCurr + dx < cqt->x0 || yCurr + dy < cqt->y0) && (xCurr + dx < cqt->x0 + (1 << cqt->log2CbSize)) && (yCurr + dy < cqt->y0 + (1 << cqt->log2CbSize)))
        {
            Snake<BlockData>::Array<32, 32, 0, 0> *snakeArrayCuCip = h;
            snakeArrayCuCip->checkPosition(before, *cqt);
            a = snakeArrayCuCip->at(xCurr + dx, yCurr + dy, h[MinCbLog2SizeY()] - 1).CuPredMode() == MODE_INTRA;
        }
        else
        {
            a = h[CuPredMode(xCurr + dx, yCurr + dy)] == MODE_INTRA;
        }
    }

    return a;
}


struct IntraReferenceAvailability
{
    template <class H>
    IntraReferenceAvailability(H &h, residual_coding rc) :
    mask(0),
    log2StepX(std::min(rc.log2TrafoSize, h[MinCbLog2SizeY()] - (h[SubWidthC()] - 1))),
    log2StepY(std::min(rc.log2TrafoSize, h[MinCbLog2SizeY()] - (h[SubHeightC()] - 1))),
    maxX(2 << (rc.log2TrafoSize - log2StepX)),
    maxY(2 << (rc.log2TrafoSize - log2StepY))
    {
        for (int i = this->maxX - 1; i >= 0; --i)
        {
            const bool a = available<Up>(h, rc.cIdx, rc.x0, rc.y0, i << log2StepX, -1);
            this->mask |= uint64_t(a) << (17 + i);
            if (a) this->xi = i;
        }

        mask |= uint64_t(available<Corner>(h, rc.cIdx, rc.x0, rc.y0, -1, -1)) << 16;

        for (int i = 0; i < this->maxY; ++i)
        {
            const bool a = available<Left>(h, rc.cIdx, rc.x0, rc.y0, -1, i << log2StepY);
            this->mask |= uint64_t(a) << (15 - i);
            if (a) this->yi = i;
        }
    }

    bool get(int dx, int dy, Left) const
    {
        assert(dx == -1);
        return !!(this->mask & (uint64_t(1) << (15 - (dy >> log2StepY))));
    }

    bool get(int dx, int dy, Corner) const
    {
        assert(dx == -1);
        assert(dx == -1);
        return !!(this->mask & (uint64_t(1) << 16));
    }

    bool get(int dx, int dy, Up) const
    {
        assert(dy == -1);
        return !!(this->mask & (uint64_t(1) << (17 + (dx >> log2StepX))));
    }

    bool any() const
    {
        return !!this->mask;
    }

    int firstX() const
    {
        assert(this->any());
        if (this->mask & 0x1ffff) return -1;
        return this->xi << log2StepX;
    }

    int firstY() const
    {
        assert(this->any());
        if (!(this->mask & 0xffff)) return -1;
        return ((this->yi + 1) << log2StepY) - 1;
    }

    const int log2StepX;
    const int log2StepY;
    const int maxX;
    const int maxY;
private:
    uint64_t mask;
    int xi, yi;
};

template <class Sample>
struct IntraReferenceSamples
{
    IntraReferenceSamples()
    {
        this->rc.cIdx = -1;
    }

    template <class ReferenceSamples>
    void substitute(const IntraReferenceAvailability &availability, ReferenceSamples samples, residual_coding rc, int bitDepth)
    {
        if (rc == this->rc) return;
        this->rc = rc;

        IntraReferenceSamples &p = *this;

        const int k = rc.log2TrafoSize;
        const int nTbS = 1 << k;

        typename ReferenceSamples::Type value = 1 << (bitDepth - 1);

        if (availability.any())
        {
            auto const firstX = availability.firstX();
            auto const firstY = availability.firstY();
            if (firstY >= 0)
            {
                assert(firstX == -1);
                value = samples(-1, firstY);
            }
            else if (firstX >= 0)
            {
                assert(firstY == -1);
                value = samples(firstX, -1);
            }
            else
            {
                assert(firstX == -1);
                assert(firstY == -1);
                value = samples(-1, -1);
            }
        }

        const int sizeX = 1 << availability.log2StepX;
        const int sizeY = 1 << availability.log2StepY;

        for (int i = availability.maxY - 1; i >= 0; --i)
        {
            const int y = i << availability.log2StepY;
            if (availability.get(-1, y, Left()))
            {
                for (int dy = sizeY - 1; dy >= 0; --dy)
                {
                    p(-1, y + dy) = value = samples(-1, y + dy);
                }
            }
            else
            {
                for (int dy = sizeY - 1; dy >= 0; --dy)
                {
                    p(-1, y + dy) = value;
                }
            }
        }

        if (availability.get(-1, -1, Corner()))
        {
            value = samples(-1, -1);
        }
        p(-1, -1) = value;

        for (int i = 0; i < availability.maxX; ++i)
        {
            const int x = i << availability.log2StepX;
            if (availability.get(x, -1, Up()))
            {
                for (int dx = 0; dx < sizeX; ++dx)
                {
                    p(x + dx, -1) = value = samples(x + dx, -1);
                }
            }
            else
            {
                for (int dx = 0; dx < sizeX; ++dx)
                {
                    p(x + dx, -1) = value;
                }
            }
        }
    }

    void filter(const IntraReferenceSamples &p, int strong_intra_smoothing_enabled_flag, int BitDepthY)
    {
        if (p.rc == this->rc) return;
        this->rc = p.rc;

        IntraReferenceSamples &pF = *this;

        const int k = rc.log2TrafoSize;
        const int nTbS = 1 << k;

        const int biIntFlag =
                (strong_intra_smoothing_enabled_flag == 1
                        && nTbS == 32
                        && std::abs(p(-1, -1) + p(nTbS * 2 - 1, -1) - 2 * p(nTbS - 1, -1)) < (1 << (BitDepthY - 5))
                        && std::abs(p(-1, -1) + p(-1, nTbS * 2 - 1) - 2 * p(-1, nTbS - 1)) < (1 << (BitDepthY - 5))) ? 1 : 0;

        if (biIntFlag == 1)
        {
            pF(-1, -1) = p(-1, -1);
            for (int y = 0; y <= 62; ++y)
            {
                pF(-1, y) = ((63 - y) * p(-1, -1) + (y + 1) * p(-1, 63) + 32) >> 6;
            }
            pF(-1, 63) = p(-1, 63);
            for (int x = 0; x <= 62; ++x)
            {
                pF(x, -1) = ((63 - x) * p(-1, -1) + (x + 1) * p(63, -1) + 32) >> 6;
            }
            pF(63, -1) = p(63, -1);
        }
        else
        {
            pF(-1, -1) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;
            for (int y = 0; y <= nTbS * 2 - 2; ++y)
            {
                pF(-1, y) = (p(-1, y + 1) + 2 * p(-1, y) + p(-1, y - 1) + 2) >> 2;
            }
            pF(-1, nTbS * 2 - 1) = p(-1, nTbS * 2 - 1);
            for (int x = 0; x <= nTbS * 2 - 2; ++x)
            {
                pF(x, -1) = (p(x - 1, -1) + 2 * p(x, -1) + p(x + 1, -1) + 2) >> 2;
            }
            pF(nTbS * 2 - 1, -1) = p(nTbS * 2 - 1, -1);
        }
    }

    Sample &operator()(int x, int y)
    {
        if (x < 0 && y < 0)
        {
        }
        else if (y < 0)
        {
            assert(y == -1);
            return this->buffer[64 + x];
        }
        else if (x < 0)
        {
            assert(x == -1);
            return this->buffer[63 - y];
        }

        assert(x == -1);
        assert(y == -1);
        return this->corner;
    }

    const Sample &operator()(int x, int y) const
    {
        if (x < 0 && y < 0)
        {
        }
        else if (y < 0)
        {
            assert(y == -1);
            return this->buffer[64 + x];
        }
        else if (x < 0)
        {
            assert(x == -1);
            return this->buffer[63 - y];
        }

        assert(x == -1);
        assert(y == -1);
        return this->corner;
    }

    std::array<Sample, 129> buffer;
    Sample corner;
    residual_coding rc;
};

#endif
