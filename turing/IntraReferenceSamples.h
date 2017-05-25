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
        if (dx > 0)
            dx *= h[SubWidthC()];
        if (dy > 0)
            dy *= h[SubHeightC()];
    }

    AvailabilityCtu *availabilityCtu = h;
    bool a = availabilityCtu->available(xCurr, yCurr, xCurr + dx, yCurr + dy, h[CtbLog2SizeY()]);

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
    template <class SubC, class H>
    static int log2Step(H &h, int log2TrafoSize)
    {
        // this is an encoder optimisation
        if (Accessible<Struct<struct StateEncode>, H>::value && !h[constrained_intra_pred_flag()])
            return log2TrafoSize;

        return std::min(log2TrafoSize, h[MinCbLog2SizeY()] - (h[SubC()] - 1));
    }

    template <class H>
    IntraReferenceAvailability(H &h, residual_coding rc) :
        mask(0),
        log2StepX(log2Step<SubWidthC>(h, rc.log2TrafoSize)),
        log2StepY(log2Step<SubHeightC>(h, rc.log2TrafoSize)),
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
        auto a = available<Corner>(h, rc.cIdx, rc.x0, rc.y0, -1, -1);

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


struct StateGrey
{
    template <class H>
    StateGrey(H h)
        :
        grey{
            1 << (h[BitDepthY()] - 1),
            1 << (h[BitDepthC()] - 1),
            1 << (h[BitDepthC()] - 1) }
    {
    }
    int const grey[3];
};


template <class Sample>
struct IntraReferenceSamples
    :
    Snake<Sample>::template Array<64, 64, 0, 0>
{
    IntraReferenceSamples()
    {
#ifdef _DEBUG
        this->rc = residual_coding(0, 0, 0, -1);
#endif
        this->resize(Turing::Rectangle{ 0, 0, 64, 64 }, 0);
    }

    IntraReferenceSamples(IntraReferenceSamples &);
    IntraReferenceSamples const &operator=(IntraReferenceSamples const &);

    // specialism used in the encoder to copy reference samples from snake
    void copy(typename Snake<Sample>::Pointer p, int nTbS)
    {
        memcpy(&(*this)(-1, 2 * nTbS - 1), &p(-1, 2 * nTbS - 1), (4 * nTbS + 1) * sizeof(Sample));
    }

    // specialism used in the decoder to copy reference samples from reconstructed plane
    void copy(Raster<Sample> r, int nTbS)
    {
        for (int y = 2 * nTbS - 1; y >= 0; --y)
            (*this)(-1, y) = r(-1, y);

        memcpy(&(*this)(-1, -1), &r(-1, -1), (2 * nTbS + 1) * sizeof(Sample));
    }

    template <class H, class ReferenceSamples>
    void substitute(H &h, ReferenceSamples samples, residual_coding rc) // review: remove rc?
    {
#ifdef _DEBUG
        assert(!(rc == this->rc));
        this->rc = rc;
#endif

        IntraReferenceAvailability const availability(h, rc);

        IntraReferenceSamples &p = *this;

        const int k = rc.log2TrafoSize;
        const int nTbS = 1 << k;

        this->copy(samples, nTbS);

        typename ReferenceSamples::Type value;

        if (availability.any())
        {
            auto const firstX = availability.firstX();
            auto const firstY = availability.firstY();
            value = p(firstX, firstY);
        }
        else
        {
            StateGrey *stateGrey = h;
            value = stateGrey->grey[rc.cIdx];
        }

        const int sizeX = 1 << availability.log2StepX;
        const int sizeY = 1 << availability.log2StepY;

        for (int i = availability.maxY - 1; i >= 0; --i)
        {
            const int y = i << availability.log2StepY;
            if (availability.get(-1, y, Left()))
                value = p(-1, y);
            else
                for (int dy = sizeY - 1; dy >= 0; --dy)
                    p(-1, y + dy) = value;
        }

        if (availability.get(-1, -1, Corner()))
            value = p(-1, -1);
        else
            p(-1, -1) = value;

        for (int i = 0; i < availability.maxX; ++i)
        {
            const int x = i << availability.log2StepX;
            if (availability.get(x, -1, Up()))
                value = p(x + sizeX - 1, -1);
            else
                for (int dx = 0; dx < sizeX; ++dx)
                    p(x + dx, -1) = value;
        }
    }

    static bool availableTopRight(unsigned x, unsigned y)
    {
        unsigned p = ~x;
        unsigned s = (y | p) - 1;
        unsigned q = (x & y & ~s) | (p & s);
        return p > q;
    }

    static bool availableBottomLeft(unsigned x, unsigned y)
    {
        unsigned yNot = ~y;
        unsigned s = (x | yNot) - 1;
        unsigned p = (~x & yNot & ~s) | (x & s);
        return p > x;
    }

    static void fill(Sample *p, Sample sample, int n)
    {
        while (n--)
            *p++ = sample;
    }

    template <class H, class ReferenceSamples>
    void substituteFast(H &h, ReferenceSamples samples, residual_coding rc)
    {
        // this optimised function makes a few assumptions:
        assert(!h[constrained_intra_pred_flag()]);
        assert(!h[tiles_enabled_flag()]);
        assert(h[SliceAddrRs()] == 0);

#ifdef _DEBUG
        assert(!(rc == this->rc));
        this->rc = rc;
#endif

        auto const log2Size = rc.log2TrafoSize + (rc.cIdx ? 1 : 0);
        auto const size = 1 << log2Size;
        auto const sizeMinus1 = size - 1;
        auto const ctbSizeY = h[CtbSizeY()];

        IntraReferenceSamples &p = *this;

        const auto nTbS = 1 << rc.log2TrafoSize;

        this->copy(samples, nTbS);

        bool const left = rc.x0 != 0;
        bool const top = rc.y0 != 0;
        bool const corner = top && left;
        bool const topRight = top && availableTopRight(rc.x0 + sizeMinus1, rc.y0 & (ctbSizeY - 1));
        bool const bottomLeft = left && availableBottomLeft(rc.x0 | ctbSizeY, rc.y0 + sizeMinus1);

        typename ReferenceSamples::Type value;


        if (left)
            value = p(0, nTbS);
        else if (top)
            value = p(0, -1);
        else
        {
            StateGrey *stateGrey = h;
            value = stateGrey->grey[rc.cIdx];
        }

        if (!bottomLeft)
        {
            auto *q = &p(0, 2 * nTbS);
            if (left)
                fill(q, value, nTbS);
            else
                fill(q, value, 2 * nTbS);
        }

        if (left)
            value = p(-1, top ? -1 : 0);

        p(-1, -1) = value;

        if (!topRight)
            if (top)
            {
                value = p(nTbS, 0);
                fill(&p(nTbS, -1), value, nTbS);
            }
            else
                fill(&p(0, -1), value, 2 * nTbS);

#ifdef _DEBUG
        IntraReferenceSamples check;
        check.substitute(h, samples, rc);
        if (!(check == *this))
        {
            std::cerr << "problem at " << rc << "\n";
            check.substitute(h, samples, rc);
            assert(false);
        }
#endif
    }

#ifdef _DEBUG
    // for debugging
    bool operator==(IntraReferenceSamples const &other) const
    {
        for (int y = (2 << rc.log2TrafoSize) - 1; y >= -1; --y)
            if ((*this)(-1, y) != other(-1, y))
                return false;
        for (int x = 0; x < 2 << rc.log2TrafoSize; ++x)
            if ((*this)(x, -1) != other(x, -1))
                return false;
        return true;
    }
#endif

    void filter(const IntraReferenceSamples &p, int strong_intra_smoothing_enabled_flag, int BitDepthY, int nTbS)
    {
        using std::abs;
#ifdef _DEBUG
        assert(!(p.rc == this->rc));
        this->rc = p.rc;
#endif

        IntraReferenceSamples &pF = *this;

        auto const biIntFlag =
            strong_intra_smoothing_enabled_flag == 1 &&
            nTbS == 32 &&
            abs(p(-1, -1) + p(nTbS * 2 - 1, -1) - 2 * p(nTbS - 1, -1)) < (1 << (BitDepthY - 5)) &&
            abs(p(-1, -1) + p(-1, nTbS * 2 - 1) - 2 * p(-1, nTbS - 1)) < (1 << (BitDepthY - 5));

        // review: are we filtering more positions than needed; are we providing enough input samples?

        if (biIntFlag)
        {
            pF(-1, -1) = p(-1, -1);

            for (int y = 0; y <= 62; ++y)
                pF(-1, y) = ((63 - y) * p(-1, -1) + (y + 1) * p(-1, 63) + 32) >> 6;

            pF(-1, 63) = p(-1, 63);

            for (int x = 0; x <= 62; ++x)
                pF(x, -1) = ((63 - x) * p(-1, -1) + (x + 1) * p(63, -1) + 32) >> 6;

            pF(63, -1) = p(63, -1);
        }
        else
        {
            pF(-1, -1) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;

            for (int y = 0; y <= nTbS * 2 - 2; ++y)
                pF(-1, y) = (p(-1, y + 1) + 2 * p(-1, y) + p(-1, y - 1) + 2) >> 2;

            pF(-1, nTbS * 2 - 1) = p(-1, nTbS * 2 - 1);

            for (int x = 0; x <= nTbS * 2 - 2; ++x)
                pF(x, -1) = (p(x - 1, -1) + 2 * p(x, -1) + p(x + 1, -1) + 2) >> 2;

            pF(nTbS * 2 - 1, -1) = p(nTbS * 2 - 1, -1);
        }
    }

#ifdef _DEBUG
    residual_coding rc;
#endif
};

#endif
