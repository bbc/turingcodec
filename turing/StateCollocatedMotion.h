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

#ifndef INCLUDED_StateCollocateMotion_h
#define INCLUDED_StateCollocateMotion_h

#pragma once

#include "BlockData.h"
#include "Global.h"
#include "Handlers.h"
#include "Picture.h"
#include "Violation.h"
#include <list>
#include <memory>
#include <array>
#include <stdint.h>


template <class H>
static int32_t DiffPicOrderCnt(H &h, int32_t pocA, int32_t pocB)
{
    int32_t const diff = pocA - pocB;
    if (diff < -(1 << 15) || diff >= (1 << 15))
    {
        h(Violation("8.3.2", "DiffPicOrderCnt(picA, picB){%1%} is not in the range of -2^15 to 2^15 - 1, inclusive ") % diff); // CondCheck 8.3.1-C
    }
    return diff;
}


typedef int32_t PicOrderCnt;

struct StateCollocatedMotion
{
    template <class Dpb>
    StateCollocatedMotion(int w, int h, const Dpb &dpb, int32_t poc)
    :
    poc(poc)
    {
        this->stride = (w + 15) / 16;
        const int height = (h + 15) / 16;
        this->entries.resize(this->stride * height);

        for (const auto &i : dpb)
        {
            this->lookup.push_back(LookupEntry());
            this->lookup.back().reference = i->reference;
            this->lookup.back().poc = (*i)[PicOrderCntVal()];
        }
    }

    const PuData &operator()(int x, int y) const
    {
        return this->entries[x/16+ this->stride * y/16];
    }

    template <class Rectangular>
    void fillRectangle(Rectangular const &rectangular, PuData const &puData)
    {
        int const x0 = xPositionOf(rectangular);
        int const y0 = yPositionOf(rectangular);
        int const x1 = x0 + widthOf(rectangular);
        int const y1 = y0 + heightOf(rectangular);

        for (int y = (y0 + 15) / 16; 16 * y < y1; ++y)
        {
            for (int x = (x0 + 15) / 16; 16 * x < x1; ++x)
            {
                this->entries[x + this->stride * y] = puData;
            }
        }
    }

    template <class Rectangular>
    void fillRectangleIntra(const Rectangular &rectangular)
    {
        PuData puData;
        puData.reset();
        fillRectangle(rectangular, puData);
    }

    std::intptr_t stride;
    std::vector<PuData> entries;

    struct LookupEntry
    {
        int poc;
        int reference;
    };

    std::vector<LookupEntry> lookup;

    int getReference(const PuData &puData, int i)
    {
        return this->lookup[puData.getDpbIndex(i)].reference;
    }

    int getPoc(const PuData &puData, int i)
    {
        return this->lookup[puData.getDpbIndex(i)].poc;
    }

    PicOrderCnt const poc;
};

#endif
