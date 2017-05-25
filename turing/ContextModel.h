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

#ifndef INCLUDED_ContextModel_h
#define INCLUDED_ContextModel_h
#pragma once

#include <cstdint>
#include <algorithm>


struct ContextModel
{
    uint8_t state;

    template <class T, class U>
    static U Clip3(T min, T max, U value)
    {
        return std::min(std::max(value, min), max);
    }

    ContextModel()
    {
    }

    ContextModel(int SliceQpY , int initValue)
    {
        const int m = (initValue / 16) * 5 - 45;
        const int n = (initValue % 16) * 8 - 16;
        const int initState = Clip3(1, 126, ((m * Clip3(0, 51, SliceQpY) ) >> 4) + n);
        const bool mostProbable = initState >= 64;
        if (mostProbable)
        {
            this->state = static_cast<uint8_t>(2 * (initState - 64) + 1);
        }
        else
        {
            this->state = static_cast<uint8_t>(2 * (63 - initState));
        }
    }

    uint8_t getState() const
    {
        return this->state >> 1;
    }

    uint8_t getMps() const
    {
        return this->state & 0x1;
    }

    void updateMps()
    {
        static const uint8_t nextState[128] =
        {
                2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
                34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
                66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
                82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
                98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113,
                114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 124, 125, 126, 127
        };
        this->state = nextState[this->state];
    }

    void updateLps()
    {
        static const uint8_t nextState[128] =
        {
                1, 0, 0, 1, 2, 3, 4, 5, 4, 5, 8, 9, 8, 9, 10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 18, 19, 22, 23, 22, 23, 24, 25,
                26, 27, 26, 27, 30, 31, 30, 31, 32, 33, 32, 33, 36, 37, 36, 37,
                38, 39, 38, 39, 42, 43, 42, 43, 44, 45, 44, 45, 46, 47, 48, 49,
                48, 49, 50, 51, 52, 53, 52, 53, 54, 55, 54, 55, 56, 57, 58, 59,
                58, 59, 60, 61, 60, 61, 60, 61, 62, 63, 64, 65, 64, 65, 66, 67,
                66, 67, 66, 67, 68, 69, 68, 69, 70, 71, 70, 71, 70, 71, 72, 73,
                72, 73, 72, 73, 74, 75, 74, 75, 74, 75, 76, 77, 76, 77, 126, 127
        };
        this->state = nextState[this->state];
    }

    static uint8_t const nextStateBoth[256];

    inline void update(int32_t mask)
    {
        int i = this->state;
        i ^= mask;
        this->state = nextStateBoth[128 + i];
    }


    bool operator==(ContextModel const& other) const
    {
        return this->state == other.state;
    }

    bool operator!=(ContextModel const& other) const
    {
        return this->state != other.state;
    }
};

#endif
