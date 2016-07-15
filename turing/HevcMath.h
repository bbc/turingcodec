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

#ifndef INCLUDED_HevcMath_h
#define INCLUDED_HevcMath_h

#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>


template <class T, class U>
U Clip3(T min, T max, U value)
{
    return std::min(std::max(value, min), max);
}


static int ceilLog2(int x)
{
    int n = 0;
    while (x > (1 << n))
    {
        assert ((1 << n) > 0);
        ++n;
    }
    return n;
}


static int ceilDiv(int quotient, int divisor)
{
    return (quotient + divisor - 1) / divisor;
}


template <class T>
static T Abs(T x)
{
    return x<T(0) ? -x : x;
}


template <class T>
static T Sign(T x)
{
    if (x > T(0)) return T(1);
    if (x < T(0)) return T(-1);
    return T(0);
}


template <class T>
static T Clip1Y(T x, int BitDepthY)
{
    return Clip3(0, (1 << BitDepthY) - 1, x);
}


template <class T>
static T Clip1C(T x, int BitDepthC)
{
    return Clip3(0, (1 << BitDepthC) - 1, x);
}

#endif
