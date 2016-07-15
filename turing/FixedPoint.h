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

#ifndef INCLUDED_FixedPoint_h
#define INCLUDED_FixedPoint_h

#pragma once

#include <stdint.h>
#include <ostream>


// Representation of fixed-point integers with parameters for base type and position of the radix point.
template <class U, int k>
struct FixedPoint
{
    template <class V, int j> friend struct FixedPoint;

    U value;

    static const int radixPosition = k;

    static constexpr FixedPoint make(U u, int shift)
    {
        return FixedPoint{ u << (k - shift) };
    }

    void set(U u, int shift)
    {
        this->value = u << (k - shift);
    }

    void set(double d)
    {
        this->value = static_cast<U>(d * (1 << k) + 0.5);
        if (d == std::numeric_limits<double>::max())
        {
            // review: remove
            this->value = std::numeric_limits<U>::max();
        }
    }

    double asDouble() const { return double(value) / double(1 << k); }
    FixedPoint operator+=(FixedPoint other) { this->value += other.value; return *this; }
    FixedPoint operator-=(FixedPoint other) { this->value -= other.value; return *this; }
    friend FixedPoint operator+(FixedPoint x, FixedPoint y) { return x += y; }
    friend FixedPoint operator-(FixedPoint x, FixedPoint y) { return x -= y; }
    friend FixedPoint operator-(FixedPoint x) { x.value = -x.value; return x; }
    FixedPoint operator>>=(int n) { this->value >>= n; return *this; }
    FixedPoint operator<<=(int n) { this->value <<= n; return *this; }
    FixedPoint operator>>(int n) { auto temp = *this; return temp >>= n; }
    FixedPoint operator<<(int n) { auto temp = *this; return temp <<= n; }
    bool operator<(FixedPoint other) const { return this->value < other.value; }
    bool operator<=(FixedPoint other) const { return this->value <= other.value; }
    bool operator>(FixedPoint other) const { return this->value > other.value; }
    bool operator>=(FixedPoint other) const { return this->value >= other.value; }
    bool operator==(FixedPoint other) const { return this->value == other.value; }
    bool operator!=(FixedPoint other) const { return this->value != other.value; }
    FixedPoint operator++() { this->value += 1 << k; return *this; }
};

template <int k> FixedPoint<int64_t, k> operator*(FixedPoint<int32_t, k> x, int32_t y) { return FixedPoint < int64_t, k>{x.value * int64_t(y)}; }
template <int k> FixedPoint<int64_t, k> operator*(int32_t x, FixedPoint<int32_t, k> y) { return FixedPoint < int64_t, k>{int64_t(x) * y.value}; }

template <class U, int k>
static std::ostream& operator<<(std::ostream& o, FixedPoint<U, k> const x)
{
    return o << x.asDouble();
}

namespace std {

    template <class U, int k>
    class numeric_limits < FixedPoint<U, k> >
    {
    public:
        static FixedPoint<U, k> max()
        {
            return FixedPoint<U, k> { std::numeric_limits<U>::max() };
        }
    };

}

#endif
