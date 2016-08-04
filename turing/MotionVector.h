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

#ifndef INCLUDED_MotionVector_h
#define INCLUDED_MotionVector_h

#pragma once

#include <iostream>
#include <stdint.h>


// Representation of a motion vector, or motion vector difference having 16-bit signed components
struct MotionVector
{
    typedef int16_t ComponentType;

    ComponentType &operator[](int i)
    {
        return i ? this->y : this->x;
    }

    const ComponentType &operator[](int i) const
    {
        return i ? this->y : this->x;
    }

    bool operator==(MotionVector const &other) const
    {
        return this->x == other.x && this->y == other.y;
    }

    bool operator!=(MotionVector const &other) const
    {
        return this->x != other.x || this->y != other.y;
    }

    MotionVector const &operator+=(MotionVector const &other)
    {
        this->x += other.x;
        this->y += other.y;
        return *this;
    }

    MotionVector const &operator+=(ComponentType i)
    {
        this->x += i;
        this->y += i;
        return *this;
    }

    template <class Other>
    MotionVector operator+(Other const &other) const
    {
        MotionVector result = *this;
        result += other;
        return result;
    }

    MotionVector const &operator-=(MotionVector const &other)
    {
        this->x -= other.x;
        this->y -= other.y;
        return *this;
    }

    MotionVector const &operator-=(ComponentType i)
    {
        this->x -= i;
        this->y -= i;
        return *this;
    }

    template <class Other>
    MotionVector operator-(Other const &other) const
    {
        MotionVector result = *this;
        result -= other;
        return result;
    }

    MotionVector const &operator*=(ComponentType i)
    {
        this->x *= i;
        this->y *= i;
        return *this;
    }

    template <class Other>
    MotionVector operator*(Other const &other) const
    {
        MotionVector result = *this;
        result *= other;
        return result;
    }

    template <class Other>
    MotionVector operator/=(Other const &other)
    {
        this->x /= other;
        this->y /= other;
        return *this;
    }

    template <class Other>
    MotionVector operator/(Other const &other) const
    {
        MotionVector result = *this;
        result /= other;
        return result;
    }

    template <class Other>
    MotionVector operator%=(Other const &other)
    {
        this->x %= other;
        this->y %= other;
        return *this;
    }

    template <class Other>
    MotionVector operator%(Other const &other) const
    {
        MotionVector result = *this;
        result %= other;
        return result;
    }

    MotionVector const &operator>>=(int n)
    {
        this->x >>= n;
        this->y >>= n;
        return *this;
    }

    MotionVector operator>>(int n) const
    {
        MotionVector result = *this;
        result >>= n;
        return result;
    }

    MotionVector const &operator<<=(int n)
    {
        this->x <<= n;
        this->y <<= n;
        return *this;
    }

    MotionVector operator<<(int n) const
    {
        MotionVector result = *this;
        result <<= n;
        return result;
    }

    MotionVector const &operator&=(ComponentType i)
    {
        this->x &= i;
        this->y &= i;
        return *this;
    }

    MotionVector operator&(ComponentType i) const
    {
        MotionVector result = *this;
        result &= i;
        return result;
    }

    MotionVector const &operator|=(ComponentType i)
    {
        this->x |= i;
        this->y |= i;
        return *this;
    }

    MotionVector operator|(ComponentType i) const
    {
        MotionVector result = *this;
        result |= i;
        return result;
    }

    MotionVector operator-() const
    {
        return MotionVector{ (int16_t)(-this->x), (int16_t)(-this->y) };
    }

    explicit operator bool() const
    {
        return this->x || this->y;
    }

    friend MotionVector operator+(int i, MotionVector mv)
    {
        mv += i;
        return mv;
    }

    friend MotionVector operator-(int i, MotionVector mv)
    {
        mv -= i;
        return mv;
    }

    friend MotionVector operator*(int i, MotionVector mv)
    {
        mv *= i;
        return mv;
    }

    friend std::ostream &operator<<(std::ostream &o, MotionVector mv)
    {
        return o << "(" << mv.x << ", " << mv.y << ")";
    }

    ComponentType x, y;
};

#endif
