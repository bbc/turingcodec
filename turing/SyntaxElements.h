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

#ifndef INCLUDED_SyntaxElements_h
#define INCLUDED_SyntaxElements_h

#pragma once

#include <cassert>


// Pure syntatic sugar - allows us to write e.g. "ae(v)" as in the standard.
static const int v = -1;


struct ue
{
    ue(int) { }
};

struct se
{
    se(int) { }
};

struct u
{
    u(int n) : n(n) { }
    int n;
};

struct uv
{
};

// Handle the "u(v)" case - written uv() in syntax here for type safety and performance
// Specialisations of this struct have static function get() to compute number of bits at runtime dependent on other parameters.
template <class V> struct NumberOfBitsUv;

// Defines the number of bits in a u(v) element V
#define NUMBER_OF_BITS_UV(V, N) \
        \
    template <> struct NumberOfBitsUv<V> \
    { \
    template <class H> static int get(V v, H &h) \
    { \
            return N; \
    } \
    };

// Defines the number of bits in a u(v) element when it depends on the value of another element
#define NUMBER_OF_BITS_MINUS1(V, MINUS1) NUMBER_OF_BITS_UV(V, h[MINUS1()] + 1)

struct i
{
    i(int n) : n(n) { }
    int n;
};

struct b
{
    b(int n) : n(n)
    {
        assert(n % 8 == 0);
    }
    int n;
};

struct f
{
    f(int n) : n(n) { }
    int n;
};

struct ae
{
    ae(int) { }
};

#endif
