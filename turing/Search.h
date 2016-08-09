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

#ifndef INCLUDED_Search_h
#define INCLUDED_Search_h

#pragma once

#include "Global.h"
#include "StateValues.h"


template <class F> struct Search :
    Syntax<F>
    {
    };

template <unsigned u> struct Mode
{
    static constexpr unsigned value = u;
};

template <unsigned u> struct Search<Mode<u>> :
    ValueConst<transquant_bypass_enabled_flag, 0>
    {
    };

template <>
struct Search<sao>
{
    template <class H> static void go(const sao &sao, H &h);
};

template <>
struct Search<coding_quadtree>
{
    template <class H> static void go(const coding_quadtree &cqt, H &h);
};

template <class Direction>
struct Search<Deleted<coding_quadtree, Direction>>
{
    template <class H> static void go(const Deleted<coding_quadtree, Direction> &cqt, H &h);
};

template <class F> struct SearchMerge2Nx2N :
    Search<F>,
    ValueConst<Neighbouring<PartMode, Current>, PART_2Nx2N>
    {
    };

#endif
