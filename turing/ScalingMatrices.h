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

#ifndef INCLUDED_ScalingMatrices_h
#define INCLUDED_ScalingMatrices_h

#pragma once


struct ScalingMatrices
{
    typedef int Type;

    void initialise(struct ScalingListState *scalingListState);

    static int matrixId(int sizeId, int cIdx, int cuPredMode);

    Type* getMatrix(int sizeId, int matrixId);

    Type &operator()(int sizeId, int matrixId, int x, int y);

private:

    // review: std::array<> here?
    Type matrix4x4[6][4*4];
    Type matrix8x8[6][8*8];
    Type matrix16x16[6][16*16];
    Type matrix32x32[2][32*32];
};

#endif
