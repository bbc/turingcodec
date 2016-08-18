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

#ifndef INCLUDED_Padding_h
#define INCLUDED_Padding_h

#pragma once

#include "Padding.h"
#include <cstring>
#include <cstring>
#include <algorithm>


namespace Padding {

template <class T>
void padImage(T *p, int width, int height, int stride, int pad)
{
    // Horizontal pad
    for (int y=0; y<height; ++y)
    {
        T *left = p + y*stride;
        std::fill(left - pad, left, *left);

        T *right = p + y*stride+width;
        std::fill(right, right + pad, *(right-1));
    }

    // Vertical pad
    const size_t size = (pad + width + pad) * sizeof(T);

    T *top = p - pad;
    for (int i=1; i<=pad; ++i)
        std::memcpy(top - i*stride, top, size);

    T *bottom = p + (height-1)*stride - pad;
    for (int i=1; i<=pad; ++i)
        std::memcpy(bottom + i*stride, bottom, size);
}


template <class T>
void padBlock(T *p, int width, int height, std::intptr_t stride, int pad, bool top, bool bottom, bool left, bool right)
{
    size_t widthOverall = width;
    T *pVerticalPad = p;

    if (left)
    {
        for (int y = 0; y < height; ++y)
        {
            T *l = p + y*stride;
            std::fill(l - pad, l, *l);
        }
        pVerticalPad -= pad;
        widthOverall += pad;
    }

    if (right)
    {
        for (int y = 0; y < height; ++y)
        {
            T *r = p + width + y*stride;
            std::fill(r, r + pad, *(r - 1));
        }

        widthOverall += pad;
    }

    if (top)
        for (int i = 1; i <= pad; ++i)
            std::memcpy(pVerticalPad - i*stride, pVerticalPad, widthOverall * sizeof(T));

    pVerticalPad += (height - 1) * stride;

    if (bottom)
        for (int i = 1; i <= pad; ++i)
            std::memcpy(pVerticalPad + i*stride, pVerticalPad, widthOverall * sizeof(T));
}

template <class Sample>
static void padPicture(Picture<Sample> &picture)
{
    const int pad = 80; // todo: be more intelligent
    Padding::padImage<Sample>(&picture[0](0, 0), picture[0].width, picture[0].height, (int)picture[0].stride, pad);
    Padding::padImage<Sample>(&picture[1](0, 0), picture[1].width, picture[1].height, (int)picture[1].stride, pad / 2);
    Padding::padImage<Sample>(&picture[2](0, 0), picture[2].width, picture[2].height, (int)picture[2].stride, pad / 2);
}

}

#endif
