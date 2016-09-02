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

#ifndef INCLUDED_Picture_h
#define INCLUDED_Picture_h

#include "MotionVector.h"
#include <iostream>
#include <stdint.h>
#include <cassert>
#include <memory>


// Pair of pointer and stride that allows access into a 2D plane of samples
template <typename Sample>
struct Raster
{
    typedef Sample Type;

    Raster() { }

    template <class U, typename Stride>
    Raster(U* u, Stride stride, int xOffset = 0, int yOffset = 0) :
        p(u + xOffset + yOffset * stride),
        stride(stride)
    {
    }

    template <class U>
    Raster(const Raster<U> &other, int xOffset = 0, int yOffset = 0) :
        Raster(other.p, other.stride, xOffset, yOffset)
    {
    }

    Sample &operator()(int x, int y)
    {
        return this->p[x + this->stride * y];
    }

    const Sample &operator()(int x, int y) const
    {
        return this->p[x + this->stride * y];
    }

    bool operator==(const Raster<Sample> &other) const
    {
        return this->p == other.p && this->stride == other.stride;
    }

    Raster<Sample> offset(int x, int y)
    {
        return Raster<Sample>(*this, x, y);
    }

    Raster const &operator +=(MotionVector const &mv)
    {
        this->p += mv[0] + mv[1] * this->stride;
        return *this;
    }

    Raster operator +(MotionVector const &mv) const
    {
        Raster temp;
        temp += mv;
        return temp;
    }

    Sample *p;
    std::intptr_t stride;
};


template <class Sample>
void fillRectangle(Raster<Sample> p, Sample value, int width, int height)
{
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            p(x, y) = value;
}


// Representation of a rectangular region of a 2D plane of samples
template <class Sample>
struct Plane :
    Raster<Sample>
{
    Plane() : width(0), height(0) { }

    template <class U>
    Plane(Raster<U> sp, int width, int height) :
        Raster<Sample>(sp),
        width(width),
        height(height)
    {
    }

    bool contiguous() const
    {
        return this->width == this->stride;
    }

    const int width;
    const int height;
};


// Three planar rectangles of samples
template <class Sample>
struct ThreePlanes
{
    typedef Sample Type;

    ThreePlanes(int subWidthC = 1, int subHeightC = 1) :
        subWidthC(subWidthC),
        subHeightC(subHeightC)
    {
        assert(subWidthC == 1 || subWidthC == 2);
        assert(subHeightC == 1 || subHeightC == 2);
    }

    // Copy constuctor with optional crop
    ThreePlanes(ThreePlanes const& src, int left = 0, int top = 0, int right = 0, int bottom = 0) :
        subWidthC(src.subWidthC),
        subHeightC(src.subHeightC)
    {
        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            if (cIdx == 1)
            {
                left /= src.subWidthC;
                right /= src.subWidthC;
                top /= src.subHeightC;
                bottom /= src.subHeightC;
            }

            const int width = src[cIdx].width - left - right;
            const int height = src[cIdx].height - top - bottom;

            Raster<Sample> sp(src[cIdx], left, top);

            this->planes[cIdx].~Plane<Sample>();
            new (&this->planes[cIdx]) Plane<Sample>(sp, width, height);
        }
    }

    ThreePlanes const &operator=(ThreePlanes const &other)
    {
        this->~ThreePlanes();
        new (this) ThreePlanes(other);
        return *this;
    }

    Plane<Sample> &operator[](int cIdx)
    {
        return this->planes[cIdx];
    }

    Plane<Sample> const &operator[](int cIdx) const
    {
        return this->planes[cIdx];
    }

    Raster<Sample> operator()(int x, int y, int cIdx)
    {
        if (cIdx)
        {
            x >>= (this->subWidthC - 1);
            y >>= (this->subHeightC - 1);
        }
        return Raster<Sample>((*this)[cIdx], x, y);
    }

    void write(std::ostream &) const;
    void read(std::istream &);

    const int subWidthC;
    const int subHeightC;

protected:
    Plane<Sample> planes[3];
};


std::ostream &operator<<(std::ostream &, ThreePlanes<uint8_t> &);
std::ostream &operator<<(std::ostream &, ThreePlanes<uint16_t> &);

std::istream &operator >> (std::istream &, ThreePlanes<uint8_t> &);
std::istream &operator >> (std::istream &, ThreePlanes<uint16_t> &);


template <class Sample>
static void writeFields(std::ostream &os, ThreePlanes<Sample> &planestop, ThreePlanes<Sample> &planesbottom)
{
    for (int cIdx = 0; cIdx < 3; ++cIdx)
        for (int y = 0; y < planesbottom[cIdx].height; ++y)
        {
            os.write(reinterpret_cast<const char *>(&planestop[cIdx](0, y)), planesbottom[cIdx].width);
            os.write(reinterpret_cast<const char *>(&planesbottom[cIdx](0, y)), planesbottom[cIdx].width);
        }
}


// Three planar rectangles of samples in allocated memory
template <class Sample>
struct Picture :
    ThreePlanes<Sample>
{
    Picture(int width, int height, int chromaFormat, int paddingX, int paddingY, int alignment);

    Picture(const Picture &) = delete;

    ~Picture();

private:
    Sample *buffers[3 /*cIdx*/];
};


struct PictureWrapper
{
    virtual ~PictureWrapper() { }
    int sampleSize = 0;
    int fieldTB = 0; /* 0 = FRAME, 1 = TOP, 2 = BOTTOM */
    int64_t pts;
};

template <typename Sample>
struct PictureWrap :
    PictureWrapper,
    Picture<Sample>
{
    using Picture<Sample>::Picture;
};

#endif
