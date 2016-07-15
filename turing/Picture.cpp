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

#include "Picture.h"
#include "Memory.h"
#include "HevcTypes.h"


template <typename Sample> void ThreePlanes<Sample>::read(std::istream &is)
{
    auto &planes = *this;
    for (int cIdx = 0; cIdx < 3; ++cIdx)
    {
        if (planes[cIdx].stride == planes[cIdx].width)
        {
            is.read(reinterpret_cast<char *>(&planes[cIdx](0, 0)), planes[cIdx].width * planes[cIdx].height * sizeof(Sample));
        }
        else for (int y = 0; y < planes[cIdx].height; ++y)
        {
            is.read(reinterpret_cast<char *>(&planes[cIdx](0, y)), planes[cIdx].width * sizeof(Sample));
        }
    }
}


template <typename Sample> void ThreePlanes<Sample>::write(std::ostream &os) const
        {
    auto &planes = *this;
    for (int cIdx = 0; cIdx < 3; ++cIdx)
    {
        if (planes[cIdx].stride == planes[cIdx].width)
        {
            os.write(reinterpret_cast<const char *>(&planes[cIdx](0, 0)), planes[cIdx].width * planes[cIdx].height * sizeof(Sample));
        }
        else for (int y = 0; y < planes[cIdx].height; ++y)
        {
            os.write(reinterpret_cast<const char *>(&planes[cIdx](0, y)), planes[cIdx].width * sizeof(Sample));
        }
    }
        }


std::ostream &operator<<(std::ostream &os, ThreePlanes<uint8_t> &planes)
{
    planes.write(os);
    return os;
}


std::ostream &operator<<(std::ostream &os, ThreePlanes<uint16_t> &planes)
{
    planes.write(os);
    return os;
}


std::istream &operator >> (std::istream &is, ThreePlanes<uint8_t> &planes)
{
    planes.read(is);
    return is;
}


std::istream &operator >> (std::istream &is, ThreePlanes<uint16_t> &planes)
{
    planes.read(is);
    return is;
}




template <class Sample>
Picture<Sample>::Picture(int width, int height, int chromaFormat, int paddingX, int paddingY, int alignment)
:
ThreePlanes<Sample>(::subWidthC(chromaFormat), ::subHeightC(chromaFormat))
{
    for (int cIdx = 0; cIdx < 3; ++cIdx)
    {
        if (cIdx == 1)
        {
            width /= this->subWidthC;
            height /= this->subHeightC;
            paddingX /= this->subWidthC;
            paddingY /= this->subHeightC;
        }

        const int n = alignment / sizeof(Sample);

        int paddingLeft = paddingX;
        int paddingRight = paddingX;

        // Additional padding may also be required between video lines
        // to ensure the first sample of every picture is aligned.
        int stride = paddingLeft + width + paddingRight;
        if (stride % n)
        {
            paddingRight += n - stride % n;
            stride = paddingLeft + width + paddingRight;
        }

        // Some extra space at the front of the allocated memory may be necessary to ensure that
        // the first sample of the picture is aligned.
        const int paddingPicture = (paddingLeft % n) ? (n - paddingLeft % n) : 0;

        Memory::allocate(this->buffers[cIdx], paddingPicture + stride * (paddingY + height + paddingY), alignment);

        Raster<Sample> sp(this->buffers[cIdx], stride, paddingPicture + paddingLeft, paddingY);

        this->planes[cIdx].~Plane<Sample>();
        new (&this->planes[cIdx]) Plane<Sample>(sp, width, height);
    }
}


template <class Sample>
Picture<Sample>::~Picture()
{
    for (int cIdx = 0; cIdx < 3; ++cIdx)
    {
        Memory::free(this->buffers[cIdx]);
    }
}


// explicit specialisations
template struct Picture<uint8_t>;
template struct Picture<uint16_t>;
template struct Picture<int16_t>;
