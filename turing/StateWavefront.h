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

#ifndef INCLUDED_StateWavefront_h
#define INCLUDED_StateWavefront_h

// State to manage wavefront threading.

#pragma once

#include "Global.h"


struct WavefrontStatus
{
    WavefrontStatus(int w=0, int h=0) :
        previousRowPosition(h + 1),
        width(w),
        height(h)
        {
            this->previousRowPosition[0] = this->width + 1;
        }

    bool done() const
    {
        return this->previousRowPosition.back() == this->width;
    }

    bool operator()(int rx, int ry) const
    {
        return this->previousRowPosition[(ry + 1)] > rx;
    }

    void set(int rx, int ry)
    {
        assert(ry <= this->height);
        assert(rx + 1 > this->previousRowPosition[(ry + 1)]);
        this->previousRowPosition[(ry + 1)] = rx + 1;
    }

private:
    int width;
    int height;
    std::vector<int> previousRowPosition;
};


struct StateWavefront
{
    template <class H>
    void resize(H &h)
    {
        this->encoded = WavefrontStatus(h[PicWidthInCtbsY()], h[PicHeightInCtbsY()]);
        this->deblocked = WavefrontStatus(h[PicWidthInCtbsY()], h[PicHeightInCtbsY()]);
        this->saoed = WavefrontStatus(h[PicWidthInCtbsY()], h[PicHeightInCtbsY()]);

        this->nSubstreamsUnfinished = h[entropy_coding_sync_enabled_flag()] ? h[PicHeightInCtbsY()] : 1;
        this->encodingDone = false;
    }

    WavefrontStatus encoded;
    WavefrontStatus deblocked;
    WavefrontStatus saoed;

    int nSubstreamsUnfinished;
    bool encodingDone;
};

#endif
