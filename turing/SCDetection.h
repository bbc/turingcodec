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

#ifndef INCLUDED_SCDetection_h
#define INCLUDED_SCDetection_h

#include <vector>
#include <memory>


#pragma once

class ShotChangeDetection {
public:
    ShotChangeDetection(const char* filename, int fileBitDepth, int width,
                        int height, int frameSize, int frameSkip, int frameNum) :
                            filename_(filename), fileBitDepth_(fileBitDepth), width_(width)
, height_(height), frameSize_(frameSize), frameSkip_(frameSkip),
frameNum_(frameNum) {}
    ~ShotChangeDetection() {}

    using FramePtr = std::shared_ptr<std::vector<unsigned char>>;
    double getLikelihood(FramePtr prev, FramePtr cur);
    void processSeq(std::vector<int>& shotChangeList);
private:
    const char* filename_;
    int fileBitDepth_;
    int width_;
    int height_;
    int frameSize_;
    int frameSkip_;
    int frameNum_;

private:
    // no copy or assigment
    ShotChangeDetection& operator= (const ShotChangeDetection&);
    ShotChangeDetection(const ShotChangeDetection&);
};

#endif
