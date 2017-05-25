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

#ifndef INCLUDED_InputQueue_h
#define INCLUDED_InputQueue_h

#pragma once

#include "Picture.h"
#include "EstimateIntraComplexity.h"
#include "AdaptiveQuantisation.h"
#include "SCDetection.h"
#include <memory>
#include <vector>
#include <set>


struct InputQueue
{
    InputQueue(int maxGopN, int maxGopM, bool fieldCoding, bool shotChange, int segmentLength, int baseQP);

    ~InputQueue();

    // Addes an input picture to the queue.
    void append(std::shared_ptr<PictureWrapper> picture, std::shared_ptr<AdaptiveQuantisation> aqInfo);

    // Called to tell the queue that no more pictures will be appended - this is the end of video input.
    void endOfInput();

    void preanalyse();

    // Returns true if endOfInput() was previously called
    bool eos() const;

    struct References
    {
        std::set<int> positive;
        std::set<int> negative;
        std::set<int> after;
    };

    // Instructions on how to encode a video frame. 
    // Review: use StateEncodePicture instead [there is 1:1 relationship]?
    struct Docket
    {
        std::shared_ptr<PictureWrapper> picture;
        std::shared_ptr<EstimateIntraComplexity> icInfo;
        std::shared_ptr<AdaptiveQuantisation> aqInfo;
        int poc;
        int absolutePoc;
        int nut;
        int sliceType;
        int qpOffset;
        double qpFactor;
        int currentGopSize;
        int sopLevel;
        int pocInSop;
        int sopId;
        References references;
        int64_t dts;
        bool isShotChange;
        int hierarchyLevel;
        int numSameHierarchyLevel;
        int intraFramePoc;
        int segmentPoc;
    };

    // Retrieves a docket: a packet of work containing a picture and instructions of how to encode it.
    std::shared_ptr<InputQueue::Docket> getDocket();

private:
    struct State;
    std::unique_ptr<State> state;
};

#endif
