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
#include <boost/program_options.hpp>
#include <memory>


struct InputQueue
{
    InputQueue(const boost::program_options::variables_map &vm);

    // Declare, but do not define, copy constructor to pevent object copying
    InputQueue(const InputQueue&);

    struct References
    {
        std::set<int> positive;
        std::set<int> negative;
        std::set<int> after;
    };

    // Instructions on how to encode a video frame. (consider renaming as Recipe<AccessUnit> ???)
    struct Docket
    {
        std::shared_ptr<PictureWrapper> picture;
        std::shared_ptr<EstimateIntraComplexity> icInfo;
        std::shared_ptr<AdaptiveQuantisation> aqInfo;
        int poc;
        int nut;
        int sliceType;
        int qpOffset;
        double qpFactor;
        int currentGopSize;
        References references;
    };

    // Addes an input picture to the queue.
    void append(std::shared_ptr<PictureWrapper> picture, std::shared_ptr<AdaptiveQuantisation> aqInfo);

    // Called to tell the queue that no more pictures will be appended - this is the end of video input.
    void endOfInput();

    // Returns true if endOfInput() was previously called
    bool eos() const;

    // Retrieves a docket: a packet of work containing a picture and instructions of how to encode it.
    std::shared_ptr<InputQueue::Docket> getDocket();
    void setShotChangeList(std::vector<int>& shotChangeList) { if (shotChangeList.size()) m_shotChangeList.swap(shotChangeList); };

    // Returns true if getDocket() would return a docket.
    bool docketAvailable() const;

protected:
    std::vector<int> m_shotChangeList;

private:
    struct State;
    std::shared_ptr<State> state;
};

#endif
