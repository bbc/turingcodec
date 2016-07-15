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

#ifndef INCLUDED_TaskEncodeSubstream_h
#define INCLUDED_TaskEncodeSubstream_h

// Substream encoder top-level task.

#pragma once

#include "StateWavefront.h"
#include "StateEncode.h"
#include "ThreadPool.h"
#include <memory>


template <class> struct Encode;

// Substream encoder top-level task.
// The object will delete itself in TaskEncodeSubstream::run() when it is finished.
template <typename Sample>
struct TaskEncodeSubstream :
    ThreadPool::Task
    {
        TaskEncodeSubstream(PointerTuple<StateEncodePicture2<Sample>, StateEncode> pointers, std::shared_ptr<StateEncodePicture> stateEncodePicture, int begin, int end);

        bool blocked() override;

        bool run() override;

    private:
        bool blockedLock();

        Handler<Encode<void>, Candidate<Sample>, StateEncodeSubstream<Sample>, StateEncodePicture2<Sample>, StateEncode> h;

        StateEncodeSubstream<Sample> stateEncodeSubstream;
        Candidate<Sample> candidate;
        const int begin;
        const int end;
        std::shared_ptr<StateEncodePicture> stateEncodePicture;
        const bool isLastSubstream;
        bool encodingDone;
    };

#endif
