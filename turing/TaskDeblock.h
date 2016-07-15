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

// Deblocking top-level task.
// Perhaps misleadingly named: this task also performs picture boundary padding (for out-of picture MVs)
// and will also handle SAO when fully implemented.

#ifndef INCLUDED_TaskDeblock_h
#define INCLUDED_TaskDeblock_h

#pragma once

#include "StateWavefront.h"
#include "ThreadPool.h"
#include <memory>


template <class H>
struct TaskDeblock :
    ThreadPool::Task
    {
        TaskDeblock(H &h, std::shared_ptr<struct StateEncodePicture> stateEncodePicture, int begin, int end);

        bool blocked() override;

        bool run() override;

    private:
        H h;
        std::shared_ptr<struct StateEncodePicture> stateEncodePicture;
        const int begin;
        const int end;
        WavefrontStatus *syncIn;
        WavefrontStatus *syncOut;
        int position; // CTU address, raster scan
    };


// Create a new deblock task. The object will delete itself in TaskEncodeOutput<H>::run() when it is finished.
template <class H>
TaskDeblock<H> &newTaskDeblock(H &h, std::shared_ptr<struct StateEncodePicture> stateEncodePicture, int begin, int end)
{
    return *new TaskDeblock<H>(h, stateEncodePicture, begin, end);
}

#endif
