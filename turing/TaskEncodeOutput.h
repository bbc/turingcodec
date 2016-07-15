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

// Encoder output task responsible for writing all bitstream.

#ifndef INCLUDED_TaskEncodeOutput_h
#define INCLUDED_TaskEncodeOutput_h

#pragma once

#include "ThreadPool.h"
#include "turing.h"


// Encoder output task responsible for writing all bitstream.
// Blocks when it needs more data from TaskEncodeSubstream.
template <class H>
struct TaskEncodeOutput :
    ThreadPool::Task
    {
        TaskEncodeOutput(H &h);

        bool blocked() override;

        bool run() override;

    private:
        H h;
    };


// Create a new substream encode task. The object will delete itself in TaskEncodeOutput<H>::run() when it is finished.
template <class H>
TaskEncodeOutput<H>& newTaskEncodeOutput(H &h)
{
    return *new TaskEncodeOutput<H>(h);
}

#endif
