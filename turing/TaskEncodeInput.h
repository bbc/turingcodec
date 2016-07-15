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

// Encoder frame input top-level task.

#ifndef INCLUDED_TaskEncodeInput_h
#define INCLUDED_TaskEncodeInput_h

#pragma once

#include "ThreadPool.h"
#include "Picture.h"
#include "InputQueue.h"
#include "StateEncode.h"


// Encoder frame input top-level task.
template <class H>
struct TaskEncodeInput :
    ThreadPool::Task
    {
        TaskEncodeInput(H &h);

        bool blocked() override;

        bool run() override;

        template <typename Sample>
        void startPictureEncode(StateEncode::Response &response, std::shared_ptr<InputQueue::Docket> docket, H &h);

    private:
        H h;
    };


// Create a new encode output task. The object will be deleted by TaskEncodeInput<H>::run() when it is finished.
template <class H>
TaskEncodeInput<H>& newTaskEncodeInput(H &h)
{
    return *new TaskEncodeInput<H>(h);
}

#endif
