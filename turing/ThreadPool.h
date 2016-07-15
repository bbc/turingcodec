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

#ifndef INCLUDED_ThreadPool_h
#define INCLUDED_ThreadPool_h

#pragma once

#include <cstddef>
#include <mutex>
#include <thread>
#include <vector>
#include <deque>
#include <condition_variable>


struct ThreadPool
{
    // Constructs and initialises the thread pool.
    ThreadPool(int numberOfWorkerThreads = 0);

    // Copying this object is not allowed.
    ThreadPool(ThreadPool &);

    // Joins all workers and frees resources.
    ~ThreadPool();

    // Returns number of active threads.
    size_t size() const;

    // Abstract base class for items of thread pool work.
    struct Task
    {
        // Can task be run?
        // Returns true if dependencies are unavailable and run() cannot start yet.
        // If false returned, then calling run() will perform some useful work.
        virtual bool blocked() = 0;

        // Executes a chunk of work, returning either when done or when blocked.
        // Do not call run() unless blocked() returns false.
        // Returns true if job becomes blocked. In this case, wait until blocked() returns false and call run again.
        // Returns false when job is finished. After the job is finished ThreadPool will not make any further calls to the Task object as it may have been deleted.
        virtual bool run() = 0;
    };

    // Enqueues a task for future excution by a pool thread.
    void add(Task &task);

    // Causes thread pool to check immediately whether waiting there is a waiting task that can be run (i.e. is not blocked).
    void nudge();

    std::mutex &mutex();
    void lock();
    void unlock();

    // Worker thread entry point.
    void worker();

private:
    std::vector<std::thread> threads;
    std::deque<Task *> backlog;
    Task *getNextTask();

    std::mutex poolMutex;

    std::condition_variable taskAvailable;
};

#endif
