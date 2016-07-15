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

#include "ThreadPool.h"


ThreadPool::ThreadPool(int n)
{
    if (n == 0) n = std::thread::hardware_concurrency();
    if (n == 0) n = 1;

    for (int i = 0; i < n; ++i)
    {
        this->threads.push_back(std::thread(&ThreadPool::worker, this));
    }
}


ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(this->poolMutex);

        for (unsigned i = 0; i < this->size(); ++i)
        {
            this->backlog.push_back(0);
        }
    }

    this->taskAvailable.notify_all();

    for (auto &thread : this->threads)
    {
        thread.join();
    }
}


void ThreadPool::add(Task &task)
{
    std::unique_lock<std::mutex> lock(this->poolMutex);

    this->backlog.push_back(&task);

    this->nudge();
}


ThreadPool::Task *ThreadPool::getNextTask()
{
    std::unique_lock<std::mutex> lock(this->poolMutex);

    while (true)
    {
        for (auto i = this->backlog.begin(); i != this->backlog.end(); ++i)
        {
            ThreadPool::Task *task = *i;
            if (!task || !task->blocked())
            {
                this->backlog.erase(i);
                return task;
            }
        }

        this->taskAvailable.wait(lock);
    }
}


void ThreadPool::worker()
{
    while (true)
    {
        Task *task = getNextTask();

        if (!task) break;

        const bool blocked = task->run();
        if (blocked)
        {
            std::unique_lock<std::mutex> lock(this->poolMutex);

            this->backlog.push_front(task);
        }
    }
}


void ThreadPool::lock()
{
    this->poolMutex.lock();
}


void ThreadPool::unlock()
{
    this->poolMutex.unlock();
}


std::mutex &ThreadPool::mutex()
{
    return this->poolMutex;
}


void ThreadPool::nudge()
{
    this->taskAvailable.notify_all();
}


size_t ThreadPool::size() const
{
    return this->threads.size();
}
