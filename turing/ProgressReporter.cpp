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

#include "ProgressReporter.h"
#include <ostream>
#include <sstream>
#include <iostream>
// review: try std::chrono here (but beware reports that timers are non compliant on MSVC)
#include <boost/chrono.hpp>

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#undef NOMINMAX
#endif


// PIMPL idiom to avoid polluting callee namespaces with <boost/chrono.hpp>
struct ProgressReporter::State
{
    typedef boost::chrono::process_user_cpu_clock UserClock;

    UserClock::time_point start;

    // Like the "lap" button on a stopwatch. Restart the timer and returns time elapsed since last time elapsed() was called.
    boost::chrono::milliseconds elapsed()
    {
        UserClock::time_point end = UserClock::now();
        const auto duration = end - this->start;
        this->start = end;
        return boost::chrono::duration_cast<boost::chrono::milliseconds>(duration);
    }
};


ProgressReporter::ProgressReporter(std::ostream *os, const char *message, const char *units, size_t total, bool multiline) :
p(new State()),
os(os),
message(message),
units(units),
total(total),
multiline(multiline)
{
#ifdef WIN32
    // Save previous console title
    GetConsoleTitleW(LPWSTR(this->originalConsoleTitle), sizeof(this->originalConsoleTitle) / sizeof(WCHAR));
#endif
}

ProgressReporter::~ProgressReporter()
{
#ifdef WIN32
    // Restore previous console title
    SetConsoleTitleW(LPWSTR(this->originalConsoleTitle));
#endif
}

void ProgressReporter::progress(size_t n)
{
    const auto milliseconds = this->p->elapsed();

    std::ostringstream ss;

    ss << this->message << ": " << n;
    if (this->total > 0) ss << "/" << this->total;
    ss << " " << this->units;

#ifdef WIN32
    SetConsoleTitleA(ss.str().c_str());
#endif

    if (this->os)
    {
        if (this->multiline)
        {
            if (n)
            {
                //				*this->os << milliseconds << "\n";
            }
        }
        else
        {
            *this->os << ss.str() << "\r";
            (*this->os).flush();
        }
    }
}

void ProgressReporter::finish()
{
    if (this->os) *this->os << "\n";
}

