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

#ifndef INCLUDED_ProgressReporter_h
#define INCLUDED_ProgressReporter_h

#include <ostream>
#include <memory>


struct ProgressReporter
{
    ProgressReporter(std::ostream *os, const char *message = "processed ", const char *units = " pictures", size_t total = 0, bool multiline=false);
    ~ProgressReporter();
    void progress(size_t n);
    void finish();
private:
    struct State;
    std::unique_ptr<State> p;
    const size_t total;
    const char *message;
    const char *units;
    std::ostream *os;
    char originalConsoleTitle[1000];
    const bool multiline;
};

#endif
