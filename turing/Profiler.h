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

#ifndef INCLUDED_Profiler_h
#define INCLUDED_Profiler_h


#include "havoc/havoc.h"
#include "Global.h"
#include <iomanip>
#include <atomic>
#include <sstream>


struct Profiler
{
    struct Timer
    {
        Timer() : cycleCount(0) { }
        void start()
        {
            this->cycleCount.fetch_sub(havoc_get_timestamp());
        }
        void stop()
        {
            this->cycleCount.fetch_add(havoc_get_timestamp());
        }
        havoc_timestamp get() const
        {
            return this->cycleCount;
        }
    private:
        std::atomic<havoc_timestamp> cycleCount;
    };

    // Scoped timer - measures CPU cycles spent in the current C++ scope
    struct Scope
    {
        Scope(Timer &timer)
        :
            timer(timer)
        {
            this->timer.start();
        }
        ~Scope()
        {
            this->timer.stop();
        }
    private:
        Timer &timer;
    };

    struct Timers
    {
        Timer encode;
        Timer searchInter;
        Timer searchMotionFullPel;
        Timer searchMotionHalfPel;
        Timer searchMotionQuarterPel;
        Timer searchIntra;
        Timer searchIntraSatd;
        Timer searchIntraRd;
        Timer searchIntraChroma;
        Timer searchTotal;
        Timer postprocess;
        Timer deblock;
        Timer pad;
        Timer write;
        Timer output;
        std::atomic<size_t> samples;

        Timers()
        :
            samples(0)
        {
        }

        void processed(size_t samples)
        {
            this->samples.fetch_add(samples);
        }

        double cyclesPerSamplesOf(const Timer &timer) const
        {
            return static_cast<double>(timer.get()) / this->samples;
        }

        void report(std::ostream &os) const
        {
            std::ostringstream o;
            o << std::setprecision(0);
            o << "Profiler report (units are \"CPU cycles per video sample\")\n";
            o << "Profile - encode:  " << this->cyclesPerSamplesOf(this->encode) << "\n";
            o << "Profile -   search:  " << this->cyclesPerSamplesOf(this->searchTotal) << "\n";
            o << "Profile -     inter:  " << this->cyclesPerSamplesOf(this->searchInter) << "\n";
            o << "Profile -       full-pel search:  " << this->cyclesPerSamplesOf(this->searchMotionFullPel) << "\n";
            o << "Profile -       half-pel search:  " << this->cyclesPerSamplesOf(this->searchMotionHalfPel) << "\n";
            o << "Profile -       quarter-pel search:  " << this->cyclesPerSamplesOf(this->searchMotionQuarterPel) << "\n";
            o << "Profile -     intra:  " << this->cyclesPerSamplesOf(this->searchIntra) << "\n";
            o << "Profile -       intra luma SATD:  " << this->cyclesPerSamplesOf(this->searchIntraSatd) << "\n";
            o << "Profile -       intra luma RDO:  " << this->cyclesPerSamplesOf(this->searchIntraRd) << "\n";
            o << "Profile -       intra chroma:  " << this->cyclesPerSamplesOf(this->searchIntraChroma) << "\n";
            o << "Profile -   write:  " << this->cyclesPerSamplesOf(this->write) << "\n";
            o << "Profile - postprocess:  " << this->cyclesPerSamplesOf(this->postprocess) << "\n";
            o << "Profile -   deblock:  " << this->cyclesPerSamplesOf(this->deblock) << "\n";
            o << "Profile -   pad:  " << this->cyclesPerSamplesOf(this->pad) << "\n";
            o << "Profile - output:  " << this->cyclesPerSamplesOf(this->output) << "\n";
            o << "\n";
            os << o.str();
        }
    };

};

#endif
