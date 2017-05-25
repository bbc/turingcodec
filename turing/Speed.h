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

#ifndef INCLUDED_Speed_h
#define INCLUDED_Speed_h

#pragma once


// Encoder speed control types and functions

#include "Global.h"


#define ENCODER_SPEED_PRESETS_XMACRO \
        X(slow) \
        X(medium) \
        X(fast)


struct Speed
{
    enum Type
    {
#define X(A) A,
        ENCODER_SPEED_PRESETS_XMACRO
#undef X
    };

    Speed(Type speed = slow) : speed(speed) { }

    operator Type() const
    {
        return this->speed;
    }

    bool trySplit(const coding_quadtree &cqt) const
    {
        return true;
    }

    bool useAmp() const
    {
        return false;
    }

    bool useSmp(int log2PartitionSize) const
    {
        if (*this <= medium)
        {
            return true;
        }
        else
        {
            return (log2PartitionSize <= 3);
        }
        return true;
    }

    bool useSmallSearchWindow() const
    {
        return *this >= fast;
    }

    bool useBiSmallSearchWindow() const
    {
        return *this >= fast;
    }

    bool doIntraInInter() const
    {
        return true;
    }

    bool useFdm() const
    {
        return true;
    }

    bool useFdam() const
    {
        return *this >= medium;
    }

    bool useEcu() const
    {
        return *this >= medium;
    }

    bool useEsd() const
    {
        return *this >= medium;
    }

    bool useCfm() const
    {
        return *this >= medium;
    }

    bool useMet() const
    {
        return *this >= medium;
    }

    bool useRqt() const
    {
        return *this <= slow;
    }

    bool useRdoq() const
    {
        return *this <= medium;
    }

    bool useRcuDepth() const
    {
        return *this >= medium;
    }

    bool useSdh() const
    {
        return *this <= medium;
    }

    bool useTSkip() const
    {
        return false;
    }

    bool useAps() const
    {
        return *this >= medium;
    }

    bool useSao() const
    {
        return *this <= medium;
    }

    bool useSaoSlow() const
    {
        return *this < medium;
    }

    int  setMaxNumMergeCand() const
    {
        return (*this <= medium) ? 5 : 2;
    }

    bool doIntraSearch() const
    {
        return true;
    }

    bool doHalfPelRefinement() const
    {
        return true;
    }

    bool doQuarterPelRefinement() const
    {
        return *this <= medium;
    }

    int nCandidatesIntraRefinement(int log2PartitionSize) const
    {
        if (*this <= slow) return 8;
        if (*this <= medium)
        {
            return log2PartitionSize > 3 ? 3 : 8;
        }
        else// if (*this <= medium)
        {
            return log2PartitionSize > 3 ? 3 : 4;
        }
        return 35;
    }

private:
    Type speed;
};

static std::ostream &operator<<(std::ostream& o, Speed::Type speed)
{
    switch (speed)
    {
        default: return o << "<invalid>";
#define X(A) case Speed::A: return o << #A;
        ENCODER_SPEED_PRESETS_XMACRO
#undef X
    };
}

#endif
