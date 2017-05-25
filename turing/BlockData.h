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

#ifndef INCLUDED_BlockData_h
#define INCLUDED_BlockData_h

#include "Global.h"
#include "Handlers.h"
#include "MotionVector.h"
#include "HevcTypes.h"
#include <type_traits>


int dpbIndexPlus1(struct StatePicture *p, int refList, int refIdx);


struct PuData
{
    struct Inter
    {
        MotionVector motionVector[2];
        int8_t dpbIndexPlus1[2];
    };

    struct Intra
    {
        int8_t predModeY;
        int8_t pcm;
    };

    struct // review: should (at least partly) be union?
    {
        Intra intra;
        Inter inter;
    };

    int8_t refIdxPlus1[2]; // review: do we really need this and dpbIndexPlus1?

    // Semantics same as mvLX and MvLX
    MotionVector &mv(int refList)
    {
        return this->inter.motionVector[refList];
    }

    // Semantics same as mvLX and MvLX
    const MotionVector &mv(int refList) const
    {
        return this->inter.motionVector[refList];
    }

    int getDpbIndex(int refList) const
    {
        return this->inter.dpbIndexPlus1[refList] - 1;
    }

    void reset()
    {
        static_assert(sizeof(PuData) <= 16, "PuData must be as small as possible");
        static_assert(std::is_pod<PuData>::value, "PuData must be POD");
        this->~PuData();
        new (this) PuData();
        this->intra.predModeY = INTRA_DC;
    }

    bool operator==(const PuData &other) const
    {
        if (this->refIdxPlus1[0] != other.refIdxPlus1[0]) 
            return false;
        if (this->refIdxPlus1[0])
            if (this->inter.motionVector[0] != other.inter.motionVector[0]) 
                return false;
        if (this->refIdxPlus1[1] != other.refIdxPlus1[1]) 
            return false;
        if (this->refIdxPlus1[1])
            if (this->inter.motionVector[1] != other.inter.motionVector[1]) 
                return false;
        return true;
    }

    bool operator!=(const PuData &other) const
    {
        return !this->operator==(other);
    }

    // actually, this really means "is available and is not intra"
    bool isAvailable() const
    {
        return this->refIdxPlus1[0] != 0 || this->refIdxPlus1[1] != 0;
    }

    void markUnavailable()
    {
        this->refIdxPlus1[0] = 0;
        this->refIdxPlus1[1] = 0;
    }

    void copyFrom(const PuData &other, int refList)
    {
        this->refIdxPlus1[refList] = other.refIdxPlus1[refList];
        this->inter.motionVector[refList] = other.inter.motionVector[refList];
        this->inter.dpbIndexPlus1[refList] = other.inter.dpbIndexPlus1[refList];
    }

    bool predFlag(int refList) const
    {
        return this->refIdxPlus1[refList] != 0;
    }

    template <class H>
    void setRefIdx(H &h, int refList, int refIdx)
    {
        this->refIdxPlus1[refList] = refIdx + 1;
        if (refIdx < 0)
        {
            this->inter.dpbIndexPlus1[refList] = 0;
        }
        else
        {
            this->inter.dpbIndexPlus1[refList] = dpbIndexPlus1(static_cast<struct StatePicture *>(h), refList, refIdx);
        }
    }

    int refIdx(int refList) const
    {
        return this->refIdxPlus1[refList] - 1;
    }
};


struct CuData
{
    int8_t CtDepth;
    int8_t skip;
};

struct BlockData :
    PuData,
    CuData,
    AccessOperators<BlockData>
{
    int CuPredMode() const
    {
        if (!this->predFlag(L0) && !this->predFlag(L1)) return MODE_INTRA;
        return this->skip ? MODE_SKIP : MODE_INTER;
    }

    bool operator==(const BlockData &other) const
    {
        if (static_cast<const PuData &>(*this) != static_cast<const PuData &>(other)) return false;
        if (this->CtDepth != other.CtDepth) return false;
        if (this->skip != other.skip) return false;
        return true;
    }

    bool operator!=(const BlockData &other) const
    {
        return !this->operator==(other);
    }

    void setup(coding_quadtree const *cqt, int CuPredMode)
    {
        BlockData &blockData = *this;
        blockData.reset();
        blockData.intra.pcm = 0;
        blockData.skip = (CuPredMode == MODE_SKIP);
        blockData.refIdxPlus1[0] = (CuPredMode == MODE_INTRA ? 0 : 88);
        blockData.CtDepth = cqt->cqtDepth;
    }

    void reset()
    {
        this->markUnavailable();
        this->CtDepth = -1;
        this->skip = 0;
        this->intra.predModeY = INTRA_DC;
    }

    bool isActuallyAvailable() const
    {
        return this->CtDepth >= 0;
    }
};


static std::ostream &operator<<(std::ostream &o, const PuData &puData)
{
    o << ((puData.predFlag(L0) || puData.predFlag(L1)) ? "INTER" : "INTRA");
    if (puData.predFlag(L0))
    {
        o << " L0 refIdx=" << puData.refIdx(L0) << " " << puData.mv(L0) << " ";
    }
    if (puData.predFlag(L1))
    {
        o << " L1 refIdx=" << puData.refIdx(L1) << " " << puData.mv(L1) << " ";
    }
    return o;
}

static std::ostream &operator<<(std::ostream &o, const CuData &cuData)
{
    o << "CtDepth=" << (int)cuData.CtDepth;
    o << " skip=" << (int)cuData.skip;
    return o;
}

static std::ostream &operator<<(std::ostream &o, const BlockData &blockData)
{
    o << static_cast<const CuData &>(blockData);
    o << " ";
    o << static_cast<const PuData &>(blockData);
    return o;
}


template <>
struct Access<cu_skip_flag, BlockData>
{
    typedef int8_t Type;
    static Type get(cu_skip_flag, BlockData &s)
    {
        return s.skip;
    }
};


template <>
struct Access<pcm_flag, BlockData>
{
    typedef int8_t Type;
    static Type get(pcm_flag, BlockData &s)
    {
        return s.intra.pcm;
    }
    static void set(pcm_flag, Type v, BlockData &s)
    {
        s.intra.pcm = v;
    }
};


template <>
struct Access<IntraPredModeY, BlockData>
{
    typedef int8_t Type;
    static Type get(IntraPredModeY, BlockData &s)
    {
        return s.intra.predModeY;
    }
    static void set(IntraPredModeY, Type v, BlockData &s)
    {
        s.intra.predModeY = v;
        assert(s.inter.dpbIndexPlus1[0] = 0);
//		assert(s.inter.dpbIndexPlus1[1] = 1); // available
    }
};


template <>
struct Access<CuPredMode, BlockData>
{
    typedef int Type;
    static Type get(CuPredMode, BlockData &s)
    {
        return s.CuPredMode();
    }
};

#endif
