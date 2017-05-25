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

#ifndef INCLUDED_CodedData_h
#define INCLUDED_CodedData_h

#pragma once

#include "Global.h"
#include "BitField.h"
#include "ScanOrder.h"
#include <stdint.h>
#include <cassert>
#include <fstream>
#include <iomanip>


// Encoder decisions serialised to a sequence of uint16_t, ready for entropy coding.

namespace CodedData {

    typedef uint16_t Type;

    struct PredictionUnit
    {
        Type *p;

        // Information about a motion vector
        union Metadata
        {
            uint8_t raw;
            BitField<uint8_t, 5, 1> mvp_lX_flag;
            BitField<uint8_t, 4, 1> predFlag;
            BitField<uint8_t, 0, 4> ref_idx_lX;
        };

        union Word0
        {
            Type raw;
            Metadata metadata[2 /* refList */];
        };

        Word0 &word0()
        {
            static_assert(sizeof(Word0) == 2, "assumptions about compilation of Word0 were wrong");
            return reinterpret_cast<Word0 &>(this->p[0]);
        }

        void init()
        {
            this->word0().raw = 0;
        }

        MotionVector &mvd(int x)
        {
            ASSERT(this->word0().metadata[x].predFlag);

            const bool predFlagL0 = !!this->word0().metadata[L0].predFlag;

            return *reinterpret_cast<MotionVector *>(&this->p[(predFlagL0 && x == L1) ? 3 : 1]);
        }

        // Returns a pointer to the first available CodedData entry after this prediction unit
        Type *next()
        {
            ASSERT(this->word0().metadata[L0].predFlag || this->word0().metadata[L1].predFlag);

            Type *p =  this->p;

            ++p; // word0

            if (this->word0().metadata[L0].predFlag) p += 2; // mvdL0

            if (this->word0().metadata[L1].predFlag) p += 2; // mvdL1

            return p;
        }
    };


    struct Bit
    {
        uint16_t * const p;
        uint16_t const mask;

        operator int() const
            {
            return (*this->p & mask) ? 1 : 0;
            }

        Bit const &operator=(int value)
        {
            *this->p &= ~this->mask;
            *this->p |= value ? this->mask : 0;
            return *this;
        }
    };


    struct SubBlock
    {
        friend struct Residual;

        inline SubBlock(Type *p) : p(p) {}

        inline void init(Type *p)
        {
            this->p = p;
        }

        inline Type *remaining()
        {
            return this->p + 3;
        }

        inline int sigCoeffFlag(int i) const
        {
            return (this->p[0] >> (15 - i)) & 1;
        }

        inline void  setSigCoeffFlag(int i)
        {
            this->p[0] |= bit(15 - i);
        }

        inline uint16_t sigCoeffFlags() const
        {
            return this->p[0];
        }

        inline bool significant() const
        {
            return !!this->p[0];
        }

        inline int coeffSignFlag(int i) const
        {
            return (this->p[2] >> (15 - i)) & 1;
        }

        inline void setCoeffSignFlag(int i)
        {
            this->p[2] |= bit(15 - i);
        }

        inline uint16_t coeffSignFlags()
        {
            uint16_t b = 0;
            for (int i = 0; i < 16; ++i)
                if (coeffSignFlag(i))
                    b |= 0x8000 >> i;
            return b;
        }

        inline int coeffGreater1Flag(int i) const
        {
            return (this->p[1] >> (15 - i)) & 1;
        }

        inline uint16_t greater1Flags() const
        {
            return this->p[1];
        }

        inline void setCoeffGreater1Flag(int i)
        {
            this->p[1] |= bit(15 - i);
        }

        inline void clearFlags()
        {
            this->p[0] = 0;
            this->p[1] = 0;
            this->p[2] = 0;
        }

        static inline uint16_t bit(int i)
        {
            return 1 << i;
        }

        inline int lastScanPos() const
        {
            unsigned long mask = this->sigCoeffFlags();
            ASSERT(mask);
            unsigned long index;
#if WIN32
            _BitScanForward(&index, mask);
#else
            index = 0;
            while (!(mask & (1 << index)))
                ++index;
#endif
            return 15 - index;
        }

    private:
        Type *p;
    };


    struct Residual
    {
        void init(SubBlock const &subBlock)
        {
            this->p = subBlock.p;
        }

        Type *p;

        Bit transformSkipFlag() const
        {
            return Bit{ this->p + 0, static_cast<uint16_t>(1) };
        }

        Bit codedSubBlockFlag(int i) const
        {
            return Bit{ this->p + 1 + (i >> 4), static_cast<uint16_t>(1 << (i & 0xf)) };
        }

        void clearSubBlockFlags(int log2TrafoSize)
        {
            //if (log2TrafoSize > 2)
            {
                this->p[1] = 0;
                if (log2TrafoSize == 5)
                {
                    this->p[2] = 0;
                    this->p[3] = 0;
                    this->p[4] = 0;
                }
            }
        }

        SubBlock initialSubBlock(int log2TrafoSize) const
        {
            switch (log2TrafoSize)
            {
                default:
                case 2:
                    //return{ this->p };
                case 3:
                case 4:
                    return{ this->p + 1 + 1};
                case 5:
                    return{ this->p + 1 + 4 };
            }
        }

        static int firstSetBit(uint16_t mask)
        {
            ASSERT(mask);
            unsigned long index;
#ifdef WIN32
            _BitScanReverse(&index, mask);
#else
            index = 15;
            while (!(mask & (1 << index))) --index;
#endif
            return index;
        }

        int lastSubBlock(int log2TrafoSize) const
        {
            if (log2TrafoSize == 5)
            {
                if (this->p[4]) return 48 + firstSetBit(this->p[4]);
                if (this->p[3]) return 32 + firstSetBit(this->p[3]);
                if (this->p[2]) return 16 + firstSetBit(this->p[2]);
            }
            if (this->p[1]) return firstSetBit(this->p[1]);
            return 0;
        }

        void check(const residual_coding &rc)
        {
            if (rc.log2TrafoSize == 3)
            {
                ASSERT((this->p[1] & 0xfff0) == 0);
            }

            if (rc.log2TrafoSize == 5)
            {
                ASSERT(this->p[1] || this->p[2] || this->p[3] || this->p[4]);
            }
            else
            {
                ASSERT(this->p[1]);
            }
        }
    };


    struct TransformTree
    {
        Type *p;

        union Word0
        {
            Type raw;

            // Depth of transform, equal to number of transform splits necessary to reach the residual. Same semantics as transform_tree::trafoDepth.
            BitField<Type, 0, 3> trafoDepth;

            BitFieldArray<Type, 3, 1, 3> cbf;
            BitField<Type, 3, 3> cbfWord;

            BitField<Type, 8, 2> blockIdx;
            BitField<Type, 10, 4> magic;

            BitField<Type, 15, 1> split_transform_flag;
        };

        void init(int trafoDepth, int blockIdx)
        {
            this->word0().raw = 0;
            this->word0().magic = 0xa;
            this->word0().blockIdx = blockIdx;
            this->word0().trafoDepth = trafoDepth;
        }

        void check(int trafoDepth, int blockIdx)
        {
            ASSERT(this->word0().magic == 0xa);
            ASSERT(this->word0().blockIdx == blockIdx);
            ASSERT(this->word0().trafoDepth == trafoDepth);
        }

        Word0 &word0() const
        {
            static_assert(sizeof(Word0) == 2, "assumptions about compilation of Word0 were wrong");
            return reinterpret_cast<Word0 &>(this->p[0]);
        }

        Residual firstResidual() const
        {
            ASSERT(!this->word0().split_transform_flag);
            return { this->p + 1 };
        }
    };


    struct CodingUnit
    {
        Type *p;

        union Word0
        {
            Type raw;
            BitField<Type, 0, 2> CtDepth;
            BitField<Type, 2, 2> CuPredMode;
            BitField<Type, 4, 3> part_mode;

            BitField<Type, 8, 3> cbfWord;
            BitFieldArray<Type, 8, 1, 3> cbf;

            BitField<Type, 12, 1> pcm_flag;

            BitField<Type, 13, 3> intra_chroma_pred_mode;
        };

        union Word1
        {
            Type raw;
            BitFieldArray<Type, 0, 4, 4> merge;
        };

        Word0 &word0()
        {
            static_assert(sizeof(Word0) == 2, "assumptions about compilation of Word0 were wrong");
            return reinterpret_cast<Word0 &>(this->p[0]);
        }

        Word1 &word1()
        {
            static_assert(sizeof(Word1) == 2, "assumptions about compilation of Word1 were wrong");
            return reinterpret_cast<Word1 &>(this->p[1]);
        }

        int8_t &IntraPredModeY(int partIdx)
        {
            assert(partIdx >= 0);
            assert(partIdx < (this->word0().part_mode ? 4 : 1));

            assert(this->word0().CuPredMode == MODE_INTRA);

            int8_t *p = reinterpret_cast<int8_t *>(&this->word1());
            return p[partIdx];
        }

        void init()
        {
            this->word0().raw = 0;
            this->word1().raw = 0;
        }

        void check(int CtDepth)
        {
            ASSERT(this->word0().CtDepth == CtDepth);
            if (this->word0().CuPredMode == 1)
            {
                // intra
                ASSERT(this->word0().part_mode == 0 || this->word0().part_mode == 1);
            }
            else
            {
                // inter or skip
                ASSERT(this->word0().intra_chroma_pred_mode == 0);
            }
        }

        TransformTree firstTransformTree()
        {
            ASSERT(this->word0().CuPredMode == MODE_INTRA);
            return{ &this->p[this->word0().part_mode ? 4 : 3] };
        }

        PredictionUnit firstPredictionUnit()
        {
            ASSERT(this->word0().CuPredMode != MODE_INTRA);
            return { &this->p[2] };
        }

        TransformTree firstTransformTreeChroma()
        {
            ASSERT(this->word0().CuPredMode == MODE_INTRA);
            TransformTree transformTree = this->firstTransformTree();
            transformTree.p += this->chromaOffset();
            return { transformTree.p };
        }

        Type &chromaOffset()
        {
            ASSERT(this->word0().CuPredMode == MODE_INTRA);
            return this->p[this->word0().part_mode ? 3 : 2];
        }
    };


    static void storeResidual(CodingUnit codedCu, Residual &residual, int16_t *coefficients, int log2TrafoSize, int scanIdx, bool cbf, TransformTree transformTree, int cIdx)
    {
        if (!cbf) 
            return;

        //std::cout << "storeResidual codedCu.p=" << codedCu.p << " residual.p=" << residual.p << "\n";

        codedCu.word0().cbf[cIdx] = 1;
        transformTree.word0().cbf[cIdx] = 1;

        CodedData::SubBlock subBlock = residual.initialSubBlock(log2TrafoSize);

        int const nCbS = 1 << log2TrafoSize;
        int const nSubBlocks = 1 << (2 * (log2TrafoSize - 2));

        Raster<int16_t> TransCoeffLevel(coefficients, nCbS);

        residual.clearSubBlockFlags(log2TrafoSize);

        auto *rasterC = rasterScanOrderC(log2TrafoSize, scanIdx);
        auto *rasterS = rasterScanOrderS(log2TrafoSize - 2, scanIdx);

        for (int i = nSubBlocks - 1; i >= 0; --i)
        {
            subBlock.clearFlags();

            CodedData::Type *p = subBlock.remaining();

            auto coefficientsSubBlock = coefficients + rasterS[i];

            for (int n = 15; n >= 0; --n)
            {
                auto value= coefficientsSubBlock[rasterC[n]];

                if (value != 0)
                {
                    subBlock.setSigCoeffFlag(n);

                    if (value < 0)
                    {
                        subBlock.setCoeffSignFlag(n);
                        value = -value;
                    }

                    if (value > 1)
                    {
                        subBlock.setCoeffGreater1Flag(n);
                        *p++ = value;
                    }
                }
            }

            if (subBlock.significant())
            {
                residual.codedSubBlockFlag(i) = 1;
                subBlock.init(p);
            }
        }

        residual.init(subBlock);
    }

} // namespace CodedData


// A set of pointers and accessors into sequential coded data.
struct StateCodedData
{
    // Maintains a pointer to the first word of the current coding_unit's coded data.
    CodedData::CodingUnit codedCu;

    // Maintains a pointer to the first word of the current prediction_unit's coded data.
    CodedData::PredictionUnit codedPu;

    // Current partition index. Review: duplicated in StateEncodeSubstream
    int partIdx;

    CodedData::TransformTree transformTree;

    // Different to transformTree for MODE_INTRA CUs. Maintained the same for inter CUs.
    CodedData::TransformTree transformTreeChroma;

    CodedData::Residual residual;

    // Pointer to first word of the current Candidate's coded data.
    CodedData::Type *codedDataBefore;

    // Pointer to the word after the end of the current Candidate's coded data.
    CodedData::Type *codedDataAfter;

    // Pointer to the parent buffer where the current Candidate's coded data would be copied.
    CodedData::Type *codedDataParent;

    CodedData::TransformTree transformTreeAncestry[5];

    void reset(CodedData::Type *buffer)
    {
        this->codedCu.p = buffer;
        this->codedDataBefore = buffer;
        this->codedDataAfter = buffer;
    }

    template <class F>
    void copyBefore(F, const StateCodedData &other, int copyCuWords = 3)
    {
        assert(this->codedDataBefore == this->codedDataAfter);

        if (std::is_same<F, coding_quadtree>::value ||
                std::is_same<F, coding_unit>::value ||
                std::is_same<F, IntraPartition>::value)
        {
            // This is a copy made prior to an operation upon a CU or CQT.

            // record parent position our data will be copied back to (after).
            this->codedDataParent = other.codedCu.p;

            this->startCu();

            if (std::is_same<F, IntraPartition>::value || copyCuWords < 0)
            {
                CodedData::Type const *src = other.codedDataBefore;
                CodedData::Type const *srcEnd = other.codedDataAfter;

                while (src != srcEnd)
                {
                    *this->codedDataAfter++ = *src++;
                }
            }

            // Copy any relevant prepopulated information.
            this->copyCuWords(other, copyCuWords);
        }
        else
        {
            ASSERT(0);
            this->codedDataParent = 0;
        }
    }

    void copyCuWords(const StateCodedData &other, int n)
    {
        auto src = other.codedCu.p;
        auto dst = this->codedCu.p;
        for (int i = 0; i < n; ++i)
        {
            *dst++ = *src++;
        }
    }

    template <class F>
    void copyAfter(F, const StateCodedData &other)
    {
        assert(other.codedDataParent >= this->codedDataBefore);
        assert(other.codedDataParent <= this->codedDataAfter);

        this->codedDataAfter = other.codedDataParent;

        {
            CodedData::Type const *src = other.codedDataBefore;
            CodedData::Type const *srcEnd = other.codedDataAfter;

            while (src != srcEnd)
            {
                *this->codedDataAfter++ = *src++;
            }
        }

        if (std::is_same<F, IntraPartition>::value)
        {
            int cbf[3];

            cbf[0] = this->codedCu.word0().cbf[0];
            cbf[1] = this->codedCu.word0().cbf[1];
            cbf[2] = this->codedCu.word0().cbf[2];

            // Copy CU data, including IntraPredModeY values
            this->copyCuWords(other, 3);

            this->codedCu.word0().cbf[0] = this->codedCu.word0().cbf[0] | cbf[0];
            this->codedCu.word0().cbf[1] = this->codedCu.word0().cbf[1] | cbf[1];
            this->codedCu.word0().cbf[2] = this->codedCu.word0().cbf[2] | cbf[2];
        }
    }

    // Initialise before writing a coding unit
    void startCu()
    {
        this->codedCu.p = this->codedDataAfter;
    }

    // Initialise before reading a first coding unit
    void startReading()
    {
        this->codedCu.p = this->codedDataBefore;
    }

    // Moves pointer to first prediction unit
    void firstPu()
    {
        this->codedPu = this->codedCu.firstPredictionUnit();
        this->partIdx = 0;
    }
};


template <class V, int cIdx>
struct AccessCbf
{
    typedef bool Type;
    static Type get(V v, StateCodedData &s)
    {
        if (s.codedCu.word0().CuPredMode != MODE_INTRA)
        {
            bool const n = !!s.transformTreeAncestry[v.trafoDepth].word0().cbf[cIdx];
            if (v.trafoDepth != s.transformTreeAncestry[v.trafoDepth].word0().trafoDepth)
            {
                assert(v.trafoDepth == s.transformTreeAncestry[v.trafoDepth].word0().trafoDepth - 1);
                if (s.transformTreeAncestry[v.trafoDepth + 1].word0().cbf[cIdx])
                    return true;
            }
            return n;
        }
        bool o;
        if (v.trafoDepth == 0)
        {
            o = !!s.codedCu.word0().cbf[cIdx];
        }
        else if (cIdx)
        {
            o = !!s.transformTreeChroma.word0().cbf[cIdx];
        }
        else
        {
            o = !!s.transformTree.word0().cbf[0];
        }
        return o;
    }
};


template <class S>
struct Access<cbf_luma, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
:
AccessCbf<cbf_luma, 0>
{
};


template <class S>
struct Access<cbf_cb, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
:
AccessCbf<cbf_cb, 1>
{
};


template <class S>
struct Access<cbf_cr, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
:
AccessCbf<cbf_cr, 2>
{
};


template <class S>
struct Access<rqt_root_cbf, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef bool Type;
    static Type get(rqt_root_cbf v, StateCodedData &s)
    {
        return s.codedCu.word0().CuPredMode == MODE_INTRA || !!s.codedCu.word0().cbfWord;
    }
};


template <class S>
struct Access<merge_flag, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef bool Type;
    static Type get(merge_flag v, StateCodedData &s)
    {
        return !!s.codedCu.word1().merge[s.partIdx];
    }
};


template <class S>
struct Access<transform_skip_flag, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef bool Type;
    static Type get(transform_skip_flag v, StateCodedData &s)
    {
        return !!s.residual.transformSkipFlag();
    }
};


// review: go explicit (and above)
template <class S>
struct Access<merge_idx, S, typename std::enable_if<std::is_base_of<StateCodedData, S>::value>::type>
{
    typedef int Type;
    static Type get(merge_idx v, StateCodedData &s)
    {
        ASSERT(s.codedCu.word1().merge[s.partIdx]);
        return s.codedCu.word1().merge[s.partIdx] - 1;
    }
};

template <typename Sample> struct Candidate;

template <>
struct Access<Neighbouring<CuPredMode, Current>, Candidate<uint8_t>>
{
    typedef int Type;
    static Type get(Neighbouring<CuPredMode, Current> v, StateCodedData &s)
    {
        return s.codedCu.word0().CuPredMode;
    }
};

template <>
struct Access<Neighbouring<CuPredMode, Current>, Candidate<uint16_t>>
{
    typedef int Type;
    static Type get(Neighbouring<CuPredMode, Current> v, StateCodedData &s)
    {
        return s.codedCu.word0().CuPredMode;
    }
};

#endif
