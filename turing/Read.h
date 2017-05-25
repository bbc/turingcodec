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

#ifndef INCLUDED_Read_h
#define INCLUDED_Read_h

#pragma once

#include "StateDecode.h"
#include "StreamReader.h"
#include "CodeCache.h"
#include "Syntax.h"
#include "SyntaxSei.h"
#include "SyntaxElements.h"
#include "StatePictures.h"
#include "GlobalState.h"
#include "CandModeList.h"
#include "LoopFilter.h"
#include "Violation.h"
#include "ContextModel.h"
#include "Cabac.h"
#include "StateSpatial.h"
#include "RangeLimits.h"
#include "Mvp.h"
#include <fstream>


 // review: shouldn't be here
template <class F>
struct Decode;


// checks just-decoded element value, may throw Abort if it's bad
template <class V, class H, class Enable = void>
struct ReadCheck
{
    static void go(V, H &, const typename ValueType<V>::Type &)
    {
    }
};


// If a value's RangeLimits are defined with base class Fatal it means we should abort parsing if limits are violated
template <class V, class H>
struct ReadCheck<V, H, typename std::enable_if<std::is_base_of<Fatal, RangeLimits<V>>::value>::type>
{
    static void go(V v, H &h, int value)
    {
        if (value < RangeLimits<V>::min(h) || value > RangeLimits<V>::max(h))
        {
            auto clause = RangeLimits<V>::clause(); // useful for debuging
            throw Abort();
        }
    }
};


#ifndef SATURN
static std::ofstream &flog()
{
    static std::ofstream ofs("cabac-decode.txt");
    return ofs;
}
#endif


struct CabacState
{
    int ivlCurrRange = -8;
    int bitsNeeded = 510;
    int ivlOffset = 0;
};


template <class F>
struct Read :
    Syntax<F>
{
};


template <>
struct Read<CabacRestart>
{
    template <class H> static void go(const CabacRestart &e, H &h)
    {
        CabacState &cabacState = *static_cast<CabacState *>(h);
        cabacState.ivlCurrRange = 510;
        cabacState.bitsNeeded = -8;
        cabacState.ivlOffset = read_bytes<int>(h, 1);
        cabacState.ivlOffset <<= 8;
        cabacState.ivlOffset |= read_bytes<int>(h, 1);
        if ((cabacState.ivlOffset >> 7) >= 510) // CondCheck 9.3.2.5-A
        {
            h(Violation("9.3.2.5", "ivlOffset{%1%} >= 510") % (cabacState.ivlOffset >> 7));
            throw Abort();
        }
    }
};


template <>
struct Read<rbsp_slice_segment_trailing_bits>
{
    template <class H> static void go(const rbsp_slice_segment_trailing_bits &s, H &h)
    {
        if (isRasl(h[nal_unit_type()]) && static_cast<StatePicturesBase *>(h)->NoRaslOutputFlag)
        {
            // associated IRAP picture has NoRaslOutputFlag equal to 1, the RASL picture is not output and may not be correctly decodable
            return;
        }

        Syntax<rbsp_slice_segment_trailing_bits>::go(s, h);
    }
};


template <>
struct Read<coding_tree_unit>
{
    template <class H> static void go(const coding_tree_unit &s, H &h);
};


template <typename T>
static T bitPos(BitReader& bitReader)
{
    T pos = 0;
    int m = 0x80;
    while (m != bitReader.state.mask)
    {
        m >>= 1;
        ++pos;
    }
    if (!std::is_integral<T>::value)
    {
        assert(0);
    }
    return pos + static_cast<T>(8 * (bitReader.state.p - bitReader.begin));
};






namespace {

    template <class Stream>
    bool canFastReadRbspData(Stream &)
    {
        return false;
    }

    template <>
    bool canFastReadRbspData<StreamReader>(StreamReader &streamReader)
    {
        return !!streamReader.codeCache;
    }

    template <class Stream>
    void fastReadRbspData(Stream &s, std::vector<uint8_t> &rbspData, std::streamsize endPositionNalUnit, std::vector<size_t> &eb3pPositions)
    {
    }

    template <>
    void fastReadRbspData<StreamReader>(StreamReader &s, std::vector<uint8_t> &rbspData, std::streamsize endPositionNalUnit, std::vector<size_t> &eb3pPositions)
    {
        s.readRbspData(rbspData, static_cast<size_t>(endPositionNalUnit), eb3pPositions);
    }

}

namespace {

    template <class Stream>
    size_t numBytesInNalUnitSimple(Stream& stream)
    {
        typename Stream::Bookmark mark(stream);

        try
        {
            if (stream.nextBytes(4) == 0x00000001)
            {
                stream.readBytes(1);
            }

            if (stream.readBytes(3) != 0x000001) return 0;
        }
        catch (ExceptionOverrun &)
        {
            return 0;
        }

        size_t pos1 = bitPos<size_t>(stream) / 8;

        try
        {
            while (stream.nextBytes(3) > 0x000001)
            {
                stream.readBytes(1);
            }
        }
        catch (ExceptionOverrun &)
        {
            return bitLen(stream) / 8 - pos1;
        }

        return bitPos<int>(stream) / 8 - pos1;
    }

    template <class Stream>
    size_t numBytesInNalUnit(Stream& stream)
    {
        return numBytesInNalUnitSimple(stream);
    }

    template <>
    size_t numBytesInNalUnit<StreamReader>(StreamReader& stream)
    {
        if (stream.codeCache)
        {
            return stream.numBytesInNalUnit();
        }
        return numBytesInNalUnitSimple(stream);
    }

    template <class Stream>
    void gotoNextStartCode(Stream &stream)
    {
        while (stream.nextBytes(3) > 0x000001)
        {
            stream.readBytes(1);
        }
    }


    template <>
    void gotoNextStartCode<StreamReader>(StreamReader& stream)
    {
        if (!stream.codeCache)
        {
            while (stream.nextBytes(3) > 0x000001)
            {
                stream.readBytes(1);
            }
            return;
        }

        stream.adjustCodeIt(true);
        while (stream.codeIt != stream.codeCache->codes.end() && stream.codeIt->byte > 1) ++stream.codeIt;

        if (stream.codeIt == stream.codeCache->codes.end())
        {
            seekStream(stream, bitLen(stream));
        }
        else
        {
            std::streamoff position = stream.codeIt->position;
            if (stream.codeIt->numberOfZeroes > 3)
            {
                position += stream.codeIt->numberOfZeroes - 3;
            }
            seekStream(stream, 8 * size_t(position));
        }
    }

}

template <typename T, class H>
T read_bits(H &h, size_t n)
{
    T returnValue = 0;
    while (n--)
    {
        returnValue <<= 1;
        returnValue |= h[Stream()].bit();
    }
    return returnValue;
}

template <typename T, class H>
struct ReadBytes<T, H>
{
    static T go(H &h, int n)
    {
        return h[Stream()].readBytes(n);
    }
};


// Fixed-length bitfield parsing
template <class V, class M>
struct ReadU
{
    template <class H> static void go(Element<V, M> e, H &h)
    {
        assert(e.m.n > 0);
        auto t = read_bits<typename ValueType<V>::Type>(h, e.m.n);
        Access<V, H>::set(e.v, t, h);
        ReadCheck<V, H>::go(e.v, h, t);
    }
};

template <class V> struct Read<Element<V, u>> : ReadU<V, u> { };
template <class V> struct Read<Element<V, b>> : ReadU<V, b> { };
template <class V> struct Read<Element<V, f>> : ReadU<V, f> { };


template <class V>
struct Read<Element<V, uv>>
{
    template <class H> static void go(Element<V, uv> e, H &h)
    {
        auto const n = NumberOfBitsUv<V>::get(e.v, h);
        auto t = read_bits<typename ValueType<V>::Type>(h, n);
        Access<V, H>::set(e.v, t, h);
        ReadCheck<V, H>::go(e.v, h, t);
    }
};


template <class V, class M>
struct ReadI
{
    template <class H> static void go(Element<V, M> e, H &h)
    {
        auto x = read_bits<typename ValueType<V>::Type>(h, e.m.n);
        // sign extend
        auto const signBit = 1 << (e.m.n - 1);
        if (x & signBit)
        {
            x |= ~(signBit - 1);
        }
        Access<V, H>::set(e.v, x, h);
        ReadCheck<V, H>::go(e.v, h, x);
    }
};

template <class V> struct Read<Element<V, i>> : ReadI<V, i> {};



template <class V>
struct ReadUe
{
    template <class H> static void go(Element<V, ue> fun, H &h)
    {
        int leadingZeroBits = 0;
        const int max = std::numeric_limits<typename ValueType<V>::Type>::digits - 1;
        while (!read_bits<int>(h, 1))
        {
            if (++leadingZeroBits > max)
            {
                h(Violation("7.3", "%1% ue(v) parsing error - too many leading zero bits") % TypeName<V>::value);
                throw Abort();
            }
        }
        typename ValueType<V>::Type t = (1 << leadingZeroBits) - 1;
        t += read_bits<typename ValueType<V>::Type>(h, leadingZeroBits);
        Access<V, H>::set(fun.v, t, h);
        ReadCheck<V, H>::go(fun.v, h, t);
    }
};


template <class V>
struct Read<Element<V, ue>> :
    ReadUe<V>
{
};


template <class V>
struct Read<Element<V, se>>
{
    template <class H> static void go(Element<V, se> fun, H &h)
    {
        int leadingZeroBits = 0;
        const int max = std::numeric_limits<typename ValueType<V>::Type>::digits - 1;
        while (!read_bits<int>(h, 1))
        {
            if (++leadingZeroBits > max)
            {
                h(Violation("7.3", "%1% se(v) parsing error - too many leading zero bits") % TypeName<V>::value);
                throw Abort();
            }
        }
        typename ValueType<V>::Type t = (1 << leadingZeroBits) - 1;
        t += read_bits<typename ValueType<V>::Type>(h, leadingZeroBits);

        bool negative = (t % 2 == 0);
        t = (t + 1) / 2;
        if (negative) t = -t;

        h[fun.v] = t;
        ReadCheck<V, H>::go(fun.v, h, t);
    }
};


template <>
struct Read<sei_message>
{
    template <class H> static void go(const sei_message &f, H &h)
    {
        const char *streamTypePrevious = static_cast<StatePicturesBase *>(h)->streamType;
        static_cast<StatePicturesBase *>(h)->streamType = StatePicturesBase::streamTypeSei();

        Syntax<sei_message>::go(f, h);

        static_cast<StatePicturesBase *>(h)->streamType = streamTypePrevious;
    }
};


template <class H, class V>
struct GetCtxInc
{
    static int f(int binInx, V)
    {
        return 0;
    }
};


template <class H>
struct GetCtxInc<H, split_cu_flag>
{
    static int f(H &h, split_cu_flag e, int binIdx)
    {
        Neighbourhood *neighbourhood = h;
        auto &snake = neighbourhood->snake;

        coding_quadtree cqt = *static_cast<coding_quadtree *>(h);

        return
            (snake.get<Left>(e.x0 - 1, e.y0, neighbourhood->MinCbLog2SizeYMinus1).CtDepth > cqt.cqtDepth ? 1 : 0) +
            (snake.get<Up>(e.x0, e.y0 - 1, neighbourhood->MinCbLog2SizeYMinus1).CtDepth > cqt.cqtDepth ? 1 : 0);
    }
};


template <class V>
struct DecodeDecision
{
    DecodeDecision(int *value, int ctxInc) :
        value(value),
        ctxInc(ctxInc)
    {
    }
    int *value;
    int ctxInc;
};

template <class V>
struct DecodeBypass
{
    DecodeBypass(int *value) :
        value(value)
    {
    }
    int *value;
};

template <class V>
struct DecodeTerminate
{
    DecodeTerminate(int *value) :
        value(value)
    {
    }
    int *value;
};


template <class V>
struct NameOf
{
    static const char *name()
    {
        static std::string buff;
        if (buff.empty())
        {
            buff = typeid(V).name();
            remove(buff, "");
            remove(buff, "struct ");
            remove(buff, "_0");
            remove(buff, "_1");
            remove(buff, "_2");
        }
        return buff.c_str();
    }
private:
    static void remove(std::string &str, const char *pattern)
    {
        const std::string::size_type pos = str.find(pattern);
        if (pos == std::string::npos) return;
        str.erase(pos, std::string(pattern).size());
    }
};


template <class V, class H>
struct GetContext
{
    static ContextModel &get(H &h, int ctxInc)
    {
        typedef Contexts ContextsType;
        ContextsType &contexts = *static_cast<ContextsType *>(h);
        return contexts.get<Context<V>::Type>(ctxInc);
    }
};


static inline int decodeDecision(BitReader *bitReader, CabacState *cabacState, ContextModel &contextModel)
{
    const int pStateIdx = contextModel.getState();
    const int valMPS = contextModel.getMps();

    const int qRangeIdx = (cabacState->ivlCurrRange >> 6) & 3;
    const int ivlLpsRange = rangeTabLPS(pStateIdx, qRangeIdx); // uiLPS
    cabacState->ivlCurrRange -= ivlLpsRange;
    const int32_t scaledRange = cabacState->ivlCurrRange << 7;
    int32_t mask = (cabacState->ivlOffset - scaledRange) >> 31;

    int binVal, numBits;
    binVal = valMPS ^ (1 + mask);
    cabacState->ivlOffset -= scaledRange & ~mask;
    cabacState->ivlCurrRange += (ivlLpsRange - cabacState->ivlCurrRange) & ~mask;
    contextModel.update(mask);

    static const uint8_t renormTable[] = {
        6, 5, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    numBits = renormTable[cabacState->ivlCurrRange >> 3];

    cabacState->ivlOffset <<= numBits;
    cabacState->ivlCurrRange <<= numBits;
    cabacState->bitsNeeded += numBits;

    if (cabacState->bitsNeeded >= 0)
    {
        cabacState->ivlOffset |= bitReader->readBytes(1) << cabacState->bitsNeeded;
        cabacState->bitsNeeded -= 8;

        //if ((cabacState->ivlOffset >> 7) >= 510) // CondCheck 9.3.2.5-A
        //{
        //     h(Violation("9.3.2.5", "ivlOffset{%1%} >= 510") % (cabacState->ivlOffset >> 7));
        //     throw Abort();
        //}
    }

    return binVal;
}


template <class V>
struct Read<DecodeDecision<V>>
{
    template <class H> static void go(DecodeDecision<V> e, H &h)
    {
        BitReader *bitReader = h;
        CabacState *cabacState = h;
        Contexts *contexts = h;
        ContextModel &contextModel = contexts->get<typename Context<V>::Type>(e.ctxInc);
        *e.value = decodeDecision(bitReader, cabacState, contextModel);
    }
};


static int decodeBypass(CabacState *cabacState, BitReader *bitReader)
{
    cabacState->ivlOffset *= 2;
    if (++cabacState->bitsNeeded >= 0)
    {
        cabacState->bitsNeeded = -8;
        cabacState->ivlOffset |= bitReader->readBytes(1);
    }
    int scaledRange = cabacState->ivlCurrRange << 7;
    int mask = (scaledRange - cabacState->ivlOffset - 1) >> 31;
    cabacState->ivlOffset -= scaledRange & mask;

    //if ((cabacState->ivlOffset >> 7) >= 510) // CondCheck 9.3.2.5-A
    //{
    // h(Violation("9.3.2.5", "ivlOffset{%1%} >= 510") % (cabacState->ivlOffset >> 7));
    // throw Abort();
    //}

    //if ((cabacState->ivlOffset >> 7) >= (scaledRange >> 7)) // CondCheck 9.3.4.3.4-A
    //{
    // h(Violation("9.3.4.3.4", "ivlOffset{%1%} >= ivlCurrRange{%2%} after DecodeBypass")
    //         % (cabacState->ivlOffset >> 7)
    //         % (scaledRange >> 7));
    // throw Abort();
    //}

    return mask & 1;
}


template <class V>
struct Read<DecodeBypass<V>>
{
    template <class H> static inline void go(DecodeBypass<V> e, H &h)
    {
        CabacState *cabacState = h;
        BitReader *bitReader = h;
        *e.value = decodeBypass(cabacState, bitReader);
    }
};


template <class V>
struct Read<DecodeTerminate<V>>
{
    template <class H> static void go(DecodeTerminate<V> e, H &h)
    {
        BitReader &reader = *static_cast<BitReader *>(h);
        CabacState &state = *static_cast<CabacState *>(h);

        state.ivlCurrRange -= 2;
        auto const scaledRange = state.ivlCurrRange << 7;

        int binVal;

        if (state.ivlOffset >= scaledRange)
        {
            binVal = 1;
            reader.rewind(-state.bitsNeeded);
        }
        else
        {
            binVal = 0;
            if (scaledRange < (256 << 7))
            {
                state.ivlCurrRange = scaledRange >> 6;
                state.ivlOffset <<= 1;

                if (++state.bitsNeeded == 0)
                {
                    state.bitsNeeded = -8;
                    state.ivlOffset |= read_bytes<int>(h, 1);

                    if ((state.ivlOffset >> 7) >= 510) // CondCheck 9.3.2.5-A
                    {
                        h(Violation("9.3.2.5", "ivlOffset{%1%} >= 510") % (state.ivlOffset >> 7));
                        throw Abort();
                    }
                }
            }
        }

        *e.value = binVal;
    }
};


template <>
struct Read<sao>
{
    template <class H> static void go(const sao &s, H &h)
    {
        h[SaoTypeIdx(0, s.rx, s.ry)] = 0;
        h[SaoTypeIdx(1, s.rx, s.ry)] = 0;
        h[SaoTypeIdx(2, s.rx, s.ry)] = 0;

        h[sao_merge_left_flag()] = 0;
        h[sao_merge_up_flag()] = 0;

        Syntax<sao>::go(s, h);
    }
};


template <>
struct Read<coding_quadtree>
{
    template <class H> static void go(const coding_quadtree &cqt, H &h)
    {
        Neighbourhood *neighbourhood = h;
        Snake<BlockData>::Cursor *cursor = h;
        cursor->relocate(neighbourhood->snake, cqt, neighbourhood->MinCbLog2SizeYMinus1);

        BitReader &reader = *static_cast<BitReader *>(h);

        h[split_cu_flag(cqt.x0, cqt.y0)] = cqt.log2CbSize > h[MinCbLog2SizeY()] ? 1 : 0;

        Syntax<coding_quadtree>::go(cqt, h);
    }
};


template <>
struct Read<coding_unit>
{
    template <class H> static void go(const coding_unit &cu, H &h)
    {
        StatePicture *statePictureBase = h;
        Neighbourhood *neighbourhood = h;
        QpState *qpState = h;
        const coding_quadtree *cqt = h;

        h[cu_transquant_bypass_flag()] = 0;
        h[part_mode()] = 0;
        h[rqt_root_cbf()] = 1;

        Snake<BlockData>::Cursor *cursor = h;
        cursor->relocate(neighbourhood->snake, cu, neighbourhood->MinCbLog2SizeYMinus1);
        BlockData &blockData = cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1);
        blockData.reset();
        blockData.intra.pcm = 0; // inferred
        blockData.skip = 0; // inferred

        qpState->preCu(cu, h);

        StateSubstream *stateSubstream = h;
        stateSubstream->partIdx = 0;

        Syntax<coding_unit>::go(cu, h);

        qpState->postCu(cu, h);

        if (h[IntraSplitFlag()])
        {
            // Commit the 4th partition to the snake.
            const int size = 1 << cqt->log2CbSize >> 1;
            cursor->commit(Turing::Rectangle{ cqt->x0 + size, cqt->y0 + size, size, size }, neighbourhood->MinCbLog2SizeYMinus1);
            //std::cout << "commit " << cqt->x0 + size << ", " << cqt->y0 + size << ", " << size << ", " << size << "\n";
        }

        commitCu(h, *cqt);
    }
};


struct ChromaCbf
{
    uint8_t values[2][4];
    transform_tree tt0;
    int ChromaArrayType;
    void reset(const transform_tree &tt, int ChromaArrayType)
    {
        this->ChromaArrayType = ChromaArrayType;
        if (tt.trafoDepth == 0) this->tt0 = tt;
        this->values[0][tt.trafoDepth] = 0;
        this->values[1][tt.trafoDepth] = 0;
    }
};


template <>
struct Read<transform_tree>
{
    template <class H> static void go(const transform_tree &tt, H &h)
    {
        const coding_quadtree *cqt = h;

        h[split_transform_flag()] = infer(split_transform_flag(tt.x0, tt.y0, tt.trafoDepth), h);
        static_cast<ChromaCbf *>(h)->reset(tt, h[ChromaArrayType()]);
        h[cbf_luma(tt.x0, tt.y0, tt.trafoDepth)] = 1;

        Syntax<transform_tree>::go(tt, h);
    }
};


template <>
struct Read<pcm_sample>
{
    template <class H> static void go(const pcm_sample &fun, H &h)
    {
        StatePicturesBase *statePicturesBase = h;
        statePicturesBase->streamType = StatePicturesBase::streamTypeRbsp();
        Syntax<pcm_sample>::go(fun, h);
        statePicturesBase->streamType = StatePicturesBase::streamTypeCabac();
        h(CabacRestart());
    }
};


static int subBlockRaster(int log2TrafoSize, int x, int y)
{
    switch (log2TrafoSize)
    {
    default:
        assert(false);
    case 2:
        return 0;
    case 3:
        return ((y & ~3) >> 1) | (x >> 2);
    case 4:
        return (y & ~3) | (x >> 2);
    case 5:
        return ((y & ~3) << 1) | (x >> 2);
    }
}


template <int log2TrafoSize>
struct ResidualCodingStateOptimised
    :
    ValueCache<scanIdx>
{
    template <class H>
    ResidualCodingStateOptimised(H &h)
        :
        ValueCache<scanIdx>(h)
    {
    }

    uint64_t codedSubBlockFlags = 0;
    int nSigCoeffs;
    ScanPos coeffPositions[17];
    int nn[17];
    int32_t absoluteCoefficients[17];
    uint16_t remaining;
};


template <int log2TrafoSize>
struct Access<coded_sub_block_flag, ResidualCodingStateOptimised<log2TrafoSize>>
{
    static int shift(coded_sub_block_flag v)
    {
        int a = 8 + v.xS - v.yS;
        assert(a > 0 && a < 16);
        return a;
        //return 8 + v.xS - v.yS;
        //return v.xS + (v.yS << (log2TrafoSize - 2));
    }
    typedef int Type;
    static Type get(coded_sub_block_flag v, ResidualCodingStateOptimised<log2TrafoSize> &s)
    {
        return (s.codedSubBlockFlags >> shift(v)) & 1;
    }
    static void set(coded_sub_block_flag v, Type x, ResidualCodingStateOptimised<log2TrafoSize> &s)
    {
        auto const mask = static_cast<uint64_t>(1) << shift(v);
        if (x)
            s.codedSubBlockFlags |= mask;
        else
            s.codedSubBlockFlags &= ~mask;
    }
};


template <int log2TrafoSize, class H>
void optimisedReadResidualCoding(H &hParent)
{
    ResidualCodingStateOptimised<log2TrafoSize> residualCodingState(hParent);
    auto h = hParent.extend(&residualCodingState);

    StateSubstream *stateSubstream = h;
    residual_coding *rc = h;
    assert(rc->log2TrafoSize == log2TrafoSize);

    if (std::is_same<typename H::Tag, Decode<void>>::value)
    {
        const int sizeC = 1 << rc->log2TrafoSize;
        for (int yC = 0; yC < sizeC; ++yC)
        {
            for (int xC = 0; xC < sizeC; ++xC)
            {
                h[TransCoeffLevel(rc->x0, rc->y0, rc->cIdx, xC, yC)] = 0;
            }
        }
    }

    stateSubstream->lastGreater1Flag = 1;
    stateSubstream->lastGreater1Ctx = -1;

    stateSubstream->cLastAbsLevel = 0;
    stateSubstream->firstCoeffAbsLevelRemainingInSubblock = true;

    h[transform_skip_flag()] = 0;
    h[explicit_rdpcm_flag()] = 0;

    if (h[transform_skip_enabled_flag()] && !h[cu_transquant_bypass_flag()] && (log2TrafoSize <= h[Log2MaxTransformSkipSize()]))
        h(transform_skip_flag(rc->x0, rc->y0, rc->cIdx), ae(v));

    if (h[current(CuPredMode(rc->x0, rc->y0))] == MODE_INTER  &&  h[explicit_rdpcm_enabled_flag()] &&
        (h[transform_skip_flag(rc->x0, rc->y0, rc->cIdx)] || h[cu_transquant_bypass_flag()]))
    {
        h(explicit_rdpcm_flag(rc->x0, rc->y0, rc->cIdx), ae(v));
        if (h[explicit_rdpcm_flag(rc->x0, rc->y0, rc->cIdx)])
            h(explicit_rdpcm_dir_flag(rc->x0, rc->y0, rc->cIdx), ae(v));
    }

    h(last_sig_coeff_x_prefix(), ae(v));
    h(last_sig_coeff_y_prefix(), ae(v));
    if (h[last_sig_coeff_x_prefix()] > 3)
        h(last_sig_coeff_x_suffix(), ae(v));
    if (h[last_sig_coeff_y_prefix()] > 3)
        h(last_sig_coeff_y_suffix(), ae(v));

    auto xx = h[LastSignificantCoeffX()];
    auto yy = h[LastSignificantCoeffY()];

    auto const lastSigCoeffRaster = ((h[LastSignificantCoeffY()] & 3) << 2) | (h[LastSignificantCoeffX()] & 3);
    auto const lastSigSubBlockRaster = subBlockRaster(log2TrafoSize, h[LastSignificantCoeffX()], h[LastSignificantCoeffY()]);

    constexpr auto log2TrafoSizeMinus2 = log2TrafoSize - 2;

    auto const lastScanPos = scanPosInverse<2>(h[scanIdx()], lastSigCoeffRaster);
    auto const lastSubBlock = scanPosInverse<log2TrafoSizeMinus2>(h[scanIdx()], lastSigSubBlockRaster);

    int escapeDataPresent = 0;

    auto const &sp = scanPos<log2TrafoSizeMinus2>().lookup[h[scanIdx()]];
    auto const &sp2 = scanPos<2>().lookup[h[scanIdx()]];
    int &i = stateSubstream->i;
    for (i = lastSubBlock; i >= 0; i--)
    {
        auto const s = sp[i];
        const auto &xS = s.x;
        const auto &yS = s.y;
        auto const offset = 8 + xS - yS;

        // read coded_sub_block_flag
        int inferSbDcSigCoeffFlag = 0;
        if (i == 0)
        {
        }
        else if (i < lastSubBlock)
        {
            h(coded_sub_block_flag(xS, yS), ae(v));
            inferSbDcSigCoeffFlag = 1;
        }

        auto const s2 = s << 2;

        if (i && !h[coded_sub_block_flag(offset - 8, 0)])
            continue;


        residualCodingState.nSigCoeffs = 0;
        residualCodingState.remaining = ~0;

        int max = 15;
        if (i == lastSubBlock)
        {
            residualCodingState.absoluteCoefficients[residualCodingState.nSigCoeffs] = 1;
            residualCodingState.nn[residualCodingState.nSigCoeffs] = lastScanPos;
            residualCodingState.coeffPositions[residualCodingState.nSigCoeffs++] = s2 + sp2[lastScanPos];
            max = lastScanPos - 1;
        }

        // read sig_coeff_flag

        for (int n = max; n >= inferSbDcSigCoeffFlag; n--)
        {
            auto const c = s2 + sp2[n];
            h(sig_coeff_flag(c.x, c.y), ae(v));
            residualCodingState.absoluteCoefficients[residualCodingState.nSigCoeffs] = h[sig_coeff_flag(c.x, c.y)];
            inferSbDcSigCoeffFlag &= !h[sig_coeff_flag(c.x, c.y)];
            int mask = -h[sig_coeff_flag(c.x, c.y)];
            //if (std::is_same<typename H::Tag, Decode<void>>::value)
            residualCodingState.nn[residualCodingState.nSigCoeffs] = n;
            residualCodingState.coeffPositions[residualCodingState.nSigCoeffs] = c;
            residualCodingState.nSigCoeffs += h[sig_coeff_flag(c.x, c.y)];
        }
        residualCodingState.absoluteCoefficients[residualCodingState.nSigCoeffs] = 1;
        //if (std::is_same<typename H::Tag, Decode<void>>::value)
        residualCodingState.nn[residualCodingState.nSigCoeffs] = 0;
        residualCodingState.coeffPositions[residualCodingState.nSigCoeffs] = s2;
        residualCodingState.nSigCoeffs += inferSbDcSigCoeffFlag;

        if (!residualCodingState.nSigCoeffs)
            continue;

        auto pnn = &residualCodingState.nn[residualCodingState.nSigCoeffs];
        auto pCoeffPositions = &residualCodingState.coeffPositions[residualCodingState.nSigCoeffs];

        // read coeff_abs_level_greater1_flag
        int firstSigScanPos = 16;
        int lastSigScanPos = -1;
        h[numGreater1Flag()] = 0;
        int lastGreater1ScanPos = -1;
        int lastGreater1m = -1;

        int mm = -residualCodingState.nSigCoeffs;
        do
        {
            int m = residualCodingState.nSigCoeffs + mm;
            auto const c = pCoeffPositions[mm];
            int n = pnn[mm];
            auto const &xC = c.x;
            auto const &yC = c.y;

            if (h[numGreater1Flag()] < 8)
            {
                h(coeff_abs_level_greater1_flag(n), ae(v));
                h[numGreater1Flag()]++;

                if (h[coeff_abs_level_greater1_flag(n)])
                {
                    residualCodingState.absoluteCoefficients[m] = 2;
                    if (lastGreater1ScanPos == -1)
                    {
                        lastGreater1ScanPos = n;
                        lastGreater1m = m;
                    }
                }
                else
                {
                    residualCodingState.remaining &= ~(1 << m);
                }
            }
            else
            {
                escapeDataPresent = 1;
            }

            if (lastSigScanPos == -1)
                lastSigScanPos = n;

            firstSigScanPos = n;

        } while (++mm);

        int signHidden = 0;
        if (h[sign_data_hiding_enabled_flag()])
        {
            if (h[cu_transquant_bypass_flag()] ||
                (h[current(CuPredMode(rc->x0, rc->y0))] == MODE_INTRA  &&
                    h[implicit_rdpcm_enabled_flag()] && h[transform_skip_flag(rc->x0, rc->y0, rc->cIdx)] &&
                    predModeIntraIs10or26(*rc, h)) ||
                h[explicit_rdpcm_flag(rc->x0, rc->y0, rc->cIdx)])
            {
            }
            else
            {
                signHidden = lastSigScanPos - firstSigScanPos > 3;
            }
        }

        // read coeff_abs_level_greater2_flag
        if (lastGreater1ScanPos != -1)
        {
            h(coeff_abs_level_greater2_flag(lastGreater1ScanPos), ae(v));
            if (h[coeff_abs_level_greater2_flag(lastGreater1ScanPos)])
            {
                residualCodingState.absoluteCoefficients[lastGreater1m] = 3;
                escapeDataPresent = 1;
            }
            else
            {
                residualCodingState.remaining &= ~(1 << lastGreater1m);
            }
        }

        if (h[cabac_bypass_alignment_enabled_flag()] && escapeDataPresent)
        {
            assert(!"alignment not yet implemented");
        }

        mm = -residualCodingState.nSigCoeffs + signHidden;
        do
        {
            int n = std::is_same<typename H::Tag, Decode<void>>::value ? pnn[mm - signHidden] : -1;
            h(coeff_sign_flag(n), ae(v));
        } while (++mm);

        // read coeff_abs_level_remaining
        mm = -residualCodingState.nSigCoeffs;
        do
        {
            int m = residualCodingState.nSigCoeffs + mm;

            if ((residualCodingState.remaining >> m) & 1)
            {
                h(coeff_abs_level_remaining(-1), ae(v));
                residualCodingState.absoluteCoefficients[m] += h[coeff_abs_level_remaining(-1)];
                stateSubstream->cLastAbsLevel = residualCodingState.absoluteCoefficients[m];
            }
        } while (++mm);

        if (std::is_same<typename H::Tag, Decode<void>>::value)
        {
            int sumAbsLevel = 0;
            mm = -residualCodingState.nSigCoeffs;
            do
            {
                int m = residualCodingState.nSigCoeffs + mm;
                auto const c = pCoeffPositions[mm];
                int n = pnn[mm];
                auto const &xC = c.x;
                auto const &yC = c.y;

                const int baseLevel = residualCodingState.absoluteCoefficients[m];
                sumAbsLevel ^= (baseLevel);

                if (signHidden)
                {
                    if (n == firstSigScanPos)
                        h[coeff_sign_flag(n)] = sumAbsLevel & 1;
                }

                h[TransCoeffLevel(rc->x0, rc->y0, rc->cIdx, xC, yC)] = baseLevel * (1 - 2 * h[coeff_sign_flag(n)]);

            } while (++mm);
        }
    }
}
template <>
struct Read<residual_coding>
{
    template <class H> static void go(const residual_coding &rc, H &hParent)
    {
        static void(*const read[6])(H &) =
        {
            0,
            0,
            optimisedReadResidualCoding<2, H>,
            optimisedReadResidualCoding<3, H>,
            optimisedReadResidualCoding<4, H>,
            optimisedReadResidualCoding<5, H>,
        };
        read[rc.log2TrafoSize](hParent);
     }
 };



template <class V>
struct DecisionLength1
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<V>(&binVal, 0));

        Access<V, H>::set(fun.v, binVal, h);
    }
};


template <>
struct Read<Element<sao_merge_left_flag, ae>>
{
    template <class H> static void go(Element<sao_merge_left_flag, ae> fun, H &h)
    {
        DecisionLength1<sao_merge_left_flag>::go(fun, h);

        if (h[fun.v])
        {
            sao *s = h;
            StateSpatial *stateSpatial = h;
            const SaoCtuData &saoCtuDataUp = stateSpatial->snakeSaoCtuData.at(s->rx - 1, s->ry, 0);
            stateSpatial->snakeSaoCtuData.commit(saoCtuDataUp, s->rx, s->ry, 0);
        }
    }
};


template <>
struct Read<Element<sao_merge_up_flag, ae>>
{
    template <class H> static void go(Element<sao_merge_up_flag, ae> fun, H &h)
    {
        DecisionLength1<sao_merge_up_flag>::go(fun, h);

        if (h[fun.v])
        {
            sao *s = h;
            StateSpatial *stateSpatial = h;
            const SaoCtuData &saoCtuDataUp = stateSpatial->snakeSaoCtuData.at(s->rx, s->ry - 1, 0);
            stateSpatial->snakeSaoCtuData.commit(saoCtuDataUp, s->rx, s->ry, 0);
        }
    }
};


template <class V>
struct ReadSaoTypeIdx
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        // TR cMax = 2, cRiceParam = 0

        int synVal = 0;
        int binVal0;
        h(DecodeDecision<V>(&binVal0, 0));

        if (binVal0)
        {
            int binVal1;
            h(DecodeBypass<V>(&binVal1));
            synVal = 1 + binVal1;
        }

        h[fun.v] = synVal;

        const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
        const int ry = h[yCtb()] >> h[CtbLog2SizeY()];

        if (std::is_same<V, sao_type_idx_luma>::value)
        {
            h[SaoTypeIdx(0, rx, ry)] = synVal;
        }
        else
        {
            h[SaoTypeIdx(1, rx, ry)] = synVal;
            h[SaoTypeIdx(2, rx, ry)] = synVal;
        }
    }
};

template <>
struct Read<Element<sao_type_idx_luma, ae>> : ReadSaoTypeIdx<sao_type_idx_luma> {};


template <>
struct Read<Element<sao_type_idx_chroma, ae>> : ReadSaoTypeIdx<sao_type_idx_chroma> {};


template <>
struct Read<Element<sao_offset_abs, ae>>
{
    template <class H> static void go(Element<sao_offset_abs, ae> fun, H &h)
    {
        const int bitDepth = fun.v.cIdx ? h[BitDepthC()] : h[BitDepthY()];
        const int cMax = (1 << (std::min(bitDepth, 10) - 5)) - 1;

        int binIdx;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            int value;
            h(DecodeBypass<sao_offset_abs>(&value));
            if (!value) break;
        }
        h[fun.v] = binIdx;
    }
};


template <class V, int fixedLength>
struct ReadFixedLengthBypass
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        int synVal = 0;

        for (int binIdx = 0; binIdx < fixedLength; ++binIdx)
        {
            synVal <<= 1;
            int binVal;
            h(DecodeBypass<V>(&binVal));
            synVal |= binVal;
        }
        h[fun.v] = synVal;
    }
};


template <>
struct Read<Element<sao_offset_sign, ae>> : ReadFixedLengthBypass<sao_offset_sign, 1> {};


template <>
struct Read<Element<sao_band_position, ae>> : ReadFixedLengthBypass<sao_band_position, 5> {};


template <>
struct Read<Element<sao_eo_class_luma, ae>> : ReadFixedLengthBypass<sao_eo_class_luma, 2>
{
    template <class H> static void go(Element<sao_eo_class_luma, ae> fun, H &h)
    {
        ReadFixedLengthBypass<sao_eo_class_luma, 2> ::go(fun, h);

        const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
        const int ry = h[yCtb()] >> h[CtbLog2SizeY()];
        h[SaoEoClass(0, rx, ry)] = h[fun.v];
    }
};


template <>
struct Read<Element<sao_eo_class_chroma, ae>> : ReadFixedLengthBypass<sao_eo_class_chroma, 2>
{
    template <class H> static void go(Element<sao_eo_class_chroma, ae> fun, H &h)
    {
        ReadFixedLengthBypass<sao_eo_class_chroma, 2> ::go(fun, h);

        const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
        const int ry = h[yCtb()] >> h[CtbLog2SizeY()];
        h[SaoEoClass(1, rx, ry)] = h[fun.v];
        h[SaoEoClass(2, rx, ry)] = h[fun.v];
    }
};


template <class V>
struct ReadTerminate
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        int value;
        h(DecodeTerminate<V>(&value));
        Access<V, H>::set(fun.v, value, h);
    }
};


template <>
struct Read<Element<end_of_slice_segment_flag, ae>>
{
    template <class H> static void go(Element<end_of_slice_segment_flag, ae> fun, H &h)
    {
        int value;
        h(DecodeTerminate<end_of_slice_segment_flag>(&value));
        h[fun.v] = value;

        if (value != 1 && h[CtbAddrInTs()] == h[PicSizeInCtbsY()] - 1)
        {
            h(Violation("8.1", "end_of_slice_segment_flag != 1 when CtbAddrInTs = PicSizeInCtbsY{%1%} - 1") % h[PicSizeInCtbsY()]);

            // Decoding has previously got out of sync.
            // Abort so that we do not attempt to decode CTUs outside of the picture boundary
            // which would cause access violations.

            throw Abort();
        }
    }
};


template <class H, class V>
Turing::Rectangle getIntraPartition(H &h, V &v)
{
    const coding_quadtree *cqt = h;
    Neighbourhood *neighbourhood = h;

    auto log2IntraPartSize = cqt->log2CbSize - h[IntraSplitFlag()];

    const Turing::Rectangle partition{ v.x0, v.y0, 1 << log2IntraPartSize, 1 << log2IntraPartSize };

    Snake<BlockData>::Cursor *cursor = h;
    cursor->relocate(neighbourhood->snake, partition, neighbourhood->MinCbLog2SizeYMinus1);

    return partition;
}


template <class H, class V>
void setIntraPartition(H &h, V &v, Turing::Rectangle &partition, int8_t mode, int pcm_flag)
{
    const coding_quadtree *cqt = h;
    Neighbourhood *neighbourhood = h;

    if (h[constrained_intra_pred_flag()] && v.x0 == cqt->x0 && v.y0 == cqt->y0)
    {
        Snake<BlockData>::Array<32, 32, 0, 0> *snakeArrayCuCip = h;
        snakeArrayCuCip->resize(*cqt, h[MinCbLog2SizeY()] - 1);
        //std::cout << *cqt << " a\n";
        //neighbourhood->snake.print(std::cout, *cqt, neighbourhood->MinCbLog2SizeYMinus1);
        snakeArrayCuCip->copyBlockFrom(neighbourhood->snake, before, *cqt, neighbourhood->MinCbLog2SizeYMinus1, 0, 0);
        //std::cout << *cqt << " b\n";
        //snakeArrayCuCip->print(std::cout, *cqt, neighbourhood->MinCbLog2SizeYMinus1);
    }

    Snake<BlockData>::Cursor *cursor = h;
    BlockData &blockData = cursor->current(0, 0, neighbourhood->MinCbLog2SizeYMinus1);
    blockData.reset();
    blockData.CtDepth = cqt->cqtDepth;
    blockData.skip = 0;
    blockData.intra.pcm = pcm_flag;
    blockData.intra.predModeY = mode;

    if (v.x0 == cqt->x0 || v.y0 == cqt->y0)
    {
        cursor->commit(partition, neighbourhood->MinCbLog2SizeYMinus1);
    }
}


template <>
struct Read < Element<pcm_flag, ae> >
{
    template <class H> static void go(Element<pcm_flag, ae> fun, H &h)
    {
        ReadTerminate<pcm_flag>::go(fun, h);

        if (h[fun.v] == 1)
        {
            Turing::Rectangle partition = getIntraPartition(h, fun.v);
            setIntraPartition(h, fun.v, partition, INTRA_DC, 1);
        }
    }
};


template <>
struct Read<Element<end_of_subset_one_bit, ae>>
{
    template <class H> static void go(Element<end_of_subset_one_bit, ae> fun, H &h)
    {
        ReadTerminate<end_of_subset_one_bit>::go(fun, h);

        if (h[fun.v] != 1)
        {
            h(Violation("7.4.9.1", "end_of_subset_one_bit{%1%} != 1") % h[fun.v]); // CondCheck 7.4.9.1-A
            throw Abort();
        }

        h(ContextsInitialize());
    }
};


template <>
struct Read<Element<split_cu_flag, ae>>
{
    template <class H> static void go(Element<split_cu_flag, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<split_cu_flag>(&binVal, GetCtxInc<H, split_cu_flag>::f(h, fun.v, 0)));
        h[fun.v] = binVal;
    }
};


template <>
struct Read<Element<cu_transquant_bypass_flag, ae>> : DecisionLength1<cu_transquant_bypass_flag> {};



template <>
struct Read<Element<cu_skip_flag, ae>>
{
    template <class H> static void go(Element<cu_skip_flag, ae> fun, H &h)
    {
        coding_quadtree *cqt = h;

        const int ctxInc =
            h[cu_skip_flag(fun.v.x0 - 1, fun.v.y0)] +
            h[cu_skip_flag(fun.v.x0, fun.v.y0 - 1)];

        int binVal;
        h(DecodeDecision<cu_skip_flag>(&binVal, ctxInc));

        Snake<BlockData>::Cursor *cursor = h;
        Neighbourhood *neighbourhood = h;

        cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1).reset();
        cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1).refIdxPlus1[0] = 1;
        cursor->current(cqt->x0, cqt->y0, neighbourhood->MinCbLog2SizeYMinus1).skip = binVal;

        assert(h[CuPredMode(fun.v.x0, fun.v.y0)] == (binVal ? MODE_SKIP : MODE_INTER));
    }
};


template <>
struct Read<Element<pred_mode_flag, ae>>
{
    template <class H> static void go(Element<pred_mode_flag, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<pred_mode_flag>(&binVal, 0));
        h[fun.v] = binVal;

        Snake<BlockData>::Cursor *cursor = h;
        BlockData &blockData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);
        blockData.refIdxPlus1[0] = binVal ? 0 : 1;
    }
};


template <>
struct Read<Element<part_mode, ae>>
{
    template <class H> static void go(Element<part_mode, ae> fun, H &h)
    {
        const coding_quadtree *cqt = h;

        int binVal0;
        h(DecodeDecision<part_mode>(&binVal0, 0));

        int synVal;
        if (binVal0)
        {
            // bins: "1"
            synVal = 0;
        }
        else if (h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA)
        {
            Neighbourhood *neighbourhood = h;
            Snake<BlockData>::Cursor *cursor = h;
            BlockData &blockData = cursor->current(0, 0, neighbourhood->MinCbLog2SizeYMinus1);

            // bins: "0"
            synVal = 1;
        }
        else
        {
            int x;
            h(DecodeDecision<part_mode>(&x, 1));
            // bins: "0x"
            if (cqt->log2CbSize > h[MinCbLog2SizeY()])
            {
                if (!h[amp_enabled_flag()])
                {
                    synVal = x ? 1 : 2;
                }
                else
                {
                    const int ctxInc = 2 + (cqt->log2CbSize == h[MinCbLog2SizeY()] ? 0 : 1);
                    int binVal2;
                    h(DecodeDecision<part_mode>(&binVal2, ctxInc));
                    if (binVal2)
                    {
                        // bins: "0x1"
                        synVal = x ? 1 : 2;
                    }
                    else
                    {
                        // bins: "0x0"
                        int binVal3;
                        h(DecodeBypass<part_mode>(&binVal3));
                        if (binVal3)
                        {
                            // bins: "0x01"
                            synVal = x ? 5 : 7;
                        }
                        else
                        {
                            // bins: "0x00"
                            synVal = x ? 4 : 6;
                        }
                    }
                }
            }
            else if (x)
            {
                // bins: "01"
                synVal = 1;
            }
            else
            {
                // bins: "00"
                if (cqt->log2CbSize == 3)
                {
                    synVal = 2;
                }
                else
                {
                    assert(cqt->log2CbSize > 3);

                    const int ctxInc = 2 + (cqt->log2CbSize == h[MinCbLog2SizeY()] ? 0 : 1);
                    int binVal2;
                    h(DecodeDecision<part_mode>(&binVal2, ctxInc));

                    synVal = binVal2 ? 2 : 3;
                }
            }
        }
        h[fun.v] = synVal;
    }
};


template <>
struct Read<Element<prev_intra_luma_pred_flag, ae>> : DecisionLength1<prev_intra_luma_pred_flag> {};


template <>
struct Read<Element<mpm_idx, ae>>
{
    template <class H> static void go(Element<mpm_idx, ae> fun, H &h)
    {
        const int cMax = 2;
        int binIdx;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            int binVal;
            h(DecodeBypass<mpm_idx>(&binVal));
            if (!binVal) break;
        }
        const int synVal = binIdx;
        h[fun.v] = synVal;

        Turing::Rectangle partition = getIntraPartition(h, fun.v);

        CandModeList candModeList;
        candModeList.populate(0, h, fun.v.x0, fun.v.y0);

        const auto mode = candModeList[synVal];

        setIntraPartition(h, fun.v, partition, candModeList[synVal], 0);
    }
};


template <>
struct Read<Element<rem_intra_luma_pred_mode, ae>>
{
    template <class H> static void go(Element<rem_intra_luma_pred_mode, ae> fun, H &h)
    {
        int synVal = 0;
        for (int binIdx = 0; binIdx < 5; ++binIdx)
        {
            synVal <<= 1;
            int binVal;
            h(DecodeBypass<rem_intra_luma_pred_mode>(&binVal));
            synVal |= binVal;
        }
        h[fun.v] = synVal;

        Turing::Rectangle partition = getIntraPartition(h, fun.v);

        CandModeList candModeList;
        candModeList.populate(0, h, fun.v.x0, fun.v.y0);

        candModeList.sort();

        auto mode = synVal;
        for (int i = 0; i < 3; ++i)
        {
            if (mode >= candModeList[i]) ++mode;
        }

        setIntraPartition(h, fun.v, partition, mode, 0);
    }
};


template <>
struct Read<Element<intra_chroma_pred_mode, ae>>
{
    template <class H> static void go(Element<intra_chroma_pred_mode, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<intra_chroma_pred_mode>(&binVal, 0));

        if (!binVal)
        {
            h[fun.v] = 4;
        }
        else
        {
            int val;
            h(DecodeBypass<intra_chroma_pred_mode>(&val));
            val <<= 1;
            int binVal;
            h(DecodeBypass<intra_chroma_pred_mode>(&binVal));
            val |= binVal;
            h[fun.v] = val;
        }
    }
};


template <>
struct Read<Element<rqt_root_cbf, ae>> : DecisionLength1<rqt_root_cbf> { };


// compute and store loop filter metadata for specified entity
template <class F> struct Process { };


template <>
struct Read<Process<prediction_unit>> : Null<Process<prediction_unit>> { };


template <>
struct Read<prediction_unit>
{
    template <class H> static void go(const prediction_unit &pu, H &h)
    {
        Snake<BlockData>::Cursor *cursor = h;
        Neighbourhood *neighbourhood = h;
        StateSubstream *stateSubstream = h;
        coding_quadtree *cqt = h;

        const int xPb = pu.x0;
        const int yPb = pu.y0;

        cursor->relocate(neighbourhood->snake, pu, h[MinCbLog2SizeY()] - 1);
        cursor->current(0, 0, h[MinCbLog2SizeY()] - 1).CtDepth = cqt->cqtDepth;

        h[merge_flag(xPb, yPb)] = h[current(CuPredMode(xPb, yPb))] == MODE_SKIP;
        h[merge_idx(xPb, yPb)] = 0;

        Syntax<prediction_unit>::go(pu, h);

        PuData &puData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);

        processPredictionUnit(h, pu, puData, stateSubstream->partIdx, h[merge_idx(xPb, yPb)]);

        h(Process<prediction_unit>());

        cursor->commit(pu, h[MinCbLog2SizeY()] - 1);

        commitPu(h, pu, puData);

        ++stateSubstream->partIdx;
    }
};



template <>
struct Read<Element<merge_idx, ae>>
{
    template <class H> static void go(Element<merge_idx, ae> fun, H &h)
    {
        const int cMax = h[MaxNumMergeCand()] - 1;

        int binIdx;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            if (binIdx == 0)
            {
                int binVal;
                h(DecodeDecision<merge_idx>(&binVal, 0));
                if (!binVal) break;
            }
            else
            {
                int binVal;
                h(DecodeBypass<merge_idx>(&binVal));
                if (!binVal) break;
            }
        }
        const int synVal = binIdx;
        h[fun.v] = synVal;
    }
};


template <>
struct Read<Element<merge_flag, ae>> : DecisionLength1<merge_flag> { };


template <>
struct Read<Element<inter_pred_idc, ae>>
{
    template <class H> static void go(Element<inter_pred_idc, ae> fun, H &h)
    {
        const prediction_unit &pu = *static_cast<prediction_unit *>(h);
        coding_quadtree cqt = *static_cast<coding_quadtree *>(h);

        int synVal;

        if (pu.nPbW + pu.nPbH != 12)
        {
            int t;
            h(DecodeDecision<inter_pred_idc>(&t, cqt.cqtDepth));
            if (t)
            {
                synVal = PRED_BI;
            }
            else
            {
                h(DecodeDecision<inter_pred_idc>(&synVal, 4));
            }
        }
        else
        {
            h(DecodeDecision<inter_pred_idc>(&synVal, 4));
        }

        h[fun.v] = synVal;
    }
};



template <>
struct Read<Element<cbf_cb, ae>>
{
    template <class H> static void go(Element<cbf_cb, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<cbf_cb>(&binVal, fun.v.trafoDepth));
        Access<cbf_cb, H>::set(fun.v, !!binVal, h);
    }
};



template <>
struct Read<Element<cbf_cr, ae>>
{
    template <class H> static void go(Element<cbf_cr, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<cbf_cr>(&binVal, fun.v.trafoDepth));
        Access<cbf_cr, H>::set(fun.v, !!binVal, h);
    }
};


template <>
struct Read<Element<cbf_luma, ae>>
{
    template <class H> static void go(Element<cbf_luma, ae> fun, H &h)
    {
        transform_tree const *tt = h;
        auto const ctxInc = tt->trafoDepth == 0 ? 1 : 0;
        int binVal;
        h(DecodeDecision<cbf_luma>(&binVal, ctxInc));
        Access<cbf_luma, H>::set(fun.v, binVal, h);
    }
};


template <>
struct Read<mvd_coding>
{
    template <class H> static void go(const mvd_coding &s, H &h)
    {
        h[abs_mvd_greater1_flag(0)] = 0;
        h[abs_mvd_greater1_flag(1)] = 0;
        h[abs_mvd_minus2(0)] = -1;
        h[abs_mvd_minus2(1)] = -1;
        h[mvd_sign_flag(0)] = 0;
        h[mvd_sign_flag(1)] = 0;

        Syntax<mvd_coding>::go(s, h);

        h[Mvd(s.refList, s.x0, s.y0)][0] = h[lMvd(0)];
        h[Mvd(s.refList, s.x0, s.y0)][1] = h[lMvd(1)];
    }
};


template <>
struct Read<Element<abs_mvd_greater0_flag, ae>> : DecisionLength1<abs_mvd_greater0_flag> {};


template <>
struct Read<Element<abs_mvd_greater1_flag, ae>> : DecisionLength1<abs_mvd_greater1_flag> {};


template <>
struct Read<Element<mvp_l0_flag, ae>> : DecisionLength1<mvp_l0_flag> {};


template <>
struct Read<Element<mvp_l1_flag, ae>> : DecisionLength1<mvp_l1_flag> {};



template <class V, class NumRefIdxActiveMinus1>
struct ReadRefIdx
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        const int cMax = h[NumRefIdxActiveMinus1()];
        int binIdx;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            int binVal;
            switch (binIdx)
            {
            case 0:
            case 1:
                h(DecodeDecision<V>(&binVal, binIdx));
                break;
            default:
                h(DecodeBypass<V>(&binVal));
                break;
            }
            if (!binVal) break;
        }
        h[fun.v] = binIdx;
    }
};

template <>
struct Read<Element<ref_idx_l0, ae>> :
    ReadRefIdx<ref_idx_l0, num_ref_idx_l0_active_minus1>
{
};


template <>
struct Read<Element<ref_idx_l1, ae>> :
    ReadRefIdx<ref_idx_l1, num_ref_idx_l1_active_minus1>
{
};


template <>
struct Read<Element<abs_mvd_minus2, ae>>
{
    template <class H> static void go(Element<abs_mvd_minus2, ae> fun, H &h)
    {
        const int k = 1;
        int count = k;
        unsigned u1 = 0, u2 = 0;
        int binVal;

        do
        {
            h(DecodeBypass<abs_mvd_minus2>(&binVal));
            u1 += binVal << count++;
        } while (binVal);

        while (--count)
        {
            u2 <<= 1;
            h(DecodeBypass<abs_mvd_minus2>(&binVal));
            u2 += binVal;
        }
        h[fun.v] = u1 + u2;
    }
};



template <>
struct Read<Element<mvd_sign_flag, ae>> : ReadFixedLengthBypass<mvd_sign_flag, 1> {};


template <>
struct Read<Element<cu_qp_delta_abs, ae>>
{
    template <class H> static void go(Element<cu_qp_delta_abs, ae> fun, H &h)
    {
        int binIdx;
        const int cMax = 5;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            int binVal;
            h(DecodeDecision<cu_qp_delta_abs>(&binVal, binIdx ? 1 : 0));
            if (!binVal) break;
        }

        const int prefix = binIdx;
        int suffix = 0;

        if (prefix > 4)
        {
            const int k = 0;
            int count = k;
            unsigned u1 = 0, u2 = 0;

            int binVal;
            do
            {
                h(DecodeBypass<cu_qp_delta_abs>(&binVal));
                u1 += binVal << count++;
            } while (binVal);

            while (--count)
            {
                u2 <<= 1;
                h(DecodeBypass<cu_qp_delta_abs>(&binVal));
                u2 += binVal;
            }
            suffix = u1 + u2;
        }
        h[fun.v] = prefix + suffix;

        h[IsCuQpDeltaCoded()] = 1;
    }
};




template <>
struct Read<Element<cu_qp_delta_sign_flag, ae>>
{
    template <class H> static void go(Element<cu_qp_delta_sign_flag, ae> fun, H &h)
    {
        h(DecodeBypass<cu_qp_delta_sign_flag>(&h[fun.v]));

        h[CuQpDeltaVal()] = h[cu_qp_delta_abs()] * (1 - 2 * h[cu_qp_delta_sign_flag()]);

        h[QpY()] = ((h[QpY()] + h[CuQpDeltaVal()] + 52 + 2 * h[QpBdOffsetY()]) % (52 + h[QpBdOffsetY()])) - h[QpBdOffsetY()];

        QpState *qpState = h;
        qpState->update();
    }
};


template <>
struct Read<Element<cu_chroma_qp_offset_flag, ae>>
{
    template <class H> static void go(Element<cu_chroma_qp_offset_flag, ae> fun, H &h)
    {
        h(DecodeDecision<cu_chroma_qp_offset_flag>(&h[fun.v], 0));
        h[IsCuChromaQpOffsetCoded()] = 1;
        if (h[fun.v])
        {
            h[CuQpOffsetCb()] = h[cb_qp_offset_list(0)];
            h[CuQpOffsetCr()] = h[cr_qp_offset_list(0)];
        }
        else
        {
            h[CuQpOffsetCb()] = 0;
            h[CuQpOffsetCr()] = 0;
        }
        QpState *qpState = h;
        qpState->update();
    }
};


template <>
struct Read<Element<cu_chroma_qp_offset_idx, ae>>
{
    template <class H> static void go(Element<cu_chroma_qp_offset_idx, ae> fun, H &h)
    {
        int const cMax = h[chroma_qp_offset_list_len_minus1()];

        int binIdx;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            int binVal;
            h(DecodeDecision<cu_chroma_qp_offset_idx>(&binVal, 0));
            if (!binVal) 
                break;
        }
        h[fun.v] = binIdx;
        h[CuQpOffsetCb()] = h[cb_qp_offset_list(h[fun.v])];
        h[CuQpOffsetCr()] = h[cr_qp_offset_list(h[fun.v])];
        
        QpState *qpState = h;
        qpState->update();
    }
};


template <>
struct Read<Element<log2_res_scale_abs_plus1, ae>>
{
    template <class H> static void go(Element<log2_res_scale_abs_plus1, ae> fun, H &h)
    {
        const int cMax = 4;
        int binIdx;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            int binVal;
            h(DecodeDecision<log2_res_scale_abs_plus1>(&binVal, 4 * fun.v.c + binIdx));
            if (!binVal) break;
        }
        h[fun.v] = binIdx;
    }
};


template <>
struct Read<Element<res_scale_sign_flag, ae>>
{
    template <class H> static void go(Element<res_scale_sign_flag, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<res_scale_sign_flag>(&binVal, fun.v.c));
        h[fun.v] = binVal;
    }
};


template <>
struct Read<Element<split_transform_flag, ae>>
{
    template <class H> static void go(Element<split_transform_flag, ae> fun, H &h)
    {
        transform_tree &tt = *static_cast<transform_tree *>(h);
        int t;
        h(DecodeDecision<split_transform_flag>(&t, 5 - tt.log2TrafoSize));
        h[fun.v] = t;
    }
};


template <>
struct Read<Element<transform_skip_flag, ae>>
{
    template <class H> static void go(Element<transform_skip_flag, ae> fun, H &h)
    {
        int synVal;
        if (fun.v.cIdx == 0)
        {
            h(DecodeDecision<transform_skip_flag_0>(&synVal, 0));
        }
        else
        {
            h(DecodeDecision<transform_skip_flag_1>(&synVal, 0));
        }
        h[fun.v] = synVal;
    }
};


template <>
struct Read<Element<explicit_rdpcm_flag, ae>>
{
    template <class H> static void go(Element<explicit_rdpcm_flag, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<explicit_rdpcm_flag>(&binVal, fun.v.cIdx ? 1 : 0));
        h[fun.v] = binVal;
    }
};


template <>
struct Read<Element<explicit_rdpcm_dir_flag, ae>>
{
    template <class H> static void go(Element<explicit_rdpcm_dir_flag, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<explicit_rdpcm_dir_flag>(&binVal, fun.v.cIdx ? 1 : 0));
        h[fun.v] = binVal;
    }
};


template <class V>
struct LastSigCoeffPrefix
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        residual_coding &rc = *static_cast<residual_coding *>(h);
        const int cMax = (rc.log2TrafoSize << 1) - 1;
        int binIdx;
        for (binIdx = 0; binIdx < cMax; ++binIdx)
        {
            int ctxOffset, ctxShift;
            if (rc.cIdx == 0)
            {
                ctxOffset = 3 * (rc.log2TrafoSize - 2) + ((rc.log2TrafoSize - 1) >> 2);
                ctxShift = (rc.log2TrafoSize + 1) >> 2;
            }
            else
            {
                ctxOffset = 15;
                ctxShift = rc.log2TrafoSize - 2;
            }
            const int ctxIdx = (binIdx >> ctxShift) + ctxOffset;
            int binVal;
            h(DecodeDecision<V>(&binVal, ctxIdx));
            if (!binVal) break;
        }
        h[fun.v] = binIdx;
    }
};


template <>
struct Read<Element<last_sig_coeff_x_prefix, ae>> :
    LastSigCoeffPrefix<last_sig_coeff_x_prefix>
{
};


template <class Prefix, class Suffix, class H>
int computeLast(H &h)
{
    const bool present = h[Prefix()] > 3;
    if (!present)
    {
        return h[Prefix()];
    }
    else
    {
        return
            (1 << ((h[Prefix()] >> 1) - 1))
            * (2 + (h[Prefix()] & 1))
            + h[Suffix()];
    }
}


template <class H>
void setLast(H &h)
{
    if (h[scanIdx()] == 2)
    {
        h[LastSignificantCoeffX()] = computeLast<last_sig_coeff_y_prefix, last_sig_coeff_y_suffix>(h);
        h[LastSignificantCoeffY()] = computeLast<last_sig_coeff_x_prefix, last_sig_coeff_x_suffix>(h);
    }
    else
    {
        h[LastSignificantCoeffX()] = computeLast<last_sig_coeff_x_prefix, last_sig_coeff_x_suffix>(h);
        h[LastSignificantCoeffY()] = computeLast<last_sig_coeff_y_prefix, last_sig_coeff_y_suffix>(h);
    }
    int v = 1;
    Access<coded_sub_block_flag, H>::set(coded_sub_block_flag(h[LastSignificantCoeffX()] >> 2, h[LastSignificantCoeffY()] >> 2), v, h);
}

template <>
struct Read<Element<last_sig_coeff_y_prefix, ae>> :
    LastSigCoeffPrefix<last_sig_coeff_y_prefix>
{
    template <class H> static void go(Element<last_sig_coeff_y_prefix, ae> fun, H &h)
    {
        LastSigCoeffPrefix<last_sig_coeff_y_prefix>::go(fun, h);
        if (h[last_sig_coeff_x_prefix()] <= 3 && h[last_sig_coeff_y_prefix()] <= 3) setLast(h);
    }
};


template <class V, class Prefix>
struct LastSigCoeffSuffix
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        const int fixedLength = (h[Prefix()] - 2) >> 1;
        int synVal = 0;
        for (int binIdx = 0; binIdx < fixedLength; ++binIdx)
        {
            synVal <<= 1;
            int binVal;
            h(DecodeBypass<V>(&binVal));
            synVal |= binVal;
        }
        h[fun.v] = synVal;
    }
};


template <>
struct Read<Element<last_sig_coeff_x_suffix, ae>> :
    LastSigCoeffSuffix<last_sig_coeff_x_suffix, last_sig_coeff_x_prefix>
{
    template <class H> static void go(Element<last_sig_coeff_x_suffix, ae> fun, H &h)
    {
        LastSigCoeffSuffix<last_sig_coeff_x_suffix, last_sig_coeff_x_prefix>::go(fun, h);
        if (h[last_sig_coeff_y_prefix()] <= 3) setLast(h);
    }
};


template <>
struct Read<Element<last_sig_coeff_y_suffix, ae>> :
    LastSigCoeffSuffix<last_sig_coeff_y_suffix, last_sig_coeff_y_prefix>
{
    template <class H> static void go(Element<last_sig_coeff_y_suffix, ae> fun, H &h)
    {
        LastSigCoeffSuffix<last_sig_coeff_y_suffix, last_sig_coeff_y_prefix>::go(fun, h);
        setLast(h);
    }
};


template <>
struct Read<Element<coded_sub_block_flag, ae>>
{
    template <class H> static void go(Element<coded_sub_block_flag, ae> fun, H &h)
    {
        int binVal;
        h(DecodeDecision<coded_sub_block_flag>(&binVal, ctxInc(h, fun.v)));
        Access<coded_sub_block_flag, H>::set(fun.v, binVal, h);
        if (binVal)
        {
            StateSubstream& stateSubstream = *static_cast<StateSubstream *>(h);
            stateSubstream.cLastAbsLevel = 0;
            stateSubstream.firstCoeffAbsLevelRemainingInSubblock = true;
        }
    }
    template <class H> static int ctxInc(H& h, coded_sub_block_flag v)
    {
        const residual_coding& rc = *static_cast<residual_coding *>(h);

        int csbfCtx = 0;
        if (v.xS < (1 << (rc.log2TrafoSize - 2)) - 1)
        {
            csbfCtx += h[coded_sub_block_flag(v.xS + 1, v.yS)];
        }
        if (v.yS < (1 << (rc.log2TrafoSize - 2)) - 1)
        {
            csbfCtx += h[coded_sub_block_flag(v.xS, v.yS + 1)];
        }
        return (rc.cIdx ? 2 : 0) + std::min(csbfCtx, 1);
    }
};



template <>
struct Read<Element<sig_coeff_flag, ae>>
{
    template <class H> static void go(Element<sig_coeff_flag, ae> fun, H &h)
    {
        residual_coding &rc = *static_cast<residual_coding *>(h);
        int binVal;
        h(DecodeDecision<sig_coeff_flag>(&binVal, ctxInc(h, rc.cIdx, fun.v.xC, fun.v.yC, h[scanIdx()], rc.log2TrafoSize)));
        h[fun.v] = binVal;

        const bool nIsZero = !((fun.v.xC & 3) || (fun.v.xC & 3));
        if (nIsZero)
        {
            StateSubstream& stateSubstream = *static_cast<StateSubstream *>(h);
            stateSubstream.cLastAbsLevel = 0;
            stateSubstream.firstCoeffAbsLevelRemainingInSubblock = true;
        }
    }
    // review: this duplicated in Binarization
    template <class H> static int ctxInc(H &h, int cIdx, int xC, int yC, int scanIdx, int log2TrafoSize)
    {
        //flog() << "ctxInc(" << cIdx << ", " << xC << ", " << yC << ", " << scanIdx << ", " << log2TrafoSize << "\n";
        int sigCtx;
        if (h[transform_skip_context_enabled_flag()] && (h[transform_skip_flag(xC, yC, cIdx)] || h[cu_transquant_bypass_flag()]))
        {
            sigCtx = (cIdx == 0) ? 42 : 16;
        }
        else if (log2TrafoSize == 2)
        {
            const int ctxIdxMap[16] =
            {
                0, 1, 4, 5,
                2, 3, 4, 5,
                6, 6, 8, 8,
                7, 7, 8, 8
            };
            sigCtx = ctxIdxMap[(yC << 2) + xC];
        }
        else if (xC + yC == 0)
        {
            sigCtx = 0;
        }
        else
        {
            const int xS = xC >> 2;
            const int yS = yC >> 2;
            int prevCsbf = 0;

            if (xS < (1 << (log2TrafoSize - 2)) - 1)
            {
                prevCsbf += h[coded_sub_block_flag(xS + 1, yS)];
            }
            if (yS < (1 << (log2TrafoSize - 2)) - 1)
            {
                prevCsbf += (h[coded_sub_block_flag(xS, yS + 1)] << 1);
            }

            const int xP = xC & 3;
            const int yP = yC & 3;

            if (prevCsbf == 0)
            {
                sigCtx = (xP + yP == 0) ? 2 : (xP + yP < 3) ? 1 : 0;
            }
            else if (prevCsbf == 1)
            {
                sigCtx = (yP == 0) ? 2 : (yP == 1) ? 1 : 0;
            }
            else if (prevCsbf == 2)
            {
                sigCtx = (xP == 0) ? 2 : (xP == 1) ? 1 : 0;
            }
            else
            {
                assert(prevCsbf == 3);
                sigCtx = 2;
            }

            if (cIdx == 0)
            {
                if (xS + yS > 0)
                {
                    sigCtx += 3;
                }

                if (log2TrafoSize == 3)
                {
                    sigCtx += (scanIdx == 0) ? 9 : 15;
                }
                else
                {
                    sigCtx += 21;
                }
            }
            else
            {
                if (log2TrafoSize == 3)
                {
                    sigCtx += 9;
                }
                else
                {
                    sigCtx += 12;
                }
            }
        }

        assert(sigCtx >= 0);
        assert(sigCtx < 54);

        if (cIdx == 0)
        {
            return sigCtx;
        }
        else
        {
            return 27 + sigCtx;
        }
    }
};


template <>
struct Read<Element<coeff_abs_level_greater1_flag, ae>>
{
    template <class H> static void go(Element<coeff_abs_level_greater1_flag, ae> fun, H &h)
    {
        StateSubstream& stateSubstream = *static_cast<StateSubstream *>(h);
        const residual_coding& rc = *static_cast<residual_coding *>(h);

        auto& ctxSet = stateSubstream.ctxSet;
        auto& lastGreater1Ctx = stateSubstream.lastGreater1Ctx;
        auto& greater1Ctx = stateSubstream.greater1Ctx;
        auto& lastGreater1Flag = stateSubstream.lastGreater1Flag;
        auto& previousGreater1Flag = stateSubstream.previousGreater1Flag;
        auto& numGreater1Flag = h[::numGreater1Flag()];
        auto& i = stateSubstream.i;

        if (numGreater1Flag == 0)
        {
            ctxSet = (i == 0 || rc.cIdx > 0) ? 0 : 2;
            if (lastGreater1Ctx < 0)
            {
                lastGreater1Ctx = 1;
            }
            else
            {
                lastGreater1Ctx = greater1Ctx;
                if (lastGreater1Ctx > 0)
                {
                    lastGreater1Flag = previousGreater1Flag;
                    lastGreater1Ctx = lastGreater1Flag ? 0 : lastGreater1Ctx + 1;
                }
            }
            if (lastGreater1Ctx == 0)
            {
                ctxSet = ctxSet + 1;
                assert(ctxSet < 4);
            }
            greater1Ctx = 1;
        }
        else
        {
            if (greater1Ctx > 0)
            {
                lastGreater1Flag = previousGreater1Flag;
            }
            greater1Ctx = lastGreater1Flag ? 0 : greater1Ctx + 1;
        }

        const int ctxInc = (ctxSet * 4) + std::min(3, greater1Ctx) + (rc.cIdx ? 16 : 0);

        int binVal;
        h(DecodeDecision<coeff_abs_level_greater1_flag>(&binVal, ctxInc));
        h[fun.v] = binVal;

        previousGreater1Flag = binVal;
    }
};


template <>
struct Read<Element<coeff_abs_level_greater2_flag, ae>>
{
    template <class H> static void go(Element<coeff_abs_level_greater2_flag, ae> fun, H &h)
    {
        StateSubstream& stateSubstream = *static_cast<StateSubstream *>(h);
        const residual_coding& rc = *static_cast<residual_coding *>(h);

        auto& ctxSet = stateSubstream.ctxSet;
        const int ctxInc = ctxSet + (rc.cIdx ? 4 : 0);
        int binVal;
        h(DecodeDecision<coeff_abs_level_greater2_flag>(&binVal, ctxInc));
        h[fun.v] = binVal;
    }
};


template <>
struct Read<Element<coeff_sign_flag, ae>> : ReadFixedLengthBypass<coeff_sign_flag, 1> {};


template <>
struct Read<Element<coeff_abs_level_remaining, ae>>
{
    template <class H> static void go(Element<coeff_abs_level_remaining, ae> fun, H &h)
    {
        StateSubstream& stateSubstream = *static_cast<StateSubstream *>(h);

        if (stateSubstream.firstCoeffAbsLevelRemainingInSubblock)
        {
            stateSubstream.cLastRiceParam = h[initRiceValue()];
        }

        int cRiceParam = stateSubstream.cLastRiceParam + (stateSubstream.cLastAbsLevel > (3 * (1 << stateSubstream.cLastRiceParam)) ? 1 : 0);
        if (!h[persistent_rice_adaptation_enabled_flag()])
        {
            cRiceParam = std::min(cRiceParam, 4);
        }
        stateSubstream.cLastRiceParam = cRiceParam;

        int synVal;

        int prefix = 0;
        int codeWord = 0;

        do
        {
            ++prefix;
            h(DecodeBypass<coeff_abs_level_remaining>(&codeWord));
        } while (codeWord);

        codeWord = 1 - codeWord;
        prefix -= codeWord;
        codeWord = 0;

        if (prefix < 3)
        {
            for (int i = 0; i < cRiceParam; ++i)
            {
                codeWord <<= 1;
                int binVal;
                h(DecodeBypass<coeff_abs_level_remaining>(&binVal));
                codeWord |= binVal;
            }
            synVal = (prefix << cRiceParam) + codeWord;
        }
        else
        {
            for (int i = 0; i < prefix - 3 + cRiceParam; ++i)
            {
                codeWord <<= 1;
                int binVal;
                h(DecodeBypass<coeff_abs_level_remaining>(&binVal));
                codeWord |= binVal;
            }
            synVal = (((1 << (prefix - 3)) + 3 - 1) << cRiceParam) + codeWord;
        }
        h[fun.v] = synVal;

        stateSubstream.cLastAbsLevel = 1 + h[coeff_abs_level_greater1_flag(fun.v.n)] + h[coeff_abs_level_greater2_flag(fun.v.n)] + synVal;

        if (stateSubstream.firstCoeffAbsLevelRemainingInSubblock)
        {
            stateSubstream.firstCoeffAbsLevelRemainingInSubblock = false;
            if (h[persistent_rice_adaptation_enabled_flag()])
            {
                const int sbType = h[::sbType()];

                if (synVal >= (3 << (h[StatCoeff(sbType)] / 4)))
                {
                    h[StatCoeff(sbType)]++;
                }
                else if (2 * synVal < (1 << (h[StatCoeff(sbType)] / 4)) && h[StatCoeff(sbType)] > 0)
                {
                    h[StatCoeff(sbType)]--;
                }
            }
        }
    }
};


struct UnexpectedData { };


// Eats all remaining data in stream
template <class F>
struct Greedy
{
    template <class H> static void go(const F &fun, H &h)
    {
        seek(h, bitLen(h[Stream()]));
    }
};


template <>
struct Read<UnexpectedData> : Greedy<UnexpectedData> { };


template <>
struct Read<Bitstream>
{
    template <class H> static void go(Bitstream, H &h)
    {
        try
        {
            while (!h[Stop()] && !endOfStream(h[Stream()]))
            {
                h(CodedVideoSequence());
            }
        }
        catch (ExceptionOverrun &)
        {
        }

        auto statePictures = &h[Concrete<StatePicturesBase>()];
        // flush any remaining pictures from DPB
        while (!h[::Stop()] && statePictures->bumpingProcess(h, true));
    }
};


template <>
struct Read<CodedVideoSequence>
{
    template <class H> static void go(CodedVideoSequence, H &h)
    {
        StatePictures *statePictures = h;

        h[Active<Sps>()].reset();
        h[Active<Vps>()].reset();

        ++statePictures->codedVideoSequenceId;
        statePictures->posEndCvs = h[::Stream()].len;

        while (!h[::Stop()] && h[::Stream()].state.pos < statePictures->posEndCvs)
        {
            auto newPicture = h[NewPicture()];
            auto hPicture = h.extend(newPicture.get());

            hPicture[Active<Vps>()] = h[Active<Vps>()];
            hPicture[Active<Sps>()] = h[Active<Sps>()];

            hPicture(AccessUnit());

            h[Active<Vps>()] = hPicture[Active<Vps>()];
            h[Active<Sps>()] = hPicture[Active<Sps>()];
        }
    }
};


struct NalUnitInfo
{
    template <class Stream>
    NalUnitInfo(Stream &stream)
    {
        gotoNextStartCode(stream);

        this->pos = bitPos<size_t>(stream) / 8;

        if (stream.readBytes(3) == 0x000000)
        {
            stream.readBytes(1);
        }

        // simple nal_unit_header() lookahead
        this->nal_unit_type = (stream.readBytes(2) >> 9) & 0x3f;

        if (isSliceSegment(this->nal_unit_type))
        {
            const bool first_slice_segment_in_pic_flag = !!(stream.byte() & 0x80);
            this->firstSliceSegment = first_slice_segment_in_pic_flag;
        }
        else
        {
            this->firstSliceSegment = false;
        }
    }
    int nal_unit_type;
    bool firstSliceSegment;
    size_t pos;
};


template <>
struct Read<AccessUnit>
{
    template <class H> static void go(AccessUnit au, H &h)
    {
        StatePicturesBase *statePicturesBase = h;
        statePicturesBase->sliceHeaderValid = false;

        size_t posNextAu = bitLen(h[::Stream()]) / 8;
        while (h[Position<Bytes<>>()] < posNextAu)
        {
            try
            {
                const size_t NumBytesInNalUnit = numBytesInNalUnit(h[::Stream()]);
                h(byte_stream_nal_unit(boost::numeric_cast<int>(NumBytesInNalUnit)));
            }
            catch (boost::numeric::bad_numeric_cast &)
            {
                h(Violation("7", "NAL unit too large"));
                h[Stop()] = 1;
                break;
            }
            catch (Abort &)
            {
                h[Stop()] = 1;
                break;
            }
            const bool firstSliceSegmentInPicture = isSliceSegment(h[nal_unit_type()]) && h[first_slice_segment_in_pic_flag()];
            if (firstSliceSegmentInPicture)
            {
                size_t startPos = 0;

                typename Access<Stream, H>::SetType::Bookmark mark(h[::Stream()]);
                try
                {
                    bool seenEos = false;
                    // loop over NAL units
                    while (true)
                    {
                        NalUnitInfo info(h[::Stream()]);

                        if (!startPos && (info.firstSliceSegment || startsNewAccessUnit(info.nal_unit_type)))
                        {
                            startPos = info.pos;
                        }

                        if (info.firstSliceSegment)
                        {
                            posNextAu = startPos;

                            if (isIdr(info.nal_unit_type) || isBla(info.nal_unit_type) || seenEos)
                            {
                                static_cast<StatePicturesBase *>(h)->posEndCvs = posNextAu;
                            }

                            break;
                        }

                        {
                            if (info.nal_unit_type == EOS_NUT) seenEos = true;
                        }
                    }
                }
                catch (ExceptionOverrun &)
                {
                }
            }
        }

        if (statePicturesBase->sliceHeaderValid)
        {
            h(PictureDone());
        }
    };
};


template <class V, int cIdx>
struct AccessCbfRead
{
    typedef bool Type;
    static Type get(V v, ChromaCbf &s)
    {
        int dy = 0;
        if (s.ChromaArrayType == 2)
        {
            const int log2Size = s.tt0.log2TrafoSize - v.trafoDepth;
            dy = (v.y0 & (1 << log2Size >> 1)) ? 1 : 0;
        }

        const uint8_t mask = 1 << cIdx;

        return !!(s.values[dy][v.trafoDepth] & mask);
    }

    static void set(V v, int value, ChromaCbf &s)
    {
        int dy = 0;
        if (s.ChromaArrayType == 2)
        {
            const int log2Size = s.tt0.log2TrafoSize - v.trafoDepth;
            dy = (v.y0 & (1 << log2Size >> 1)) ? 1 : 0;
        }

        const uint8_t mask = 1 << cIdx;

        if (value)
        {
            s.values[dy][v.trafoDepth] |= mask;
        }
        else
        {
            s.values[dy][v.trafoDepth] &= ~mask;
        }
    }
};

template <class S>
struct Access<cbf_cb, S, typename std::enable_if<std::is_base_of<ChromaCbf, S>::value>::type> :
    AccessCbfRead<cbf_cb, 1>
{
};

template <class S>
struct Access<cbf_cr, S, typename std::enable_if<std::is_base_of<ChromaCbf, S>::value>::type> :
    AccessCbfRead<cbf_cr, 2>
{
};



struct RbspState :
    StateSei,
    BitReader,
    CabacState,
    Contexts,
    ChromaCbf,
    ValueHolder<cbf_luma>,
    ValueHolder<rqt_root_cbf>
{
};

template <class S>
struct Access<Stream, S, typename std::enable_if<std::is_base_of<RbspState, S>::value>::type>
{
    typedef BitReader SetType;
    typedef BitReader &Type;
    static Type get(Stream, S &s)
    {
        return static_cast<RbspState &>(s);
    }
};


template <>
struct Read<nal_unit>
{
    template <class H> static void go(const nal_unit &e, H &h2);
};

template <>
struct Read<video_parameter_set_rbsp>
{
    template <class H> static void go(const video_parameter_set_rbsp &fun, H &h0)
    {
        std::shared_ptr<Vps> vps(new Vps());
        auto hh = h0.extend(&*vps);
        auto hhh = hh.extend(&vps->ptl);
        auto h = hhh.extend(&vps->hrdArray);

        try
        {
            Syntax<video_parameter_set_rbsp>::go(fun, h);

            if (!h[vps_sub_layer_ordering_info_present_flag()])
                for (int i = 0; i < h[vps_max_sub_layers_minus1()]; ++i)
                {
                    h[vps_max_dec_pic_buffering_minus1(i)] = h[vps_max_dec_pic_buffering_minus1(h[vps_max_sub_layers_minus1()])];
                    h[vps_max_num_reorder_pics(i)] = h[vps_max_num_reorder_pics(h[vps_max_sub_layers_minus1()])];
                    h[vps_max_latency_increase_plus1(i)] = h[vps_max_latency_increase_plus1(h[vps_max_sub_layers_minus1()])];
                }

            h0[Table<Vps>()][h[vps_video_parameter_set_id()]] = vps;
        }
        catch (Abort &)
        {
            // Seek to end of RBSP so as not to trigger a further error
            seek(h0, bitLen(h0[::Stream()]));
        }
    }
};


template <>
struct Read<hrd_parameters>
{
    template <class H> static void go(const hrd_parameters &fun, H &h)
    {
        HrdArray &hrdArray = *static_cast<HrdArray *>(h);

        hrdArray.hrd.push_back(Hrd());

        auto h2 = h.extend(&hrdArray.hrd.back());

        for (int i = 0; i <= fun.maxNumSubLayersMinus1; i++)
        {
            h2[fixed_pic_rate_within_cvs_flag(i)] = 1;
        }

        Syntax<hrd_parameters>::go(fun, h2);
    }
};


template <>
struct Read<sub_layer_hrd_parameters>
{
    template <class H> static void go(const sub_layer_hrd_parameters &fun, H &h)
    {
        Hrd &hrd = *static_cast<Hrd *>(h);

        hrd.sublayers.push_back(Hrd::SubLayer());

        auto h2 = h.extend(&hrd.sublayers.back());

        Syntax<sub_layer_hrd_parameters>::go(fun, h2);
    }
};


template <>
struct Read<seq_parameter_set_rbsp>
{
    template <class H> static void go(const seq_parameter_set_rbsp &fun, H &h1)
    {
        std::shared_ptr<Sps> sps(new Sps());
        auto h2 = h1.extend(&*sps);
        auto h3 = h2.extend(&sps->scalingListState);
        auto h4 = h3.extend(&sps->ptl);
        auto h = h4.extend(&sps->hrdArray);

        try
        {
            Syntax<seq_parameter_set_rbsp>::go(fun, h);

            if (!h[sps_sub_layer_ordering_info_present_flag()])
                for (int i = 0; i < h[sps_max_sub_layers_minus1()]; ++i)
                {
                    h[sps_max_dec_pic_buffering_minus1(i)] = h[sps_max_dec_pic_buffering_minus1(h[sps_max_sub_layers_minus1()])];
                    h[sps_max_num_reorder_pics(i)] = h[sps_max_num_reorder_pics(h[sps_max_sub_layers_minus1()])];
                    h[sps_max_latency_increase_plus1(i)] = h[sps_max_latency_increase_plus1(h[sps_max_sub_layers_minus1()])];
                }

            h[Table<Sps>()][h[sps_seq_parameter_set_id()]] = sps;
        }
        catch (Abort &)
        {
            // Seek to end of RBSP so as not to trigger a further error
            seek(h, bitLen(h[::Stream()]));
        }
    }
};


template <>
struct Read<pic_parameter_set_rbsp>
{
    template <class H> static void go(const pic_parameter_set_rbsp &fun, H &h2)
    {
        std::shared_ptr<Pps> pps(new Pps());

        auto h1 = h2.extend(&*pps);
        auto h = h1.extend(&pps->scalingListState);

        h[uniform_spacing_flag()] = 1;
        h[loop_filter_across_tiles_enabled_flag()] = 1;

        try
        {
            Syntax<pic_parameter_set_rbsp>::go(fun, h);

            h[Table<Pps>()][h[pps_pic_parameter_set_id()]] = pps;
        }
        catch (Abort &)
        {
            // Seek to end of RBSP so as not to trigger a further error
            seek(h, bitLen(h[::Stream()]));
        }
    }
};


template <>
struct Read<end_of_seq_rbsp>
{
    template <class H> static void go(const end_of_seq_rbsp &fun, H &h)
    {
        Syntax<end_of_seq_rbsp>::go(fun, h);
        StatePicturesBase *statePicturesBase = h;
        statePicturesBase->justSeenEndOfSeq = true;
    }
};


template <>
struct Read<short_term_ref_pic_set>
{
    template <class H> static void go(const short_term_ref_pic_set &fun, H &h1)
    {
        Strps strps;
        auto h = h1.extend(&strps);

        h[delta_idx_minus1()] = 0;
        h[inter_ref_pic_set_prediction_flag()] = 0;

        Syntax<short_term_ref_pic_set>::go(fun, h);

        const int idx = fun.stRpsIdx;
        const int RIdx = idx - (h[delta_idx_minus1()] + 1);

        if (h[inter_ref_pic_set_prediction_flag()])
        {
            {
                int i = 0;
                int j;
                for (j = h[NumPositivePics(RIdx)] - 1; j >= 0; j--)
                {
                    int dPoc = h[DeltaPocS1(RIdx, j)] + h[DeltaRPS()];

                    if (dPoc < 0 && h[use_delta_flag(h[NumNegativePics(RIdx)] + j)])
                    {
                        h[DeltaPocS0(idx, i)] = dPoc;
                        h[UsedByCurrPicS0(idx, i++)] = !!h[used_by_curr_pic_flag(h[NumNegativePics(RIdx)] + j)];
                    }
                }
                if (h[DeltaRPS()] < 0 && h[use_delta_flag(h[NumDeltaPocs(RIdx)])])
                {
                    h[DeltaPocS0(idx, i)] = h[DeltaRPS()];
                    h[UsedByCurrPicS0(idx, i++)] = !!h[used_by_curr_pic_flag(h[NumDeltaPocs(RIdx)])];
                }
                for (j = 0; j < h[NumNegativePics(RIdx)]; j++)
                {
                    const int dPoc = h[DeltaPocS0(RIdx, j)] + h[DeltaRPS()];
                    if (dPoc < 0 && h[use_delta_flag(j)])
                    {
                        h[DeltaPocS0(idx, i)] = dPoc;
                        h[UsedByCurrPicS0(idx, i++)] = !!h[used_by_curr_pic_flag(j)];
                    }
                }
                h[NumNegativePics(idx)] = i;
            }
            {
                int i = 0;
                int j;
                for (j = h[NumNegativePics(RIdx)] - 1; j >= 0; j--)
                {
                    const int dPoc = h[DeltaPocS0(RIdx, j)] + h[DeltaRPS()];

                    if (dPoc > 0 && h[use_delta_flag(j)])
                    {
                        h[DeltaPocS1(idx, i)] = dPoc;
                        h[UsedByCurrPicS1(idx, i++)] = !!h[used_by_curr_pic_flag(j)];
                    }
                }
                if (h[DeltaRPS()] > 0 && h[use_delta_flag(h[NumDeltaPocs(RIdx)])])
                {
                    h[DeltaPocS1(idx, i)] = h[DeltaRPS()];
                    h[UsedByCurrPicS1(idx, i++)] = !!h[used_by_curr_pic_flag(h[NumDeltaPocs(RIdx)])];
                }
                for (j = 0; j < h[NumPositivePics(RIdx)]; j++)
                {
                    const int dPoc = h[DeltaPocS1(RIdx, j)] + h[DeltaRPS()];
                    if (dPoc > 0 && h[use_delta_flag(h[NumNegativePics(RIdx)] + j)])
                    {
                        h[DeltaPocS1(idx, i)] = dPoc;
                        h[UsedByCurrPicS1(idx, i++)] = !!h[used_by_curr_pic_flag(h[NumNegativePics(RIdx)] + j)];
                    }
                }

                h[NumPositivePics(idx)] = i;
            }
        }
        else
        {
            h[NumNegativePics(idx)] = h[num_negative_pics()];
            h[NumPositivePics(idx)] = h[num_positive_pics()];

            for (int i = 0; i < h[NumNegativePics(idx)]; ++i)
            {
                h[UsedByCurrPicS0(idx, i)] = !!h[used_by_curr_pic_s0_flag(i)];
                if (i == 0)
                {
                    h[DeltaPocS0(idx, i)] = -(h[delta_poc_s0_minus1(i)] + 1);
                }
                else
                {
                    h[DeltaPocS0(idx, i)] = h[DeltaPocS0(idx, i - 1)] - (h[delta_poc_s0_minus1(i)] + 1);
                }
            }
            for (int i = 0; i < h[NumPositivePics(idx)]; ++i)
            {
                h[UsedByCurrPicS1(idx, i)] = !!h[used_by_curr_pic_s1_flag(i)];
                if (i == 0)
                {
                    h[DeltaPocS1(idx, i)] = +(h[delta_poc_s1_minus1(i)] + 1);
                }
                else
                {
                    h[DeltaPocS1(idx, i)] = h[DeltaPocS1(idx, i - 1)] + (h[delta_poc_s1_minus1(i)] + 1);
                }
            }
        }
    }
};


template <>
struct Read<slice_segment_header>
{
    template <class H> static void go(const slice_segment_header &v, H &h)
    {
        // reset those values in slice segment header that are not dependent on slice header
        *static_cast<SliceSegmentHeaderIndependent *>(h) = SliceSegmentHeaderIndependent();

        Syntax<slice_segment_header>::go(v, h);

        StatePictures *statePictures = h;
        statePictures->sliceHeaderDone(h);

        h[CtbAddrInTs()] = h[CtbAddrRsToTs(h[slice_segment_address()])];
    }
};


template <>
struct Read<Element<num_ref_idx_active_override_flag, u>>
{
    template <class H> static void go(Element<num_ref_idx_active_override_flag, u> e, H &h)
    {
        ReadU<num_ref_idx_active_override_flag, u>::go(e, h);

        if (!h[e.v]) // not overridden
        {
            h[num_ref_idx_l0_active_minus1()] = h[num_ref_idx_l0_default_active_minus1()];
            h[num_ref_idx_l1_active_minus1()] = h[num_ref_idx_l1_default_active_minus1()];
        }
    }
};


template <>
struct Read<slice_segment_data>
{
    template <class H> static void go(const slice_segment_data &fun, H &h)
    {
        auto statePictures = &h[Concrete<StatePicturesBase>()];
        statePictures->sliceHeaderValid = true; // move this into sliceHeaderDone() ?

        if (h[first_slice_segment_in_pic_flag()])
        {
            StatePicturesBase *statePicturesBase = h;
            // review: desirable not to have to copy this information here (current decoded metadata)
            statePicturesBase->sliceNalUnitType = h[nal_unit_type()];
            statePicturesBase->sliceType = h[slice_type()];
            statePicturesBase->idVps = h[sps_video_parameter_set_id()];
            statePicturesBase->idSps = h[pps_seq_parameter_set_id()];
            statePicturesBase->idPps = h[slice_pic_parameter_set_id()];
            statePicturesBase->picWidth = h[pic_width_in_luma_samples()];
            statePicturesBase->picHeight = h[pic_height_in_luma_samples()];
            //statePicturesBase->sliceCount[0] = 0;
            //statePicturesBase->sliceCount[1] = 0;
            //statePicturesBase->sliceCount[2] = 0;

            static_cast<ScalingMatrices *>(h)->initialise(getScalingListState(h));
        }

        h[CtbAddrInTs()] = h[CtbAddrRsToTs(h[slice_segment_address()])];

        // handle the case where first slice does not have first_slice_in_pic_flag set
        // review
        if (h[ContainerOf<CtbAddrRsToTs>()].empty()) 
            throw Abort();

        StatePicturesBase *statePicturesBase = h;
        statePicturesBase->streamType = StatePicturesBase::streamTypeCabac();

        StateSubstream stateSubstream(h);
        auto h2 = h.extend(&stateSubstream);

        h2[CtbAddrInTs()] = h[CtbAddrInTs()];
        h2[CtbAddrInRs()] = h[CtbAddrTsToRs(h[CtbAddrInTs()])];
        h2[QpY()] = h[QpY()];

        if (h[first_slice_segment_in_pic_flag()])
        {
            StateSpatial *stateSpatial = h2;
            Neighbourhood *neighbourhood = h2;
            stateSpatial->resize(h2, *neighbourhood);
        }

        CabacState *cabacState = h2;
        cabacState->ivlCurrRange = 510;

        if (isRasl(h2[nal_unit_type()]) && statePicturesBase->NoRaslOutputFlag)
        {
            // associated IRAP picture has NoRaslOutputFlag equal to 1, the RASL picture is not output and may not be correctly decodable
            Greedy<slice_segment_data>::go(fun, h2);
        }
        else
            Syntax<slice_segment_data>::go(fun, h2);

        if (h2[dependent_slice_segments_enabled_flag()])
            h2(ContextsSave(tablesDs));

        h[CtbAddrInTs()] = h2[CtbAddrInTs()];
        h[QpY()] = h2[QpY()];

        statePicturesBase->streamType = StatePicturesBase::streamTypeRbsp();
    }
};


template <>
struct Read<Element<slice_type, ue>>
{
    template <class H> static void go(Element<slice_type, ue> v, H &h)
    {
        // at this point, we know that this is not a dependent slice segment header so we can reset the dependent values
        *static_cast<SliceSegmentHeaderDependent *>(h) = SliceSegmentHeaderDependent();

        h[pic_output_flag()] = 1;
        h[collocated_from_l0_flag()] = 1;
        h[slice_beta_offset_div2()] = h[pps_beta_offset_div2()];
        h[slice_tc_offset_div2()] = h[pps_tc_offset_div2()];
        h[slice_loop_filter_across_slices_enabled_flag()] = h[pps_loop_filter_across_slices_enabled_flag()];

        ReadUe<slice_type>::go(v, h);
    }
};


template <>
struct Read<Element<slice_pic_parameter_set_id, ue>>
{
    template <class H> static void go(Element<slice_pic_parameter_set_id, ue> e, H &h)
    {
        ReadUe<slice_pic_parameter_set_id>::go(e, h);

        auto ppsEntry = h[Table<Pps>()].find(h[slice_pic_parameter_set_id()]);
        if (ppsEntry == h[Table<Pps>()].end())
        {
            h(Violation("7.4.2.4.2", "invalid slice_pic_parameter_set_id{%1%}") % h[slice_pic_parameter_set_id()]); // CondCheck 7.4.2.4.2-A
            throw Abort();
        }
        h[Active<Pps>()] = ppsEntry->second;

        auto spsEntry = h[Table<Sps>()].find(h[pps_seq_parameter_set_id()]);
        if (spsEntry == h[Table<Sps>()].end())
        {
            h(Violation("7.4.2.4.2", "invalid pps_seq_parameter_set_id{%1%}") % h[pps_seq_parameter_set_id()]); // CondCheck 7.4.2.4.2-D
            throw Abort();
        }
        h[Active<Sps>()] = spsEntry->second;

        auto vpsEntry = h[Table<Vps>()].find(h[sps_video_parameter_set_id()]);
        if (vpsEntry == h[Table<Vps>()].end())
        {
            h(Violation("7.4.2.4.2", "invalid sps_video_parameter_set_id{%1%}") % h[sps_video_parameter_set_id()]); // CondCheck 7.4.2.4.2-G
            throw Abort();
        }
        h[Active<Vps>()] = vpsEntry->second;
    }
};

template <>
struct Read<Element<pcm_sample_luma, uv>>
{
    template <class H> static void go(Element<pcm_sample_luma, uv> fun, H &h)
    {
        const auto value = read_bits<typename ValueType<pcm_sample_luma>::Type>(h, h[PcmBitDepthY()]);
        h[fun.v] = value;
        BitReader &reader = *static_cast<BitReader *>(h);
        //if (reader.log) *reader.log << "pcm_sample_luma\tu(" << h[PcmBitDepthY()] << ")\t" << std::hex << value << std::dec << "\n";
    }
};

template <>
struct Read<Element<pcm_sample_chroma, uv>>
{
    template <class H> static void go(Element<pcm_sample_chroma, uv> fun, H &h)
    {
        const auto value = read_bits<typename ValueType<pcm_sample_chroma>::Type>(h, h[PcmBitDepthC()]);
        h[fun.v] = value;
        BitReader &reader = *static_cast<BitReader *>(h);
        //if (reader.log) *reader.log << "pcm_sample_chroma\tu(" << h[PcmBitDepthC()] << ")\t" << std::hex << value << std::dec << "\n";
    }
};


template <>
struct Read<Element<cabac_zero_word, f>>
{
    template <class H> static void go(Element<cabac_zero_word, f> e, H &h)
    {
        try
        {
            ReadU<cabac_zero_word, f>::go(e, h);
        }
        catch (ExceptionOverrun &)
        {
            // Catch this exception to prevent non-conforming cabac_zero_word causing problems.
        }
    }
};

#endif
