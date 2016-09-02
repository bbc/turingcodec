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

#ifndef INCLUDED_Read_hh
#define INCLUDED_Read_hh

#pragma once

#include "Read.h"
#include "SyntaxSei.h"


template <typename T, class H>
struct NextBits<T, H>
{
    static T go(H &h, int n, bool aligned = false)
    {
        auto& stream = h[Stream()];
        typename Access<Stream, H>::SetType::Bookmark mark(stream);
        assert(bitPos<int>(stream) <= (int)bitLen(stream));
        T returnValue = 0;
        if (aligned)
        {
            assert(n % 8 == 0);
            assert(bitPos<int>(stream) % 8 == 0);
            n /= 8;
            while (n--)
            {
                if (bitPos<int>(stream) == bitLen(stream))
                {
                    return 0;
                }
                returnValue <<= 8;
                returnValue |= stream.byte();
            }
        }
        else
        {
            while (n--)
            {
                if (bitPos<size_t>(stream) >= bitLen(stream))
                {
                    return 0;
                }
                returnValue <<= 1;
                returnValue |= stream.bit();
            }
        }
        return returnValue;
    }
};


template <class H>
void Read<nal_unit>::go(const nal_unit &e, H &h2)
{
    if (e.NumBytesInNalUnit == 0)
        return;

    StateRbspData *stateRbspData = h2;

    auto &rbspData = stateRbspData->rbspData;

    rbspData.clear();
    rbspData.reserve(e.NumBytesInNalUnit);

    const std::streamsize endPositionNalUnit = bitPos<std::streamsize>(h2[::Stream()]) / 8 + e.NumBytesInNalUnit;

    if (canFastReadRbspData(h2[::Stream()]))
    {
        Syntax<nal_unit>::go(nal_unit(2), h2);
        rbspData.resize(e.NumBytesInNalUnit);
        fastReadRbspData(h2[::Stream()], rbspData, endPositionNalUnit, stateRbspData->eb3pPositions);
    }
    else
    {
        rbspData.clear();
        rbspData.reserve(e.NumBytesInNalUnit);
        Syntax<nal_unit>::go(e, h2);
    }

    RbspState rbspState;
    auto h =  h2.extend(&rbspState);

    h[::Stream()] = BitReader(rbspData.data(), rbspData.data() + rbspData.size());

    static_cast<StatePicturesBase *>(h)->streamType = StatePicturesBase::streamTypeRbsp();
    try
    {
        switch (h[nal_unit_type()])
        {
#define X(TYPE, NAME, DESCRIPTION, ELEMENT, VCLMSG) \
        case TYPE: h(ELEMENT()); break;

            NAL_UNIT_TYPES
#undef X
            default:
                assert(false);
                break;
        }
    }
    catch (Abort &)
    {
        h(Violation("", "NALU parsing aborted"));
    }
    catch (ExceptionOverrun &)
    {
        h(Violation("7", "parsing overrun (trying to read more RBSP data than exists)"));
    }

    if (h[::Stream()].state.p != h[::Stream()].end)
    {
        h(UnexpectedData());
    }

    static_cast<StatePicturesBase *>(h)->streamType = StatePicturesBase::streamTypeNal();
}


template <class H>
void Read<coding_tree_unit>::go(const coding_tree_unit &ctu, H &h)
{
    preCtu(h);
    Syntax<coding_tree_unit>::go(ctu, h);
    postCtu(h);
}

#endif
