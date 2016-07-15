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

// NAL syntax functions.
// In general these are more-or-less copy/pasted from the HEVC standard text and should not be modified.

#ifndef INCLUDED_SyntaxNal_hpp
#define INCLUDED_SyntaxNal_hpp

#pragma once

#include "Syntax.h"
#include "HevcTypes.h"
#include "Global.h"


// Annex B - Byte stream format

template <class H>
void Syntax<byte_stream_nal_unit>::go(const byte_stream_nal_unit &fun, H &h)
{
    while (next_bits<int>(h, 24, true) != 0x000001 && next_bits<int>(h, 32, true) != 0x00000001)
        h(leading_zero_8bits() /* equal to 0x00 */, f(8));

    if (next_bits<int>(h, 24, true) != 0x000001)
        h(zero_byte() /* equal to 0x00 */, f(8));

    h(start_code_prefix_one_3bytes() /* equal to 0x000001 */, f(24));

    h(nal_unit(fun.NumBytesInNalUnit));

    while (h[more_data_in_byte_stream()] && next_bits<int>(h, 24, true) != 0x000001 && next_bits<int>(h, 32, true) != 0x00000001)
        h(trailing_zero_8bits() /* equal to 0x00 */, f(8));
}


template <class H>
void  Syntax<nal_unit>::go(const nal_unit &fun, H &h)
{
    h(nal_unit_header());
    int NumBytesInRbsp = 0;
    for (int i = 2; i < fun.NumBytesInNalUnit; i++)
        if (i + 2 < fun.NumBytesInNalUnit && next_bits<int>(h, 24, true) == 0x000003)
        {
            h(rbsp_byte(NumBytesInRbsp++), b(8));
            h(rbsp_byte(NumBytesInRbsp++), b(8));
            i += 2;
            h(emulation_prevention_three_byte(), /* equal to 0x03 */ f(8));
        }
        else
            h(rbsp_byte(NumBytesInRbsp++), b(8));
}


template <class H>
void Syntax<nal_unit_header>::go(const nal_unit_header &, H &h)
{
    h(forbidden_zero_bit(), f(1));
    h(nal_unit_type(), u(6));
    h(nuh_layer_id(), u(6));
    h(nuh_temporal_id_plus1(), u(3));
}

#endif
