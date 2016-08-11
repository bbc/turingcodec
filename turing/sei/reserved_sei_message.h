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

#pragma once

#include "../SyntaxSei.h"


struct reserved_sei_message_payload_byte { };


struct ReservedSeiMessage :
    ValueHolder<reserved_sei_message_payload_byte>
    {
    };


template <>
struct Syntax<reserved_sei_message>
{
    template <class H> static void go(reserved_sei_message f, H &h)
    {
        for (int i = 0; i < f.payloadSize; i++)
            h(reserved_sei_message_payload_byte(), b(8));
    }
};


template <class H> void Read<reserved_sei_message>::go(reserved_sei_message f, H &h)
{
    ReservedSeiMessage reservedSeiMessage;
    auto hh = h.extend(&reservedSeiMessage);
    Syntax<reserved_sei_message>::go(f, hh);
}
