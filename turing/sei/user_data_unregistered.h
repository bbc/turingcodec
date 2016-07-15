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
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/string_generator.hpp>

DEFINE_STRUCT_ARITY_1(user_data_payload_byte, i);
template<> struct ValueType<user_data_payload_byte> { typedef std::uint8_t Type; };
template <> struct IsValueArray<user_data_payload_byte> : std::true_type {};
template <> struct ValueHolder<user_data_payload_byte>
{
    std::vector<uint8_t> value;
    ValueType<picture_md5>::Type &get(user_data_payload_byte data)
    {
        if(data.i == value.size())
        {
            // Trying to access the end of the vector, add an element
            value.push_back(0);
        }
        assert(data.i < value.size());
        return value[data.i];
    }
};
DEFINE_STRUCT_ARITY_1(uuid_iso_iec_11578, i);
template <> struct ValueType<uuid_iso_iec_11578> { typedef std::uint8_t Type; };
template <> struct IsValueArray<uuid_iso_iec_11578> : std::true_type {};
template <> struct ValueHolder<uuid_iso_iec_11578>
{
    std::array<uint8_t, 16> value;
    ValueType<picture_md5>::Type &get(uuid_iso_iec_11578 data)
    {
        assert(data.i < 16);
        return value[data.i];
    }
};

struct UserDataUnregistered :
    ValueHolder<uuid_iso_iec_11578>,
    ValueHolder<user_data_payload_byte>
    {
    };

template <>
struct Syntax<user_data_unregistered>
{
    template <class H> static void go(const user_data_unregistered &fun, H &h)
    {

        // Review: find a better internal representation to have a 1:1 mapping with the HEVC spec.
        //h(uuid_iso_iec_11578(), u(128));
        for(int i = 0; i < 16; i++)
        {
            h(uuid_iso_iec_11578(i), u(8));
        }
        for (int i = 16; i < fun.payloadSize; i++)
        {
            h(user_data_payload_byte(i-16), b(8));
        }
    }
};



#ifdef EXPLICIT_INSTANTIATION
    EXPLICIT_INSTANTIATION(user_data_unregistered)
#endif
