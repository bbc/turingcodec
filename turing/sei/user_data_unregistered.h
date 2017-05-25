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
#include "../StateParameterSets.h"
#include "../Read.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/string_generator.hpp>



struct uuid_iso_iec_11578 {};

struct user_data_payload_byte {};

template <> struct ValueType<uuid_iso_iec_11578>
{
    typedef boost::uuids::uuid Type;
};


struct UserDataUnregistered
    :
    ValueHolder<uuid_iso_iec_11578>,
    ValueHolder<user_data_payload_byte>
{
};


template <> struct Syntax<user_data_unregistered>
{
    template <class H> static void go(user_data_unregistered fun, H &h)
    {
        h(uuid_iso_iec_11578(), u(128));
        for (int i = 16; i < fun.payloadSize; i++)
        {
            h(user_data_payload_byte(), b(8));
        }
    }
};


template <class H> void Read<user_data_unregistered>::go(user_data_unregistered f, H &h)
{
    UserDataUnregistered userDataRegistered;
    auto h3 = h.extend(&userDataRegistered);
    Syntax<user_data_unregistered>::go(f, h3);
}


template <>
struct Read<Element<uuid_iso_iec_11578, u>>
{
    template <class H> static void go(Element<uuid_iso_iec_11578, u> e, H &h)
    {
        boost::uuids::uuid &uuid = h[e.v];
        assert(e.m.n == 8 * uuid.size());
        for (auto &byte : uuid)
            byte = read_bytes<int>(h, 1);

        StateParameterSets *stateParameterSets = h;
        stateParameterSets->unregisteredUserData[uuid].reset(new std::vector<uint8_t>());
    }
};


template <>
struct Read<Element<user_data_payload_byte, b>>
{
    template <class H> static void go(Element<user_data_payload_byte, b> e, H &h)
    {
        ReadU<user_data_payload_byte, b>::go(e, h);

        auto uuid = h[uuid_iso_iec_11578()];

        StateParameterSets *stateParameterSets = h;
        stateParameterSets->unregisteredUserData[uuid]->push_back(h[e.v]);
    }
};


struct StateWriteUserDataUnregistered
    :
    ValueHolder<uuid_iso_iec_11578>
{
    char const *p;
    char const *end;
};


template <> struct Write<user_data_unregistered>
{
    template <class H> static void go(user_data_unregistered fun, H &h)
    {
        StateWriteUserDataUnregistered *stateWriteUserDataUnregistered = h;
        fun.payloadSize = 16 + static_cast<int>(stateWriteUserDataUnregistered->end - stateWriteUserDataUnregistered->p);
        Syntax<user_data_unregistered>::go(fun, h);
    }
};


template <>
struct Write<Element<uuid_iso_iec_11578, u>>
{
    template <class H> static void go(Element<uuid_iso_iec_11578, u> e, H &h)
    {
        BitWriter *bitWriter = h;

        boost::uuids::uuid uuid = h[e.v];
        assert(e.m.n == 8 * uuid.size());

        for (auto byte : uuid)
            bitWriter->writeBits(8, byte);
    }
};


template <> struct Write<Element<user_data_payload_byte, b>>
{
    template <class H> static void go(Element<user_data_payload_byte, b> e, H &h)
    {
        StateWriteUserDataUnregistered *stateWriteUserDataUnregistered = h;
        BitWriter *bitWriter = h;
        bitWriter->writeBits(8, *stateWriteUserDataUnregistered->p++);
    }
};
