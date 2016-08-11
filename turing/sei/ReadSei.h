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

#ifndef INCLUDED_ReadSei_H
#define INCLUDED_ReadSei_H

#include "buffering_period.h"
#include "layers_not_present.h"
#include "user_data_unregistered.h"
#include "time_code.h"
#include "overlay_info.h"
#include "../Read.hpp"
#include "../Violation.h"


// buffering_period

// bp_seq_parameter_set_id specifies the SPS to use during parsing buffering_period()
// Set that SPS active
template <> struct Read<Element<bp_seq_parameter_set_id, ue>>
{
    template <class H> static void go(Element<bp_seq_parameter_set_id, ue> const &f, H &h)
    {
        ReadUe<bp_seq_parameter_set_id>::go(f, h);

        auto sps = h[Table<Sps>()][h[bp_seq_parameter_set_id()]];

        if (!sps)
        {
            h(Violation("", "no SPS associated with bp_seq_parameter_set_id{%1%}") % h[bp_seq_parameter_set_id()]);
            throw Abort();
        }

        h[Active<Sps>()] = sps;

        Hrd **hrd = h;
        *hrd = getHrd(h);

        if (!*hrd)
        {
            h(Violation("", "no HRD parameters found"));
            throw Abort();
        }
    }
};


// layers_not_present

template <class H> void Read<layers_not_present>::go(layers_not_present  f, H &h)
{
    LayersNotPreset s;
    auto h3 = h.extend(&s);
    Syntax<layers_not_present>::go(f, h3);
}


template <> struct Read<Element<lnp_sei_active_vps_id, u>>
{
    template <class H> static void go(Element<lnp_sei_active_vps_id, u> f, H &h)
    {
        ReadU<lnp_sei_active_vps_id, u>::go(f, h);

        auto found = h[Table<Vps>()].find(h[f.v]);

        if (found == h[Table<Vps>()].end())
        {
            h(Violation("F.14.3.3", "lnp_sei_active_vps_id does not indicate a valid VPS ID"));
            throw Abort();
        }
        else
        {
            h[Active<Vps>()] = found->second;
        }
    }
};


// user_data_unregistered

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
        {
            byte = read_bytes<int>(h, 1);
        }

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


// time_code

template <> struct Read<Element<time_offset_value, i>>
{
    template <class H> static void go(Element<time_offset_value, i> f, H &h)
    {
        f.m.n = h[time_offset_length(f.v.i)];
        ReadI<time_offset_value, i>::go(f, h);
    }
};


// overlay_info

template <class H> void Read<overlay_info>::go(overlay_info f, H &h)
{
    OverlayInfo overlayInfo;
    auto hh = h.extend(&overlayInfo);
    Syntax<overlay_info>::go(f, hh);
}


template <>
struct Read<Element<overlay_element_label_min, u>>
{
    template <class H> static void go(Element<overlay_element_label_min, u> fun, H &h)
    {
        const int nBits = h[overlay_element_label_value_length_minus8()] + 8;
        h[fun.v] = read_bits<typename ValueType<overlay_element_label_min>::Type>(h, nBits);
    }
};

template <>
struct Read<Element<overlay_element_label_max, u>>
{
    template <class H> static void go(Element<overlay_element_label_max, u> fun, H &h)
    {
        const int nBits = h[overlay_element_label_value_length_minus8()] + 8;
        h[fun.v] = read_bits<typename ValueType<overlay_element_label_max>::Type>(h, nBits);
    }
};


template <class V>
struct Read<Element<V, st>>
{
    template <class H> static void go(Element<V, st> fun, H &h)
    {
        while (true)
        {
            char c = read_bytes<char>(h, 1);
            if (!c) break;
            h[fun.v] += c;
        }
    }
};


#endif
