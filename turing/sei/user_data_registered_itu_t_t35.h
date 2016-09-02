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

struct itu_t_t35_country_code { };
struct itu_t_t35_country_code_extension_byte { };
struct itu_t_t35_payload_byte { };

struct UserDataRegistered :
    ValueHolder<itu_t_t35_country_code>,
    ValueHolder<itu_t_t35_country_code_extension_byte>,
    ValueHolder<itu_t_t35_payload_byte>
    {
    };

template <>
struct Syntax<user_data_registered_itu_t_t35>
{
    template <class H> static void go(user_data_registered_itu_t_t35 fun, H &h)
    {
        h(itu_t_t35_country_code(), b(8));

        int i;
        if (h[itu_t_t35_country_code()] != 0xFF)
        {
            i = 1;
        }
        else
        {
            h(itu_t_t35_country_code_extension_byte(), b(8));
            i = 2;
        }

        do
        {
            h(itu_t_t35_payload_byte(), b(8));
            i++;
        } while (i < fun.payloadSize);
    }
};

template <class H> void Read<user_data_registered_itu_t_t35>::go(user_data_registered_itu_t_t35  f, H &h)
{
    UserDataRegistered userDataRegistered;
    auto h3 = h.extend(&userDataRegistered);

    Syntax<user_data_registered_itu_t_t35>::go(f, h3);
}


