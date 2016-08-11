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


struct preferred_transfer_characteristics { };

template <>
struct Syntax<alternative_transfer_characteristics>
{
    template <class H> static void go(alternative_transfer_characteristics fun, H &h)
    {
        h(preferred_transfer_characteristics(), u(8));
    }
};

struct AlternativeTransferCharacteristics :
    ValueHolder<preferred_transfer_characteristics>
    {
    };

template <class H> void Read<alternative_transfer_characteristics>::go(alternative_transfer_characteristics f, H &h)
{
    AlternativeTransferCharacteristics atc;
    auto h2 = h.extend(&atc);

    Syntax<alternative_transfer_characteristics>::go(f, h2);
}
