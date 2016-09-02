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


struct snapshot_id { };

struct PictureSnapshot :
    ValueHolder<snapshot_id>
    {
    };


template <>
struct Syntax<picture_snapshot>
{
    template <class H> static void go(picture_snapshot fun, H &h)
    {
        h(snapshot_id(), ue(v));
    }
};


template <class H> void Read<picture_snapshot>::go(picture_snapshot f, H &h)
{
    PictureSnapshot pictureSnapshot;
    auto h3 = h.extend(&pictureSnapshot);

    Syntax<picture_snapshot>::go(f, h3);
}
