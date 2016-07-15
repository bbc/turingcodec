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

// review: is this header still necessary?

#ifndef INCLUDED_DecodedPicture_h
#define INCLUDED_DecodedPicture_h

#pragma once

#include "StateCollocatedMotion.h"
#include "Violation.h"
#include <memory>


template <class Pic>
static bool LongTermRefPic(RplEntry<Pic> picX)
{
    return picX.reference == LONG_TERM;
}



template <class P>
using DecodedPictureBuffer = std::vector<std::shared_ptr<P>>;


template <class P>
static int positionOfPictureInDpb(const DecodedPictureBuffer<P> &dpb, const P *dp)
{
    assert(dp == dpb[dp->n].get());
    return dp->n;
}


template <class P>
static const P *getDpbPictureByIndex(const DecodedPictureBuffer<P> &dpb, int dpbIndex)
{
    return dpb[dpbIndex].get();
}

#endif
