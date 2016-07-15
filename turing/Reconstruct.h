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

// Routines to reconstruct output video regions according to coding settings.

#ifndef INCLUDED_Reconstruct_h
#define INCLUDED_Reconstruct_h

#pragma once

#include "Global.h"


// performs luma intra prediction with split_transform_flag = 0 and returns SATD
template <class H>
int32_t predictIntraLuma(transform_tree const &tt, H &h);


// performs luma intra prediction, codes residual and computes SSD
template <class H>
int32_t reconstructIntraLuma(IntraPartition const &intraPartition, H &h);


// performs chroma intra prediction, codes residual and computes SSD
template <class H>
int32_t reconstructIntraChroma(transform_tree const &tt, H &h);


// performs inter prediction, codes residual and computes SSD
template <class H>
int32_t reconstructInter(transform_tree const &tt, H &h);

#endif
