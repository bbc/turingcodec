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

#ifndef INCLUDED_Sao_h
#define INCLUDED_Sao_h

#pragma once

#include <stdint.h>
#include <stddef.h>


template <typename Sample>
void sao_filter_band(Sample *__restrict dst, intptr_t dst_stride, const Sample *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict offset_table, int bitDepth);

template <typename Sample>
void sao_filter_edge(Sample *__restrict dst, intptr_t dst_stride, const Sample *__restrict src, intptr_t src_stride, int w, int h, const int16_t *__restrict SaoOffsetVal, int eoClass, int bitDepth);

#endif
