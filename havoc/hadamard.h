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


/* Functions for decoding HEVC residual: inverse transform and add the result to predictor */


#ifndef INCLUDED_hadamard_h
#define INCLUDED_hadamard_h

#include "havoc.h"


template <typename Sample>
using havoc_hadamard_satd = int(Sample const *srcA, intptr_t stride_srcA, Sample const *srcB, intptr_t stride_srcB);


template <typename Sample>
struct havoc_table_hadamard_satd
{
    havoc_hadamard_satd<Sample> *satd[3];
};


template <typename Sample>
havoc_hadamard_satd<Sample>** havoc_get_hadamard_satd(havoc_table_hadamard_satd<Sample> *table, int log2TrafoSize)
{
    return &table->satd[log2TrafoSize - 1];
}


template <typename Sample>
void havoc_populate_hadamard_satd(havoc_table_hadamard_satd<Sample> *table, havoc_code code);


void havoc_test_hadamard_satd(int *error_count, havoc_instruction_set mask);


#endif
