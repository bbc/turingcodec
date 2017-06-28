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



/* Residual calculation */


#ifndef INCLUDED_res_h
#define INCLUDED_res_h

#include "havoc.h"


template <typename Sample>
using havoc_residual = void(int16_t *res, intptr_t stride_res, Sample const *srcA, intptr_t stride_srcA, Sample const *srcB, intptr_t stride_srcB, int w, int h);


template <typename Sample>
struct havoc_table_residual
{
    havoc_residual<Sample> *res[5];
};


template <typename Sample>
static havoc_residual<Sample>** havoc_get_residual(havoc_table_residual<Sample> *table, int log2TrafoSize)
{

    return &table->res[log2TrafoSize - 2];
}


template <typename Sample>
void havoc_populate_residual(havoc_table_residual<Sample> *table, havoc_code code);


void havoc_test_residual(int *error_count, havoc_instruction_set mask);


#endif
