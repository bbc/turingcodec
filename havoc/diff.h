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


#ifndef INCLUDED_diff_h
#define INCLUDED_diff_h

#include "havoc.h"


#ifdef __cplusplus
extern "C"
{
#endif


    /* Linear SSD (Sum of Squared Differences) */
    typedef int havoc_ssd_linear(const uint8_t *src0, const uint8_t *src1, int size);

    havoc_ssd_linear* havoc_get_ssd_linear(int size, havoc_code code);

    havoc_test_function havoc_test_ssd_linear;


#ifdef __cplusplus
}
#endif


#endif
