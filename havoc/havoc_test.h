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


#ifndef INCLUDED_havoc_test_h
#define INCLUDED_havoc_test_h

#include "havoc.h"

#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif


    static const char *havoc_instruction_set_as_text(havoc_instruction_set set)
    {
#define X(value, name, description) if (set == (1 << value)) return #name;
        HAVOC_INSTRUCTION_SET_XMACRO
#undef X
        return 0;
    }

    typedef void havoc_bound_invoke(void *bound, int n);

    typedef int havoc_bound_mismatch(void *boundRef, void *boundTest);

    typedef int havoc_bound_get(void *bound, havoc_code code);


    int havoc_count_average_cycles(
            void *boundRef, void *boundTest,
            havoc_bound_invoke *f,
            havoc_bound_mismatch *m,
            double *first_result,
            havoc_instruction_set set,
            int iterations);


    int havoc_test(
            void *ref,
            void *test,
            havoc_bound_get *get,
            havoc_bound_invoke *invoke,
            havoc_bound_mismatch *mismatch,
            havoc_instruction_set mask,
            int iterations);


#ifdef __cplusplus
}
#endif


#endif
