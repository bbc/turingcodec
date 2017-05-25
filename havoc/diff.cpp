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

#include "diff.h"
#include "havoc_test.h"
#include "Jit.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>


static int havoc_ssd_linear_c_ref(const uint8_t *p1, const uint8_t *p2, int n)
{
    int sum = 0;
    for (int i = 0; i<n; ++i)
    {
        const int diff = p1[i] - p2[i];
        sum += diff * diff;
    }
    return sum;
}


havoc_ssd_linear * havoc_get_ssd_linear(int size, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);

    havoc_ssd_linear *f = 0;


    if (buffer.isa & (HAVOC_C_REF | HAVOC_C_OPT)) f = havoc_ssd_linear_c_ref;

    // review:
    //if (mask & HAVOC_AVX) f = havoc_ssd_linear_avx;

    return f;
}


#define BLOCK_SIZE 0x200


typedef struct
{
    HAVOC_ALIGN(32, uint8_t, data[2][BLOCK_SIZE]);
    havoc_ssd_linear *f;
    int ssd;
}
bound_ssd_linear;


int get_ssd_linear(void *p, havoc_code code)
{
    auto &buffer = *reinterpret_cast<Jit::Buffer *>(code.implementation);
    bound_ssd_linear *s = (bound_ssd_linear *)p;

    s->f = havoc_get_ssd_linear(BLOCK_SIZE, code);

    if (buffer.isa == HAVOC_C_REF) printf("\t%d:", BLOCK_SIZE);

    return !!s->f;
}


void invoke_ssd_linear(void *p, int n)
{
    bound_ssd_linear *s = (bound_ssd_linear *)p;
    while (n--)
    {
        s->ssd = s->f(s->data[0], s->data[1], BLOCK_SIZE);
    }
}


int mismatch_ssd_linear(void *boundRef, void *boundTest)
{
    bound_ssd_linear *ref = (bound_ssd_linear *)boundRef;
    bound_ssd_linear *test = (bound_ssd_linear *)boundTest;

    return  ref->ssd != test->ssd;
}


// review: not called
void havoc_test_ssd_linear(int *error_count, havoc_instruction_set mask)
{
    printf("\nhavoc_ssd_linear - Linear Sum of Square Differences\n");

    bound_ssd_linear b[2];

    for (int n = 0; n < 2; ++n)
    {
        for (int x = 0; x < BLOCK_SIZE; ++x)
        {
            b[0].data[n][x] = rand();
        }
    }

    b[1] = b[0];

    *error_count += havoc_test(&b[0], &b[1], get_ssd_linear, invoke_ssd_linear, mismatch_ssd_linear, mask, 10);
}
