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

#include "havoc_test.h"


int do_measure_speed = 1;


int havoc_test_run(
    void *boundRef, void *boundTest,
    havoc_bound_invoke *f,
    havoc_bound_mismatch *m,
    double *first_result,
    havoc_instruction_set set,
    int iterations)
{
    havoc_timestamp sum = 0;
    int warmup = 100;
    int count = 0;
    if (!do_measure_speed)
    {
        // not measuring time, just checking for mismatch
        f(boundTest, 1);
    }
    else while (count < iterations)
    {
        const havoc_timestamp start = havoc_get_timestamp();
        f(boundTest, 40);
        const havoc_timestamp duration = havoc_get_timestamp() - start;

        if (warmup == 0)
        {
            if (8 * duration * count < 7 * 4 * sum)
            {
                // duration lower than mean
                sum = 0;
                count = 0;
                warmup = 10;
            }
            else 	if (7 * duration * count <= 8 * 4 * sum)
            {
                // duration close to mean
                count += 40;
                sum += duration;
            }
            else
            {
                // duration higher than mean
                warmup = 10;
            }
        }
        else
        {
            --warmup;
        }
    }

    printf(" %s", havoc_instruction_set_as_text(set));

    if (do_measure_speed)
    {
        const int average = (int)((sum + count / 2) / count);

        printf(":%d", average);
        if (*first_result != 0.0)
        {
            printf("(x%.2f)", *first_result / average);
        }
        else
        {
            *first_result = average;
        }
    }

    if (boundRef)
    {
        f(boundRef, 1);
        if (m(boundRef, boundTest))
        {
            printf("-MISMATCH");
            return 1;
        }
    }

    return 0;
}


int havoc_test(
    void *ref,
    void *test,
    havoc_bound_get *get,
    havoc_bound_invoke *invoke,
    havoc_bound_mismatch *mismatch,
    havoc_instruction_set mask,
    int iterations)
{
    int error_count = 0;

    havoc_code code_ref;
    code_ref = havoc_new_code(HAVOC_C_REF, 1000000);

    if (get(ref, code_ref))
    {
        invoke(ref, 1);

        double first_result = 0.0;
        havoc_instruction_set set;
        for (set = HAVOC_C_OPT; set; set <<= 1)
        {
            havoc_code code_test = havoc_new_code(set & mask, 100000000);
            if (get(test, code_test))
            {
                error_count += havoc_test_run(ref, test, invoke, mismatch, &first_result, set, iterations);
            }
            havoc_delete_code(code_test);
        }
        printf("\n");
    }

    havoc_delete_code(code_ref);

    return error_count;
}
