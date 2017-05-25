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
#include <map>

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


int havoc_test_run(
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

namespace havoc {

template <class Bound>
int runTest(
    Bound b[2], 
    double *first_result,
    havoc_instruction_set set,
    int iterations)
{
    bool const do_measure_speed = !!iterations;

    havoc_timestamp sum = 0;
    int warmup = 100;
    int count = 0;

    if (!do_measure_speed)
        b[1].invoke(1); // not measuring time, just checking for mismatch
    else 
        while (count < iterations)
        {
            const havoc_timestamp start = havoc_get_timestamp();
            b[1].invoke(4);
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
                    count += 4;
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
            printf("(x%.2f)", *first_result / average);
        else
            *first_result = average;
    }

    b[0].invoke(1);
    if (b[0].mismatch(b[1]))
    {
        printf("-MISMATCH");
        return 1;
    }

    return 0;
}

template <class Table>
struct Tester
{
    std::map<havoc_instruction_set, havoc_code> codes;
    std::map<havoc_instruction_set, Table> tables;

    Tester(havoc_instruction_set mask, int codeBytes = 1000000)
    {
#define X(I, S, D) \
        this->codes[HAVOC_ ## S] = havoc_new_code(HAVOC_ ## S, codeBytes); \
        this->tables[HAVOC_ ## S].populate(this->codes[HAVOC_ ## S]); \

        HAVOC_INSTRUCTION_SET_XMACRO
#undef X
    }

    ~Tester()
    {
#define X(I, S, D) \
        havoc_delete_code(this->codes[HAVOC_ ## S]); \

        HAVOC_INSTRUCTION_SET_XMACRO
#undef X
    }

    template <class Bound>
    int test(
        Bound b[2],
        havoc_instruction_set mask,
        int iterations)
    {
        int error_count = 0;

        if (b[0].get(this->tables[HAVOC_C_REF]))
        {
            b[0].print();

            b[0].invoke(1);

            double firstResult = 0.0;
            havoc_instruction_set sets[] = {
#define X(I, S, D) HAVOC_ ## S,
                HAVOC_INSTRUCTION_SET_XMACRO
#undef X
            };

            for (auto set : sets)
                if (set & mask)
                {
                    if (b[1].get(this->tables[set]) && (firstResult == 0 || b[1] != b[0]))
                        error_count += havoc::runTest(b, &firstResult, set, iterations);
                }

            printf("\n");
        }

        return error_count;
    }
};



}

#endif

#endif
