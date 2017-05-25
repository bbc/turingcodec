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

#include "pred_inter.h"
#include "pred_intra.h"
#include "transform.h"
#include "sad.h"
#include "ssd.h"
#include "diff.h"
#include "quantize.h"
#include "hadamard.h"
#include "havoc.h"
#include "Jit.h"
#include <stdint.h>
#include <type_traits>
#include <string>
#ifdef WIN32
#include <Windows.h>
#endif
#ifdef _MSC_VER
#include <intrin.h>
#endif


#ifdef __GNUC__

static void __cpuidex(int cpuInfo[4], int function_id, int subfunction_id)
{
    __asm__ __volatile__("cpuid" :
    "=a" ((cpuInfo)[0]),
        "=b" ((cpuInfo)[1]),
        "=c" ((cpuInfo)[2]),
        "=d" ((cpuInfo)[3])
        :
        "0" (function_id),
        "2" (subfunction_id));
}


static uint64_t _xgetbv(uint32_t index)
{
    uint32_t eax, edx;
    __asm__ __volatile__("xgetbv" :
    "=a" (eax),
        "=d" (edx)
        :
        "c" (index));

    return ((uint64_t)edx << 32) | eax;
}


#endif


static int bit_is_set(int value, int n)
{
    return value & (1 << n);
}


havoc_instruction_set havoc_instruction_set_support()
{
    std::underlying_type<havoc_instruction_set>::type mask = HAVOC_C_REF | HAVOC_C_OPT;

    enum { eax = 0, ebx = 1, ecx = 2, edx = 3 };

    int cpuInfo[4]; // eax ... edx

    __cpuidex(cpuInfo, 0, 0);

    const int max_standard_level = cpuInfo[0];

    if (max_standard_level == 0) return (havoc_instruction_set)mask;

    __cpuidex(cpuInfo, 1, 0);

    if (bit_is_set(cpuInfo[edx], 26)) mask |= HAVOC_SSE2;

    if (bit_is_set(cpuInfo[ecx], 1)) mask |= HAVOC_SSE3;
    if (bit_is_set(cpuInfo[ecx], 9)) mask |= HAVOC_SSSE3;
    if (bit_is_set(cpuInfo[ecx], 19)) mask |= HAVOC_SSE41;
    if (bit_is_set(cpuInfo[ecx], 20)) mask |= HAVOC_SSE42;
    if (bit_is_set(cpuInfo[ecx], 23)) mask |= HAVOC_POPCNT;

    if (bit_is_set(cpuInfo[ecx], 28) && bit_is_set(cpuInfo[ecx], 27))
    {
        uint64_t xcr0 = _xgetbv(0);

        if ((xcr0 & 0x2) && (xcr0 & 0x4))
        {
            mask |= HAVOC_AVX;
            if (max_standard_level >= 7)
            {
                __cpuidex(cpuInfo, 7, 0);

                if (bit_is_set(cpuInfo[ebx], 5))
                    mask |= HAVOC_AVX2;
            }
        }
    }

    __cpuidex(cpuInfo, 0x80000000, 0);
    uint32_t cap = static_cast<uint32_t>(cpuInfo[eax]);
    if (cap >= 0x80000001)
    {
        __cpuidex(cpuInfo, 0x80000001, 0);
        if (cpuInfo[ecx] & 0x00000020)
            mask |= HAVOC_LZCNT;
    }

    return (havoc_instruction_set)mask;
}


void havoc_print_instruction_set_support(FILE *f, havoc_instruction_set mask)
{
    f = stdout;
    fprintf(f, "havoc processor instruction set support:\n");
#define X(value, name, description) fprintf(f, "[%c] " #name " (" description ")\n", ((1 << value) & mask) ? 'x' : ' ');
    HAVOC_INSTRUCTION_SET_XMACRO
#undef X
        fprintf(f, "\n");
}


havoc_code havoc_new_code(havoc_instruction_set mask, int size)
{
    havoc_code code;
    code.implementation = new Jit::Buffer(size, mask);
    return code;
}


void havoc_delete_code(havoc_code code)
{
    delete (Jit::Buffer *)code.implementation;
}


extern "C" int do_measure_speed;


int havoc_main(int argc, const char *argv[])
{
    std::cout << "HAVOC [Handcoded Assembly for VideO Codecs] self test";
    std::cout << "Performs the following checks:\n";
    std::cout << "\t1) Retrieve information on the machine's instruction set\n";
    std::cout << "\t2) Run some basic computations of SSD, transform, etc. and check the results for errors\n";
    std::cout << "";

    if (argc == 2)
    {
        if (argv[1] != std::string("--no-measure-speed"))
        {
            std::cout << "Usage: " << argv[0] << "[--no-measure-speed]\n";
            return -1;
        }
        do_measure_speed = 0;
    }

    havoc_instruction_set mask = havoc_instruction_set_support();

#ifdef WIN32
    if (!SetProcessAffinityMask(GetCurrentProcess(), 1))
        std::cout << "** SetProcessAffinityMask() failed **\n\n";
#endif

    havoc_print_instruction_set_support(stdout, mask);

    int error_count = 0;

    havoc::testSubtractBi<uint8_t>(&error_count, mask);
    havoc::testSubtractBi<uint16_t>(&error_count, mask);
    havoc_test_sad_multiref(&error_count, mask);
    havoc_test_sad(&error_count, mask);
    havoc_test_ssd(&error_count, mask);
    havoc::intra::test<uint8_t>(&error_count, mask);
    havoc::intra::test<uint16_t>(&error_count, mask);
    havoc_test_hadamard_satd(&error_count, mask);
    havoc_test_quantize_inverse(&error_count, mask);
    havoc_test_quantize(&error_count, mask);
    havoc_test_quantize_reconstruct(&error_count, mask);
    havoc_test_pred_uni(&error_count, mask);
    havoc_test_pred_bi(&error_count, mask);
    havoc::test_inverse_transform_add<uint8_t>(&error_count, mask);
    havoc::test_inverse_transform_add<uint16_t>(&error_count, mask);
    havoc::test_transform(&error_count, mask);

    printf("\n");
    printf("havoc self test: %d errors\n", error_count);

    return error_count;
}
