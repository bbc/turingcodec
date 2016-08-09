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


#ifndef INCLUDED_havoc_prediction_inter_h
#define INCLUDED_havoc_prediction_inter_h

// HEVC inter prediction

// Note that these functions may write data to the right of the destination block.


#include "havoc.h"


// HEVC uni prediction
template <typename Sample>
using HavocPredUni = void (Sample *dst, intptr_t stride_dst, Sample const *ref, intptr_t stride_ref, int nPbW, int nPbH, int xFrac, int yFrac, int bitDepth);
typedef HavocPredUni<uint8_t> havoc_pred_uni_8to8;
typedef HavocPredUni<uint16_t> havoc_pred_uni_16to16;

template <typename Sample>
struct HavocTablePredUni
{
    HavocPredUni<Sample>* p[3][2][17][2][2];
};


template <typename Sample>
static HavocPredUni<Sample>** havocGetPredUni(HavocTablePredUni<Sample> *table, int taps, int w, int h, int xFrac, int yFrac, int bitDepth)
{
    return &table->p[bitDepth - 8][taps / 4 - 1][(w + taps - 1) / taps][xFrac ? 1 : 0][yFrac ? 1 : 0];
}


template <typename Sample>
void havocPopulatePredUni(HavocTablePredUni<Sample> *table, havoc_code code);


havoc_test_function havoc_test_pred_uni;


// HEVC bi prediction

template <typename Sample>
using HavocPredBi = void(Sample *dst0, intptr_t stride_dst, const Sample *ref0, const Sample *ref1, intptr_t stride_ref, int nPbW, int nPbH, int xFrac0, int yFrac0, int xFrac1, int yFrac1, int bitDepth);

template <typename Sample>
struct HavocTablePredBi
{
    HavocPredBi<Sample> * p[3][2][9][2];
};

template <typename Sample>
static HavocPredBi<Sample>** havocGetPredBi(HavocTablePredBi<Sample> *table, int taps, int w, int h, int xFracA, int yFracA, int xFracB, int yFracB, int bitDepth)
{
    const int frac = xFracA || yFracA || xFracB || yFracB;
    return &table->p[bitDepth - 8][taps / 4 - 1][(w + 2 * taps - 1) / (2 * taps)][frac];
}

template <typename Sample>
void havocPopulatePredBi(HavocTablePredBi<Sample> *table, havoc_code code);

havoc_test_function havoc_test_pred_bi;

// plan eventually to use havoc namespace for entire library
namespace havoc {

    template <typename Sample>
    using SubtractBi = void(Sample *dst0, intptr_t stride_dst, const Sample *ref0, intptr_t stride_ref, const Sample *src, intptr_t stride_src, int nPbW, int nPbH, int bitDepth);

    template <typename Sample>
    struct TableSubtractBi
    {
        SubtractBi<Sample>* p;

        SubtractBi<Sample>*& get()
            {
            return this->p;
            }
    };

    template <typename Sample>
    void populateSubtractBi(TableSubtractBi<Sample> *table, havoc_code code, int bitDepth=0);

    template <typename Sample>
    void testSubtractBi(int *error_count, havoc_instruction_set mask);

}

#endif
