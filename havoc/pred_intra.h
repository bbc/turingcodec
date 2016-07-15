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

#ifndef INCLUDED_havoc_prediction_intra_h
#define INCLUDED_havoc_prediction_intra_h

// HEVC intra prediction


#include "havoc.h"


#ifdef __cplusplus
extern "C"
{
#endif

    // HEVC intra prediction

    typedef int havoc_pred_intra_packed;

    // Function prototype compatible with that of f265
    typedef void havoc_pred_intra(uint8_t *dst, const uint8_t *neighbours, int intraPredMode, havoc_pred_intra_packed packed);

    // Packing compatible with f265
    static __inline havoc_pred_intra_packed havoc_pred_intra_pack(int cIdx, int log2CbSize)
    {
        const int bit_depth = 8;
        const int edge_flag = (cIdx == 0 && log2CbSize < 5) ? 1 : 0;
        return (bit_depth << 8) | (log2CbSize << 1) | edge_flag;
    }

    typedef struct
    {
        havoc_pred_intra *p[4 /* log2CbSize - 2 */][35 /* intraPredModeY */];
        havoc_pred_intra *p_dc_filter[3 /* log2CbSize - 2 */];
    }
    havoc_table_pred_intra;

    static __inline havoc_pred_intra** havoc_get_pred_intra(havoc_table_pred_intra *table, int intraPredMode, havoc_pred_intra_packed packed)
    {
        const int edge_flag = packed & 0x1;
        const int log2CbSize = (packed >> 1) & 0x7f;
        if (intraPredMode == 1 && edge_flag)
        {
            return &table->p_dc_filter[log2CbSize - 2];
        }
        else
        {
            return &table->p[log2CbSize - 2][intraPredMode];
        }
    }

    void havoc_populate_pred_intra(havoc_table_pred_intra *table, havoc_code code);

    havoc_test_function havoc_test_pred_intra;

#ifdef __cplusplus
}
#endif

#endif
