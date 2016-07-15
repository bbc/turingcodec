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

#include "ScalingMatrices.h"
#include "ScanOrder.h"
#include "HevcTypes.h"
#include "StateParameterSets.h"
#include "Handlers.h"
#include "Syntax.h"

#include <algorithm>
#include <cassert>


namespace {

    int coefNum(int sizeId)
    {
        return std::min(64, (1 << (4 + (sizeId << 1))));
    }

    int defaultScalingList(int sizeId, int matrixId, int i)
    {
        if (sizeId == 0) return 16;

        const bool intra = matrixId < (sizeId == 3 ? 1 : 3);

        if (intra)
        {
            const int defaults[64] = {
                    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 16, 17, 16, 17, 18,
                    17, 18, 18, 17, 18, 21, 19, 20, 21, 20, 19, 21, 24, 22, 22, 24,
                    24, 22, 22, 24, 25, 25, 27, 30, 27, 25, 25, 29, 31, 35, 35, 31,
                    29, 36, 41, 44, 41, 36, 47, 54, 54, 47, 65, 70, 65, 88, 88, 115 };
            return defaults[i];
        }
        else
        {
            const int defaults[64] = {
                    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18,
                    18, 18, 18, 18, 18, 20, 20, 20, 20, 20, 20, 20, 24, 24, 24, 24,
                    24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 28, 28, 28, 28, 28,
                    28, 33, 33, 33, 33, 33, 41, 41, 41, 41, 54, 54, 54, 71, 71, 91 };
            return defaults[i];
        }
    }

}


void ScalingMatrices::initialise(ScalingListState *scalingListState)
{
    const bool useDefaultScalingLists = !scalingListState;
    ScalingListState defaultScalingListState;

    auto &h = useDefaultScalingLists ? defaultScalingListState : *scalingListState;

    for (int sizeId = 0; sizeId <= 3; ++sizeId)
    {
        const int nMatrices = ((sizeId == 3) ? 2 : 6);
        for (int matrixId = 0; matrixId < nMatrices; ++matrixId)
        {
            if (useDefaultScalingLists)
            {
                h[scaling_list_pred_mode_flag(sizeId, matrixId)] = 0;
                h[scaling_list_pred_matrix_id_delta(sizeId, matrixId)] = 0;
            }

            if (h[scaling_list_pred_mode_flag(sizeId, matrixId)] == 0)
            {
                // the values of the scaling list are the same as the values of a reference scaling list

                if (h[scaling_list_pred_matrix_id_delta(sizeId, matrixId)] == 0)
                {
                    // the scaling list is inferred from the default scaling list
                    for (int i = 0; i < coefNum(sizeId); ++i)
                    {
                        h[ScalingList(sizeId, matrixId, i)] = defaultScalingList(sizeId, matrixId, i);
                    }

                    if (sizeId >= 2)
                    {
                        // not in standard
                        h[scaling_list_dc_coef_minus8(sizeId - 2, matrixId)] = h[ScalingList(sizeId, matrixId, 0)] - 8;
                    }
                }
                else
                {
                    // the scaling list is inferred from the reference scaling list as follows
                    const int refMatrixId = matrixId - h[scaling_list_pred_matrix_id_delta(sizeId, matrixId)];
                    for (int i = 0; i < coefNum(sizeId); ++i)
                    {
                        h[ScalingList(sizeId, matrixId, i)] = h[ScalingList(sizeId, refMatrixId, i)];
                    }

                    if (sizeId >= 2)
                    {
                        // not in standard
                        h[scaling_list_dc_coef_minus8(sizeId - 2, matrixId)] = h[scaling_list_dc_coef_minus8(sizeId - 2, refMatrixId)];
                    }
                }

            }
            for (int i = 0; i < coefNum(sizeId); ++i)
            {
                assert(h[ScalingList(sizeId, matrixId, i)] >= 0);
                assert(h[ScalingList(sizeId, matrixId, i)] <= 255);
            }
        }
    }

    // The elements of the quantization matrix of size 4x4, ScalingFactor[ 0 ][ matrixId ][ ][ ], are derived as follows:
    for (int matrixId = 0; matrixId <= 5; ++matrixId)
    {
        for (int i = 0; i <= 15; ++i)
        {
            int x = ScanOrder(2, 0, i, 0);
            int y = ScanOrder(2, 0, i, 1);
            (*this)(0, matrixId, x, y) = h[ScalingList(0, matrixId, i)];
        }
    }

    // The elements of the quantization matrix of size 8x8, ScalingFactor[ 1 ][ matrixId ][ ][ ], are derived as follows:
    for (int matrixId = 0; matrixId <= 5; ++matrixId)
    {
        for (int i = 0; i <= 63; ++i)
        {
            int x = ScanOrder(3, 0, i, 0);
            int y = ScanOrder(3, 0, i, 1);
            (*this)(1, matrixId, x, y) = h[ScalingList(1, matrixId, i)];
        }
    }

    // The elements of the quantization matrix of size 16x16, ScalingFactor[ 2 ][ matrixId ][ ][ ], are derived as follows:
    for (int matrixId = 0; matrixId <= 5; ++matrixId)
    {
        for (int i = 0; i <= 63; ++i)
        {
            int x = ScanOrder(3, 0, i, 0);
            int y = ScanOrder(3, 0, i, 1);
            for (int j = 0; j <= 1; ++j)
            {
                for (int k = 0; k <= 1; ++k)
                {
                    (*this)(2, matrixId, x * 2 + k, y * 2 + j) = h[ScalingList(2, matrixId, i)];
                }
            }
        }
        (*this)(2, matrixId, 0, 0) = h[scaling_list_dc_coef_minus8(0, matrixId)] + 8;
    }

    // The elements of the quantization matrix of size 32x32, ScalingFactor[ 3 ][ matrixId ][ ][ ], are derived as follows:
    for (int matrixId = 0; matrixId <= 1; ++matrixId)
    {
        for (int i = 0; i <= 63; ++i)
        {
            int x = ScanOrder(3, 0, i, 0);
            int y = ScanOrder(3, 0, i, 1);
            for (int j = 0; j <= 3; ++j)
            {
                for (int k = 0; k <= 3; ++k)
                {
                    (*this)(3, matrixId, x * 4 + k, y * 4 + j) = h[ScalingList(3, matrixId, i)];
                }
            }
        }
        (*this)(3, matrixId, 0, 0) = h[scaling_list_dc_coef_minus8(1, matrixId)] + 8;
    }
}

int ScalingMatrices::matrixId(int sizeId, int cIdx, int cuPredMode)
{

    if (sizeId < 3)
    {
        return (cuPredMode == MODE_INTRA ? 0 : 3) + cIdx;
    }
    else
    {
        assert(cIdx == 0);
        return cuPredMode == MODE_INTRA ? 0 : 1;
    }
}

ScalingMatrices::Type* ScalingMatrices::getMatrix(int sizeId, int matrixId)
{
    switch (sizeId)
    {
        default:
            assert(!"invalid sizeId");
        case 0: return this->matrix4x4[matrixId];
        case 1: return this->matrix8x8[matrixId];
        case 2: return this->matrix16x16[matrixId];
        case 3: return this->matrix32x32[matrixId];
    }
}

ScalingMatrices::Type &ScalingMatrices::operator()(int sizeId, int matrixId, int x, int y)
{
    return this->getMatrix(sizeId, matrixId)[x + (y << (2 + sizeId))];
}



