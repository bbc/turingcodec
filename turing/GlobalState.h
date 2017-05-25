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

// State structures for top level state common to parsing, writing, decoding and encoding.
// Review: rename / reorganise?

#ifndef INCLUDED_GlobalState_h
#define INCLUDED_GlobalState_h

#pragma once

#include "IntraReferenceSamples.h"
#include "StateParameterSets.h"
#include "Cabac.h"
#include "StatePictures.h"
#include "HevcTypes.h"
#include "Snake.h"
#include "StateSpatial.h"
#include "ScanOrder.h"
#include "ScalingMatrices.h"
#include "StateParameterSets.h"
#include "CandModeList.h"
#include "Profiler.h"
#include "StateFunctionTables.h"
#include "QpState.h"
#include <cstdint>
#include <type_traits>
#include <vector>
#include <map>
#include <memory>
#include <algorithm>



// Decoder or encoder state that persists for the whole video sequence
struct StateSequence :
    StateParameterSets,
    ValueHolder<Active<Vps>>,
    ValueHolder<Active<Sps>>,
    StatePictures,
    Profiler::Timers
{
};


namespace Mvp {

struct Predictors
{
    PuData merge[5];
    MotionVector mvp[15 /* refIdx */][2 /* refList */][2 /* mvp_flag */];
};

}



struct StateSubstream :
    AccessOperators<StateSubstream>,
    ValueCache<CtbLog2SizeY>,
    ValueCache<MinCbLog2SizeY>,
    ValueCache<PicWidthInCtbsY>,
    ValueCache<PicHeightInCtbsY>,
    ValueHolder<end_of_slice_segment_flag>,
    ValueHolder<end_of_subset_one_bit>,
    ValueHolder<CtbAddrInRs>,
    ValueHolder<CtbAddrInTs>,
    ValueHolder<xCtb>,
    ValueHolder<yCtb>,
    AvailabilityCtu,
    QpState,
    Snake<BlockData>::Array<32, 32, 0, 0>,
    CandModeList,
    Mvp::Predictors,
    coding_tree_unit,
    ValueHolder<pcm_sample_luma>,
    ValueHolder<pcm_sample_chroma>,
    ValueHolder<pcm_alignment_one_bit>,
    ValueHolder<pcm_alignment_zero_bit>,
    sao,
    ValueHolder<sao_merge_left_flag>,
    ValueHolder<sao_merge_up_flag>,
    ValueHolder<sao_type_idx_luma>,
    ValueHolder<sao_type_idx_chroma>,
    ValueHolder<sao_eo_class_luma>,
    ValueHolder<sao_eo_class_chroma>,
    coding_quadtree,
    ValueHolder<split_cu_flag>,
    ValueHolder<cu_transquant_bypass_flag>,
    ValueHolder<MaxTrafoDepth>,
    ValueHolder<pred_mode_flag>,
    ValueHolder<part_mode>,
    ValueHolder<split_transform_flag>,
    ValueHolder<IsCuQpDeltaCoded>,
    ValueHolder<CuQpDeltaVal>,
    ValueHolder<IsCuChromaQpOffsetCoded>,
    ValueHolder<cu_qp_delta_abs>,
    ValueHolder<cu_qp_delta_sign_flag>,
    ValueHolder<cu_chroma_qp_offset_flag>,
    ValueHolder<cu_chroma_qp_offset_idx>,
    prediction_unit,
    ValueHolder<merge_flag>,
    ValueHolder<merge_idx>,
    ValueHolder<inter_pred_idc>,
    transform_tree,
    transform_unit,
    ValueHolder<transform_skip_flag>,
    residual_coding,
    ValueHolder<explicit_rdpcm_flag>,
    ValueHolder<explicit_rdpcm_dir_flag>,
    ValueHolder<last_sig_coeff_x_prefix>,
    ValueHolder<last_sig_coeff_y_prefix>,
    ValueHolder<last_sig_coeff_x_suffix>,
    ValueHolder<last_sig_coeff_y_suffix>,
    ValueHolder<LastSignificantCoeffX>,
    ValueHolder<LastSignificantCoeffY>,
    ValueHolder<coded_sub_block_flag>,
    ValueHolder<numGreater1Flag>,
    ValueHolder<sig_coeff_flag>,
    ValueHolder<coeff_abs_level_greater1_flag>,
    ValueHolder<coeff_abs_level_greater2_flag>,
    ValueHolder<coeff_sign_flag>,
    ValueHolder<coeff_abs_level_remaining>,
    ValueHolder<TransCoeffLevel>,
    cross_comp_pred,
    ValueHolder<log2_res_scale_abs_plus1>,
    ValueHolder<res_scale_sign_flag>,
    mvd_coding,
    ValueHolder<ref_idx_l0>,
    ValueHolder<ref_idx_l1>,
    ValueHolder<abs_mvd_greater0_flag>,
    ValueHolder<abs_mvd_greater1_flag>,
    ValueHolder<abs_mvd_minus2>,
    ValueHolder<mvd_sign_flag>,
    ValueHolder<mvp_l0_flag>,
    ValueHolder<mvp_l1_flag>,
    ValueHolder<Mvd>,
    StateGrey
{
    using AccessOperators<StateSubstream>::operator[];

    template <class H>
    StateSubstream(H &h)
        :
        stateSpatial(h),
        QpState(h),
        ValueCache<CtbLog2SizeY>(h),
        ValueCache<MinCbLog2SizeY>(h),
        ValueCache<PicWidthInCtbsY>(h),
        ValueCache<PicHeightInCtbsY>(h),
        chromaArrayType(h[::ChromaArrayType()]),
        StateGrey(h)
    {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                prev_intra_luma_pred_flag[i][j] = 0;
    }

    Profiler::Timers timers;

    int rem_intra_luma_pred_mode[2][2];
    int mpm_idx[2][2];
    int prev_intra_luma_pred_flag[2][2];
    int intra_chroma_pred_mode[2][2];
    int chromaArrayType;
    int ctxSet; // used in Read.h Call<Element<coeff_abs_level_greater1_flag, ae>,  H>
    int greater1Ctx; // used in Read.h Call<Element<coeff_abs_level_greater1_flag, ae>,  H>
    int previousGreater1Flag; // used in Read.h Call<Element<coeff_abs_level_greater1_flag, ae>,  H>
    int lastGreater1Flag; // used in Read.h Call<residual_coding, H>
    int lastGreater1Ctx; // used in Read.h Call<residual_coding, H>
    int cLastAbsLevel; // used in Read.h Call<residual_coding, H>
    int cLastRiceParam; // used in Read.h Call<residual_coding, H>
    bool firstCoeffAbsLevelRemainingInSubblock;
    int i; // used in Syntax<residual_coding, H>
    StateSpatial *stateSpatial;
    int partIdx;
};


template <class V, class Derived>
struct Access<V, Derived, typename std::enable_if<std::is_base_of<Neighbourhood, Derived>::value && Accessible<V, BlockData>::value>::type>
{
    typedef const typename Access<V, BlockData>::Type Type;
    static Type get(V v, Neighbourhood &s)
    {
        Snake<BlockData>::Cursor &cursor = s;
        BlockData &blockData = const_cast<BlockData&>(cursor.offset<Anywhere>(v.x0 - cursor.x0, v.y0 - cursor.y0, s.MinCbLog2SizeYMinus1));
        return Access<V, BlockData>::get(v, blockData);
    }
    static void set(V v, Type value, Neighbourhood &s)
    {
        Snake<BlockData>::Cursor &cursor = s;
        BlockData &blockData = const_cast<BlockData&>(cursor.offset<Anywhere>(v.x0 - cursor.x0, v.y0 - cursor.y0, s.MinCbLog2SizeYMinus1));
        return Access<V, BlockData>::set(v, value, blockData);
    }
};


template <class Direction, class V, class Derived>
struct Access<Neighbouring<V, Direction>, Derived, typename std::enable_if<std::is_base_of<Neighbourhood, Derived>::value && Accessible<V, BlockData>::value>::type>
{
    typedef const typename Access<V, BlockData>::Type Type;
    static Type get(Neighbouring<V, Direction> v, Neighbourhood &s)
    {
        Snake<BlockData>::Cursor &cursor = s;
        BlockData &blockData = const_cast<BlockData&>(cursor.offset<Direction>(v.dx, v.dy, s.MinCbLog2SizeYMinus1));
        return Access<V, BlockData> ::get(v.v, blockData);
    }
};


template <typename Sample>
struct StateEncodeSubstream;


template <class S> struct Access<rem_intra_luma_pred_mode, S, typename std::enable_if<std::is_base_of<StateSubstream, S>::value>::type>
{
    typedef int& Type;
    static Type get(rem_intra_luma_pred_mode v, StateSubstream &s)
    {
        const coding_quadtree &cqt = s;
        return s.rem_intra_luma_pred_mode[v.x0 == cqt.x0 ? 0 : 1][v.y0 == cqt.y0 ? 0 : 1];
    }
    static void set(rem_intra_luma_pred_mode v, int i, StateSubstream &s)
    {
        get(v, s) = i;
    }
};


template <class S> struct Access<mpm_idx, S, typename std::enable_if<std::is_base_of<StateSubstream, S>::value>::type>
{
    typedef int &Type;
    static Type get(mpm_idx v, StateSubstream &s)
    {
        const coding_quadtree &cqt = s;
        return s.mpm_idx[v.x0 == cqt.x0 ? 0 : 1][v.y0 == cqt.y0 ? 0 : 1];
    }
    static void set(mpm_idx v, int i, StateSubstream &s)
    {
        get(v, s) = i;
    }
};


template <class S> struct Access<prev_intra_luma_pred_flag, S, typename std::enable_if<std::is_base_of<StateSubstream, S>::value>::type>
{
    typedef int &Type;
    static Type get(prev_intra_luma_pred_flag v, StateSubstream &s)
    {
        const coding_quadtree &cqt = s;
        return s.prev_intra_luma_pred_flag[v.x0 == cqt.x0 ? 0 : 1][v.y0 == cqt.y0 ? 0 : 1];
    }
    static void set(prev_intra_luma_pred_flag v, int i, StateSubstream &s)
    {
        get(v, s) = i;
    }
};


template <class S> struct Access<intra_chroma_pred_mode, S, typename std::enable_if<std::is_base_of<StateSubstream, S>::value>::type>
{
    typedef int &Type;
    static Type get(intra_chroma_pred_mode v, StateSubstream &s)
    {
        if (s.chromaArrayType == 3)
        {
            if (s[part_mode()] == 1 /* PART_NxN */)
            {
                coding_quadtree const &cqt = s;
                if (v.x0 != cqt.x0) assert((v.x0 - cqt.x0) >= (1 << cqt.log2CbSize >> 1));
                if (v.y0 != cqt.y0) assert((v.y0 - cqt.y0) >= (1 << cqt.log2CbSize >> 1));
                return s.intra_chroma_pred_mode[v.x0 == cqt.x0 ? 0 : 1][v.y0 == cqt.y0 ? 0 : 1];
            }
        }

        return s.intra_chroma_pred_mode[0][0];
    }
};

#endif
