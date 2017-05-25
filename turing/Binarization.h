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

// CABAC binarization for HEVC bitstream writing.

#ifndef INCLUDED_Binarization_h
#define INCLUDED_Binarization_h

#pragma once

#include "Global.h"
#include "SyntaxElements.h"
#include "CandModeList.h"
#include "StateEncode.h"
#include <type_traits>
#include <cassert>


template <class V>
DEFINE_STRUCT_ARITY_2(EncodeDecision, binVal, ctxInc);

template <class V>
DEFINE_STRUCT_ARITY_1(EncodeTerminate, binVal);

template <class V>
DEFINE_STRUCT_ARITY_1(EncodeBypass, binVal);


struct FinishCabac;


template <class F> struct Write;


template <class V>
struct Write<Element<V, ae>>
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        ASSERT(!"not yet implemented");
    }
};


// slice_segment_data( )

template <>
struct Write<Element<end_of_slice_segment_flag, ae>>
{
    template <class H> static void go(Element<end_of_slice_segment_flag, ae> fun, H &h)
    {
        const int ctbAddrEnd = h[PicSizeInCtbsY()];

        int binVal = (h[CtbAddrInTs()] == ctbAddrEnd - 1) ? 1 : 0;

        Access<end_of_slice_segment_flag, H>::set(fun.v, binVal, h);

        h(EncodeTerminate<end_of_slice_segment_flag>(binVal));

        if (binVal)
        {
            h(FinishCabac(0));
        }
    }
};

template <>
struct Write<Element<end_of_subset_one_bit, ae>>
{
    template <class H> static void go(Element<end_of_subset_one_bit, ae> fun, H &h)
    {
        h(EncodeTerminate<end_of_subset_one_bit>(1));
        h(FinishCabac(0));
    }
};


// sao ()



// coding_quadtree()

template <>
struct Write<Element<split_cu_flag, ae>>
{
    template <class H> static void go(Element<split_cu_flag, ae> fun, H &h)
    {
        coding_quadtree cqt = *static_cast<coding_quadtree *>(h);

        const int binVal = h[fun.v];

        split_cu_flag e = fun.v;

        Neighbourhood *neighbourhood = h;

        const int ctxInc =
                (neighbourhood->snake.get<Left>(e.x0 - 1, e.y0, neighbourhood->MinCbLog2SizeYMinus1).CtDepth > cqt.cqtDepth ? 1 : 0) +
                (neighbourhood->snake.get<Up>(e.x0, e.y0 - 1, neighbourhood->MinCbLog2SizeYMinus1).CtDepth > cqt.cqtDepth ? 1 : 0);

        //std::cout << "\tsplit_cu_flag\t" << binVal << "\t" << ctxInc << "\n";

        h(EncodeDecision<split_cu_flag>(binVal, ctxInc));
    }
};

// sao( )
template <>
struct Write<Element<sao_type_idx_luma, ae>>
{
    template <class H> static void go(Element<sao_type_idx_luma, ae> fun, H &h)
    {
        //coding_quadtree cqt = *static_cast<coding_quadtree *>(h);
        const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
        const int ry = h[yCtb()] >> h[CtbLog2SizeY()];

        const int synVal = h[SaoTypeIdx(0, rx, ry)];
        const int binVal0 = (synVal == 0 ? 0 : 1);
        h(EncodeDecision<sao_type_idx_luma>(binVal0, 0));
        if (binVal0 != 0)
        {
            const int binVal1 = synVal - 1;
            h(EncodeBypass<sao_type_idx_luma>(binVal1));
        }
    }
};

template <>
struct Write<Element<sao_type_idx_chroma, ae>>
{
    template <class H> static void go(Element<sao_type_idx_chroma, ae> fun, H &h)
    {
        //coding_quadtree cqt = *static_cast<coding_quadtree *>(h);
        const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
        const int ry = h[yCtb()] >> h[CtbLog2SizeY()];

        const int synVal = h[SaoTypeIdx(1, rx, ry)];
        const int binVal0 = (synVal == 0 ? 0 : 1);
        h(EncodeDecision<sao_type_idx_chroma>(binVal0, 0));
        if (binVal0 != 0)
        {
            const int binVal1 = synVal - 1;
            h(EncodeBypass<sao_type_idx_chroma>(binVal1));
        }
    }
};

template <>
struct Write<Element<sao_offset_abs, ae>>
{
    template <class H> static void go(Element<sao_offset_abs, ae> fun, H &h)
    {
        const int bitDepth = fun.v.cIdx ? h[BitDepthC()] : h[BitDepthY()];
        const int cMax = (1 << (std::min(bitDepth, 10) - 5)) - 1;

        int synVal = h[fun.v];
        int binVal = synVal;
        for (int loop = 0; loop < cMax; loop++)
        {
            int currBin = (binVal == 0 ? 0 : 1);
            h(EncodeBypass<sao_offset_abs>(currBin));
            if (!binVal) break;
            binVal -= 1;
        }
    }
};

template <>
struct Write<Element<sao_eo_class_luma, ae>> // FL, cMax=3 | Bypass
{
    template <class H> static void go(Element<sao_eo_class_luma, ae> fun, H &h)
    {
        const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
        const int ry = h[yCtb()] >> h[CtbLog2SizeY()];
        const int synVal = h[SaoEoClass(0, rx, ry)];// h[fun.v]; //GetValueForWrite<rem_intra_luma_pred_mode, H>::get(fun.v, h);
        int fl = 2;
        for (int binIdx = 0; binIdx < fl; binIdx++)
        {
            int binValue = 1 << (fl - 1 - binIdx);
            binValue &= synVal;
            binValue >>= (fl - 1 - binIdx);
            h(EncodeBypass<sao_eo_class_luma>(binValue));
        }

        //int binVal0 = synVal <= 2 ? 1 : 0;
        //h(EncodeBypass<sao_eo_class_luma>(binVal0));
        //int binVal1 = (synVal == 3 || synVal == 1)? 1 : 0;
        //h(EncodeBypass<sao_eo_class_luma>(binVal1));
    }
};

template <>
struct Write<Element<sao_eo_class_chroma, ae>> // FL, cMax=3 | Bypass
{
    template <class H> static void go(Element<sao_eo_class_chroma, ae> fun, H &h)
    {
        const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
        const int ry = h[yCtb()] >> h[CtbLog2SizeY()];
        const int synVal = h[SaoEoClass(1, rx, ry)];// h[fun.v]; //GetValueForWrite<rem_intra_luma_pred_mode, H>::get(fun.v, h);
        int fl = 2;
        for (int binIdx = 0; binIdx < fl; binIdx++)
        {
            int binValue = 1 << (fl - 1 - binIdx);
            binValue &= synVal;
            binValue >>= (fl - 1 - binIdx);
            h(EncodeBypass<sao_eo_class_chroma>(binValue));
        }
        //int binVal0 = synVal <= 2 ? 1 : 0;
        //h(EncodeBypass<sao_eo_class_chroma>(binVal0));
        //int binVal1 = (synVal == 3 || synVal == 1) ? 1 : 0;
        //h(EncodeBypass<sao_eo_class_chroma>(binVal1));
    }
};

template <>
struct Write<Element<sao_band_position, ae>> // FL, cMax=31 | Bypass
{
    template <class H> static void go(Element<sao_band_position, ae> fun, H &h)
    {
        const int synVal = h[fun.v];//GetValueForWrite<sao_band_position, H>::get(fun.v, h);

        h(EncodeBypass<sao_band_position>((synVal & 0x10) ? 1 : 0));
        h(EncodeBypass<sao_band_position>((synVal & 0x08) ? 1 : 0));
        h(EncodeBypass<sao_band_position>((synVal & 0x04) ? 1 : 0));
        h(EncodeBypass<sao_band_position>((synVal & 0x02) ? 1 : 0));
        h(EncodeBypass<sao_band_position>((synVal & 0x01) ? 1 : 0));
    }
};


template <>
struct Write<Element<sao_offset_sign, ae>> // FL, cMax=1 | Bypass
{
    template <class H> static void go(Element<sao_offset_sign, ae> fun, H &h)
    {
        const int synVal = h[fun.v]; //GetValueForWrite<rem_intra_luma_pred_mode, H>::get(fun.v, h);
        int binVal0 = synVal == 1 ? 1 : 0;
        h(EncodeBypass<sao_offset_sign>(binVal0));
    }
};

template <>
struct Write<Element<sao_merge_left_flag, ae>> // FL, cMax=1 | Bypass
{
    template <class H> static void go(Element<sao_merge_left_flag, ae> fun, H &h)
    {
        const int synVal = h[fun.v];
        int binVal0 = synVal == 1 ? 1 : 0;
        h(EncodeDecision<sao_merge_left_flag>(binVal0, 0));
    }
};

template <>
struct Write<Element<sao_merge_up_flag, ae>> // FL, cMax=1 | Bypass
{
    template <class H> static void go(Element<sao_merge_up_flag, ae> fun, H &h)
    {
        const int synVal = h[fun.v];
        int binVal0 = synVal == 1 ? 1 : 0;
        h(EncodeDecision<sao_merge_up_flag>(binVal0, 0));
    }
};

// coding_unit( )

template <>
struct Write<Element<cu_skip_flag, ae>>
{
    template <class H> static void go(Element<cu_skip_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int ctxInc =
                h[cu_skip_flag(fun.v.x0 - 1, fun.v.y0)] +
                h[cu_skip_flag(fun.v.x0, fun.v.y0 - 1)];
        const int synVal = h[cu_skip_flag(fun.v.x0, fun.v.y0)];
        h(EncodeDecision<cu_skip_flag>(synVal, ctxInc));
    }
};

template <>
struct Write<Element<pred_mode_flag, ae>>
{
    template <class H> static void go(Element<pred_mode_flag, ae> fun, H &h)
    {
        const coding_quadtree *cqt = h;

        const int synVal = h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTRA ? 1 : 0;

        // FL cMax=1
        h(EncodeDecision<pred_mode_flag>(synVal, 0));
    }
};

template <>
struct Write<Element<part_mode, ae>>
{
    typedef ValueType<PartMode>::Type Type;
    typedef std::underlying_type<ValueType<PartMode>::Type>::type UnderlyingType;

    static UnderlyingType binVal1(Type value)
    {
        // convenient branch-free way to derive value of second bin
        return (~UnderlyingType(value) >> 1) & 1;
    }

    static UnderlyingType binVal2(Type value)
    {
        // convenient branch-free way to derive value of third bin
        return (~UnderlyingType(value) >> 2) & 1;
    }

    static UnderlyingType binVal3(Type value)
    {
        // convenient branch-free way to derive value of fourth bin
        return UnderlyingType(value) & 1;
    }

    template <class H> static void go(Element<part_mode, ae> fun, H &h)
    {
        const auto partMode = h[PartMode()];

        if (partMode == PART_2Nx2N)
        {
            h(EncodeDecision<part_mode>(1, 0));
        }
        else
        {
            const coding_quadtree *cqt = h;

            h(EncodeDecision<part_mode>(0, 0));

            if (h[current(CuPredMode(cqt->x0, cqt->y0))] == MODE_INTER)
            {
                h(EncodeDecision<part_mode>(binVal1(partMode), 1));

                if (cqt->log2CbSize == h[MinCbLog2SizeY()])
                {
                    if (cqt->log2CbSize > 3 && partMode != PART_2NxN)
                    {
                        h(EncodeDecision<part_mode>(binVal2(partMode), 2));
                    }
                }
                else if (h[amp_enabled_flag()])
                {
                    assert(cqt->log2CbSize > h[MinCbLog2SizeY()]);
                    assert(partMode != PART_NxN);

                    h(EncodeDecision<part_mode>(binVal2(partMode), 3));

                    if (isAmp(partMode))
                    {
                        h(EncodeBypass<part_mode>(binVal3(partMode)));
                    }
                }
            }
        }
    }
};

template <>
struct Write<Element<pcm_flag, ae>>
{
    template <class H> static void go(Element<pcm_flag, ae> fun, H &h)
    {
        const int binVal = h[fun.v];
        h(EncodeTerminate<pcm_flag>(binVal));
        if (binVal)
        {
            h(FinishCabac(0));
        }
    }
};


template <class F> struct EstimateRateLuma;


template <>
struct Write<Element<prev_intra_luma_pred_flag, ae>>
{
    template <class H> static void go(Element<prev_intra_luma_pred_flag, ae> fun, H &h)
    {
        // review: these const?
        StateCodedData *stateCodedData = h;
        CandModeList *candModeList = h;
        coding_quadtree const *cqt = h;

        int partIdx = 0;
        if (fun.v.x0 != cqt->x0) partIdx += 1;
        if (fun.v.y0 != cqt->y0) partIdx += 2;
        auto const value = stateCodedData->codedCu.IntraPredModeY(partIdx);

        if (!std::is_same<typename H::Tag, EstimateRateLuma<void>>::value)
        {
            Neighbourhood *neighbourhood = h;

            const int size = 1 << cqt->log2CbSize >> (h[PartMode()] == PART_2Nx2N ? 0 : 1);

            Turing::Rectangle partition{ fun.v.x0, fun.v.y0, size, size };

            neighbourhood->relocate(neighbourhood->snake, partition, neighbourhood->MinCbLog2SizeYMinus1);

            candModeList->populate(0, h, fun.v.x0, fun.v.y0);

            BlockData blockData;
            blockData.reset();
            blockData.skip = 0;
            blockData.intra.pcm = 0;
            blockData.intra.predModeY = value;
            coding_quadtree *cqt = h;
            blockData.CtDepth = cqt->cqtDepth;
            neighbourhood->snake.commitRectangle(Turing::Rectangle{ fun.v.x0, fun.v.y0, size, size }, blockData, neighbourhood->MinCbLog2SizeYMinus1);
        }

        int n = 0;
        int x;
        for (x=0; x<3; ++x)
        {
            if (value > (*candModeList)[x])
            {
                ++n;
            }
            else if (value == (*candModeList)[x])
            {
                h[prev_intra_luma_pred_flag(fun.v.x0, fun.v.y0)] = 1;
                h[mpm_idx(fun.v.x0, fun.v.y0)] = x;
                h(EncodeDecision<prev_intra_luma_pred_flag>(1, 0));
                return;
            }
        }
        h[rem_intra_luma_pred_mode(fun.v.x0, fun.v.y0)] = value - n;
        h[prev_intra_luma_pred_flag(fun.v.x0, fun.v.y0)] = 0;
        h(EncodeDecision<prev_intra_luma_pred_flag>(0, 0));
    }
};

template <>
struct Write<Element<mpm_idx, ae>>
{
    template <class H> static void go(Element<mpm_idx, ae> fun, H &h)
    {
        const int synVal = h[fun.v];
        assert(synVal >= 0);

        // Truncated Rice (TR) binarization process
        const int cMax = 2;
        const int cRiceParam = 0;

        const int prefixVal = synVal  >>  cRiceParam;

        // Prefix
        if (prefixVal < (cMax  >>  cRiceParam))
        {
            for (int binIdx = 0; binIdx < prefixVal; ++binIdx)
            {
                h(EncodeBypass<mpm_idx>(1));
            }
            h(EncodeBypass<mpm_idx>(0));
        }
        else
        {
            for (int binIdx = 0; binIdx < (cMax  >>  cRiceParam); ++binIdx)
            {
                h(EncodeBypass<mpm_idx>(1));
            }
        }

        assert(synVal <= cMax || !"should be no suffix, i.e. synVal < 3");
    }
};

template <>
struct Write<Element<rem_intra_luma_pred_mode, ae>> // FL, cMax=31 | Bypass
{
    template <class H> static void go(Element<rem_intra_luma_pred_mode, ae> fun, H &h)
    {
        const int synVal = h[fun.v];

        h(EncodeBypass<rem_intra_luma_pred_mode>((synVal & 0x10) ? 1 : 0));
        h(EncodeBypass<rem_intra_luma_pred_mode>((synVal & 0x08) ? 1 : 0));
        h(EncodeBypass<rem_intra_luma_pred_mode>((synVal & 0x04) ? 1 : 0));
        h(EncodeBypass<rem_intra_luma_pred_mode>((synVal & 0x02) ? 1 : 0));
        h(EncodeBypass<rem_intra_luma_pred_mode>((synVal & 0x01) ? 1 : 0));
    }
};

template <>
struct Write<Element<intra_chroma_pred_mode, ae>>
{
    template <class H> static void go(Element<intra_chroma_pred_mode, ae> fun, H &h)
    {
        const int synVal = h[intra_chroma_pred_mode()];

        if (synVal == 4)
        {
            h(EncodeDecision<intra_chroma_pred_mode>(0, 0));
        }
        else
        {
            h(EncodeDecision<intra_chroma_pred_mode>(1, 0));
            h(EncodeBypass<intra_chroma_pred_mode>(synVal >> 1));
            h(EncodeBypass<intra_chroma_pred_mode>(synVal & 1));
        }
    }
};

template <>
struct Write<Element<rqt_root_cbf, ae>>
{
    template <class H> static void go(Element<rqt_root_cbf, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<rqt_root_cbf>(synVal, 0));
    }
};


// prediction_unit( )

template <>
struct Write<Element<merge_flag, ae>>
{
    template <class H> static void go(Element<merge_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<merge_flag>(synVal, 0));
    }
};

template <>
struct Write<Element<merge_idx, ae>>
{
    template <class H> static void go(Element<merge_idx, ae> fun, H &h)
    {
        // TR	cMax = MaxNumMergeCand - 1, cRiceParam = 0
        // no suffix
        const int cMax = h[MaxNumMergeCand()] - 1;

        const int synVal = h[fun.v];
        const int prefixVal = synVal;

        int binIdx = 0;
        for ( ; binIdx < prefixVal; ++binIdx)
        {
            encodeBin(binIdx, 1, h);
        }
        if (prefixVal < cMax)
        {
            encodeBin(binIdx, 0, h);
        }
    }
    template <class H> static void encodeBin(int binIdx, int binVal, H &h)
    {
        if (binIdx == 0)
        {
            h(EncodeDecision<merge_idx>(binVal, 0));
        }
        else
        {
            h(EncodeBypass<merge_idx>(binVal));
        }
    }
};


template <>
struct Write<Element<inter_pred_idc, ae>>
{
    template <class H> static void go(Element<inter_pred_idc, ae> fun, H &h)
    {
        const int synVal = h[fun.v];
        const prediction_unit pu = *static_cast<prediction_unit *>(h);
        const coding_quadtree cqt = *static_cast<coding_quadtree *>(h);

        if (pu.nPbW + pu.nPbH != 12)
        {
            if (synVal == PRED_BI)
            {
                h(EncodeDecision<inter_pred_idc>(1, cqt.cqtDepth));
            }
            else
            {
                h(EncodeDecision<inter_pred_idc>(0, cqt.cqtDepth));
                h(EncodeDecision<inter_pred_idc>(synVal, 4));
            }
        }
        else
        {
            assert(synVal != PRED_BI);
            h(EncodeDecision<inter_pred_idc>(synVal, 4));
        }
    }
};


// transform_tree( )

template <>
struct Write<Element<split_transform_flag, ae>>
{
    template <class H> static void go(Element<split_transform_flag, ae> fun, H &h)
    {
        transform_tree &tt = *static_cast<transform_tree *>(h);
        const auto ctxInc = 5 - tt.log2TrafoSize;

        StateCodedData *stateCodedData = h;

        // review: have handler read h[split_transform_flag()] direct from coded data during write
        h[fun.v] = stateCodedData->transformTree.word0().split_transform_flag;

        const int synVal = h[fun.v];

        h(EncodeDecision<split_transform_flag>(synVal, ctxInc));
    }
};


template <>
struct Write<Element<cbf_luma, ae>>
{
    template <class H> static void go(Element<cbf_luma, ae> fun, H &h)
    {
        const int binVal = h[fun.v];

        StateEncode &stateEncode = *static_cast<StateEncode *>(h);

        transform_tree &tt = *static_cast<transform_tree *>(h);
        const int ctxInc = tt.trafoDepth == 0 ? 1 : 0;

        h(EncodeDecision<cbf_luma>(binVal, ctxInc));
    }
};

template <class V>
struct WriteCbfC
{
    template <class H> static void go(Element<V, ae> e, H &h)
    {
        auto const binVal = h[e.v];
        auto const ctxInc = e.v.trafoDepth;
        h(EncodeDecision<V>(binVal, ctxInc));
    }
};

template <> struct Write<Element<cbf_cb, ae>> : WriteCbfC<cbf_cb> { };

template <> struct Write<Element<cbf_cr, ae>> : WriteCbfC<cbf_cr> { };



template <>
struct Write<Element<cu_qp_delta_abs, ae>>
{
    template <class H> static void go(Element<cu_qp_delta_abs, ae> fun, H &h)
    {
        int valueToWrite    = h[fun.v];
        int absValue        = abs(valueToWrite);
        int signFlagToWrite = h[cu_qp_delta_sign_flag()];

        int prefixVal = std::min<int>(absValue, 5);

        // Clause 9.3.3.2 for truncated Rice binarisation process, no suffix needed since cRiceParam=0
        int cMax = 5;
        int cRiceParam = 0;

        // Prefix
        int trPrefixVal = prefixVal >> cRiceParam;
        if(trPrefixVal < (cMax >> cRiceParam))
        {
            int binIdx = 0;
            for(; binIdx < trPrefixVal; binIdx++)
                h(EncodeDecision<cu_qp_delta_abs>(1, binIdx ? 1 : 0));
            h(EncodeDecision<cu_qp_delta_abs>(0, binIdx ? 1 : 0));
        }
        else
        {
            h(EncodeDecision<cu_qp_delta_abs>(1, 0));
            for(int binIdx = 1; binIdx < (cMax >> cRiceParam); binIdx++)
                h(EncodeDecision<cu_qp_delta_abs>(1, 1));
        }

        // Suffix for 9.3.3.9: EGk with k = 0
        if(prefixVal > 4)
        {
            int suffixVal = absValue - 5;
            int k = 0;
            int stopLoop = 0;
            do
            {
                if(suffixVal >= (1 << k))
                {
                    h(EncodeBypass<cu_qp_delta_abs>(1));
                    suffixVal -= (1 << k);
                    k++;
                }
                else
                {
                    h(EncodeBypass<cu_qp_delta_abs>(0));
                    while(k--)
                        h(EncodeBypass<cu_qp_delta_abs>((suffixVal >> k) & 1));
                    stopLoop = 1;
                }
            }while(!stopLoop);
        }
        h[IsCuQpDeltaCoded()] = 1;
        h[CuQpDeltaVal()] = absValue * (1 - 2 * signFlagToWrite);
        int qpY = static_cast<QpState* >(h)->getQp(0) - h[QpBdOffsetY()];
        qpY = Clip3(0, 51, qpY);
        static_cast<QpState* >(h)->setCodedQp(qpY);
    }
};

template <>
struct Write<Element<cu_qp_delta_sign_flag, ae>>
{
    template <class H> static void go(Element<cu_qp_delta_sign_flag, ae> fun, H &h)
    {
        int sign = h[fun.v];
        h(EncodeBypass<cu_qp_delta_sign_flag>(sign));
    }
};


// mvd_coding( )

template <>
struct Write<Element<abs_mvd_greater0_flag, ae>>
{
    template <class H> static void go(Element<abs_mvd_greater0_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<abs_mvd_greater0_flag>(synVal, 0));
    }
};

template <>
struct Write<Element<abs_mvd_greater1_flag, ae>>
{
    template <class H> static void go(Element<abs_mvd_greater1_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<abs_mvd_greater1_flag>(synVal, 0));
    }
};

template <>
struct Write<Element<abs_mvd_minus2, ae>>
{
    template <class H> static void go(Element<abs_mvd_minus2, ae> fun, H &h)
    {
        unsigned synVal = h[fun.v];

        // EGk
        int k = 1;
        int absV = synVal;
        int stopLoop = 0;
        do
        {
            if (absV >= (1 << k))
            {
                h(EncodeBypass<abs_mvd_minus2>(1));
                absV = absV - (1 << k);
                k++;
            }
            else
            {
                h(EncodeBypass<abs_mvd_minus2>(0));
                while (k--)
                {
                    h(EncodeBypass<abs_mvd_minus2>((absV >> k) & 1));
                }
                stopLoop = 1;
            }
        }
        while (!stopLoop);
    }
};

template <>
struct Write<Element<mvd_sign_flag, ae>>
{
    template <class H> static void go(Element<mvd_sign_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeBypass<mvd_sign_flag>(synVal));
    }
};

template <>
struct Write<Element<mvp_l0_flag, ae>>
{
    template <class H> static void go(Element<mvp_l0_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<mvp_l0_flag>(synVal, 0));
    }
};

template <>
struct Write<Element<mvp_l1_flag, ae>>
{
    template <class H> static void go(Element<mvp_l1_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<mvp_l1_flag>(synVal, 0));
    }
};


// residual_coding( )

template <>
struct Write<Element<transform_skip_flag, ae>>
{
    template <class H> static void go(Element<transform_skip_flag, ae> fun, H &h)
    {
        auto const synVal = h[fun.v];
        if (fun.v.cIdx == 0)
        {
            h(EncodeDecision<transform_skip_flag_0>(synVal, 0));
        }
        else
        {
            h(EncodeDecision<transform_skip_flag_1>(synVal, 0));
        }
    }
};


template <class V>
struct WriteLastSigPrefix
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        const int synVal = h[fun.v];

        const residual_coding& rc = *static_cast<residual_coding *>(h);

        // Truncated Rice (TR) binarization process
        const int cMax = ( rc.log2TrafoSize  <<  1 ) - 1;
        const int cRiceParam = 0;

        const int prefixVal = synVal  >>  cRiceParam;

        // Prefix
        if (prefixVal < (cMax  >>  cRiceParam))
        {
            int binIdx = 0;
            for ( ; binIdx < prefixVal; ++binIdx)
            {
                h(EncodeDecision<V>(1, ctxInc(binIdx, rc.cIdx, rc.log2TrafoSize)));
            }
            h(EncodeDecision<V>(0, ctxInc(binIdx, rc.cIdx, rc.log2TrafoSize)));
        }
        else
        {
            for (int binIdx = 0; binIdx < (cMax  >>  cRiceParam); ++binIdx)
            {
                h(EncodeDecision<V>(1, ctxInc(binIdx, rc.cIdx, rc.log2TrafoSize)));
            }
        }

        if (synVal > cMax)
        {
            // Suffix
            assert(0);
        }
    }

    static int ctxInc(int binIdx, int cIdx, int log2TrafoSize)
    {
        const int ctxOffset = cIdx ? 15 : ( 3 * ( log2TrafoSize - 2 ) + ( ( log2TrafoSize - 1 )  >>  2 ) );
        const int ctxShift = cIdx ? ( log2TrafoSize - 2 ) : ( ( log2TrafoSize + 1 )  >>  2 );
        return ( binIdx  >>  ctxShift ) + ctxOffset;
    }
};

template <>
struct Write<Element<last_sig_coeff_x_prefix, ae>>
:
WriteLastSigPrefix<last_sig_coeff_x_prefix>
{
};

template <>
struct Write<Element<last_sig_coeff_y_prefix, ae>>
:
WriteLastSigPrefix<last_sig_coeff_y_prefix>
{
};

template <class V, class Prefix>
struct WriteLastSigSuffix
{
    template <class H> static void go(Element<V, ae> fun, H &h)
    {
        const int synVal = h[fun.v];
        const int prefixVal = h[Prefix()];
        // FL
        const int cMax = ( 1  <<  ( ( prefixVal  >>  1 ) - 1 ) ) - 1;
        int n = ( ( prefixVal  >>  1 ) - 1 );
        while (n--)
        {
            h(EncodeBypass<V>( (synVal >> n) & 1));
        }
    }
};

template <>
struct Write<Element<last_sig_coeff_x_suffix, ae>>
:
WriteLastSigSuffix<last_sig_coeff_x_suffix, last_sig_coeff_x_prefix>
{
};

template <>
struct Write<Element<last_sig_coeff_y_suffix, ae>>
:
WriteLastSigSuffix<last_sig_coeff_y_suffix, last_sig_coeff_y_prefix>
{
};

// These two integrated into EncodeResidual
template <> struct Write<Element<coeff_abs_level_greater1_flag, ae>>;
template <> struct Write<Element<coeff_abs_level_greater2_flag, ae>>;

template <>
struct Write<Element<coeff_sign_flag, ae>>
{
    template <class H> static void go(Element<coeff_sign_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeBypass<coeff_sign_flag>(synVal));
    }
};

template <>
struct Write<Element<cu_transquant_bypass_flag, ae>>
{
    template <class H> static void go(Element<cu_transquant_bypass_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<cu_transquant_bypass_flag>(synVal, 0));
    }
};

template <>
struct Write<Element<sig_coeff_flag, ae>>
{
    template <class H> static void go(Element<sig_coeff_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        const residual_coding& rc = *static_cast<residual_coding *>(h);
        h(EncodeDecision<sig_coeff_flag>(synVal, ctxInc(h, rc.cIdx, fun.v.xC, fun.v.yC, rc.log2TrafoSize)));
    }

    template <class H> static int ctxInc(H &h, int cIdx, int xC, int yC, int log2TrafoSize)
    {
        int sigCtx;
        if (log2TrafoSize == 2)
        {
            const int ctxIdxMap[16] =
            {
                    0, 1, 4, 5,
                    2, 3, 4, 5,
                    6, 6, 8, 8,
                    7, 7, 8, 8
            };
            sigCtx = ctxIdxMap[(yC << 2) + xC];
        }
        else if (xC + yC == 0)
        {
            sigCtx = 0;
        }
        else
        {
            const int xS = xC >> 2;
            const int yS = yC >> 2;
            int prevCsbf = 0;

            if (xS < (1 << (log2TrafoSize - 2)) - 1)
            {
                prevCsbf += h[coded_sub_block_flag(xS + 1, yS)];
            }
            if (yS < (1 << (log2TrafoSize - 2)) - 1)
            {
                prevCsbf += (h[coded_sub_block_flag(xS, yS + 1)] << 1);
            }

            const int xP = xC & 3;
            const int yP = yC & 3;

            if (prevCsbf == 0)
            {
                sigCtx = (xP + yP == 0) ? 2 : (xP + yP < 3) ? 1 : 0;
            }
            else if (prevCsbf == 1)
            {
                sigCtx = (yP == 0) ? 2 : (yP == 1) ? 1 : 0;
            }
            else if (prevCsbf == 2)
            {
                sigCtx = (xP == 0) ? 2 : (xP == 1) ? 1 : 0;
            }
            else
            {
                assert(prevCsbf == 3);
                sigCtx = 2;
            }

            if (cIdx == 0)
            {
                if (xS + yS > 0)
                {
                    sigCtx += 3;
                }

                if (log2TrafoSize == 3)
                {
                    sigCtx += (h[scanIdx()] == 0) ? 9 : 15;
                }
                else
                {
                    sigCtx += 21;
                }
            }
            else
            {
                if (log2TrafoSize == 3)
                {
                    sigCtx += 9;
                }
                else
                {
                    sigCtx += 12;
                }
            }
        }

        assert(sigCtx >= 0);
        assert(sigCtx < 54);

        if (cIdx == 0)
        {
            return sigCtx;
        }
        else
        {
            return 27 + sigCtx;
        }
    }
};



template <>
struct Write<Element<coded_sub_block_flag, ae>>
{
    template <class H> static void go(Element<coded_sub_block_flag, ae> fun, H &h)
    {
        // FL cMax=1
        const int synVal = h[fun.v];
        h(EncodeDecision<coded_sub_block_flag>(synVal, ctxInc(h, fun.v)));
    }

    template <class H> static int ctxInc(H& h, coded_sub_block_flag v)
    {
        StateEncode& stateEncode = *static_cast<StateEncode *>(h);
        const residual_coding& rc = *static_cast<residual_coding *>(h);

        int csbfCtx = 0;
        if (v.xS < ( 1  <<  ( rc.log2TrafoSize - 2 ) ) - 1)
        {
            csbfCtx += h[ coded_sub_block_flag( v.xS + 1, v.yS )];
        }
        if (v.yS < ( 1  <<  ( rc.log2TrafoSize - 2 ) ) - 1)
        {
            csbfCtx += h[ coded_sub_block_flag( v.xS , v.yS + 1)];
        }
        return (rc.cIdx?2:0) + std::min(csbfCtx, 1);
    }
};

template <>
struct Write<Element<coeff_abs_level_remaining, ae>>
{
    template <class H> static void go(Element<coeff_abs_level_remaining, ae> fun, H &h);

    template <class H> static void write(int synVal, int cAbsLevel, int &cRiceParam, H &h)
    {
        const int cMax = 4 << cRiceParam;

        // coeff_abs_level_remaining prefix
        const int prefixVal = std::min(cMax, synVal);

        bool fourOnes = false;
        {
            // TR
            const int synValTr = prefixVal;
            const int prefixValTr = synValTr >> cRiceParam;

            // TR prefix
            if (prefixValTr < (cMax >> cRiceParam))
            {
                for (int binIdx = 0; binIdx < prefixValTr; ++binIdx)
                {
                    h(EncodeBypass<coeff_abs_level_remaining>(1));
                }
                h(EncodeBypass<coeff_abs_level_remaining>(0));
            }
            else
            {
                for (int binIdx = 0; binIdx < (cMax >> cRiceParam); ++binIdx)
                {
                    h(EncodeBypass<coeff_abs_level_remaining>(1));
                }

                assert((cMax >> cRiceParam) == 4);
                fourOnes = true;
            }

            // TR suffix
            if (cMax > synValTr)
            {
                assert(!fourOnes);
                const int suffixValTr = synValTr - ((prefixValTr) << cRiceParam);
                int shift = cRiceParam;
                while (shift--)
                {
                    h(EncodeBypass<coeff_abs_level_remaining>((suffixValTr >> shift) & 1));
                }
            }
        }

        if (fourOnes)
        {
            // coeff_abs_level_remaining suffix
            const int suffixVal = synVal - cMax;

            // EGk
            int k = cRiceParam + 1;
            int absV = std::abs(suffixVal);
            int stopLoop = 0;
            do
            {
                if (absV >= (1 << k))
                {
                    h(EncodeBypass<coeff_abs_level_remaining>(1));
                    absV = absV - (1 << k);
                    k++;
                }
                else
                {
                    h(EncodeBypass<coeff_abs_level_remaining>(0));
                    while (k--)
                    {
                        h(EncodeBypass<coeff_abs_level_remaining>((absV >> k) & 1));
                    }
                    stopLoop = 1;
                }
            } while (!stopLoop);
        }

        cRiceParam = std::min(cRiceParam + (cAbsLevel > (3 * (1 << cRiceParam)) ? 1 : 0), 4);
    }
};


template <class F> struct EstimateRate;


template <>
struct EstimateRate<Element<coeff_abs_level_remaining, ae>>
{
    // this function no longer used -  nowinlined into residual encoding
    template <class H> static void go(Element<coeff_abs_level_remaining, ae> fun, H &h);

    template <class H> static inline void write(int synVal, int cAbsLevel, int &cRiceParam, H &h)
    {
        StateEstimateRate *stateEstimateRate = h;
        
        int bits = cRiceParam + 1 + 3;
        auto const a = (synVal >> cRiceParam) - 3;

        // TR prefix
        if (a < 0)
        {
            bits += a;
        }
        else
        {
            unsigned long index;
#ifdef WIN32
            _BitScanReverse(&index, a + 1);
#elif __GNUC__
            index = __builtin_clz(a + 1) ^ 31;
#else
            index = 15;
            while (!((a + 1) & (1 << index)))
                --index;
#endif

            bits += 2 * index;
        }


        cRiceParam = std::min(cRiceParam + (cAbsLevel > (3 * (1 << cRiceParam)) ? 1 : 0), 4);

        stateEstimateRate->rate += Cost::make(bits, 0);
    }
};

#endif
