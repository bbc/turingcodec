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

// CTU syntax functions.
// In general these are more-or-less copy/pasted from the HEVC standard text and should not be modified.


#ifndef INCLUDED_SyntaxCtu_hpp
#define INCLUDED_SyntaxCtu_hpp

#include "Syntax.h"
#include "SyntaxElements.h"
#include "HevcTypes.h"
#include "Global.h"
#include "GlobalState.h"


template <class H>
void Syntax<coding_tree_unit>::go(const coding_tree_unit&, H &h)
{
    h[xCtb()] = (h[CtbAddrInRs()] % h[PicWidthInCtbsY()]) << h[CtbLog2SizeY()];
    h[yCtb()] = (h[CtbAddrInRs()] / h[PicWidthInCtbsY()]) << h[CtbLog2SizeY()];
    if (h[slice_sao_luma_flag()] || h[slice_sao_chroma_flag()])
        h(sao(h[xCtb()] >> h[CtbLog2SizeY()], h[yCtb()] >> h[CtbLog2SizeY()]));
    h(coding_quadtree(h[xCtb()], h[yCtb()], h[CtbLog2SizeY()], 0));
}


template <class H>
void Syntax<sao>::go(const sao &e, H &h)
{
    if (e.rx > 0)
    {
        const bool leftCtbInSliceSeg = h[CtbAddrInRs()] > h[SliceAddrRs()];
        const bool leftCtbInTile = h[TileId(h[CtbAddrInTs()])] == h[TileId(h[CtbAddrRsToTs(h[CtbAddrInRs()] - 1)])];
        if (leftCtbInSliceSeg && leftCtbInTile)
            h(sao_merge_left_flag(), ae(v));
    }
    if (e.ry > 0 && !h[sao_merge_left_flag()])
    {
        const bool upCtbInSliceSeg = (h[CtbAddrInRs()] - h[PicWidthInCtbsY()]) >= h[SliceAddrRs()];
        const bool upCtbInTile = h[TileId(h[CtbAddrInTs()])] == h[TileId(h[CtbAddrRsToTs(h[CtbAddrInRs()] - h[PicWidthInCtbsY()])])];
        if (upCtbInSliceSeg && upCtbInTile)
            h(sao_merge_up_flag(), ae(v));
    }
    if (!h[sao_merge_up_flag()] && !h[sao_merge_left_flag()])
        for (int cIdx = 0; cIdx < 3; cIdx++)
            if ((h[slice_sao_luma_flag()] && cIdx == 0) || (h[slice_sao_chroma_flag()] && cIdx > 0))
            {
                if (cIdx == 0)
                    h(sao_type_idx_luma(), ae(v));
                else if (cIdx == 1)
                    h(sao_type_idx_chroma(), ae(v));
                if (h[SaoTypeIdx(cIdx, e.rx, e.ry)] != 0)
                {
                    for (int i = 0; i < 4; i++)
                        h(sao_offset_abs(cIdx, e.rx, e.ry, i), ae(v));
                    if (h[SaoTypeIdx(cIdx, e.rx, e.ry)] == 1)
                    {
                        for (int i = 0; i < 4; i++)
                            if (h[sao_offset_abs(cIdx, e.rx, e.ry, i)] != 0)
                                h(sao_offset_sign(cIdx, e.rx, e.ry, i), ae(v));
                        h(sao_band_position(cIdx, e.rx, e.ry), ae(v));
                    }
                    else
                    {
                        if (cIdx == 0)
                            h(sao_eo_class_luma(), ae(v));
                        if (cIdx == 1)
                            h(sao_eo_class_chroma(), ae(v));
                    }
                }
            }
}


template <class H>
void Syntax<coding_quadtree>::go(const coding_quadtree &e, H &h)
{
    if (e.x0 + (1 << e.log2CbSize) <= h[pic_width_in_luma_samples()] &&
            e.y0 + (1 << e.log2CbSize) <= h[pic_height_in_luma_samples()] &&
            e.log2CbSize > h[MinCbLog2SizeY()])
    {
        h(split_cu_flag(e.x0, e.y0), ae(v));
    }
    if (h[cu_qp_delta_enabled_flag()] && e.log2CbSize >= h[Log2MinCuQpDeltaSize()])
    {
        h[IsCuQpDeltaCoded()] = 0;
        h[CuQpDeltaVal()] = 0;
    }
    if (h[cu_chroma_qp_offset_enabled_flag()] && e.log2CbSize >= h[Log2MinCuChromaQpOffsetSize()])
    {
        h[IsCuChromaQpOffsetCoded()] = 0;
    }
    if (h[split_cu_flag(e.x0, e.y0)])
    {
        int x1 = e.x0 + (1 << (e.log2CbSize - 1));
        int y1 = e.y0 + (1 << (e.log2CbSize - 1));

        h(coding_quadtree(e.x0, e.y0, e.log2CbSize - 1, e.cqtDepth + 1));

        if (x1 < h[pic_width_in_luma_samples()])
            h(coding_quadtree(x1, e.y0, e.log2CbSize - 1, e.cqtDepth + 1));
        else
            h(Deleted<coding_quadtree, Right>(x1, e.y0, e.log2CbSize - 1, e.cqtDepth + 1));

        if (y1 < h[pic_height_in_luma_samples()])
            h(coding_quadtree(e.x0, y1, e.log2CbSize - 1, e.cqtDepth + 1));
        else
            h(Deleted<coding_quadtree, Right>(e.x0, y1, e.log2CbSize - 1, e.cqtDepth + 1));

        if (x1 < h[pic_width_in_luma_samples()] && y1 < h[pic_height_in_luma_samples()])
            h(coding_quadtree(x1, y1, e.log2CbSize - 1, e.cqtDepth + 1));
        else
            if (x1 < h[pic_width_in_luma_samples()])
                h(Deleted<coding_quadtree, Down>(x1, y1, e.log2CbSize - 1, e.cqtDepth + 1));
            else
                h(Deleted<coding_quadtree, Right>(x1, y1, e.log2CbSize - 1, e.cqtDepth + 1));
    }
    else
        h(coding_unit(e.x0, e.y0, e.log2CbSize));
}


template <class H>
void Syntax<coding_unit>::go(const coding_unit &cu, H &h)
{
    if (h[transquant_bypass_enabled_flag()])
        h(cu_transquant_bypass_flag(), ae(v));

    if (h[slice_type()] != I)
        h(cu_skip_flag(cu.x0, cu.y0), ae(v));

    const int nCbS = (1 << cu.log2CbSize);

    if (h[current(cu_skip_flag(cu.x0, cu.y0))])
        h(prediction_unit(cu.x0, cu.y0, nCbS, nCbS));
    else
    {
        if (h[slice_type()] != I)
            h(pred_mode_flag(), ae(v));

        if (h[current(CuPredMode(cu.x0, cu.y0))] != MODE_INTRA || cu.log2CbSize == h[MinCbLog2SizeY()])
            h(part_mode(), ae(v));

        if (h[current(CuPredMode(cu.x0, cu.y0))] == MODE_INTRA)
        {
            if (h[PartMode()] == PART_2Nx2N && h[pcm_enabled_flag()] &&
                    cu.log2CbSize >= h[Log2MinIpcmCbSizeY()] &&
                    cu.log2CbSize <= h[Log2MaxIpcmCbSizeY()])
                h(pcm_flag(cu.x0, cu.y0), ae(v));

            if (h[current(pcm_flag(cu.x0, cu.y0))])
            {
                h(pcm_alignment_one_bit(), f(1)); // https://hevc.hhi.fraunhofer.de/trac/hevc/ticket/1013#comment:3
                while (!h[byte_aligned()])
                    h(pcm_alignment_zero_bit(), f(1));
                h(pcm_sample(cu.x0, cu.y0, cu.log2CbSize));
            }
            else
            {
                const int pbOffset = (h[PartMode()] == PART_NxN) ? (nCbS / 2) : nCbS;

                for (int j = 0; j < nCbS; j = j + pbOffset)
                    for (int i = 0; i < nCbS; i = i + pbOffset)
                        h(prev_intra_luma_pred_flag(cu.x0 + i, cu.y0 + j), ae(v));

                for (int j = 0; j < nCbS; j = j + pbOffset)
                    for (int i = 0; i < nCbS; i = i + pbOffset)
                        if (h[prev_intra_luma_pred_flag(cu.x0 + i, cu.y0 + j)])
                            h(mpm_idx(cu.x0 + i, cu.y0 + j), ae(v));
                        else
                            h(rem_intra_luma_pred_mode(cu.x0 + i, cu.y0 + j), ae(v));

                if (h[ChromaArrayType()] == 3)
                {
                    for (int j = 0; j < nCbS; j = j + pbOffset)
                        for (int i = 0; i < nCbS; i = i + pbOffset)
                            h(intra_chroma_pred_mode(cu.x0 + i, cu.y0 + j), ae(v));
                }
                else if (h[ChromaArrayType()] != 0)
                {
                    h(intra_chroma_pred_mode(cu.x0, cu.y0), ae(v));
                }
            }
        }
        else
        {
            if (h[PartMode()] == PART_2Nx2N)
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS, nCbS));
            }
            else if (h[PartMode()] == PART_2NxN)
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS, nCbS / 2));
                h(prediction_unit(cu.x0, cu.y0 + (nCbS / 2), nCbS, nCbS / 2));
            }
            else if (h[PartMode()] == PART_Nx2N)
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS / 2, nCbS));
                h(prediction_unit(cu.x0 + (nCbS / 2), cu.y0, nCbS / 2, nCbS));
            }
            else if (h[PartMode()] == PART_2NxnU)
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS, nCbS / 4));
                h(prediction_unit(cu.x0, cu.y0 + (nCbS / 4), nCbS, nCbS * 3 / 4));
            }
            else if (h[PartMode()] == PART_2NxnD)
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS, nCbS * 3 / 4));
                h(prediction_unit(cu.x0, cu.y0 + (nCbS * 3 / 4), nCbS, nCbS / 4));
            }
            else if (h[PartMode()] == PART_nLx2N)
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS / 4, nCbS));
                h(prediction_unit(cu.x0 + (nCbS / 4), cu.y0, nCbS * 3 / 4, nCbS));
            }
            else if (h[PartMode()] == PART_nRx2N)
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS * 3 / 4, nCbS));
                h(prediction_unit(cu.x0 + (nCbS * 3 / 4), cu.y0, nCbS / 4, nCbS));
            }
            else /* PART_NxN */
            {
                h(prediction_unit(cu.x0, cu.y0, nCbS / 2, nCbS / 2));
                h(prediction_unit(cu.x0 + (nCbS / 2), cu.y0, nCbS / 2, nCbS / 2));
                h(prediction_unit(cu.x0, cu.y0 + (nCbS / 2), nCbS / 2, nCbS / 2));
                h(prediction_unit(cu.x0 + (nCbS / 2), cu.y0 + (nCbS / 2), nCbS / 2, nCbS / 2));
            }
        }

        if (!h[current(pcm_flag(cu.x0, cu.y0))])
        {
            if (h[current(CuPredMode(cu.x0, cu.y0))] != MODE_INTRA &&
                    !(h[PartMode()] == PART_2Nx2N && h[merge_flag(cu.x0, cu.y0)]))
            {
                h(rqt_root_cbf(), ae(v));
            }

            h[MaxTrafoDepth()] = (h[current(CuPredMode(cu.x0, cu.y0))] == MODE_INTRA)
                        ? h[max_transform_hierarchy_depth_intra()] + h[IntraSplitFlag()]
                                                                       : h[max_transform_hierarchy_depth_inter()];

            h(IfCbf<rqt_root_cbf, transform_tree>{ rqt_root_cbf(), transform_tree(cu.x0, cu.y0, cu.x0, cu.y0, cu.log2CbSize, 0, 0) });
        }
    }
}


template <class H>
void Syntax<prediction_unit>::go(const prediction_unit &e, H &h)
{
    if (h[cu_skip_flag(e.x0, e.y0)])
    {
        if (h[MaxNumMergeCand()] > 1)
            h(merge_idx(e.x0, e.y0), ae(v));
    }
    else
    {
        /* MODE_INTER */
        h(merge_flag(e.x0, e.y0), ae(v));
        if (h[merge_flag(e.x0, e.y0)])
        {
            if (h[MaxNumMergeCand()] > 1)
                h(merge_idx(e.x0, e.y0), ae(v));
        }
        else
        {
            if (h[slice_type()] == B)
                h(inter_pred_idc(e.x0, e.y0), ae(v));

            if (h[inter_pred_idc(e.x0, e.y0)] != PRED_L1)
            {
                if (h[num_ref_idx_l0_active_minus1()] > 0)
                    h(ref_idx_l0(e.x0, e.y0), ae(v));

                h(mvd_coding(e.x0, e.y0, 0));
                h(mvp_l0_flag(e.x0, e.y0), ae(v));
            }

            if (h[inter_pred_idc(e.x0, e.y0)] != PRED_L0)
            {
                if (h[num_ref_idx_l1_active_minus1()] > 0)
                    h(ref_idx_l1(e.x0, e.y0), ae(v));

                if (h[mvd_l1_zero_flag()] && h[inter_pred_idc(e.x0, e.y0)] == PRED_BI)
                {
                    h[Mvd(L1, e.x0, e.y0)] = MotionVector{ 0, 0 };
                }
                else
                    h(mvd_coding(e.x0, e.y0, 1));

                h(mvp_l1_flag(e.x0, e.y0), ae(v));
            }
        }
    }
}


template <class H>
void Syntax<pcm_sample>::go(const pcm_sample &e, H &h)
{
    for (int i = 0; i < 1 << (e.log2CbSize << 1); i++)
        h(pcm_sample_luma(i), uv());

    if (h[ChromaArrayType()] != 0)
        for (int i = 0; i < ((2 << (e.log2CbSize << 1)) / (h[SubWidthC()] * h[SubHeightC()])); i++)
            h(pcm_sample_chroma(i), uv());
}


template <class H>
void Syntax<transform_tree>::go(const transform_tree &tt, H &h)
{
    if (tt.log2TrafoSize <= h[MaxTbLog2SizeY()] &&
            tt.log2TrafoSize > h[MinTbLog2SizeY()] &&
            tt.trafoDepth < h[MaxTrafoDepth()] && !(h[IntraSplitFlag()] && (tt.trafoDepth == 0)))
    {
        h(split_transform_flag(tt.x0, tt.y0, tt.trafoDepth), ae(v));
    }

    if ((tt.log2TrafoSize > 2 && h[ChromaArrayType()] != 0) || h[ChromaArrayType()] == 3 )
    {
        if (tt.trafoDepth == 0 || h[cbf_cb(tt.xBase, tt.yBase, tt.trafoDepth - 1)])
        {
            h(cbf_cb(tt.x0, tt.y0, tt.trafoDepth), ae(v));

            if (h[ChromaArrayType()] == 2 && (!h[split_transform_flag(tt.x0, tt.y0, tt.trafoDepth)] || tt.log2TrafoSize == 3))
                h(cbf_cb(tt.x0, tt.y0 + (1 << (tt.log2TrafoSize - 1)), tt.trafoDepth), ae(v));
        }
        if (tt.trafoDepth == 0 || h[cbf_cr(tt.xBase, tt.yBase, tt.trafoDepth - 1)])
        {
            h(cbf_cr(tt.x0, tt.y0, tt.trafoDepth), ae(v));

            if (h[ChromaArrayType()] == 2 && (!h[split_transform_flag(tt.x0, tt.y0, tt.trafoDepth)] || tt.log2TrafoSize == 3))
                h(cbf_cr(tt.x0, tt.y0 + (1 << (tt.log2TrafoSize - 1)), tt.trafoDepth), ae(v));
        }
    }

    if (h[split_transform_flag(tt.x0, tt.y0, tt.trafoDepth)])
    {
        const int x1 = tt.x0 + (1 << (tt.log2TrafoSize - 1));
        const int y1 = tt.y0 + (1 << (tt.log2TrafoSize - 1));
        h(transform_tree(tt.x0, tt.y0, tt.x0, tt.y0, tt.log2TrafoSize - 1, tt.trafoDepth + 1, 0));
        h(transform_tree(x1, tt.y0, tt.x0, tt.y0, tt.log2TrafoSize - 1, tt.trafoDepth + 1, 1));
        h(transform_tree(tt.x0, y1, tt.x0, tt.y0, tt.log2TrafoSize - 1, tt.trafoDepth + 1, 2));
        h(transform_tree(x1, y1, tt.x0, tt.y0, tt.log2TrafoSize - 1, tt.trafoDepth + 1, 3));
    }
    else
    {
        if (h[current(CuPredMode(tt.x0, tt.y0))] == MODE_INTRA || tt.trafoDepth != 0 ||
                h[cbf_cb(tt.x0, tt.y0, tt.trafoDepth)] || h[cbf_cr(tt.x0, tt.y0, tt.trafoDepth)] ||
                (h[ChromaArrayType()] == 2 &&
                        (h[cbf_cb(tt.x0, tt.y0 + (1 << (tt.log2TrafoSize - 1)), tt.trafoDepth)] ||
                                h[cbf_cr(tt.x0, tt.y0 + (1 << (tt.log2TrafoSize - 1)), tt.trafoDepth)])))
        {
            h(cbf_luma(tt.x0, tt.y0, tt.trafoDepth), ae(v));
        }

        h(transform_unit(tt.x0, tt.y0, tt.xBase, tt.yBase, tt.log2TrafoSize, tt.trafoDepth, tt.blkIdx));
    }
}


template <class H>
void Syntax<mvd_coding>::go(const mvd_coding &e, H &h)
{
    h(abs_mvd_greater0_flag(0), ae(v));
    h(abs_mvd_greater0_flag(1), ae(v));

    if (h[abs_mvd_greater0_flag(0)])
        h(abs_mvd_greater1_flag(0), ae(v));
    if (h[abs_mvd_greater0_flag(1)])
        h(abs_mvd_greater1_flag(1), ae(v));

    if (h[abs_mvd_greater0_flag(0)])
    {
        if (h[abs_mvd_greater1_flag(0)])
            h(abs_mvd_minus2(0), ae(v));
        h(mvd_sign_flag(0), ae(v));
    }
    if (h[abs_mvd_greater0_flag(1)])
    {
        if (h[abs_mvd_greater1_flag(1)])
            h(abs_mvd_minus2(1), ae(v));
        h(mvd_sign_flag(1), ae(v));
    }
}


template <class> struct Decode;


template <class H>
void Syntax<transform_unit>::go(const transform_unit &tu, H &h)
{
    const int log2TrafoSizeC = std::max(2, tu.log2TrafoSize - (h[ChromaArrayType()] == 3 ? 0 : 1));
    const int cbfDepthC = tu.trafoDepth - (h[ChromaArrayType()] != 3 && tu.log2TrafoSize == 2 ? 1 : 0);

    const int xC = (h[ChromaArrayType()] != 3 && tu.log2TrafoSize == 2) ? tu.xBase : tu.x0;
    const int yC = (h[ChromaArrayType()] != 3 && tu.log2TrafoSize == 2) ? tu.yBase : tu.y0;

    const int cbfLuma = h[cbf_luma(tu.x0, tu.y0, tu.trafoDepth)];
    const bool cbfChroma =
            h[cbf_cb(xC, yC, cbfDepthC)] ||
            h[cbf_cr(xC, yC, cbfDepthC)] ||
            (h[ChromaArrayType()] == 2 &&
                    (h[cbf_cb(xC, yC + (1 << log2TrafoSizeC), cbfDepthC)] ||
                            h[cbf_cr(xC, yC + (1 << log2TrafoSizeC), cbfDepthC)]));

    if ((cbfLuma || cbfChroma) && h[cu_qp_delta_enabled_flag()] && !h[IsCuQpDeltaCoded()])
    {
        h(cu_qp_delta_abs(), ae(v));
        if (h[cu_qp_delta_abs()])
            h(cu_qp_delta_sign_flag(), ae(v));
        // review: move these specifics out of the syntax function
        if (!std::is_same<typename H::Tag, Decode<void>>::value && static_cast<QpState *>(h)->getCanWrite())
        {
            int rowQgModulo = (tu.yBase & (h[CtbSizeY()] - 1)) >> 3;
            int colQgModulo = (tu.xBase & (h[CtbSizeY()] - 1)) >> 3;
            int qgSize = std::max<int>(1, 1 << ((tu.log2TrafoSize + tu.trafoDepth - 3) << 1));
            int qpValue = static_cast<QpState *>(h)->getQp(0) - h[QpBdOffsetY()];
            static_cast<QpState *>(h)->setQpInternal(rowQgModulo, colQgModulo, qgSize, qpValue);
        }
    }

    if (h[cu_chroma_qp_offset_enabled_flag()] && cbfChroma  &&
            !h[cu_transquant_bypass_flag()] && !h[IsCuChromaQpOffsetCoded()])
    {
        h(cu_chroma_qp_offset_flag(), ae(v));
        if (h[cu_chroma_qp_offset_flag()]  &&  h[chroma_qp_offset_list_len_minus1()] > 0)
            h(cu_chroma_qp_offset_idx(), ae(v));
    }

    h(IfCbf<cbf_luma, residual_coding>{
        cbf_luma(tu.x0, tu.y0, tu.trafoDepth),
                residual_coding(tu.x0, tu.y0, tu.log2TrafoSize, 0) });

    if (tu.log2TrafoSize > 2 || h[ChromaArrayType()] == 3)
    {
        if (h[cross_component_prediction_enabled_flag()] && cbfLuma &&
                (h[CuPredMode(tu.x0, tu.y0)] == MODE_INTER ||
                        h[intra_chroma_pred_mode(tu.x0, tu.y0)] == 4))
        {
            h(cross_comp_pred(tu.x0, tu.y0, 0));
        }

        for (int tIdx = 0; tIdx < (h[ChromaArrayType()] == 2 ? 2 : 1); tIdx++)
        {
            h(IfCbf<cbf_cb, residual_coding>{
                cbf_cb(tu.x0, tu.y0 + (tIdx << log2TrafoSizeC), tu.trafoDepth),
                        residual_coding(tu.x0, tu.y0 + (tIdx << log2TrafoSizeC), log2TrafoSizeC, 1) });
        }

        if (h[cross_component_prediction_enabled_flag()] && cbfLuma &&
                (h[CuPredMode(tu.x0, tu.y0)] == MODE_INTER ||
                        h[intra_chroma_pred_mode(tu.x0, tu.y0)] == 4))
        {
            h(cross_comp_pred(tu.x0, tu.y0, 1));
        }

        for (int tIdx = 0; tIdx < (h[ChromaArrayType()] == 2 ? 2 : 1); tIdx++)
        {
            h(IfCbf<cbf_cr, residual_coding>{
                cbf_cr(tu.x0, tu.y0 + (tIdx << log2TrafoSizeC), tu.trafoDepth),
                        residual_coding(tu.x0, tu.y0 + (tIdx << log2TrafoSizeC), log2TrafoSizeC, 2) });
        }
    }
    else if (tu.blkIdx == 3)
    {
        for (int tIdx = 0; tIdx < (h[ChromaArrayType()] == 2 ? 2 : 1); tIdx++)
        {
            h(IfCbf<cbf_cb, residual_coding>{
                cbf_cb(tu.xBase, tu.yBase + (tIdx << log2TrafoSizeC), tu.trafoDepth - 1),
                        residual_coding(tu.xBase, tu.yBase + (tIdx << log2TrafoSizeC), log2TrafoSizeC, 1) });
        }

        for (int tIdx = 0; tIdx < (h[ChromaArrayType()] == 2 ? 2 : 1); tIdx++)
        {
            h(IfCbf<cbf_cr, residual_coding>{
                cbf_cr(tu.xBase, tu.yBase + (tIdx << log2TrafoSizeC), tu.trafoDepth - 1),
                        residual_coding(tu.xBase, tu.yBase + (tIdx << log2TrafoSizeC), log2TrafoSizeC, 2) });
        }
    }
}


template <class H>
bool predModeIntraIs10or26(residual_coding e, H &h)
{
    const int predModeIntra = e.cIdx ? h[IntraPredModeC(e.x0, e.y0)] : h[IntraPredModeY(e.x0, e.y0)];
    return predModeIntra == 10 || predModeIntra == 26;
}


template <class H>
void Syntax<residual_coding>::go(const residual_coding &rc, H &h)
{
    if (h[transform_skip_enabled_flag()] && !h[cu_transquant_bypass_flag()] && (rc.log2TrafoSize <= h[Log2MaxTransformSkipSize()]))
        h(transform_skip_flag(rc.x0, rc.y0, rc.cIdx), ae(v));

    if (h[current(CuPredMode(rc.x0, rc.y0))] == MODE_INTER  &&  h[explicit_rdpcm_enabled_flag()] &&
            (h[transform_skip_flag(rc.x0, rc.y0, rc.cIdx)] || h[cu_transquant_bypass_flag()]))
    {
        h(explicit_rdpcm_flag(rc.x0, rc.y0, rc.cIdx), ae(v));
        if (h[explicit_rdpcm_flag(rc.x0, rc.y0, rc.cIdx)])
            h(explicit_rdpcm_dir_flag(rc.x0, rc.y0, rc.cIdx), ae(v));
    }

    h(last_sig_coeff_x_prefix(), ae(v));
    h(last_sig_coeff_y_prefix(), ae(v));
    if (h[last_sig_coeff_x_prefix()] > 3)
        h(last_sig_coeff_x_suffix(), ae(v));
    if (h[last_sig_coeff_y_prefix()] > 3)
        h(last_sig_coeff_y_suffix(), ae(v));

    int xC, yC;
    int lastScanPos = 16;
    int lastSubBlock = (1 << (rc.log2TrafoSize - 2)) * (1 << (rc.log2TrafoSize - 2)) - 1;
    int escapeDataPresent = 0;

    do
    {
        if (lastScanPos == 0)
        {
            lastScanPos = 16;
            lastSubBlock--;
        }
        lastScanPos--;
        const int xS = ScanOrder(rc.log2TrafoSize - 2, h[scanIdx()], lastSubBlock, 0);
        const int yS = ScanOrder(rc.log2TrafoSize - 2, h[scanIdx()], lastSubBlock, 1);
        xC = (xS << 2) + ScanOrder(2, h[scanIdx()], lastScanPos, 0);
        yC = (yS << 2) + ScanOrder(2, h[scanIdx()], lastScanPos, 1);
    } while ((xC != h[LastSignificantCoeffX()]) || (yC != h[LastSignificantCoeffY()]));

    StateSubstream &stateSubstream = *static_cast<StateSubstream *>(h);
    int &i = stateSubstream.i;
    for (i = lastSubBlock; i >= 0; i--)
    {
        const int xS = ScanOrder(rc.log2TrafoSize - 2, h[scanIdx()], i, 0);
        const int yS = ScanOrder(rc.log2TrafoSize - 2, h[scanIdx()], i, 1);

        int inferSbDcSigCoeffFlag = 0;
        if ((i < lastSubBlock) && (i > 0))
        {
            h(coded_sub_block_flag(xS, yS), ae(v));
            inferSbDcSigCoeffFlag = 1;
        }
        for (int n = (i == lastSubBlock) ? lastScanPos - 1 : 15; n >= 0; n--)
        {
            const int xC = (xS << 2) + ScanOrder(2, h[scanIdx()], n, 0);
            const int yC = (yS << 2) + ScanOrder(2, h[scanIdx()], n, 1);
            if (h[coded_sub_block_flag(xS, yS)] && (n > 0 || !inferSbDcSigCoeffFlag))
            {
                h(sig_coeff_flag(xC, yC), ae(v));
                if (h[sig_coeff_flag(xC, yC)])
                    inferSbDcSigCoeffFlag = 0;
            }
        }

        int firstSigScanPos = 16;
        int lastSigScanPos = -1;
        h[numGreater1Flag()] = 0;
        int lastGreater1ScanPos = -1;
        for (int n = 15; n >= 0; n--)
        {
            const int xC = (xS << 2) + ScanOrder(2, h[scanIdx()], n, 0);
            const int yC = (yS << 2) + ScanOrder(2, h[scanIdx()], n, 1);
            if (h[sig_coeff_flag(xC, yC)])
            {
                if (h[numGreater1Flag()] < 8)
                {
                    h(coeff_abs_level_greater1_flag(n), ae(v));
                    h[numGreater1Flag()]++;

                    if (h[coeff_abs_level_greater1_flag(n)] && lastGreater1ScanPos == -1)
                        lastGreater1ScanPos = n;
                    else if (h[coeff_abs_level_greater1_flag(n)])
                        escapeDataPresent = 1;
                }
                else
                {
                    escapeDataPresent = 1;
                }

                if (lastSigScanPos == -1)
                    lastSigScanPos = n;

                firstSigScanPos = n;
            }
        }

        int signHidden;

        if (h[cu_transquant_bypass_flag()] ||
                (h[current(CuPredMode(rc.x0, rc.y0))] == MODE_INTRA  &&
                        h[implicit_rdpcm_enabled_flag()] && h[transform_skip_flag(rc.x0, rc.y0, rc.cIdx)] &&
                        predModeIntraIs10or26(rc, h)) ||
                        h[explicit_rdpcm_flag(rc.x0, rc.y0, rc.cIdx)])
        {
            signHidden = 0;
        }
        else
        {
            signHidden = lastSigScanPos - firstSigScanPos > 3;
        }

        if (lastGreater1ScanPos != -1)
        {
            h(coeff_abs_level_greater2_flag(lastGreater1ScanPos), ae(v));
            if (h[coeff_abs_level_greater2_flag(lastGreater1ScanPos)])
                escapeDataPresent = 1;
        }

        if (h[cabac_bypass_alignment_enabled_flag()] && escapeDataPresent)
        {
            assert(!"alignment not yet implemented");
        }

        for (int n = 15; n >= 0; n--)
        {
            const int xC = (xS << 2) + ScanOrder(2, h[scanIdx()], n, 0);
            const int yC = (yS << 2) + ScanOrder(2, h[scanIdx()], n, 1);
            if (h[sig_coeff_flag(xC, yC)] && (!h[sign_data_hiding_enabled_flag()] || !signHidden || (n != firstSigScanPos)))
                h(coeff_sign_flag(n), ae(v));
        }

        int numSigCoeff = 0;
        int sumAbsLevel = 0;
        for (int n = 15; n >= 0; n--)
        {
            const int xC = (xS << 2) + ScanOrder(2, h[scanIdx()], n, 0);
            const int yC = (yS << 2) + ScanOrder(2, h[scanIdx()], n, 1);
            if (h[sig_coeff_flag(xC, yC)])
            {
                const int baseLevel = 1 + h[coeff_abs_level_greater1_flag(n)] + h[coeff_abs_level_greater2_flag(n)];
                if (baseLevel == ((numSigCoeff < 8) ? ((n == lastGreater1ScanPos) ? 3 : 2) : 1))
                    h(coeff_abs_level_remaining(n), ae(v));
                h[TransCoeffLevel(rc.x0, rc.y0, rc.cIdx, xC, yC)] =
                        (h[coeff_abs_level_remaining(n)] + baseLevel) * (1 - 2 * h[coeff_sign_flag(n)]);
                if (h[sign_data_hiding_enabled_flag()] && signHidden)
                {
                    sumAbsLevel += (h[coeff_abs_level_remaining(n)] + baseLevel);
                    if ((n == firstSigScanPos) && ((sumAbsLevel % 2) == 1))
                        h[TransCoeffLevel(rc.x0, rc.y0, rc.cIdx, xC, yC)] = -h[TransCoeffLevel(rc.x0, rc.y0, rc.cIdx, xC, yC)];
                }

                numSigCoeff++;
            }
        }
    }
}


template <>
struct Syntax<cross_comp_pred>
{
    template <class H> static void go(const cross_comp_pred &ccp, H &h)
    {
        h(log2_res_scale_abs_plus1(ccp.c), ae(v));
        if (h[log2_res_scale_abs_plus1(ccp.c)] != 0)
            h(res_scale_sign_flag(ccp.c), ae(v));
    }
};


// The following are not true HEVC constructs and exist purely for convenience of implementation.

template <>
struct Syntax<IntraPartitionPrediction>
{
    template <class H> static void go(const IntraPartitionPrediction &e, H &h)
    {
        const int i = (e.blkIdx & 1) << (e.log2CbSize - 1);
        const int j = (e.blkIdx >> 1) << (e.log2CbSize - 1);

        h(prev_intra_luma_pred_flag(e.x0 + i, e.y0 + j), ae(v));

        if (h[prev_intra_luma_pred_flag(e.x0 + i, e.y0 + j)])
            h(mpm_idx(e.x0 + i, e.y0 + j), ae(v));
        else
            h(rem_intra_luma_pred_mode(e.x0 + i, e.y0 + j), ae(v));
    }
};


template <>
struct Syntax<IntraPartition>
{
    template <class H> static void go(const IntraPartition &e, H &h)
    {
        h(IntraPartitionPrediction(e.x0, e.y0, e.log2CbSize, e.blkIdx));

        const int i = (e.blkIdx & 1) << (e.log2CbSize - 1);
        const int j = (e.blkIdx >> 1) << (e.log2CbSize - 1);

        const int IntraSplitFlag = h[::IntraSplitFlag()];

        assert(e.split == IntraSplitFlag);

        h[MaxTrafoDepth()] = h[max_transform_hierarchy_depth_intra()] + IntraSplitFlag;

        h(transform_tree(e.x0 + i, e.y0 + j, e.x0, e.y0, e.log2CbSize - e.split, IntraSplitFlag, e.blkIdx));
    }
};


#endif
