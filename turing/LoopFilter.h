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

#ifndef INCLUDED_LoopFilter_h
#define INCLUDED_LoopFilter_h

#pragma once

// Implementation of HEVC deblocking and SAO filters


#include "Global.h"
#include "HevcTypes.h"
#include "StateSpatial.h"
#include "sao.h"
#include "havoc/pred_inter.h"
#include "QpState.h"
#include <stdint.h>
#include <vector>
#include <vector>
#include <string.h>


static const int EDGE_VER = 0; // vertical edge
static const int EDGE_HOR = 1; // horizontal edge

// review: these should not be necessary, LoopFilter should be generic
template <class> struct Write;
template <class> struct Decode;

namespace LoopFilter
{
    // Information about an 8x8 picture region
    struct Block
    {
        int8_t data;
        uint8_t packedBs;

        template <class H>
        static int8_t packData(H &h, coding_unit cu, int qpy)
        {
            const bool disable =
                    (h[pcm_loop_filter_disabled_flag()] && h[current(pcm_flag(cu.x0, cu.y0))])
                    || h[cu_transquant_bypass_flag()];

            return  (qpy << 1) | (disable ? 1 : 0);
        }

        int8_t getQp() const
        {
            return this->data >> 1;
        }
        bool enabled() const
        {
            return !(this->data & 1);
        }
        int getBs(int edgeType, int position) const
        {
            assert(edgeType < 2 && position < 2);
            return 0x3 & (this->packedBs >> (4 * edgeType + 2 * position));
        }
        void setBs(int edgeType, int position, int bS)
        {
            assert(edgeType < 2 && position < 2 && bS < 4);
            this->packedBs &= ~(0x3 << (4 * edgeType + 2 * position));
            this->packedBs |= (bS << (4 * edgeType + 2 * position));
        }
        void increaseBs(int edgeType, int position, int minBs)
        {
            assert(edgeType < 2 && position < 2 && minBs < 4);
            if (minBs > this->getBs(edgeType, position)) this->setBs(edgeType, position, minBs);
        }
    };

    struct Ctu
    {
        int tc_offset_div2;
        int beta_offset_div2;
        int top, bottom, left, right;
        bool topLeft, bottomRight;
        bool topRight, bottomLeft;

        struct Plane
        {
            int SaoTypeIdx;
            union Union
            {
                int eoClass;
                int saoLeftClass;
            };
            Union u;
            int16_t SaoOffsetVal[5];
        };

        Plane planes[3];

        template <class H>
        void set(H &h)
        {
            this->tc_offset_div2 = h[slice_tc_offset_div2()];
            this->beta_offset_div2 = h[slice_beta_offset_div2()];

            const int rx = h[xCtb()] >> h[CtbLog2SizeY()];
            const int ry = h[yCtb()] >> h[CtbLog2SizeY()];

            for (int cIdx = 0; cIdx < 3; ++cIdx)
            {
                if (cIdx ? !!(h[slice_sao_chroma_flag()]) : !!(h[slice_sao_luma_flag()]))
                {
                    this->planes[cIdx].SaoTypeIdx = h[SaoTypeIdx(cIdx, rx, ry)];
                }
                else
                {
                    this->planes[cIdx].SaoTypeIdx = 0;
                }

                int offsetSign[4];

                if (this->planes[cIdx].SaoTypeIdx == 2)
                {
                    offsetSign[2] = -1;
                    offsetSign[3] = -1;
                    offsetSign[0] = 1;
                    offsetSign[1] = 1;

                    this->planes[cIdx].u.eoClass = h[SaoEoClass(cIdx, rx, ry)];
                }
                else if (this->planes[cIdx].SaoTypeIdx == 1)
                {
                    for (int i = 0; i < 4; ++i)
                    {
                        offsetSign[i] = h[sao_offset_sign(cIdx, rx, ry, i)] ? -1 : 1;
                    }

                    this->planes[cIdx].u.saoLeftClass = h[sao_band_position(cIdx, rx, ry)];
                }

                for (int i = 0; i < 4; ++i)
                {
                    const int bitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];

                    this->planes[cIdx].SaoOffsetVal[i + 1] = offsetSign[i] * h[sao_offset_abs(cIdx, rx, ry, i)] << (bitDepth - std::min(bitDepth, 10));
                }
            }
        }
    };

    template <typename Sample, int edgeType>
    struct BlockEdge;

    template <typename Sample>
    struct BlockEdge<Sample, EDGE_VER>
    {

        BlockEdge(Raster<Sample> s, Raster<Block> blocks, int bitDepth, int position = 0) :
        blockP(blocks(-1, 0)),
        blockQ(blocks(0, 0)),
        bitDepth(bitDepth),
        s(s.offset(0, 4*position))
        {
        }

        Sample &p(int i, int k)
        {
            return this->s(-i - 1, k);
        }
        Sample &q(int i, int k)
        {
            return this->s(i, k);
        }
        Block const &blockP, &blockQ;
        int const bitDepth;
        Raster<Sample> s;
    };

    template <typename Sample>
    struct BlockEdge<Sample, EDGE_HOR>
    {
        BlockEdge(Raster<Sample> s, Raster<Block> blocks, int bitDepth, int position = 0) :
        blockP(blocks(0, -1)),
        blockQ(blocks(0, 0)),
        bitDepth(bitDepth),
        s(s.offset(4*position, 0))
        {
        }

        Sample &p(int i, int k)
        {
            return this->s(k, -i - 1);
        }
        Sample &q(int i, int k)
        {
            return this->s(k, i);
        }
        Block const &blockP, &blockQ;
        int const bitDepth;
        Raster<Sample> s;
    };

    static int betaTable(int Q)
    {
        static const int lookup[52] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64 };
        return lookup[Q];
    }

    static int tCTable(int Q)
    {
        static const int lookup[54] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10, 11, 13, 14, 16, 18, 20, 22, 24 };
        return lookup[Q];
    }

    template <typename Sample, int edgeType>
    struct LumaBlockEdge :
        BlockEdge<Sample, edgeType>
        {
            template <class H>
            LumaBlockEdge(H &h, Raster<Sample> recSamplesL, Raster<Block> blocks, const Ctu &ctu, int position) :
            BlockEdge<Sample, edgeType>(recSamplesL, blocks, h[BitDepthY()], position)
            {
                this->dE = 0;
                this->dEp = 0;
                this->dEq = 0;

                const int bS = this->blockQ.getBs(edgeType, position);

                if (bS != 0)
                {
                    const auto QpP = this->blockP.getQp();
                    const auto QpQ = this->blockQ.getQp();

                    const auto qPL = (QpQ + QpP + 1) >> 1;

                    {
                        const int Q = Clip3(0, 51, qPL + (ctu.beta_offset_div2 << 1));

                        beta = betaTable(Q);
                        beta *= (1 << (h[BitDepthY()] - 8));
                    }

                    {
                        const int Q = Clip3(0, 53, qPL + 2 * (bS - 1) + (ctu.tc_offset_div2 << 1));

                        tC = tCTable(Q);
                        tC *= (1 << (h[BitDepthY()] - 8));
                    }

                    const int dp0 = std::abs(this->p(2, 0) - 2 * this->p(1, 0) + this->p(0, 0));
                    const int dp3 = std::abs(this->p(2, 3) - 2 * this->p(1, 3) + this->p(0, 3));
                    const int dq0 = std::abs(this->q(2, 0) - 2 * this->q(1, 0) + this->q(0, 0));
                    const int dq3 = std::abs(this->q(2, 3) - 2 * this->q(1, 3) + this->q(0, 3));

                    const int dpq0 = dp0 + dq0;
                    const int dpq3 = dp3 + dq3;
                    const int dp = dp0 + dp3;
                    const int dq = dq0 + dq3;
                    const int d = dpq0 + dpq3;

                    if (d < beta)
                    {
                        const int dSam0 = dSam(this->p(0, 0), this->p(3, 0), this->q(0, 0), this->q(3, 0), 2 * dpq0, beta, tC);
                        const int dSam3 = dSam(this->p(0, 3), this->p(3, 3), this->q(0, 3), this->q(3, 3), 2 * dpq3, beta, tC);

                        this->dE = 1;

                        if (dSam0 == 1 && dSam3 == 1) this->dE = 2;

                        if (dp < ((beta + (beta >> 1)) >> 3)) this->dEp = 1;
                        if (dq < ((beta + (beta >> 1)) >> 3)) this->dEq = 1;
                    }
                }
            }

            static int dSam(Sample p0, Sample p3, Sample q0, Sample q3, int dpq, int beta, int tC)
            {
                if (dpq < (beta >> 2)
                        && std::abs(p3 - p0) + std::abs(q0 - q3) < (beta >> 3)
                        && std::abs(p0 - q0) < ((5 * tC + 1) >> 1))
                {
                    return 1;
                }
                else
                {
                    return 0;
                }
            }

            void filter()
            {
                for (int k = 0; k < 4; ++k)
                {
                    const int p0 = this->p(0, k);
                    const int p1 = this->p(1, k);
                    const int p2 = this->p(2, k);
                    const int p3 = this->p(3, k);
                    const int q0 = this->q(0, k);
                    const int q1 = this->q(1, k);
                    const int q2 = this->q(2, k);
                    const int q3 = this->q(3, k);

                    if (dE == 2)
                    {
                        // strong filtering
                        if (this->blockP.enabled())
                        {
                            this->p(0, k) = Clip3(p0 - 2 * tC, p0 + 2 * tC, (p2 + 2 * p1 + 2 * p0 + 2 * q0 + q1 + 4) >> 3);
                            this->p(1, k) = Clip3(p1 - 2 * tC, p1 + 2 * tC, (p2 + p1 + p0 + q0 + 2) >> 2);
                            this->p(2, k) = Clip3(p2 - 2 * tC, p2 + 2 * tC, (2 * p3 + 3 * p2 + p1 + p0 + q0 + 4) >> 3);
                        }
                        if (this->blockQ.enabled())
                        {
                            this->q(0, k) = Clip3(q0 - 2 * tC, q0 + 2 * tC, (p1 + 2 * p0 + 2 * q0 + 2 * q1 + q2 + 4) >> 3);
                            this->q(1, k) = Clip3(q1 - 2 * tC, q1 + 2 * tC, (p0 + q0 + q1 + q2 + 2) >> 2);
                            this->q(2, k) = Clip3(q2 - 2 * tC, q2 + 2 * tC, (p0 + q0 + q1 + 3 * q2 + 2 * q3 + 4) >> 3);
                        }
                    }
                    else if (dE == 1)
                    {
                        int delta = (9 * (q0 - p0) - 3 * (q1 - p1) + 8) >> 4;
                        if (std::abs(delta) < tC * 10)
                        {
                            delta = Clip3(-tC, tC, delta);
                            if (this->blockP.enabled()) this->p(0, k) = Clip1Y(p0 + delta, this->bitDepth);
                            if (this->blockQ.enabled()) this->q(0, k) = Clip1Y(q0 - delta, this->bitDepth);
                            if (dEp == 1 && this->blockP.enabled())
                            {
                                int deltap = Clip3(-(tC >> 1), tC >> 1, (((p2 + p0 + 1) >> 1) - p1 + delta) >> 1);
                                this->p(1, k) = Clip1Y(p1 + deltap, this->bitDepth);
                            }
                            if (dEq == 1 && this->blockQ.enabled())
                            {
                                int deltaq = Clip3(-(tC >> 1), tC >> 1, (((q2 + q0 + 1) >> 1) - q1 - delta) >> 1);
                                this->q(1, k) = Clip1Y(q1 + deltaq, this->bitDepth);
                            }
                        }
                    }
                }
            }
            int dE, dEp, dEq;
            int beta, tC;
        };

    template <typename Sample, int edgeType>
    struct ChromaBlockEdge :
        BlockEdge<Sample, edgeType>
        {
            template <class H>
            ChromaBlockEdge(H &h, Raster<Sample> s, Raster<Block> blocks, Ctu ctu, int cQpPicOffset, int position) :
            BlockEdge<Sample, edgeType>(s, blocks, h[BitDepthC()], position)
            {
                if (this->blockQ.getBs(edgeType, 0) == 2)
                {
                    const int QpP = this->blockP.getQp();
                    const int QpQ = this->blockQ.getQp();

                    const int qPi = ((QpQ + QpP + 1) >> 1) + cQpPicOffset;
                    const int QpC = ::QpC(qPi);

                    tC = tCTable(Clip3(0, 53, QpC + 2 + (ctu.tc_offset_div2 << 1)));
                    tC *= (1 << (this->bitDepth - 8));
                }
            }

            void filter()
            {
                if (this->blockQ.getBs(edgeType, 0) == 2)
                {
                    for (int k = 0; k < 4; ++k)
                    {
                        const int p0 = this->p(0, k);
                        const int p1 = this->p(1, k);
                        const int q0 = this->q(0, k);
                        const int q1 = this->q(1, k);

                        const int delta = Clip3(-tC, tC, ((((q0 - p0) << 2) + p1 - q1 + 4) >> 3));

                        if (this->blockP.enabled()) this->p(0, k) = Clip1C(p0 + delta, this->bitDepth);
                        if (this->blockQ.enabled()) this->q(0, k) = Clip1C(q0 - delta, this->bitDepth);
                    }
                }
            }

            int tC;
        };

    static bool sameMotion(const MotionVector &mvA, const int8_t &dpbIndexA, const MotionVector &mvB, const int8_t &dpbIndexB)
    {
        if (dpbIndexA != dpbIndexB) return false;
        if (dpbIndexA < 0) return true;
        if (std::abs(mvA[0] - mvB[0]) >= 4) return false;
        if (std::abs(mvA[1] - mvB[1]) >= 4) return false;
        return true;
    }

    static bool sameMotion(const PuData &puDataA, const PuData &puDataB)
    {
        const bool sameL0AL0B = sameMotion(puDataA.mv(L0), puDataA.getDpbIndex(L0), puDataB.mv(L0), puDataB.getDpbIndex(L0));
        const bool sameL1AL1B = sameMotion(puDataA.mv(L1), puDataA.getDpbIndex(L1), puDataB.mv(L1), puDataB.getDpbIndex(L1));
        if (sameL0AL0B && sameL1AL1B) return true;

        const bool sameL0AL1B = sameMotion(puDataA.mv(L0), puDataA.getDpbIndex(L0), puDataB.mv(L1), puDataB.getDpbIndex(L1));
        const bool sameL1AL0B = sameMotion(puDataA.mv(L1), puDataA.getDpbIndex(L1), puDataB.mv(L0), puDataB.getDpbIndex(L0));
        if (sameL0AL1B && sameL1AL0B) return true;

        return false;
    }

    // Review: rename process*() functions - "process" is perhaps misleading: we don't apply filtering in these functions, just capture information for use in later filtering?
    struct Picture
    {
        std::vector<Block> blocks;
        std::vector<Ctu> ctus;
        std::intptr_t stride;

        Block &blockAt(int x, int y)
        {
            return this->blocks[this->stride * y + x];
        }

        template <class H>
        Picture(H &h) :
        stride((h[PicWidthInCtbsY()] << h[CtbLog2SizeY()] >> 3) + 1)
        {
            const int height = (h[PicHeightInCtbsY()] << h[CtbLog2SizeY()] >> 3) + 1;
            this->blocks.resize(this->stride * height);
            this->ctus.resize(h[PicSizeInCtbsY()]);
        }

        template <class H>
        bool neighbourCtuAvailable(H &h, int dx, int dy)
        {
            const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
            const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];

            if (rx + dx < 0) return false;
            if (rx + dx >= h[PicWidthInCtbsY()]) return false;
            if (ry + dy < 0) return false;
            if (ry + dy >= h[PicHeightInCtbsY()]) return false;

            const int neighbourCtbAddrInRs = h[CtbAddrInRs()] + dx + dy * h[PicWidthInCtbsY()];
            const int neighbourCtbAddrInTs = h[CtbAddrRsToTs(neighbourCtbAddrInRs)];

            if (h[TileId(h[CtbAddrInTs()])] != h[TileId(neighbourCtbAddrInTs)] && !h[loop_filter_across_tiles_enabled_flag()])
            {
                return false;
            }

            const int sliceAddrTs = h[CtbAddrRsToTs(h[SliceAddrRs()])];

            if (neighbourCtbAddrInTs < sliceAddrTs && !h[slice_loop_filter_across_slices_enabled_flag()])
            {
                return false;
            }

            return true;
        }

        template <class H>
        void processCtu(H &h, coding_tree_unit e)
        {
            Ctu &ctu = this->ctus[h[CtbAddrInRs()]];
            ctu.set(h);

            ctu.left = 0;
            ctu.top = 0;
            ctu.right = h[pic_width_in_luma_samples()];
            ctu.bottom = h[pic_height_in_luma_samples()];

            const bool filterLeftEdge = neighbourCtuAvailable(h, -1, 0);
            if (!filterLeftEdge)
            {
                for (int y = h[yCtb()] / 8; y < (h[yCtb()] + h[CtbSizeY()]) / 8; ++y)
                {
                    this->blockAt(h[xCtb()] / 8, y).setBs(EDGE_VER, 0, 0);
                    this->blockAt(h[xCtb()] / 8, y).setBs(EDGE_VER, 1, 0);
                }
                ctu.left = h[xCtb()];
                if (h[xCtb()] > 0)
                {
                    Ctu &ctuLeft = this->ctus[h[CtbAddrInRs()] - 1];
                    ctuLeft.right = h[xCtb()];
                }
            }

            const bool filterTopEdge = neighbourCtuAvailable(h, 0, -1);
            if (!filterTopEdge)
            {
                for (int x = h[xCtb()] / 8; x < (h[xCtb()] + h[CtbSizeY()]) / 8; ++x)
                {
                    this->blockAt(x, h[yCtb()] / 8).setBs(EDGE_HOR, 0, 0);
                    this->blockAt(x, h[yCtb()] / 8).setBs(EDGE_HOR, 1, 0);
                }
                ctu.top = h[yCtb()];
                if (h[yCtb()] > 0)
                {
                    Ctu &ctuAbove = this->ctus[h[CtbAddrInRs()] - h[PicWidthInCtbsY()]];
                    ctuAbove.bottom = h[yCtb()];
                }
            }

            const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
            const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];

            ctu.topLeft = neighbourCtuAvailable(h, -1, -1);
            if (rx && ry)
            {
                Ctu &ctuTopLeft = this->ctus[h[CtbAddrInRs()] - 1 - h[PicWidthInCtbsY()]];
                ctuTopLeft.bottomRight = ctu.topLeft;
            }

            ctu.topRight = neighbourCtuAvailable(h, 1, -1);
            if (rx < h[PicWidthInCtbsY()] - 1 && ry)
            {
                Ctu &ctuTopRight = this->ctus[h[CtbAddrInRs()] + 1 - h[PicWidthInCtbsY()]];
                ctuTopRight.bottomLeft = ctu.topRight;
            }

            ctu.bottomRight = false;
            ctu.bottomLeft = false;
        }

        template <class H>
        void processCu(H &h, coding_unit cu)
        {
            if (h[current(pcm_flag(cu.x0, cu.y0))])
            {
                const int bS = 2;
                const int sizeCu = 1 << cu.log2CbSize;

                if (cu.x0 % 8 == 0)
                {
                    // left edge of PCM coding_unit
                    for (int y = (cu.y0 + 3) / 4; y < (cu.y0 + sizeCu) / 4; ++y)
                    {
                        const int x = cu.x0 / 8;
                        this->blockAt(x, y / 2).increaseBs(EDGE_VER, y % 2, bS);
                    }
                }
                if ((cu.x0 + sizeCu) % 8 == 0)
                {
                    // right edge of PCM coding_unit
                    for (int y = (cu.y0 + 3) / 4; y < (cu.y0 + sizeCu) / 4; ++y)
                    {
                        const int x = (cu.x0 + sizeCu) / 8;
                        this->blockAt(x, y / 2).increaseBs(EDGE_VER, y % 2, bS);
                    }
                }
                if (cu.y0 % 8 == 0)
                {
                    // top edge of PCM coding_unit
                    for (int x = (cu.x0 + 3) / 4; x < (cu.x0 + sizeCu) / 4; ++x)
                    {
                        const int y = cu.y0 / 8;
                        this->blockAt(x / 2, y).increaseBs(EDGE_HOR, x % 2, bS);
                    }
                }
                if ((cu.y0 + sizeCu) % 8 == 0)
                {
                    // bottom edge of PCM coding_unit
                    for (int x = (cu.x0 + 3) / 4; x < (cu.x0 + sizeCu) / 4; ++x)
                    {
                        const int y = (cu.y0 + sizeCu) / 8;
                        this->blockAt(x / 2, y).increaseBs(EDGE_HOR, x % 2, bS);
                    }
                }
            }

            int8_t data;

            const int size = (1 << cu.log2CbSize);
            QpState *qpState = h;

            for (int y = cu.y0 / 8; y < (cu.y0 + size) / 8; ++y)
            {
                for (int x = cu.x0 / 8; x < (cu.x0 + size) / 8; ++x)
                {
                    int qpy = h[QpY()];

                    if(!std::is_same<typename H::Tag, Decode<void>>::value && h[cu_qp_delta_enabled_flag()])
                    {
                        int rowModulo = ((y<<3) & qpState->getMaskCtb()) >> 3;
                        int colModulo = ((x<<3) & qpState->getMaskCtb()) >> 3;
                        qpy = qpState->getQpInternal(rowModulo, colModulo);
                    }

                    data = Block::packData(h, cu, qpy);
                    this->blockAt(x, y).data = data;
                }
            }
        }

        template <class H>
        void processPu(H &h, const prediction_unit &pu)
        {
            static const int bS = 1;

            if (!h[slice_deblocking_filter_disabled_flag()])
            {
                Snake<BlockData>::Cursor *cursor = h;
                const PuData &puData = cursor->current(0, 0, h[MinCbLog2SizeY()] - 1);

                if (pu.x0 % 8 == 0)
                {
                    // left edge of prediction unit
                    for (int y = pu.y0; y < pu.y0 + pu.nPbH; y += 4)
                    {
                        const PuData &puDataLeft = cursor->offset<Left>(-1, y - pu.y0, h[MinCbLog2SizeY()] - 1);
                        if (!sameMotion(puDataLeft, puData))
                        {
                            this->blockAt(pu.x0 / 8, y / 8).increaseBs(EDGE_VER, (y / 4) % 2, bS);
                        }
                    }
                }
                if (pu.y0 % 8 == 0)
                {
                    // top edge of prediction unit
                    for (int x = pu.x0; x < pu.x0 + pu.nPbW; x += 4)
                    {
                        const PuData &puDataAbove = cursor->offset<Up>(x - pu.x0, -1, h[MinCbLog2SizeY()] - 1);
                        if (!sameMotion(puDataAbove, puData))
                        {
                            this->blockAt(x / 8, pu.y0 / 8).increaseBs(EDGE_HOR, (x / 4) % 2, bS);
                        }
                    }
                }
            }
        }

        template <class H>
        void processTu(H &h, transform_unit tu)
        {
            if (h[current(CuPredMode(tu.x0, tu.y0))] == MODE_INTRA && !h[slice_deblocking_filter_disabled_flag()])
            {
                const int sizeTu = 1 << tu.log2TrafoSize;
                const int bS = 2;
                if (tu.x0 % 8 == 0)
                {
                    // left edge of luma transform block coded with intra prediction mode
                    for (int y = (tu.y0 + 3) / 4; y < (tu.y0 + sizeTu) / 4; ++y)
                    {
                        const int x = tu.x0 / 8;
                        this->blockAt(x, y / 2).increaseBs(EDGE_VER, y % 2, bS);
                    }
                }
                if ((tu.x0 + sizeTu) % 8 == 0)
                {
                    // right edge of luma transform block coded with intra prediction mode
                    for (int y = (tu.y0 + 3) / 4; y < (tu.y0 + sizeTu) / 4; ++y)
                    {
                        const int x = (tu.x0 + sizeTu) / 8;
                        this->blockAt(x, y / 2).increaseBs(EDGE_VER, y % 2, bS);
                    }
                }
                if (tu.y0 % 8 == 0)
                {
                    // top edge of luma transform block coded with intra prediction mode
                    for (int x = (tu.x0 + 3) / 4; x < (tu.x0 + sizeTu) / 4; ++x)
                    {
                        const int y = tu.y0 / 8;
                        this->blockAt(x / 2, y).increaseBs(EDGE_HOR, x % 2, bS);
                    }
                }
                if ((tu.y0 + sizeTu) % 8 == 0)
                {
                    // bottom edge of luma transform block coded with intra prediction mode
                    for (int x = (tu.x0 + 3) / 4; x < (tu.x0 + sizeTu) / 4; ++x)
                    {
                        const int y = (tu.y0 + sizeTu) / 8;
                        this->blockAt(x / 2, y).increaseBs(EDGE_HOR, x % 2, bS);
                    }
                }
            }
        }

        template <class H>
        void processRc(H &h, residual_coding rc)
        {
            if (rc.cIdx == 0 && h[current(CuPredMode(rc.x0, rc.y0))] != MODE_INTRA && !h[slice_deblocking_filter_disabled_flag()])
            {
                const int sizeRc = (1 << rc.log2TrafoSize);

                const int bS = 1;

                if (rc.x0 % 8 == 0)
                {
                    // left edge of luma transform block which contains one or more non-zero transform coefficient levels
                    for (int y = (rc.y0 + 3) / 4; y < (rc.y0 + sizeRc) / 4; ++y)
                    {
                        const int x = rc.x0 / 8;
                        this->blockAt(x, y / 2).increaseBs(EDGE_VER, y % 2, bS);
                    }
                }
                if ((rc.x0 + sizeRc) % 8 == 0)
                {
                    // right edge of luma transform block which contains one or more non-zero transform coefficient levels
                    for (int y = (rc.y0 + 3) / 4; y < (rc.y0 + sizeRc) / 4; ++y)
                    {
                        const int x = (rc.x0 + sizeRc) / 8;
                        this->blockAt(x, y / 2).increaseBs(EDGE_VER, y % 2, bS);
                    }
                }
                if (rc.y0 % 8 == 0)
                {
                    // top edge of luma transform block which contains one or more non-zero transform coefficient levels
                    for (int x = (rc.x0 + 3) / 4; x < (rc.x0 + sizeRc) / 4; ++x)
                    {
                        const int y = rc.y0 / 8;
                        this->blockAt(x / 2, y).increaseBs(EDGE_HOR, x % 2, bS);
                    }
                }
                if ((rc.y0 + sizeRc) % 8 == 0)
                {
                    // bottom edge of luma transform block which contains one or more non-zero transform coefficient levels
                    for (int x = (rc.x0 + 3) / 4; x < (rc.x0 + sizeRc) / 4; ++x)
                    {
                        const int y = (rc.y0 + sizeRc) / 8;
                        this->blockAt(x / 2, y).increaseBs(EDGE_HOR, x % 2, bS);
                    }
                }
            }
        }

        template <int edgeType, class H, class Pic>
        void deblock(H &h, Pic recPictureL, Pic recPictureCb, Pic recPictureCr, int xBegin, int yBegin, int xEnd, int yEnd)
        {
            typedef typename Pic::Type Sample;

            for (int y = yBegin / 8; y < yEnd / 8; ++y)
            {
                for (int x = xBegin / 8; x < xEnd / 8; ++x)
                {
                    int poc = h[PicOrderCntVal()];
                    const Ctu &ctu = this->ctus[h[PicWidthInCtbsY()] * (y << 3 >> h[CtbLog2SizeY()]) + (x << 3 >> h[CtbLog2SizeY()])];

                    Raster<Block> block(&this->blockAt(x, y), this->stride);

                    for (int position = 0; position < 2; ++position)
                    {
                        LumaBlockEdge<Sample, edgeType> lumaBlockEdge(h, recPictureL.offset(8 * x, 8 * y), block, ctu, position);
                        lumaBlockEdge.filter();
                    }

                    const bool doFilterChroma =
                            (edgeType == EDGE_HOR && y % h[SubHeightC()] == 0) ||
                            (edgeType == EDGE_VER && x % h[SubWidthC()] == 0);

                    if (doFilterChroma)
                    {
                        auto const positions = 2 / (edgeType == EDGE_HOR ? h[SubWidthC()] : h[SubHeightC()]);

                        for (int position = 0; position < positions; ++position)
                        {
                            ChromaBlockEdge<Sample, edgeType> chromaBlockEdgeCb(h, recPictureCb.offset(8 * x / h[SubWidthC()], 8 * y / h[SubHeightC()]), block, ctu, h[pps_cb_qp_offset()], position);
                            chromaBlockEdgeCb.filter();
                            ChromaBlockEdge<Sample, edgeType> chromaBlockEdgeCr(h, recPictureCr.offset(8 * x / h[SubWidthC()], 8 * y / h[SubHeightC()]), block, ctu, h[pps_cr_qp_offset()], position);
                            chromaBlockEdgeCr.filter();
                        }
                    }
                }
            }
        }

        template <typename Sample, class H>
        void applySao2(H &h, ThreePlanes<Sample> &saoPicture, ThreePlanes<Sample> &recPicture)
        {
            for (int ry = 0; ry < h[PicHeightInCtbsY()]; ++ry)
            {
                for (int rx = 0; rx < h[PicWidthInCtbsY()]; ++rx)
                {
                    const int nCtbS = 1 << h[CtbLog2SizeY()];
                    filterBlockSao<Sample>(h, saoPicture[0], recPicture[0], rx, ry, nCtbS, nCtbS, 0);
                    filterBlockSao<Sample>(h, saoPicture[1], recPicture[1], rx, ry, nCtbS / h[SubWidthC()], nCtbS / h[SubHeightC()], 1);
                    filterBlockSao<Sample>(h, saoPicture[2], recPicture[2], rx, ry, nCtbS / h[SubWidthC()], nCtbS / h[SubHeightC()], 2);
                }
            }
        }

        template <typename Sample, class H>
        void applySaoCTU(H &h, int rx, int ry)
        {
            const int nCtbS = 1 << h[CtbLog2SizeY()];
            StateReconstructedPicture<Sample> *reconstructedPicture = h;
            auto &recPicture = *reconstructedPicture->picture;
            auto &saoPicture = *reconstructedPicture->saoPicture;

            if (h[slice_sao_luma_flag()])
            {
                filterBlockSao<Sample>(h, recPicture[0], saoPicture[0], rx, ry, nCtbS, nCtbS, 0);
            }
            if (h[slice_sao_chroma_flag()])
            {
                filterBlockSao<Sample>(h, recPicture[1], saoPicture[1], rx, ry, nCtbS / h[SubWidthC()], nCtbS / h[SubHeightC()], 1);
                filterBlockSao<Sample>(h, recPicture[2], saoPicture[2], rx, ry, nCtbS / h[SubWidthC()], nCtbS / h[SubHeightC()], 2);
            }
        }

        template <class H>
        void applySao(H &h, Raster<uint8_t> recPictureL, Raster<uint8_t> recPictureCb, Raster<uint8_t> recPictureCr)
        {
            std::vector<uint8_t> temporaryBuffer[3];
            temporaryBuffer[0].resize(h[PicSizeInCtbsY()] * h[CtbSizeY()] * h[CtbSizeY()]);
            temporaryBuffer[1].resize(h[PicSizeInCtbsY()] * h[CtbSizeY()] * h[CtbSizeY()] / 2);
            temporaryBuffer[2].resize(h[PicSizeInCtbsY()] * h[CtbSizeY()] * h[CtbSizeY()] / 2);

            Raster<uint8_t> saoPictureL(&temporaryBuffer[0][0], h[PicWidthInCtbsY()] << h[CtbLog2SizeY()]);
            Raster<uint8_t> saoPictureCb(&temporaryBuffer[1][0], h[PicWidthInCtbsY()] << h[CtbLog2SizeY()]);
            Raster<uint8_t> saoPictureCr(&temporaryBuffer[2][0], h[PicWidthInCtbsY()] << h[CtbLog2SizeY()]);

            Raster<uint8_t> recPicture[3] = { recPictureL, recPictureCb, recPictureCr };
            Raster<uint8_t> saoPicture[3] = { saoPictureL, saoPictureCb, saoPictureCr };

            applySao2(h, saoPicture, recPicture);

            for (int y = 0; y < h[pic_height_in_luma_samples()]; ++y)
            {
                for (int x = 0; x < h[pic_width_in_luma_samples()]; ++x)
                {
                    recPictureL(x, y) = saoPictureL(x, y);
                }
            }

            for (int y = 0; y < h[pic_height_in_luma_samples()] / 2; ++y)
            {
                for (int x = 0; x < h[pic_width_in_luma_samples()] / 2; ++x)
                {
                    recPictureCb(x, y) = saoPictureCb(x, y);
                    recPictureCr(x, y) = saoPictureCr(x, y);
                }
            }
        }

        // Undoes SAO and deblocking in regions where PCM and trans/quant bypass disables the loop filter by copying unfiltered reconstruction blocks into the filtered picture.
        template <typename Sample, class H>
        void restoreUnfilteredRegions(H &h, Raster<Sample> saoPicture, Raster<Sample> recPicture, int rx, int ry, int nCtbSw, int nCtbSh, int cIdx)
        {
            const int xCtb = rx * nCtbSw;
            const int yCtb = ry * nCtbSh;

            const int nW = cIdx ? 8 / h[SubWidthC()] : 8;
            const int nH = cIdx ? 8 / h[SubHeightC()] : 8;

            for (int j = 0; j < nCtbSh; j += nH)
            {
                for (int i = 0; i < nCtbSw; i += nW)
                {
                    const int xSi = xCtb + i;
                    const int ySj = yCtb + j;

                    const int xYi = (cIdx == 0) ? xSi : (xSi * h[SubWidthC()]);
                    const int yYj = (cIdx == 0) ? ySj : (ySj * h[SubHeightC()]);

                    if (!this->blockAt(xYi / 8, yYj / 8).enabled())
                    {
                        for (int dy = 0; dy < nH; ++dy)
                        {
                            memcpy(&saoPicture(xSi, ySj + dy), &recPicture(xSi, ySj + dy), nW * sizeof(Sample));
                        }
                    }
                }
            }
        }

        template <typename Sample>
        static void copySample(Raster<Sample> saoPicture, Raster<Sample> recPicture, int x, int y)
        {
            saoPicture(x, y) = recPicture(x, y);
        }

        template <typename Sample, class H>
        void filterBlockSao(H &h, Raster<Sample> saoPicture, Raster<Sample> recPicture, int rx, int ry, int nCtbSw, int nCtbSh, int cIdx)
        {
            auto const bitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];

            const Ctu &ctu = this->ctus[h[PicWidthInCtbsY()] * ry + rx];

            const int xCtb = rx * nCtbSw;
            const int yCtb = ry * nCtbSh;

            // review: these look wrong for chroma
            const int width = std::min(nCtbSw, h[pic_width_in_luma_samples()] - xCtb);
            const int height = std::min(nCtbSh, h[pic_height_in_luma_samples()] - yCtb);

            const auto SaoTypeIdx = ctu.planes[cIdx].SaoTypeIdx;

            const int SaoTypeIdxOff = 0;
            const int SaoTypeIdxBand = 1;
            const int SaoTypeIdxEdge = 2;

            if (SaoTypeIdx == SaoTypeIdxEdge)
            {
                const auto eoClass = ctu.planes[cIdx].u.eoClass;

                sao_filter_edge<Sample>(&saoPicture(xCtb, yCtb), saoPicture.stride, &recPicture(xCtb, yCtb), recPicture.stride, width, height, ctu.planes[cIdx].SaoOffsetVal, eoClass, bitDepth);

                restoreUnfilteredRegions(h, saoPicture, recPicture, rx, ry, nCtbSw, nCtbSh, cIdx);

                int top = ctu.top / (cIdx ? h[SubHeightC()] : 1);
                int left = ctu.left / (cIdx ? h[SubWidthC()] : 1);
                int right = ctu.right / (cIdx ? h[SubWidthC()] : 1);
                int bottom = ctu.bottom / (cIdx ? h[SubHeightC()] : 1);

                const bool availableL = left < xCtb;
                const bool availableR = right > xCtb + nCtbSw;
                const bool availableT = top < yCtb;
                const bool availableB = bottom > yCtb + nCtbSh;
                const bool availableTL = ctu.topLeft;
                const bool availableTR = ctu.topRight;
                const bool availableBL = ctu.bottomLeft;
                const bool availableBR = ctu.bottomRight;

                int undoT = 0;
                int undoL = 0;
                int undoR = 0;
                int undoB = 0;

                if (eoClass == 2)
                {
                    // Diagonal filter (\)

                    if (!availableTL)
                    {
                        ++undoT;
                        ++undoL;
                    }

                    if (!availableBR)
                    {
                        ++undoR;
                        ++undoB;
                    }
                }

                if (eoClass != 1)
                {
                    // All modes except vertical filter (|)
                    if (!availableL) undoL = nCtbSh;
                    if (!availableR) undoR = nCtbSh;
                }

                if (eoClass != 0)
                {
                    // All modes except horizontal filter (-)
                    if (!availableT) undoT = nCtbSw;
                    if (!availableB) undoB = nCtbSw;
                }

                if (eoClass == 3)
                {
                    // Diagonal filter (/)

                    if (availableTR)
                    {
                        --undoT;
                        --undoR;
                    }

                    if (availableBL)
                    {
                        --undoL;
                        --undoB;
                    }
                }

                if (right > xCtb + nCtbSw) right = xCtb + nCtbSw;
                if (bottom > yCtb + nCtbSh) bottom = yCtb + nCtbSh;

                for (int x = 0; x < undoT; ++x) copySample<Sample>(saoPicture, recPicture, xCtb + x, yCtb);
                for (int y = 0; y < undoL; ++y) copySample<Sample>(saoPicture, recPicture, xCtb, yCtb + y);
                for (int y = nCtbSh - undoR; y < nCtbSh; ++y) copySample<Sample>(saoPicture, recPicture, right - 1, yCtb + y);
                for (int x = nCtbSw - undoB; x < nCtbSw; ++x) copySample<Sample>(saoPicture, recPicture, xCtb + x, bottom - 1);

            }
            else if (SaoTypeIdx == SaoTypeIdxBand)
            {
                int16_t offset_table[32];

                for (int k = 0; k < 32; ++k)
                {
                    offset_table[k] = ctu.planes[cIdx].SaoOffsetVal[0];
                }

                const auto saoLeftClass = ctu.planes[cIdx].u.saoLeftClass;

                for (int k = 0; k < 4; k++)
                {
                    offset_table[(k + saoLeftClass) & 31] = ctu.planes[cIdx].SaoOffsetVal[k + 1];
                }

                sao_filter_band<Sample>(&saoPicture(xCtb, yCtb), saoPicture.stride, &recPicture(xCtb, yCtb), recPicture.stride, width, height, offset_table, bitDepth);

                restoreUnfilteredRegions<Sample>(h, saoPicture, recPicture, rx, ry, nCtbSw, nCtbSh, cIdx);
            }
            else
            {
                assert(SaoTypeIdx == SaoTypeIdxOff);

                HavocTablePredUni<Sample> *table = h;
                auto *f = *havocGetPredUni(table, 8, nCtbSw, nCtbSh, 0, 0, bitDepth);

                f(&saoPicture(xCtb, yCtb), saoPicture.stride, &recPicture(xCtb, yCtb), recPicture.stride, nCtbSw, nCtbSh, 0, 0, bitDepth);
            }
        }
    };
}

#endif
