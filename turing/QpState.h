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
 /*
  * QpState.h
  *
  *  Created on: 10 Jun 2016
  *      Author: Matteo Naccari
  *      Helper class for quantisation processing
  */

#ifndef Included_QpState_h
#define Included_QpState_h

struct QpState :
    AccessOperators<QpState>,
    ValueHolder<QpY>,
    ValueHolder<CuQpOffsetCb>,
    ValueHolder<CuQpOffsetCr>
{
    template <class H>
    QpState(H &h) :
        maskMinCuQpDeltaSize(h[MinCuQpDeltaSize()] - 1),
        maskCtbLog2SizeY(h[CtbSizeY()] - 1),
        cb_qp_offset(h[pps_cb_qp_offset()] + h[slice_cb_qp_offset()]),
        cr_qp_offset(h[pps_cr_qp_offset()] + h[slice_cr_qp_offset()]),
        QpBdOffsetY(h[::QpBdOffsetY()]),
        QpBdOffsetC(h[::QpBdOffsetC()]),
        ChromaArrayType(h[::ChromaArrayType()]),
        isIntraSlice(h[slice_type()] == I)
    {
    }

    // must be called whenever h[QpY()] has been changed, for example when QP delta is parsed
    // must be called also if h[cu_chroma_qp_offset_idx()] has been changed
    void update()
    {
        auto &QpY = this->ValueHolder<::QpY>::get(::QpY());

        const int QP_Y = QpY + this->QpBdOffsetY; // QP_Y is QP'Y

        int qPiCb = Clip3(-this->QpBdOffsetC, 57, QpY + this->cb_qp_offset + (*this)[CuQpOffsetCb()]);
        int qPiCr = Clip3(-this->QpBdOffsetC, 57, QpY + this->cr_qp_offset + (*this)[CuQpOffsetCr()]);

        int qPCb, qPCr;

        if (this->ChromaArrayType == 1)
        {
            qPCb = QpC(qPiCb);
            qPCr = QpC(qPiCr);
        }
        else
        {
            qPCb = qPiCb;
            qPCr = qPiCr;
        }

        auto const Qp_Cb = qPCb + this->QpBdOffsetC;
        auto const Qp_Cr = qPCr + this->QpBdOffsetC;

        this->table[0].Qp = QP_Y;
        this->table[1].Qp = Qp_Cb;
        this->table[2].Qp = Qp_Cr;

        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            auto const qp = this->table[cIdx].Qp;

            static const int levelScale[6] = { 40, 45, 51, 57, 64, 72 };
            this->table[cIdx].scale = levelScale[qp % 6] << (qp / 6);

            static const int scaleLookup[6] = { 26214, 23302, 20560, 18396, 16384, 14564 };
            this->table[cIdx].quantiseScale = scaleLookup[qp % 6];

            this->table[cIdx].offsetQuantiseShifted = (isIntraSlice ? 171 : 85) << 7;

            auto const bitDepth = (cIdx ? this->QpBdOffsetC : this->QpBdOffsetY) / 6 + 8;
            this->table[cIdx].quantizeShift = 29 - bitDepth + qp / 6;
        }
    }

    int getQp(int cIdx) const
    {
        return this->table[cIdx].Qp;
    }

    int getScale(int cIdx) const
    {
        return this->table[cIdx].scale;
    }

    int getQuantiseScale(int cIdx) const
    {
        return this->table[cIdx].quantiseScale;
    }

    template <class H>
    void preCu(coding_unit cu, H &h)
    {
        auto &QpY = this->ValueHolder<::QpY>::get(::QpY());

        const bool isLeftmostCuInQG(!(cu.x0 & this->maskMinCuQpDeltaSize));
        const bool isTopmostCuInQG(!(cu.y0 & this->maskMinCuQpDeltaSize));
        const int x((cu.x0 & this->maskCtbLog2SizeY) / h[MinCuQpDeltaSize()]);
        const int y((cu.y0 & this->maskCtbLog2SizeY) / h[MinCuQpDeltaSize()]);

        if (isLeftmostCuInQG && isTopmostCuInQG)
        {
            const int qPY_PREV = QpY;
            const int qPY_A = x ? this->qpLeft[y] : qPY_PREV;
            const int qPY_B = y ? this->qpAbove[x] : qPY_PREV;

            const int qPY_PRED = (qPY_A + qPY_B + 1) >> 1;
            QpY = qPY_PRED;

            this->update();
        }
    }

    template <class H>
    void postCu(coding_unit cu, H &h)
    {
        const bool isLeftmostCuInQG(!(cu.x0 & this->maskMinCuQpDeltaSize));
        const bool isTopmostCuInQG(!(cu.y0 & this->maskMinCuQpDeltaSize));
        const int x((cu.x0 & this->maskCtbLog2SizeY) >> h[Log2MinCuQpDeltaSize()]);
        const int y((cu.y0 & this->maskCtbLog2SizeY) >> h[Log2MinCuQpDeltaSize()]);
        const bool isRightmostCuInQG = !((cu.x0 + (1 << cu.log2CbSize)) & this->maskMinCuQpDeltaSize);
        const int n = std::max(1, 1 << cu.log2CbSize >> h[Log2MinCuQpDeltaSize()]);
        if (isRightmostCuInQG && isTopmostCuInQG)
        {
            for (int i = 0; i < n; ++i)
            {
                this->qpLeft[y + i] = this->ValueHolder<::QpY>::get(::QpY());
            }
        }
        const bool isBottomostCuInQG = !((cu.y0 + (1 << cu.log2CbSize)) & this->maskMinCuQpDeltaSize);
        if (isBottomostCuInQG && isLeftmostCuInQG)
        {
            for (int i = 0; i < n; ++i)
            {
                this->qpAbove[x + i] = this->ValueHolder<::QpY>::get(::QpY());
            }
        }
    }

    template <class H>
    int getQpYPred(coding_unit cu, H &h)
    {
        // review: divisions here

        const int x(((cu.x0 & this->maskCtbLog2SizeY) / h[MinCuQpDeltaSize()]) << (h[Log2MinCuQpDeltaSize()] - 3));
        const int y(((cu.y0 & this->maskCtbLog2SizeY) / h[MinCuQpDeltaSize()]) << (h[Log2MinCuQpDeltaSize()] - 3));

        // Derive qPy_prev according to Clause 8.6.1
        int qPY_PREV, qPY_A, qPY_B, qPY_PRED;

        qPY_PREV = getLastCodedQp(cu, h);

        const bool availableA = !!((cu.x0 & this->maskCtbLog2SizeY) / h[MinCuQpDeltaSize()]);
        const bool availableB = !!((cu.y0 & this->maskCtbLog2SizeY) / h[MinCuQpDeltaSize()]);

        qPY_A = availableA ? this->currentQpInternal[rasterToZscan[y][x - 1]] : qPY_PREV;
        qPY_B = availableB ? this->currentQpInternal[rasterToZscan[y - 1][x]] : qPY_PREV;

        qPY_PRED = (qPY_A + qPY_B + 1) >> 1;

        return qPY_PRED;
    }

    void setCanWrite(bool flag) { canWriteQp = flag; }
    bool getCanWrite() { return canWriteQp; }
    void setCodedQp(int qp) { codedQp = qp; }
    int  getCodedQp() { return codedQp; }
    void setQpInternal(int row, int col, int size, int value, int flagValue = 1)
    {
        int idx = rasterToZscan[row][col];
        assert(0 <= row && row < 8);
        assert(0 <= col && col < 8);
        assert(idx + size <= 64);
        memset(currentQpInternal + idx, value, size);
        memset(currentValid + idx, flagValue, size);
    }
    int getQpInternal(int row, int col)
    {
        assert(0 <= row && row < 8);
        assert(0 <= col && col < 8);
        int idx = rasterToZscan[row][col];
        return currentQpInternal[idx];
    }

    template<class H>
    void initInternalMemory(H &h)
    {
        currentQpInternal = qpInternal[0];
        previousQpInternal = qpInternal[1];
        currentValid = validBLock[0];
        previousValid = validBLock[1];
        memset(currentValid, 0, 64);
        memset(previousValid, 0, 64);
        int QpY = h[SliceQpY()];
        memset(currentQpInternal, QpY, 64);
        memset(previousQpInternal, QpY, 64);
    }

    void swapInternalMemory()
    {
        uint8_t *temp;
        temp = previousQpInternal;
        previousQpInternal = currentQpInternal;
        currentQpInternal = temp;
        temp = previousValid;
        previousValid = currentValid;
        currentValid = temp;
    }

    template<class H>
    int getLastCodedQp(coding_unit &cu, H &h)
    {
        int rowQg = (cu.y0 - (cu.y0 & this->maskCtbLog2SizeY));
        int colQg = (cu.x0 - (cu.x0 & this->maskCtbLog2SizeY));
        int rowModulo = ((cu.y0 & this->maskCtbLog2SizeY) >> h[Log2MinCuQpDeltaSize()]) << (h[Log2MinCuQpDeltaSize()] - 3);
        int colModulo = ((cu.x0 & this->maskCtbLog2SizeY) >> h[Log2MinCuQpDeltaSize()]) << (h[Log2MinCuQpDeltaSize()] - 3);
        int zScanIdx = rasterToZscan[rowModulo][colModulo];
        int rsIdxQg = rowQg * (h[pic_width_in_luma_samples()] >> h[Log2MinCuQpDeltaSize()]) + colQg;

        int lastPartIdx = getLastPartIdx(zScanIdx, currentValid);

        bool sliceQpCondition = false;

        bool firstQgInCtb = !rowModulo && !colModulo;

        if (firstQgInCtb)
        {
            // Check conditions for Step #1 of Clause 8.6.1
            sliceQpCondition |= (rsIdxQg == h[slice_segment_address()]);
            sliceQpCondition |= (rsIdxQg == h[CtbAddrRsToTs()]);
            sliceQpCondition |= ((h[CtbAddrInRs()] % h[PicWidthInCtbsY()] == 0) && h[entropy_coding_sync_enabled_flag()]);
        }

        if (sliceQpCondition)
            return h[SliceQpY()];

        if (lastPartIdx >= 0)
        {
            return currentQpInternal[lastPartIdx];
        }
        else
        {
            lastPartIdx = getLastPartIdx(64, previousValid);
            assert(lastPartIdx >= 0);
            return previousQpInternal[lastPartIdx];
        }
    }

    void setQpValue(int value)
    {
        auto &QpY = this->ValueHolder<::QpY>::get(::QpY());
        QpY = value;
        this->update();
    }

    const int getMaskCtb() { return maskCtbLog2SizeY; }

    struct Entry
    {
        int Qp;
        int scale;
        int quantiseScale;
        int offsetQuantiseShifted;
        int quantizeShift; // review: consistent spelling of quantize (as HEVC standard doc - check)
    };

    Entry const &lookup(int cIdx) const
    {
        return this->table[cIdx];
    }

private:
    bool const isIntraSlice;
    const int QpBdOffsetY;
    const int QpBdOffsetC;
    const int maskMinCuQpDeltaSize;
    const int maskCtbLog2SizeY;
    const int cr_qp_offset;
    const int cb_qp_offset;
    const int ChromaArrayType;
    int qpLeft[8];
    int qpAbove[8];
    Entry table[3];
    bool canWriteQp;
    int codedQp;
    int rasterToZscan[8][8] = { {0,  1,  4,  5,  16, 17, 20, 21},
            {2,  3,  6,  7,  18, 19, 22, 23},
            {8,  9,  12, 13, 24, 25, 28, 30},
            {10, 11, 14, 15, 26, 27, 29, 31},
            {32, 33, 36, 37, 48, 49, 52, 53},
            {34, 35, 38, 39, 50, 51, 54, 55},
            {40, 41, 44, 45, 56, 57, 60, 61},
            {42, 43, 46, 47, 58, 59, 62, 63} };
    uint8_t qpInternal[2][64];
    uint8_t validBLock[2][64];
    uint8_t *currentQpInternal;
    uint8_t *previousQpInternal;
    uint8_t *currentValid;
    uint8_t *previousValid;

    int getLastPartIdx(int zScanIdx, uint8_t *flag)
    {
        assert(0 <= zScanIdx && zScanIdx <= 64);
        int idx = zScanIdx - 1;
        while (idx >= 0 && flag[idx] == 0)
        {
            idx--;
        }
        return idx;
    }
};


#endif /* QPSTATE_H_ */
