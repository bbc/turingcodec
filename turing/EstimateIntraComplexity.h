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

#ifndef INCLUDED_EstimateIntraComplexity_h
#define INCLUDED_EstimateIntraComplexity_h

#include "Picture.h"
#include "StatePicture.h"

struct EstimateIntraComplexity
{
private:
    int                       *m_satdArray;
    int                       m_satdSum;
    int                       m_picHeightInBlocks;
    int                       m_picWidthInBlocks;
    bool                      m_isIntra;
public:
    EstimateIntraComplexity() : m_satdArray(NULL), m_satdSum(0), m_picHeightInBlocks(0), m_picWidthInBlocks(0), m_isIntra(false) {}

    void allocate(int picHeightInLumaSamples, int picWidthInLumaSamples)
    {
        m_isIntra = true;
        m_picHeightInBlocks = picHeightInLumaSamples >> 3;
        m_picWidthInBlocks = picWidthInLumaSamples >> 3;
        const int n = m_picHeightInBlocks * m_picWidthInBlocks;
        m_satdArray = new int[n];
        m_satdSum = 0;
    }

    ~EstimateIntraComplexity()
    {
        if (m_satdArray)
            delete[] m_satdArray;
    }


    template <typename Sample>
    static int computeSatd8x8(const Sample *p, intptr_t stride)
    {
        int k, i, j, jj, sad = 0;
        int diff[64], m1[8][8], m2[8][8], m3[8][8];

        for (k = 0; k < 64; k += 8)
        {
            diff[k + 0] = p[0];
            diff[k + 1] = p[1];
            diff[k + 2] = p[2];
            diff[k + 3] = p[3];
            diff[k + 4] = p[4];
            diff[k + 5] = p[5];
            diff[k + 6] = p[6];
            diff[k + 7] = p[7];
            p += stride;
        }

        //horizontal
        for (j = 0; j < 8; j++)
        {
            jj = j << 3;
            m2[j][0] = diff[jj] + diff[jj + 4];
            m2[j][1] = diff[jj + 1] + diff[jj + 5];
            m2[j][2] = diff[jj + 2] + diff[jj + 6];
            m2[j][3] = diff[jj + 3] + diff[jj + 7];
            m2[j][4] = diff[jj] - diff[jj + 4];
            m2[j][5] = diff[jj + 1] - diff[jj + 5];
            m2[j][6] = diff[jj + 2] - diff[jj + 6];
            m2[j][7] = diff[jj + 3] - diff[jj + 7];

            m1[j][0] = m2[j][0] + m2[j][2];
            m1[j][1] = m2[j][1] + m2[j][3];
            m1[j][2] = m2[j][0] - m2[j][2];
            m1[j][3] = m2[j][1] - m2[j][3];
            m1[j][4] = m2[j][4] + m2[j][6];
            m1[j][5] = m2[j][5] + m2[j][7];
            m1[j][6] = m2[j][4] - m2[j][6];
            m1[j][7] = m2[j][5] - m2[j][7];

            m2[j][0] = m1[j][0] + m1[j][1];
            m2[j][1] = m1[j][0] - m1[j][1];
            m2[j][2] = m1[j][2] + m1[j][3];
            m2[j][3] = m1[j][2] - m1[j][3];
            m2[j][4] = m1[j][4] + m1[j][5];
            m2[j][5] = m1[j][4] - m1[j][5];
            m2[j][6] = m1[j][6] + m1[j][7];
            m2[j][7] = m1[j][6] - m1[j][7];
        }

        //vertical
        for (i = 0; i < 8; i++)
        {
            m3[0][i] = m2[0][i] + m2[4][i];
            m3[1][i] = m2[1][i] + m2[5][i];
            m3[2][i] = m2[2][i] + m2[6][i];
            m3[3][i] = m2[3][i] + m2[7][i];
            m3[4][i] = m2[0][i] - m2[4][i];
            m3[5][i] = m2[1][i] - m2[5][i];
            m3[6][i] = m2[2][i] - m2[6][i];
            m3[7][i] = m2[3][i] - m2[7][i];

            m1[0][i] = m3[0][i] + m3[2][i];
            m1[1][i] = m3[1][i] + m3[3][i];
            m1[2][i] = m3[0][i] - m3[2][i];
            m1[3][i] = m3[1][i] - m3[3][i];
            m1[4][i] = m3[4][i] + m3[6][i];
            m1[5][i] = m3[5][i] + m3[7][i];
            m1[6][i] = m3[4][i] - m3[6][i];
            m1[7][i] = m3[5][i] - m3[7][i];

            m2[0][i] = m1[0][i] + m1[1][i];
            m2[1][i] = m1[0][i] - m1[1][i];
            m2[2][i] = m1[2][i] + m1[3][i];
            m2[3][i] = m1[2][i] - m1[3][i];
            m2[4][i] = m1[4][i] + m1[5][i];
            m2[5][i] = m1[4][i] - m1[5][i];
            m2[6][i] = m1[6][i] + m1[7][i];
            m2[7][i] = m1[6][i] - m1[7][i];
        }

        for (i = 0; i < 8; i++)
        {
            for (j = 0; j < 8; j++)
            {
                sad += abs(m2[i][j]);
            }
        }
        sad -= abs(m2[0][0]);

        sad = ((sad + 2) >> 2);
        if (sizeof(Sample) == 2)
            sad >>= 2;

        return sad;
    }

    template<typename Sample>
    void preAnalysis(std::shared_ptr<PictureWrapper> picture)
    {
        auto &pictureInput = static_cast<PictureWrap<Sample> &>(*picture);
        int k = 0;
        for (int j = 0; j < m_picHeightInBlocks; ++j)
        {
            for (int i = 0; i < m_picWidthInBlocks; ++i)
            {
                int rx = i << 3;
                int ry = j << 3;
                Sample *p = pictureInput(rx, ry, 0 /*LUMA*/).p;
                auto stride = pictureInput(rx, ry, 0 /*LUMA*/).stride;
                m_satdArray[k] = computeSatd8x8<Sample>(p, stride);
                m_satdSum += m_satdArray[k++];
            }
        }
    }

    int getSatdSum() { return m_isIntra ? m_satdSum : -1; }
    int getSatd(int rx, int ry) //in block coordinates
    {
        if (!m_isIntra) return -1;
        return m_satdArray[rx + ry * m_picWidthInBlocks];
    }
    int getSatdCtu(int rowCtu, int colCtu)  //in ctu coordinates
    {
        if (!m_isIntra)
            return -1;

        const int rowBlock = rowCtu << 3;
        const int colBlock = colCtu << 3;
        const int height = rowBlock + 8 > m_picHeightInBlocks ? m_picHeightInBlocks - rowBlock : 8;
        const int width = colBlock + 8 > m_picWidthInBlocks ? m_picWidthInBlocks - colBlock : 8;
        int satdSum = 0;

        for (int r = rowBlock; r < rowBlock + height; r++)
        {
            for (int c = colBlock; c < colBlock + width; c++)
            {
                satdSum += m_satdArray[r*m_picWidthInBlocks + c];
            }
        }
        return satdSum;
    }
};

#endif
