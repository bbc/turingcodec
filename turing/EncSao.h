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

#pragma once

//#include "SyntaxCtu.hpp"
#include <algorithm>
///#include "Search.h"
//#include "StateSpatial.h"
#include "Global.h"
//#include "StatePicture.h"
//#include "StateEncode.h"
//#include "Cost.h"

class EncSao
{
public : 
    EncSao() : m_distLuma(-1), m_distChroma(-1) { }

    static int sign3(int x)
    {
        return ((x > 0) ? 1 : ((x == 0) ? 0 : -1));
    }
    
    inline double roundSao(int bitDepth, double x)
    {
        if(bitDepth == 8)
            return ((x) >= 0 ? ((int)((x)+0.5)) : ((int)((x)-0.5)));
        else
            return ((x)>0) ? (int)(((int)(x)+(1 << (bitDepth - 8 - 1))) / (1 << (bitDepth - 8))) : ((int)(((int)(x)-(1 << (bitDepth - 8 - 1))) / (1 << (bitDepth - 8))));
    }
    int64_t estSaoDist(int64_t num_in_cat, int64_t off, int64_t diff, int shift)
    {
        int64_t dist = num_in_cat * off* off - 2 *off *diff;
        if(shift == 0)
            return dist;

        if (dist >= 0 )
           return (dist >> (2 * shift));
        else
            return -1*((-dist) >> (2 * shift));
    }
    
  
    template <typename Sample>
    static int band_offset_chroma_stats(const Sample *srcUP, const Sample *srcVP, intptr_t stride_src, const Sample *refUP, const Sample *refVP, intptr_t stride_ref, int64_t* E_long, int64_t* num_in_band_long, int height, int width, int shift)
    {
        const Sample *srcU = srcUP;
        const Sample *refU = refUP;
        const Sample *srcV = srcVP;
        const Sample *refV = refVP;

        int hist_max = 0, band_start = 0;
        int64_t max_cumul = 0, curr_cumul = 0;
        int band;

        //start from second row:
        refU += stride_ref;
        srcU += stride_src;
        refV += stride_ref;
        srcV += stride_src;
        for (int y = 1; y < (height - 1); ++y)
        {
            for (int x = 1; x < (width - 1); x++)
            {
                band = (refU[x] >> (3 + shift));
                num_in_band_long[band]++;
                E_long[band] += (srcU[x] - refU[x]);
                band = (refV[x] >> (3 + shift));
                num_in_band_long[band]++;
                E_long[band] += (srcV[x] - refV[x]);
            }
            refU += stride_ref;
            srcU += stride_src;
            refV += stride_ref;
            srcV += stride_src;
        }

        for (int b = 0; b<29; b++)
        {
            curr_cumul = num_in_band_long[b] + num_in_band_long[b + 1] + num_in_band_long[b + 2] + num_in_band_long[b + 3];
            if (curr_cumul > max_cumul)
            {
                max_cumul = curr_cumul;
                hist_max = static_cast<int>(num_in_band_long[b]);
                band_start = b;
            }
        }
        band_start = band_start + 1;
        band_start = (band_start < 2 ? 2 : band_start);
        return band_start;
    }
  
    template <typename Sample>
    static int band_offset_luma_stats(const Sample *srcP, intptr_t stride_src, const Sample *refP, intptr_t stride_ref, int64_t* E_long, int64_t* num_in_band_long, int height, int width, int shift)
    {
        const Sample *src = srcP;
        const Sample *ref = refP;

        int hist_max = 0, band_start = 0;
        int64_t max_cumul = 0, curr_cumul = 0;
        int band;

        //start from second row:
        ref += stride_ref;
        src += stride_src;
        for (int y = 1; y < (height - 1); ++y)
        {
            for (int x = 1; x < (width - 1); x++)
            {
                band = (ref[x] >> (3 + shift));
                num_in_band_long[band]++;
                E_long[band] += (src[x] - ref[x]);
            }
            ref += stride_ref;
            src += stride_src;
        }

        for (int b = 0; b<29; b++)
        {
            curr_cumul = num_in_band_long[b] + num_in_band_long[b + 1] + num_in_band_long[b + 2] + num_in_band_long[b + 3];
            if (curr_cumul > max_cumul)
            {
                max_cumul = curr_cumul;
                hist_max = static_cast<int>(num_in_band_long[b]);
                band_start = b;
            }
        }
        band_start = band_start + 1;
        band_start = (band_start < 2 ? 2 : band_start);
        return band_start;
    }

    template <typename Sample>
    static void edge_offset_stats_class0(const Sample *srcP, intptr_t stride_src, const Sample *refP, intptr_t stride_ref, int64_t* E, int64_t* num_in_cat, int height, int width)
    {
        const Sample *src = srcP;
        const Sample *ref = refP;
        int signA, signB, edgeIdx, cat;
        int edgeIdx2category[5] = { 1,2,0,3,4 };

        for (int cat = 0; cat<5; cat++)
        {
            E[cat] = 0;
            num_in_cat[cat] = 0;
        }
        //start from second row:
        ref += stride_ref;
        src += stride_src;
        for (int y = 1; y < (height - 1); ++y)
        {
            signA = sign3((int)(ref[1] - ref[0]));  
            signB = sign3((int)(ref[1] - ref[2]));
            edgeIdx = 2 + signA + signB; 
            cat = edgeIdx2category[edgeIdx];
            num_in_cat[cat]++;    
            E[cat] += (src[1] - ref[1]);
            for (int x = 1; x < (width - 1); x++)
            {
                signA = -signB;
                signB = sign3((int)(ref[x] - ref[x + 1]));
                edgeIdx = 2 + signA + signB;
                cat = edgeIdx2category[edgeIdx];
                num_in_cat[cat]++;
                E[cat] += (src[x] - ref[x]);
            }
            ref += stride_ref;
            src += stride_src;
        }
    }

    template <typename Sample>
    static void edge_offset_stats_class1(const Sample *srcP, intptr_t stride_src, const Sample *refP, intptr_t stride_ref, int64_t* E, int64_t* num_in_cat, int height, int width)
    {
        const Sample *src = srcP;
        const Sample *ref = refP;
        int signA, signB, edgeIdx, cat;
        int edgeIdx2category[5] = { 1,2,0,3,4 };

        for (int cat = 0; cat<5; cat++)
        {
            E[cat] = 0;
            num_in_cat[cat] = 0;
        }
        //start from second row:
        ref += stride_ref;
        src += stride_src;
        for (int y = 1; y < (height - 1); ++y)
        {
            for (int x = 1; x < (width - 1); x++)
            {
                signA = sign3((int)(ref[x] - ref[x - stride_ref]));
                signB = sign3((int)(ref[x] - ref[x + stride_ref]));
                edgeIdx = 2 + signA + signB;
                cat = edgeIdx2category[edgeIdx];
                num_in_cat[cat]++;
                E[cat] += (src[x] - ref[x]);
            }
            ref += stride_ref;
            src += stride_src;
        }
    }

    template <typename Sample>
    static void edge_offset_stats_class2(const Sample *srcP, intptr_t stride_src, const Sample *refP, intptr_t stride_ref, int64_t* E, int64_t* num_in_cat, int height, int width)
    {
        const Sample *src = srcP;
        const Sample *ref = refP;
        int signA, signB, edgeIdx, cat;
        int edgeIdx2category[5] = { 1,2,0,3,4 };

        for (int cat = 0; cat<5; cat++)
        {
            E[cat] = 0;
            num_in_cat[cat] = 0;
        }
        //start from second row:
        ref += stride_ref;
        src += stride_src;
        for (int y = 1; y < (height - 1); ++y)
        {
            for (int x = 1; x < (width - 1); x++)
            {
                signA = sign3((int)(ref[x] - ref[x - stride_ref - 1]));
                signB = sign3((int)(ref[x] - ref[x + stride_ref + 1]));
                edgeIdx = 2 + signA + signB;
                cat = edgeIdx2category[edgeIdx];
                num_in_cat[cat]++;
                E[cat] += (src[x] - ref[x]);
            }
            ref += stride_ref;
            src += stride_src;
        }
    }

    template <typename Sample>
    static void edge_offset_stats_class3(const Sample *srcP, intptr_t stride_src, const Sample *refP, intptr_t stride_ref, int64_t* E, int64_t* num_in_cat, int height, int width)
    {
        const Sample *src = srcP;
        const Sample *ref = refP;
        int signA, signB, edgeIdx, cat;
        int edgeIdx2category[5] = { 1,2,0,3,4 };

        for (int cat = 0; cat<5; cat++)
        {
            E[cat] = 0;
            num_in_cat[cat] = 0;
        }
        //start from second row:
        ref += stride_ref;
        src += stride_src;

        for (int y = 1; y < (height - 1); ++y)
        {
            for (int x = 1; x < (width - 1); x++)
            {
                signA = sign3((int)(ref[x] - ref[x - stride_ref + 1]));
                signB = sign3((int)(ref[x] - ref[x + stride_ref - 1]));
                edgeIdx = 2 + signA + signB;
                cat = edgeIdx2category[edgeIdx];
                num_in_cat[cat]++;
                E[cat] += (src[x] - ref[x]);
            }
            ref += stride_ref;
            src += stride_src;
        }
    }

    template <typename Sample, class H>
    void saoRdEstimateLuma(H &h, PictureWrapper &orgPicWrapper, StateReconstructedPicture<Sample> *recPic)
    {
        StateEncode* stateEncode = h;
        const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
        const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];
        int xCtb = rx << h[CtbLog2SizeY()];
        int yCtb = ry << h[CtbLog2SizeY()];
        int xEnd = std::min(((rx + 1) << h[CtbLog2SizeY()]), h[pic_width_in_luma_samples()]);
        int yEnd = std::min(((ry + 1) << h[CtbLog2SizeY()]), h[pic_height_in_luma_samples()]);
        int width = xEnd - xCtb;
        int height = yEnd - yCtb;

        auto &orgPicture = static_cast<PictureWrap<Sample> &>(orgPicWrapper);
        Raster<Sample> sourceSamples = orgPicture(xCtb, yCtb, 0);

        Raster<Sample> recSamples = ((ThreePlanes<Sample>)(*(recPic->picture)))(xCtb, yCtb, 0);
        if(stateEncode->saoslow)
            recSamples = ((ThreePlanes<Sample>)(*(recPic->deblockPicture)))(xCtb, yCtb, 0);

        int shift = h[BitDepthY()] - 8;

        int bestTypeIdx, bestBandPosition, bestClass;
        int64_t best_offset[5] = { 0,0,0,0,0 };

        int64_t offset[5] = { 0,0,0,0,0 };
   
        int64_t num_in_cat[5] = { 0,0,0,0,0 };
    
        int64_t E[5] = { 0,0,0,0,0 };
        int64_t num_in_band_long[32];
        int64_t E_long[32];

        int64_t off_start, off_end, temp_rate, deltaD, sign;
        int startBand, endBand;
        double curr_delta_J, deltaJ, totDeltaJ, totDeltaJclass;
        Lambda lam = getReciprocalLambda(h);
        const double lambda = 1 / (lam.asDouble());

        bestTypeIdx = 0;
        totDeltaJ = 0.0;

        edge_offset_stats_class0(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, E, num_in_cat, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthY()], static_cast<double>(abs(E[cat])) / num_in_cat[cat]))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) -1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat],sign*off_start,E[cat],shift);//((num_in_cat[cat] * off_start* off_start - 2 * sign*off_start * E[cat]) >> (2 * shift));
            temp_rate = (off_start + 1);
            deltaJ = deltaD + lambda*temp_rate;
            offset[cat] = (sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);// ((num_in_cat[cat] * off* off - 2 * sign*off * E[cat]) >> (2 * shift));
                temp_rate = (off + 1);
                curr_delta_J = deltaD + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = (sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 0;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }

        edge_offset_stats_class1(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, E, num_in_cat, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthY()], static_cast<double>((abs(E[cat])) / num_in_cat[cat])))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat], sign*off_start, E[cat], shift);
            temp_rate = (off_start + 1);
            deltaJ = deltaD + lambda*temp_rate;
            offset[cat] = (sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);
                temp_rate = (off + 1);
                curr_delta_J = deltaD + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = (sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 1;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }

        edge_offset_stats_class2(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, E, num_in_cat, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthY()], static_cast<double>(abs(E[cat])) / num_in_cat[cat]))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat], sign*off_start, E[cat], shift);
            temp_rate = (off_start + 1);
            deltaJ = deltaD + lambda*temp_rate;
            offset[cat] = (sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);
                temp_rate = (off + 1);
                curr_delta_J = deltaD + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = (sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 2;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }

        edge_offset_stats_class3(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, E, num_in_cat, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthY()], static_cast<double>(abs(E[cat])) / num_in_cat[cat]))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat], sign*off_start, E[cat], shift);
            temp_rate = (off_start + 1);
            deltaJ = deltaD + lambda*temp_rate;
            offset[cat] = (sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);
                temp_rate = (off + 1);
                curr_delta_J = deltaD + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = (sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 3;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }


        for (int b = 0; b<32; b++)
        {
            E_long[b] = 0;
            num_in_band_long[b] = 0;
        }
        startBand = band_offset_luma_stats(sourceSamples.p, sourceSamples.stride, recSamples.p, recSamples.stride, E_long, num_in_band_long, height, width, shift);
        endBand = std::min(0,abs(startBand - 4));
        for (int bandPosition = startBand; bandPosition >= endBand; bandPosition--)
        {
            totDeltaJclass = 0.0;
            for (int band = 1; band < 5; band++)
            {
                offset[band] = (num_in_band_long[band + bandPosition - 1] == 0 ? 0 : (int)roundSao(h[BitDepthY()], static_cast<double>(abs(E_long[band + bandPosition - 1])) / num_in_band_long[band + bandPosition - 1]));
                sign = (E_long[band + bandPosition - 1] >= 0 ? 1 : -1);
                offset[band] = sign* ((abs(offset[band]) < ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1)) ? abs(offset[band]) : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1));
                deltaD = estSaoDist(num_in_band_long[band + bandPosition - 1], offset[band], E_long[band + bandPosition - 1],shift);
                temp_rate = (abs(offset[band]) + 2);
                totDeltaJclass += deltaD + lambda*temp_rate;
            }
            if (totDeltaJclass < totDeltaJ)
            {
                totDeltaJ = totDeltaJclass;
                bestTypeIdx = 1;
                bestBandPosition = bandPosition;
                for (int cat = 1; cat<5; cat++)
                {
                    best_offset[cat] = offset[cat];
                }
            }

        }

        if (abs(best_offset[1]) + abs(best_offset[2]) + abs(best_offset[3]) + abs(best_offset[4]) == 0)
            bestTypeIdx = 0;

        {
            h[SaoTypeIdx(0, rx, ry)] = bestTypeIdx;
            if (bestTypeIdx == 1)
            {
                h[sao_band_position(0, rx, ry)] = bestBandPosition;
            }
            else if (bestTypeIdx == 2)
            {
                h[SaoEoClass(0, rx, ry)] = bestClass;
            }

            for (int cat = 1; cat<5; cat++)
            {
                h[sao_offset_abs(0, rx, ry, cat - 1)] = static_cast<int>(abs(best_offset[cat]));
                if (bestTypeIdx == 1)
                    h[sao_offset_sign(0, rx, ry, cat - 1)] = (best_offset[cat] < 0);
            }
        }
    }

    template <typename Sample, class H>
    void saoRdEstimateChroma(H &h, PictureWrapper &orgPicWrapper, StateReconstructedPicture<Sample> *recPic)
    {
        StateEncode* stateEncode = h;
        const int rx = h[CtbAddrInRs()] % h[PicWidthInCtbsY()];
        const int ry = h[CtbAddrInRs()] / h[PicWidthInCtbsY()];
        int xCtb = (rx << h[CtbLog2SizeY()]) >> 1;
        int yCtb = (ry << h[CtbLog2SizeY()]) >> 1;
        int xEnd = (std::min(((rx + 1) << h[CtbLog2SizeY()]), h[pic_width_in_luma_samples()]) >> 1);
        int yEnd = (std::min(((ry + 1) << h[CtbLog2SizeY()]), h[pic_height_in_luma_samples()]) >> 1);
        int width = xEnd - xCtb;
        int height = yEnd - yCtb;

        auto &orgPicture = static_cast<PictureWrap<Sample> &>(orgPicWrapper);
        Raster<Sample> sourceSamplesU = orgPicture(xCtb, yCtb, 1);
        Raster<Sample> sourceSamplesV = orgPicture(xCtb, yCtb, 2);

        Raster<Sample> recSamplesU = ((ThreePlanes<Sample>)(*(recPic->picture)))(xCtb, yCtb, 1);
        Raster<Sample> recSamplesV = ((ThreePlanes<Sample>)(*(recPic->picture)))(xCtb, yCtb, 2);
        if(stateEncode->saoslow)
        {
            recSamplesU = ((ThreePlanes<Sample>)(*(recPic->deblockPicture)))(xCtb, yCtb, 1);
            recSamplesV = ((ThreePlanes<Sample>)(*(recPic->deblockPicture)))(xCtb, yCtb, 2);
        }

        int bestTypeIdx, bestBandPosition, bestClass;
        int best_offset[5] = { 0,0,0,0,0 };

        int shift = h[BitDepthC()] - 8;

        int offset[5] = { 0,0,0,0,0 };

        int64_t num_in_cat[5] = { 0,0,0,0,0 };
        int64_t num_in_catU[5] = { 0,0,0,0,0 };
        int64_t num_in_catV[5] = { 0,0,0,0,0 };

        int64_t E[5] = { 0,0,0,0,0 };
        int64_t EU[5] = { 0,0,0,0,0 };
        int64_t EV[5] = { 0,0,0,0,0 };

        int64_t num_in_band_long[32];
        int64_t E_long[32];

        int64_t off_start, off_end, temp_rate, deltaD, sign;
        int startBand, endBand;
        double curr_delta_J, deltaJ, totDeltaJ, totDeltaJclass;
        Lambda lam = getReciprocalLambda(h);
        const double lambda = 1 / (lam.asDouble());
        int distScale = 4;

        bestTypeIdx = 0;
        totDeltaJ = 0.0;

        edge_offset_stats_class0(sourceSamplesU.p, sourceSamplesU.stride, recSamplesU.p, recSamplesU.stride, EU, num_in_catU, height, width);
        edge_offset_stats_class0(sourceSamplesV.p, sourceSamplesV.stride, recSamplesV.p, recSamplesV.stride, EV, num_in_catV, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            E[cat] = EU[cat] + EV[cat];
            num_in_cat[cat] = num_in_catU[cat] + num_in_catV[cat];
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthC()], static_cast<double>(abs(E[cat])) / num_in_cat[cat]))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthC()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat], sign*off_start, E[cat], shift);
            temp_rate = (off_start + 1);
            deltaJ = deltaD * distScale + lambda*temp_rate;
            offset[cat] = static_cast<int>(sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);
                temp_rate = (off + 1);
                curr_delta_J = deltaD * distScale + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = static_cast<int>(sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 0;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }

        edge_offset_stats_class1(sourceSamplesU.p, sourceSamplesU.stride, recSamplesU.p, recSamplesU.stride, EU, num_in_catU, height, width);
        edge_offset_stats_class1(sourceSamplesV.p, sourceSamplesV.stride, recSamplesV.p, recSamplesV.stride, EV, num_in_catV, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            E[cat] = EU[cat] + EV[cat];
            num_in_cat[cat] = num_in_catU[cat] + num_in_catV[cat];
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthC()], static_cast<double>(abs(E[cat])) / num_in_cat[cat]))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthC()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat], sign*off_start, E[cat], shift);
            temp_rate = (off_start + 1);
            deltaJ = deltaD * distScale + lambda*temp_rate;
            offset[cat] = static_cast<int>(sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);
                temp_rate = (off + 1);
                curr_delta_J = deltaD * distScale + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = static_cast<int>(sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 1;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }

        edge_offset_stats_class2(sourceSamplesU.p, sourceSamplesU.stride, recSamplesU.p, recSamplesU.stride, EU, num_in_catU, height, width);
        edge_offset_stats_class2(sourceSamplesV.p, sourceSamplesV.stride, recSamplesV.p, recSamplesV.stride, EV, num_in_catV, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            E[cat] = EU[cat] + EV[cat];
            num_in_cat[cat] = num_in_catU[cat] + num_in_catV[cat];
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthC()], static_cast<double>(abs(E[cat])) / num_in_cat[cat]))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthC()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat], sign*off_start, E[cat], shift);
            temp_rate = (off_start + 1);
            deltaJ = deltaD * distScale + lambda*temp_rate;
            offset[cat] = static_cast<int>(sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);
                temp_rate = (off + 1);
                curr_delta_J = deltaD * distScale + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = static_cast<int>(sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 2;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }

        edge_offset_stats_class3(sourceSamplesU.p, sourceSamplesU.stride, recSamplesU.p, recSamplesU.stride, EU, num_in_catU, height, width);
        edge_offset_stats_class3(sourceSamplesV.p, sourceSamplesV.stride, recSamplesV.p, recSamplesV.stride, EV, num_in_catV, height, width);
        totDeltaJclass = 0.0;
        for (int cat = 1; cat<5; cat++)
        {
            E[cat] = EU[cat] + EV[cat];
            num_in_cat[cat] = num_in_catU[cat] + num_in_catV[cat];
            sign = (cat <= 2 ? 1 : -1);
            off_start = abs(num_in_cat[cat] == 0 ? 0 : ((int)roundSao(h[BitDepthC()], static_cast<double>(abs(E[cat])) / num_in_cat[cat]))) + 1;
            off_start = off_start < ((1 << (std::min(h[BitDepthC()], 10) - 5)) - 1) ? off_start : ((1 << (std::min(h[BitDepthY()], 10) - 5)) - 1);
            off_end = std::min(0, abs((int)off_start - 2));
            deltaD = estSaoDist(num_in_cat[cat], sign*off_start, E[cat], shift);
            temp_rate = (off_start + 1);
            deltaJ = deltaD * distScale + lambda*temp_rate;
            offset[cat] = static_cast<int>(sign*off_start);
            for (int64_t off = off_start - 1; off >= off_end; off--)
            {
                deltaD = estSaoDist(num_in_cat[cat], sign*off, E[cat], shift);
                temp_rate = (off + 1);
                curr_delta_J = deltaD * distScale + lambda*temp_rate;
                if (curr_delta_J < deltaJ)
                {
                    deltaJ = curr_delta_J;
                    offset[cat] = static_cast<int>(sign*off);
                }
            }
            totDeltaJclass += deltaJ;
        }
        if (totDeltaJclass < totDeltaJ)
        {
            totDeltaJ = totDeltaJclass;
            bestTypeIdx = 2;
            bestClass = 3;
            for (int cat = 1; cat<5; cat++)
            {
                best_offset[cat] = offset[cat];
            }
        }

        for (int b = 0; b<32; b++)
        {
            E_long[b] = 0;
            num_in_band_long[b] = 0;
        }
        startBand = band_offset_chroma_stats(sourceSamplesU.p, sourceSamplesV.p, sourceSamplesU.stride, recSamplesU.p, recSamplesV.p, recSamplesU.stride, E_long, num_in_band_long, height, width, shift);
        endBand = std::min(0, abs(startBand - 4));
        for (int bandPosition = startBand; bandPosition >= endBand; bandPosition--)
        {
            totDeltaJclass = 0.0;
            for (int band = 1; band < 5; band++)
            {
                offset[band] = (num_in_band_long[band + bandPosition - 1] == 0 ? 0 : (int)roundSao(h[BitDepthC()], static_cast<double>(abs(E_long[band + bandPosition - 1])) / num_in_band_long[band + bandPosition - 1]));
                sign = (E_long[band + bandPosition - 1] >= 0 ? 1 : -1);
                offset[band] = static_cast<int>(sign* ((abs(offset[band]) < ((1 << (std::min(h[BitDepthC()], 10) - 5)) - 1)) ? abs(offset[band]) : ((1 << (std::min(h[BitDepthC()], 10) - 5)) - 1)));
                deltaD = estSaoDist(num_in_band_long[band + bandPosition - 1], offset[band], E_long[band + bandPosition - 1], shift);
                temp_rate = (abs(offset[band]) + 2);
                totDeltaJclass += deltaD * distScale + lambda*temp_rate;
            }
            if (totDeltaJclass < totDeltaJ)
            {
                totDeltaJ = totDeltaJclass;
                bestTypeIdx = 1;
                bestBandPosition = bandPosition;
                for (int cat = 1; cat<5; cat++)
                {
                    best_offset[cat] = static_cast<int>(offset[cat]);
                }
            }

        }

        if (abs(best_offset[1]) + abs(best_offset[2]) + abs(best_offset[3]) + abs(best_offset[4]) == 0)
            bestTypeIdx = 0;

        {
            h[SaoTypeIdx(1, rx, ry)] = bestTypeIdx;
            h[SaoTypeIdx(2, rx, ry)] = bestTypeIdx;
            if (bestTypeIdx == 1)
            {
                h[sao_band_position(1, rx, ry)] = bestBandPosition;
                h[sao_band_position(2, rx, ry)] = bestBandPosition;
            }
            else if (bestTypeIdx == 2)
            {
                h[SaoEoClass(1, rx, ry)] = (bestClass == -1 ? 0 : bestClass);
                h[SaoEoClass(2, rx, ry)] = (bestClass == -1 ? 0 : bestClass);
            }

            for (int cat = 1; cat<5; cat++)
            {
                h[sao_offset_abs(1, rx, ry, cat - 1)] = abs(best_offset[cat]);
                h[sao_offset_abs(2, rx, ry, cat - 1)] = abs(best_offset[cat]);
                if (bestTypeIdx == 1)
                {
                    h[sao_offset_sign(1, rx, ry, cat - 1)] = (best_offset[cat] < 0);
                    h[sao_offset_sign(2, rx, ry, cat - 1)] = (best_offset[cat] < 0);
                }
            }
        }
    }

    template <typename Sample>
    static uint32_t ssd(const Sample *pA, intptr_t strideA, const Sample *pB, intptr_t strideB, int w, int h)
    {
        uint32_t ssd = 0;
        for (int y = 0; y < h; ++y)
        {
            for (int x = 0; x < w; ++x)
            {
                const int diff = pA[x + y * strideA] - pB[x + y * strideB];
                ssd += diff * diff;
            }
        }
        if (sizeof(Sample) == 2)
            ssd >>= 4;
        return ssd;
    }

    template <typename Sample, class H>
    int  computeSaoDistortion(H &h, PictureWrapper &orgPicWrapper, StateReconstructedPicture<Sample> *currPic, int rx, int ry)
    {
        StateEncode* stateEncode = h;
        LoopFilter::Ctu *ctu = new LoopFilter::Ctu();
        int xCtb = rx << h[CtbLog2SizeY()];
        int yCtb = ry << h[CtbLog2SizeY()];
        int xEnd = std::min(((rx + 1) << h[CtbLog2SizeY()]), h[pic_width_in_luma_samples()]);
        int yEnd = std::min(((ry + 1) << h[CtbLog2SizeY()]), h[pic_height_in_luma_samples()]);
        int width = xEnd - xCtb;
        int height = yEnd - yCtb;

        int distortion = 0, bestTypeIdx;

        auto &orgPicture = static_cast<PictureWrap<Sample> &>(orgPicWrapper);
        Raster<Sample> sourceSamplesY = orgPicture(xCtb, yCtb, 0);
        Raster<Sample> sourceSamplesCb = orgPicture(xCtb, yCtb, 1);
        Raster<Sample> sourceSamplesCr = orgPicture(xCtb, yCtb, 2);

        auto &recPic = *currPic->picture;
        auto &saoPic = *currPic->saoPicture;
        Raster<Sample> recPictureY = recPic[0];
        Raster<Sample> saoPictureY = saoPic[0];
        Raster<Sample> recPictureCb = recPic[1];
        Raster<Sample> saoPictureCb = saoPic[1];
        Raster<Sample> recPictureCr = recPic[2];
        Raster<Sample> saoPictureCr = saoPic[2];
        Raster<Sample> recSamplesY = ((ThreePlanes<Sample>)(*(currPic->picture)))(xCtb, yCtb, 0);
        Raster<Sample> saoSamplesY = ((ThreePlanes<Sample>)(*(currPic->saoPicture)))(xCtb, yCtb, 0);
        Raster<Sample> recSamplesCb = ((ThreePlanes<Sample>)(*(currPic->picture)))(xCtb, yCtb, 1);
        Raster<Sample> saoSamplesCb = ((ThreePlanes<Sample>)(*(currPic->saoPicture)))(xCtb, yCtb, 1);
        Raster<Sample> recSamplesCr = ((ThreePlanes<Sample>)(*(currPic->picture)))(xCtb, yCtb, 2);
        Raster<Sample> saoSamplesCr = ((ThreePlanes<Sample>)(*(currPic->saoPicture)))(xCtb, yCtb, 2);

        if (stateEncode->saoslow)
        {
            auto &debPic = *currPic->deblockPicture;
            recPictureY = debPic[0];
            recPictureCb = debPic[1];
            recPictureCr = debPic[2];
            recSamplesY = ((ThreePlanes<Sample>)(*(currPic->deblockPicture)))(xCtb, yCtb, 0);
            recSamplesCb = ((ThreePlanes<Sample>)(*(currPic->deblockPicture)))(xCtb, yCtb, 1);
            recSamplesCr = ((ThreePlanes<Sample>)(*(currPic->deblockPicture)))(xCtb, yCtb, 2);
        }

        //luma:
        bestTypeIdx = h[SaoTypeIdx(0, rx, ry)];
        if (bestTypeIdx == 1)
        {
            int16_t offset_table[32];
            ctu->set(h);
            for (int k = 0; k < 32; ++k)
            {
                offset_table[k] = ctu->planes[0].SaoOffsetVal[0];
            }

            const auto saoLeftClass = ctu->planes[0].u.saoLeftClass;

            for (int k = 0; k < 4; k++)
            {
                offset_table[(k + saoLeftClass) & 31] = ctu->planes[0].SaoOffsetVal[k + 1];
            }
            sao_filter_band<Sample>(&saoPictureY(xCtb, yCtb), saoPictureY.stride, &recPictureY(xCtb, yCtb), recPictureY.stride, width, height, offset_table, h[BitDepthY()]);
        }
        else if (bestTypeIdx == 2)
        {
            ctu->set(h);
            const auto eoClass = ctu->planes[0].u.eoClass;
            sao_filter_edge<Sample>(&saoPictureY(xCtb, yCtb), saoPictureY.stride, &recPictureY(xCtb, yCtb), recPictureY.stride, width, height, ctu->planes[0].SaoOffsetVal, eoClass, h[BitDepthY()]);

        }

        if (bestTypeIdx == 0)
        {
            if(m_distLuma == -1)
                m_distLuma = ssd(sourceSamplesY.p, sourceSamplesY.stride, recSamplesY.p, recSamplesY.stride, width, height);
            distortion += m_distLuma;
        }
        else
        {
            distortion += ssd(sourceSamplesY.p, sourceSamplesY.stride, saoSamplesY.p, saoSamplesY.stride, width, height);
        }

        //chroma:
        int distScale = 4;
        xCtb >>=1;
        yCtb >>= 1;
        width >>= 1;
        height >>= 1;
        bestTypeIdx = h[SaoTypeIdx(1, rx, ry)];
        if (bestTypeIdx == 1)
        {
            int16_t offset_table[32];
            ctu->set(h);
            for (int k = 0; k < 32; ++k)
            {
                offset_table[k] = ctu->planes[1].SaoOffsetVal[0];
            }

            const auto saoLeftClass = ctu->planes[1].u.saoLeftClass;

            for (int k = 0; k < 4; k++)
            {
                offset_table[(k + saoLeftClass) & 31] = ctu->planes[1].SaoOffsetVal[k + 1];
            }
            sao_filter_band<Sample>(&saoPictureCb(xCtb, yCtb), saoPictureCb.stride, &recPictureCb(xCtb, yCtb), recPictureCb.stride, width, height, offset_table, h[BitDepthY()]);
            sao_filter_band<Sample>(&saoPictureCr(xCtb, yCtb), saoPictureCr.stride, &recPictureCr(xCtb, yCtb), recPictureCr.stride, width, height, offset_table, h[BitDepthY()]);
        }
        else if (bestTypeIdx == 2)
        {
            ctu->set(h);
            const auto eoClass = ctu->planes[1].u.eoClass;
            sao_filter_edge<Sample>(&saoPictureCb(xCtb, yCtb), saoPictureCb.stride, &recPictureCb(xCtb, yCtb), recPictureCb.stride, width, height, ctu->planes[1].SaoOffsetVal, eoClass, h[BitDepthY()]);
            sao_filter_edge<Sample>(&saoPictureCr(xCtb, yCtb), saoPictureCr.stride, &recPictureCr(xCtb, yCtb), recPictureCr.stride, width, height, ctu->planes[1].SaoOffsetVal, eoClass, h[BitDepthY()]);
        }

        if (bestTypeIdx == 0)
        {
            if (m_distChroma == -1)
            {
                m_distChroma = (ssd(sourceSamplesCb.p, sourceSamplesCb.stride, recSamplesCb.p, recSamplesCb.stride, width, height) * distScale);
                m_distChroma += (ssd(sourceSamplesCr.p, sourceSamplesCr.stride, recSamplesCr.p, recSamplesCr.stride, width, height) * distScale);
            }
            distortion += m_distChroma;
        }
        else
        {
            distortion += (ssd(sourceSamplesCb.p, sourceSamplesCb.stride, saoSamplesCb.p, saoSamplesCb.stride, width, height) * distScale);
            distortion += (ssd(sourceSamplesCr.p, sourceSamplesCr.stride, saoSamplesCr.p, saoSamplesCr.stride, width, height) * distScale);
        }
        return distortion;
    }

    template <typename Sample, class H>
    void rdSao(H &h, PictureWrapper& orgPic, StateReconstructedPicture<Sample> *recPic, int rx, int ry)
    {
        StateEncode* stateEncode = h;
        Candidate<Sample> *candidate = h;
        auto contextsCostBefore = *static_cast<ContextsAndCost *>(candidate);
        auto h2 = h.template change<Search<void>>();
        sao* sample_adaptive_offset = h;
        Cost bestCost;
        int distortion = 0;
    
        bool allowMergeUp = (ry != 0);
        bool allowMergeLeft = (rx != 0);

        //if sao_slow_mode is used, SAO decisions are applied after deblocking. Deblocking is applied here:
        if (stateEncode->saoslow)
        {
            int xCtb = rx << h[CtbLog2SizeY()];
            int yCtb = ry << h[CtbLog2SizeY()];
            Picture<Sample> *picture = recPic->deblockPicture.get();
            auto &debSamplesL = (*picture)[0];
            auto &debSamplesU = (*picture)[1];
            auto &debSamplesV = (*picture)[2];


            for (int cIdx = 0; cIdx < 3; cIdx++)
            {
                int xBegin = rx << h[CtbLog2SizeY()];
                int yBegin = ry << h[CtbLog2SizeY()];
                int xEnd = std::min(((rx + 1) << h[CtbLog2SizeY()]), h[pic_width_in_luma_samples()]);
                int yEnd = std::min(((ry + 1) << h[CtbLog2SizeY()]), h[pic_height_in_luma_samples()]);
                if (cIdx != 0)
                {
                    xBegin >>= 1;
                    yBegin >>= 1;
                    xEnd >>= 1;
                    yEnd >>= 1;
                }
                for (int y = yBegin; y < yEnd; ++y)
                {
                    for (int x = xBegin; x < xEnd; ++x)
                    {
                        (*(recPic->deblockPicture))[cIdx](x, y) = (*(recPic->picture))[cIdx](x, y);
                    }
                }
            }

            if (!h[slice_deblocking_filter_disabled_flag()])
            {
                int xBegin = rx << h[CtbLog2SizeY()];
                int yBegin = ry << h[CtbLog2SizeY()];
                if (rx) xBegin += 8;
                if (ry) yBegin += 8;
                int xEnd = std::min(((rx + 1) << h[CtbLog2SizeY()]), h[pic_width_in_luma_samples()]);
                int yEnd = std::min(((ry + 1) << h[CtbLog2SizeY()]), h[pic_height_in_luma_samples()]);

                static_cast<StatePicture *>(h)->loopFilterPicture->deblock<EDGE_VER>(h, debSamplesL, debSamplesU, debSamplesV, xBegin, yBegin, xEnd, yEnd);
            }
            if (!h[slice_deblocking_filter_disabled_flag()])
            {
                int xBegin = rx << h[CtbLog2SizeY()];
                int yBegin = ry << h[CtbLog2SizeY()];
                if (ry) yBegin += 8;
                int xEnd = std::min(((rx + 1) << h[CtbLog2SizeY()]) - 8, h[pic_width_in_luma_samples()]);
                int yEnd = std::min(((ry + 1) << h[CtbLog2SizeY()]), h[pic_height_in_luma_samples()]);

                static_cast<StatePicture *>(h)->loopFilterPicture->deblock<EDGE_HOR>(h, debSamplesL, debSamplesU, debSamplesV, xBegin, yBegin, xEnd, yEnd);
            }
        }


        //RD decisions based on rate and distortion estimates:
        if (h[slice_sao_luma_flag()])
            saoRdEstimateLuma<Sample>(h, orgPic, recPic);

        if (h[slice_sao_chroma_flag()])
            saoRdEstimateChroma<Sample>(h, orgPic, recPic);

        //Full RD cost computations:
        {
            h[sao_merge_left_flag()] = 0;
            h[sao_merge_up_flag()] = 0;
            distortion = computeSaoDistortion<Sample>(h, orgPic, recPic, rx,ry);
            *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;

            Search<sao>::go(*sample_adaptive_offset, h2);
            auto contextsCostSao = *static_cast<ContextsAndCost *>(candidate);
            contextsCostSao.lambdaDistortion += distortion * getReciprocalLambda(h);
            bestCost = contextsCostSao.cost2();

            StateSpatial *stateSpatial = h;
            sao *s = h;
            SaoCtuData saoBestCtuData = stateSpatial->snakeSaoCtuData.at(s->rx, s->ry, 0);
            int bestIsMerge = 0;

            //Cost of NO SAO:
            if(h[SaoTypeIdx(0, rx, ry)] != 0 || h[SaoTypeIdx(1, rx, ry)] != 0)
            {
                h[SaoTypeIdx(0, rx, ry)] = 0;
                h[SaoTypeIdx(1, rx, ry)] = 0;
                h[SaoTypeIdx(2, rx, ry)] = 0;
                distortion = computeSaoDistortion<Sample>(h, orgPic, recPic, rx, ry);
                *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;

                Search<sao>::go(*sample_adaptive_offset, h2);
                auto contextsCostSao = *static_cast<ContextsAndCost *>(candidate);
                contextsCostSao.lambdaDistortion += distortion * getReciprocalLambda(h);
                if (contextsCostSao.cost2() < bestCost)
                {
                    saoBestCtuData = stateSpatial->snakeSaoCtuData.at(s->rx, s->ry, 0);
                    bestCost = contextsCostSao.cost2();
                }
            }

            if (allowMergeUp)
            {
                const SaoCtuData saoCtuDataMerge = stateSpatial->snakeSaoCtuData.at(s->rx, s->ry - 1, 0);
                stateSpatial->snakeSaoCtuData.commit(saoCtuDataMerge, s->rx, s->ry, 0);
                distortion = computeSaoDistortion<Sample>(h, orgPic, recPic, rx, ry);
                *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;
                h[sao_merge_up_flag()] = 1;
                h[sao_merge_left_flag()] = 0;

                Search<sao>::go(*sample_adaptive_offset, h2);
                auto contextsCostSao = *static_cast<ContextsAndCost *>(candidate);
                contextsCostSao.lambdaDistortion += distortion * getReciprocalLambda(h);
                if (contextsCostSao.cost2() < bestCost)
                {
                    bestIsMerge = 1;
                    bestCost = contextsCostSao.cost2();
                }
            }

            if (allowMergeLeft)
            {
                const SaoCtuData saoCtuDataMerge = stateSpatial->snakeSaoCtuData.at(s->rx - 1, s->ry, 0);
                stateSpatial->snakeSaoCtuData.commit(saoCtuDataMerge, s->rx, s->ry, 0);
                distortion = computeSaoDistortion<Sample>(h, orgPic, recPic, rx, ry);
                *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;
                h[sao_merge_left_flag()] = 1;
                h[sao_merge_up_flag()] = 0;

                Search<sao>::go(*sample_adaptive_offset, h2);
                auto contextsCostSao = *static_cast<ContextsAndCost *>(candidate);
                contextsCostSao.lambdaDistortion += distortion * getReciprocalLambda(h);
                if (contextsCostSao.cost2() < bestCost)
                {
                    bestIsMerge = 2;
                }
            }

            if(bestIsMerge==0)
            {
                h[sao_merge_left_flag()] = 0;
                h[sao_merge_up_flag()] = 0;
                stateSpatial->snakeSaoCtuData.commit(saoBestCtuData, s->rx, s->ry, 0);
            }
            else if (bestIsMerge==1)
            {
                h[sao_merge_left_flag()] = 0;
                h[sao_merge_up_flag()] = 1;
                const SaoCtuData saoCtuDataMerge = stateSpatial->snakeSaoCtuData.at(s->rx, s->ry - 1, 0);
                stateSpatial->snakeSaoCtuData.commit(saoCtuDataMerge, s->rx, s->ry, 0);
            }
            else if (bestIsMerge==2)
            {
                h[sao_merge_left_flag()] = 1;
                h[sao_merge_up_flag()] = 0;
                const SaoCtuData saoCtuDataMerge = stateSpatial->snakeSaoCtuData.at(s->rx - 1, s->ry, 0);
                stateSpatial->snakeSaoCtuData.commit(saoCtuDataMerge, s->rx, s->ry, 0);
            } 
        } 
    *static_cast<ContextsAndCost *>(candidate) = contextsCostBefore;
    }

    int m_distLuma;
    int m_distChroma;
};