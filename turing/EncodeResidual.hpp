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

#ifndef INCLUDED_EncodeResidual_hpp
#define INCLUDED_EncodeResidual_hpp

#include "EncodeResidual.h"
#include "EstimateRate.h"
#include "Write.h"
#include<type_traits>


template <class, class> struct Invoke;

template <template <typename> class Verb, class F>
struct Invoke<Verb<void>, F> : Verb<F> {};


template <bool havePopcnt, bool is4x4, class H>
void EncodeResidual::inner(H &h)
{
    residual_coding *rc = h;
    auto const log2TrafoSize = rc->log2TrafoSize;

    using Sample = typename SampleType<H>::Type;

    StateEncodeSubstreamBase *stateEncodeSubstreamBase = h;
    StateReconstructionCache<Sample> *stateReconstructionCache = h;

    StateCodedData *stateCodedData = h;
    CodedData::Residual &residual = stateCodedData->residual;
    CodedData::SubBlock subBlock = residual.initialSubBlock(log2TrafoSize);

    int const scanIdx = h[::scanIdx()];

    // Find last subblock and last significant coefficient in that block
    auto const lastScanPos = subBlock.lastScanPos();
    auto const lastSubBlock = residual.lastSubBlock(log2TrafoSize);

    auto const rasterS = rasterScanOrder(log2TrafoSize - 2, scanIdx);
    auto const rasterC = rasterScanOrder(2, scanIdx);

    {
        int const xS = is4x4 ? 0 : (rasterS[lastSubBlock] & ((1 << (log2TrafoSize - 2)) - 1));
        int const yS = is4x4 ? 0 : (rasterS[lastSubBlock] >> (log2TrafoSize - 2));
        int const xC = (xS << 2) + (rasterC[lastScanPos] & 3);
        int const yC = (yS << 2) + (rasterC[lastScanPos] >> 2);
        setLastSignificantCoeff<last_sig_coeff_x_prefix, last_sig_coeff_x_suffix>(h, *rc, (scanIdx == 2) ? yC : xC);
        setLastSignificantCoeff<last_sig_coeff_y_prefix, last_sig_coeff_y_suffix>(h, *rc, (scanIdx == 2) ? xC : yC);
    }

    //Encode transform_skip_flag
    if (h[transform_skip_enabled_flag()] && !h[cu_transquant_bypass_flag()] && (log2TrafoSize <= h[Log2MaxTransformSkipSize()]))
        h(transform_skip_flag(rc->x0, rc->y0, rc->cIdx), ae(v));

    h(last_sig_coeff_x_prefix(), ae(v));
    h(last_sig_coeff_y_prefix(), ae(v));
    if (!is4x4)
    {
        if (h[last_sig_coeff_x_prefix()] > 3)
            h(last_sig_coeff_x_suffix(), ae(v));
        if (h[last_sig_coeff_y_prefix()] > 3)
            h(last_sig_coeff_y_suffix(), ae(v));
    }

    int lastGreater1Flag = 1;
    int greater1Ctx = -1;
    uint16_t snake = 0;

    // loop over sub blocks
    int i = lastSubBlock;
    do
    {
        int const xS = is4x4 ? 0 : (rasterS[i] & ((1 << (log2TrafoSize - 2)) - 1));
        int const yS = is4x4 ? 0 : (rasterS[i] >> (log2TrafoSize - 2));

        bool const codedSubBlockFlag = !!residual.codedSubBlockFlag(i);

        uint16_t const sigCoeffFlagsActual = codedSubBlockFlag ? subBlock.sigCoeffFlags() : 0;

        auto const neighbouringCsbFlags = (snake >> (7 + yS - xS)) & 3;
        uint16_t const bits = 3 << (7 + yS - xS);
        snake &= ~bits;
        if (codedSubBlockFlag)
            snake |= bits;

        int inferSbDcSigCoeffFlag = 0;
        if ((i != lastSubBlock) && (i != 0))
        {
            h(EncodeDecision<coded_sub_block_flag>(codedSubBlockFlag, (rc->cIdx ? 2 : 0) + (neighbouringCsbFlags ? 1 : 0)));
            inferSbDcSigCoeffFlag = 1;
        }

        if (codedSubBlockFlag || i == 0)
        {
            CtxIncSigCoeffFlagCalculator<is4x4> ctxIncSigCoeffFlagCalculator(scanIdx, i == 0, rc->cIdx, neighbouringCsbFlags, rc->log2TrafoSize);

            uint16_t sigCoeffFlags = sigCoeffFlagsActual;

            if ((sigCoeffFlags & 0x7fff) != 0)
                inferSbDcSigCoeffFlag = 0;
            uint16_t mask = inferSbDcSigCoeffFlag ? 0x7fff : 0xffff;

            bool const isLastSubBlock = (i == lastSubBlock);
            auto *scan = rasterC + (isLastSubBlock ? lastScanPos - 1 : 15);
            if (isLastSubBlock)
            {
                mask >>= 16 - lastScanPos;
                sigCoeffFlags >>= 16 - lastScanPos;
            }

            // loop for sig_coeff_flag
            for (; mask; mask >>= 1)
            {
                auto const ctxInc = ctxIncSigCoeffFlagCalculator.ctxInc(*scan--);
                h(EncodeDecision<sig_coeff_flag>(sigCoeffFlags & 1, ctxInc));
                sigCoeffFlags >>= 1;
            }

            // coeff_abs_level_greater1_flag
            int numGreater1Flag = 0;
            uint16_t hideSignMask = 0;
            int greater2 = -1;
#ifdef WIN32
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif
            CodedData::Type * RESTRICT pAbsCoeff = subBlock.remaining();

            int const ctxSet
                = ((i != 0 && rc->cIdx == 0) ? 2 : 0)
                | (((greater1Ctx > 0 && lastGreater1Flag) || greater1Ctx == 0) ? 1 : 0);

            greater1Ctx = 1;

            int const ctxIncGreater1Base = (ctxSet * 4) + (rc->cIdx ? 16 : 0);

            uint16_t greater1Flags = subBlock.greater1Flags();
            sigCoeffFlags = sigCoeffFlagsActual;
            do
            {
                if (sigCoeffFlags & 1)
                {
                    int greater1Flag = greater1Flags & 1;
                    h(EncodeDecision<coeff_abs_level_greater1_flag>(greater1Flag, ctxIncGreater1Base + greater1Ctx));
                    if (greater1Ctx > 0)
                        lastGreater1Flag = greater1Flag;
                    if (greater1Flag && greater2 < 0)
                    {
                        greater2 = *pAbsCoeff > 2 ? 1 : 0;
                        // if estimating rate (i.e. not actually writing bitstream), we can process coeff_abs_level_greater2_flag now
                        if (std::is_same<typename H::Tag, EstimateRate<void>>::value)
                            h(EncodeDecision<coeff_abs_level_greater2_flag>(*pAbsCoeff > 2 ? 1 : 0, ctxSet + (rc->cIdx ? 4 : 0)));
                    }
                    if (++numGreater1Flag == 8)
                        break;
                    pAbsCoeff += greater1Flag;
                    if (lastGreater1Flag)
                        greater1Ctx = 0;
                    else if (greater1Ctx < 3)
                        ++greater1Ctx;
                }
                greater1Flags >>= 1;
                sigCoeffFlags >>= 1;
            } while (sigCoeffFlags);

            // coeff_abs_level_greater2_flag
            if (!std::is_same<typename H::Tag, EstimateRate<void>>::value && greater2 >= 0)
                h(EncodeDecision<coeff_abs_level_greater2_flag>(greater2, ctxSet + (rc->cIdx ? 4 : 0)));

            // coeff_sign_flag
            uint16_t signFlagPresent = sigCoeffFlags = sigCoeffFlagsActual;
            bool signHidden = false;
            if (!h[cu_transquant_bypass_flag()] && h[sign_data_hiding_enabled_flag()] && (sigCoeffFlags & 0x1fff))
            {
                // derivation of firstSigScanPos and lastSigScanPos benefits from x86 clz and ctz
#if defined(_MSC_VER) || defined(__GNUC__)
                if (0)
                {
                    unsigned long forward, reverse;
#if defined(_MSC_VER) 
                    _BitScanForward(&forward, sigCoeffFlags);
                    _BitScanReverse(&reverse, sigCoeffFlags);
#else
                    forward = __builtin_ctz(sigCoeffFlags);
                    reverse = __builtin_clz(sigCoeffFlags) ^ 31;
#endif
                    signHidden = reverse - forward > 3;
                    if (std::is_same<typename H::Tag, Write<void>>::value)
                        if (signHidden)
                            signFlagPresent &= ~(1 << reverse);
                }
                else
#endif
                {
                    int firstSigScanPos = 0;
                    while (!(sigCoeffFlags & (1 << (15 - firstSigScanPos))))
                        ++firstSigScanPos;
                    int lastSigScanPos = 15;
                    while (!(sigCoeffFlags & (1 << (15 - lastSigScanPos))))
                        --lastSigScanPos;
                    signHidden = lastSigScanPos - firstSigScanPos > 3;
                    if (std::is_same<typename H::Tag, Write<void>>::value)
                        if (signHidden)
                            signFlagPresent &= ~(0x8000 >> firstSigScanPos);
                }
            }

            if (std::is_same<typename H::Tag, EstimateRate<void>>::value)
            {
                StateEstimateRate *stateEstimateRate = h;
                int bits;
#if defined (_MSC_VER) || defined( __GCC__)
                if (havePopcnt)
#if defined (_MSC_VER)
                    bits = __popcnt16(signFlagPresent);
#else
                    bits = __builtin_popcount(signFlagPresent);
#endif
                else
#endif
                {
                    uint16_t x = signFlagPresent;
                    x = (x & 0x5555) + ((x >> 1) & 0x5555);
                    x = (x & 0x3333) + ((x >> 2) & 0x3333);
                    x = (x & 0x0707) + ((x >> 4) & 0x0707);
                    bits = (x & 0xf) + (x >> 8);
                }
                bits -= signHidden;
                stateEstimateRate->rate += Cost::make(bits, 0);
            }
            else
            {
                auto signFlags = subBlock.coeffSignFlags();
                do
                {
                    if (signFlagPresent & 1)
                        h(EncodeBypass<coeff_sign_flag>(signFlags & 1));
                    signFlags >>= 1;
                    signFlagPresent >>= 1;
                } while (signFlagPresent);
            }

            // coeff_abs_level_remaining
            pAbsCoeff = subBlock.remaining();
            int cRiceParam = 0;
            int countdown1 = 0b100111;
            int countdown2 = 0b011110;
            sigCoeffFlags = sigCoeffFlagsActual;
            greater1Flags = subBlock.greater1Flags();
            do
            {
                if (sigCoeffFlags & 1)
                {
                    auto greater1 = greater1Flags & 1;
                    int absCoeff = greater1 ? *pAbsCoeff : 1;

                    int base = countdown1 >> 4;
                    base += (base + countdown2) >> 5;

                    --countdown1;
                    countdown2 -= greater1;
                    pAbsCoeff += greater1;

                    int const remaining = absCoeff - base;
                    if (remaining >= 0)
                        Invoke<typename H::Tag, Element<coeff_abs_level_remaining, ae>>::template write(remaining, absCoeff, cRiceParam, h);
                }
                greater1Flags >>= 1;
                sigCoeffFlags >>= 1;
            } while (sigCoeffFlags);

            if (!codedSubBlockFlag)
                break;

            // move to next subblock
            subBlock.init(pAbsCoeff);
        }
    } while (--i >= 0);

    // move to next residual_coding
    residual.init(subBlock);
}

#endif