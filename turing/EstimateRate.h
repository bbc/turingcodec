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

// Bitrate estimation


#ifndef INCLUDED_EstimateRate_h
#define INCLUDED_EstimateRate_h

#include "Write.h"
#include "Binarization.h"
#include "StateEncode.h"


// EstimateRate estimates bitrates to StateEstimateRate and *does not* update contexts' states.
template <class F> struct EstimateRateNoContextUpdate : Write<F> { };
template <> struct EstimateRateNoContextUpdate<void>;


// EstimateRate estimates bitrate to StateEstimateRate and *does* update contexts' states.
template <class F> struct EstimateRate : EstimateRateNoContextUpdate<F> { };
template <> struct EstimateRate<void>;
// MN: to be fixed
template <> struct EstimateRate<Element<cu_qp_delta_abs, ae>> : Null<Element<cu_qp_delta_abs, ae>>{};
template <> struct EstimateRate<Element<cu_qp_delta_sign_flag, ae>> : Null<Element<cu_qp_delta_sign_flag, ae>>{};


template <class V>
struct EstimateRateBin
{
    ContextModel contextModel;

    template <class H> EstimateRateBin(H &h, int ctxInc)
            {
        Contexts *contexts = h;
        this->contextModel = contexts->template get<typename Context<V>::Type>(ctxInc);
            }

    Cost rate(int binVal)
    {
        return entropyEstimate.table[this->contextModel.state ^ binVal];
    }
};

template <class V>
struct EstimateRateNoContextUpdate<EncodeDecision<V>>
{
    template <class H> static void go(EncodeDecision<V> encodeDecision, H &h)
    {
        Contexts *contexts = h;

        ContextModel &contextModel = contexts->template get<typename Context<V>::Type>(encodeDecision.ctxInc);

        StateEstimateRate *stateEstimateRate = h;
        Cost rate = entropyEstimate.table[contextModel.state ^ encodeDecision.binVal];
        stateEstimateRate->rate += rate;
    }
};


template <class V>
struct EstimateRate<EncodeDecision<V>>
{
    template <class H> static void go(EncodeDecision<V> encodeDecision, H &h)
    {
        StateEstimateRate *stateEstimateRate = h;
        Contexts *contexts = h;
        ContextModel &contextModel = contexts->template get<typename Context<V>::Type>(encodeDecision.ctxInc);

        auto rate = measureEncodeDecision(contextModel, encodeDecision.binVal);
        stateEstimateRate->rate += rate;
    }
};


template <class V>
struct EstimateRate<EncodeTerminate<V>>
{
    template <class H> static void go(EncodeTerminate<V> encodeTerminate, H &h)
    {
        // ignore termination bins during rate estimation
    }
};


template <class V>
struct EstimateRateNoContextUpdate<EncodeBypass<V>>
{
    template <class H> static void go(EncodeBypass<V>, H &h)
    {
        StateEstimateRate *stateEstimateRate = h;
        ++stateEstimateRate->rate;
    }
};


template <class F> struct EstimateRateLuma : EstimateRate<F> { };
template <> struct EstimateRateLuma<Element<cbf_cb, ae>> : Null<Element<cbf_cb, ae>>{};
template <> struct EstimateRateLuma<Element<cbf_cr, ae>> : Null<Element<cbf_cr, ae>>{};
template <> struct EstimateRateLuma<IfCbf<cbf_cb, residual_coding>> : Null<IfCbf<cbf_cb, residual_coding>>{};
template <> struct EstimateRateLuma<IfCbf<cbf_cr, residual_coding>> : Null<IfCbf<cbf_cr, residual_coding>>{};
template <> struct EstimateRateLuma<Element<intra_chroma_pred_mode, ae>> : Null<Element<intra_chroma_pred_mode, ae>>{};
// MN: to be fixed
//template <> struct EstimateRate<Element<cu_qp_delta_abs, ae>> : Null<Element<cu_qp_delta_abs, ae>>{};
//template <> struct EstimateRate<Element<cu_qp_delta_sign_flag, ae>> : Null<Element<cu_qp_delta_sign_flag, ae>>{};

template <class F> struct EstimateRateChroma : EstimateRate<F> { };
template <> struct EstimateRateChroma<Element<cbf_luma, ae>> : Null<Element<cbf_luma, ae>>{};
template <> struct EstimateRateChroma<IfCbf<cbf_luma, residual_coding>> : Null<IfCbf<cbf_luma, residual_coding>>{};
template <> struct EstimateRateChroma<Element<split_transform_flag, ae>> : Null<Element<split_transform_flag, ae>>{};
// MN: to be fixed
//template <> struct EstimateRate<Element<cu_qp_delta_abs, ae>> : Null<Element<cu_qp_delta_abs, ae>>{};
//template <> struct EstimateRate<Element<cu_qp_delta_sign_flag, ae>> : Null<Element<cu_qp_delta_sign_flag, ae>>{};


template <class H, class F>
Cost estimateRateNoContextUpdate(H &h, const F &f)
{
    StateEstimateRate *stateEstimateRate = h;
    Cost const before2 = stateEstimateRate->rate;
    auto m = h.template change<EstimateRateNoContextUpdate<void>>();
    m(f);
    Cost rate = stateEstimateRate->rate - before2;
    stateEstimateRate->rate = before2;
    return rate;
}


#endif
