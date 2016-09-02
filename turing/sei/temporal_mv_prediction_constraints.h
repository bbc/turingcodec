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

#include "../SyntaxSei.h"


struct prev_pics_not_used_flag { };
struct no_intra_layer_col_pic_flag { };


template <>
struct Syntax<temporal_mv_prediction_constraints>
{
    template <class H> static void go(temporal_mv_prediction_constraints const &fun, H &h)
    {
        h(prev_pics_not_used_flag(), u(1));
        h(no_intra_layer_col_pic_flag(), u(1));
    }
};


struct TemporalMvPredictionConstraints :
    ValueHolder<prev_pics_not_used_flag>,
    ValueHolder<no_intra_layer_col_pic_flag>
    {
    };


template <class H> void Read<temporal_mv_prediction_constraints>::go(temporal_mv_prediction_constraints f, H &h)
{
    TemporalMvPredictionConstraints temporalMvPredictionConstraints;
    auto hh = h.extend(&temporalMvPredictionConstraints);
    Syntax<temporal_mv_prediction_constraints>::go(f, hh);
}
