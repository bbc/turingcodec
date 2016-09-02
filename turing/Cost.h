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

// Rate and distortion measurement and estimation

#ifndef INCLUDED_Cost_h
#define INCLUDED_Cost_h

#include "Global.h"
#include "FixedPoint.h"
#include <stdint.h>


typedef struct SumOfAbsoluteTransformedDifferences Satd;

typedef FixedPoint<int64_t, 16> Cost;
typedef FixedPoint<int32_t, 16> Lambda;


namespace RateDistortion {

    struct Metric
    {
        typedef Cost Type;

        Metric()
        {
            this->chromaComponent.set(0, 0);
            this->value = std::numeric_limits<Type>::max();
            this->bits = std::numeric_limits<Type>::max();
            this->distortionY = 0;
            this->distortionC = 0;
        }

        Metric(double bits, int32_t sseLuma, int32_t sseChroma, Lambda lambda)
        {
            this->chromaComponent = Type(lambda * sseChroma); // should optimise to imull ??
            this->bits.set(bits);
            this->value = this->bits + lambda * sseLuma + this->chromaComponent;
            this->distortionY = sseLuma;
            this->distortionC = sseChroma;
        }
        void operator+=(const Metric &other)
            {
            this->chromaComponent += other.chromaComponent;
            this->value += other.value;
            this->bits += other.bits;
            this->distortionY += other.distortionY;
            this->distortionC += other.distortionC;
            }
        static Metric max() { return Metric(); }
        bool operator<(Metric b) const { return this->value < b.value; }

        double getValue() const { return this->value.asDouble(); }
        double getChromaComponent() const { return this->chromaComponent.asDouble(); }
        double getBits() const { return this->bits.asDouble(); }
        double getDistortionY() const { return this->distortionY; }
        double getDistortionC() const { return this->distortionC; }

        Type value;
        Type chromaComponent;
        Type bits;
        int32_t distortionY;
        int32_t distortionC;
    };

}



#endif
