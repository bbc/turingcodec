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

#ifndef INCLUDED_Primitive_h
#define INCLUDED_Primitive_h

#include <cassert>


namespace Primitive {


    template <class T>
    struct FunctionType;

    template <class T>
    struct FunctionNumber
    {
        static const int value = 1;
    };

    template <class T, int i>
    struct Unoptimized;

    struct InstructionSetSupport
    {
        bool sse2() const { return true; }
    };

    template <class T, int n>
    struct Initialize;

    template <class T>
    struct Function
    {
        Function(InstructionSetSupport s)
        {
            Initialize<T, FunctionNumber<T>::value>::init(s, this->f);
        }
        typename FunctionType<T>::Type *f[FunctionNumber<T>::value];
    };

    template <class T, class Functions>
    typename FunctionType<T>::Type &get(Functions &functions, int i=0)
    {
        assert(i >= 0 && i < FunctionNumber<T>::value);
        return *functions.Function<T>::f[i];
    }

    static InstructionSetSupport instructionSetSupport;

} // namespace Primitive

#endif
