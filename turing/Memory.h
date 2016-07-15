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

#ifndef INCLUDED_Memory_h
#define INCLUDED_Memory_h

#pragma once

// Memory handling and alignment


#ifdef WIN32

#define ALIGN(n, T, v) \
        __declspec(align(n)) T v

#endif

#ifdef __GNUC__

#define ALIGN(n, T, v) \
        T v __attribute__((aligned(n)))

#endif


namespace Memory {

    template <class T>
    void allocate(T *&t, size_t size, std::intptr_t alignment = 32)
    {
        size *= sizeof(T);

        void *p;

#ifdef WIN32
        p = _aligned_malloc(size, alignment);
#else
        if (posix_memalign((void **)&p, alignment, size)) p = 0;
#endif

        assert(p);

        t = static_cast<T *>(p);
    }

    template <class T>
    void free(T *t)
    {
#ifdef WIN32
        _aligned_free(t);
#else
        ::free(t);
#endif
    }

}

#endif
