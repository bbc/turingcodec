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

// Cache of reconstructed CTU samples. Reconstruction is stored as a list of "pieces", not as a single contiguous planar array.

#ifndef INCLUDED_ReconstructionCache_h
#define INCLUDED_ReconstructionCache_h

#pragma once

#include "Memory.h"
#include <vector>
#include <cstddef>
#include <cinttypes>
#include <cassert>
#include <memory>


// Cache of square blocks of different sizes.
template <typename Sample>
struct ReconstructionCache
{
    void allocateBuffer(int log2SizeMin, int log2SizeMax, int n)
    {
        this->log2SizeMin = log2SizeMin;

        assert(log2SizeMax - log2SizeMin <= 6 - 2);

        int log2Step = 2 * (log2SizeMax - log2SizeMin);
        this->step = 1 << log2Step;

        Memory::allocate(this->buffer, n << 2 * log2SizeMax);

        const int size = (n << 2 * (log2SizeMax - log2SizeMin)) / (8 * sizeof(*this->busyMask));

        Memory::allocate(this->busyMask, size);
        this->busyMaskEnd = this->busyMask + size;

        this->pieces.resize(1ull << (2 * (log2SizeMax - log2SizeMin)));
    }

    void freeBuffer()
    {
        Memory::free(this->buffer);
        Memory::free(this->busyMask);
    }

    void reset()
    {
        std::fill(this->busyMask, this->busyMaskEnd, 0);
        this->piece = &this->pieces.front();
    }

    int16_t findAvailableBlock(int log2Size, int z)
    {
        z >>= 2 * this->log2SizeMin;
        z &= (this->step - 1);
        uint8_t *p = reinterpret_cast<uint8_t *>(this->busyMask);

        switch (log2Size - this->log2SizeMin)
        {
            default:
                ASSERT(0);
            case 0: // 1x1
            {
                const uint8_t mask = 1 << (z & 7);
                while (p[z >> 3] & mask) z += this->step;
                p[z >> 3] |= mask;
                ASSERT(&p[z >> 3] < reinterpret_cast<uint8_t *>(this->busyMaskEnd));
                return z;
            }
            case 1: // 2x2
            {
                ASSERT(z % 4 == 0);
                const uint8_t mask = 0x0f << (z & 7);
                while (p[z >> 3] & mask) z += this->step;
                p[z >> 3] |= mask;
                ASSERT(&p[z >> 3] < reinterpret_cast<uint8_t *>(this->busyMaskEnd));
                return z;
            }
            case 2: // 4x4
            {
                ASSERT(z % 16 == 0);
                while (reinterpret_cast<uint16_t *>(&p[z >> 3])[0]) z += this->step;
                reinterpret_cast<uint16_t *>(&p[z >> 3])[0] = 0xffff;
                ASSERT(&p[z >> 3] < reinterpret_cast<uint8_t *>(this->busyMaskEnd));
                return z;
            }
            case 3: // 8x8
            {
                ASSERT(z % 64 == 0);
                while (reinterpret_cast<uint64_t *>(&p[z >> 3])[0]) z += this->step;
                reinterpret_cast<uint64_t *>(&p[z >> 3])[0] = ~uint64_t(0);
                ASSERT(&p[z >> 3] < reinterpret_cast<uint8_t *>(this->busyMaskEnd));
                return z;
            }
            case 4: // 16x16
            {
                ASSERT(z % 256 == 0);
                while (
                        reinterpret_cast<uint64_t *>(&p[z >> 3])[0] |
                        reinterpret_cast<uint64_t *>(&p[z >> 3])[1] |
                        reinterpret_cast<uint64_t *>(&p[z >> 3])[2] |
                        reinterpret_cast<uint64_t *>(&p[z >> 3])[3])	z += this->step;

                reinterpret_cast<uint64_t *>(&p[z >> 3])[0] = ~uint64_t(0);
                reinterpret_cast<uint64_t *>(&p[z >> 3])[1] = ~uint64_t(0);
                reinterpret_cast<uint64_t *>(&p[z >> 3])[2] = ~uint64_t(0);
                reinterpret_cast<uint64_t *>(&p[z >> 3])[3] = ~uint64_t(0);

                ASSERT(&p[z >> 3] < reinterpret_cast<uint8_t *>(this->busyMaskEnd));
                return z;
            }
        }
    }

    // Object describing the size and location in candidate cache of a square block of reconstructed samples
    // A list of Piece objects is generated during encoding and used to reconstruct the output image after each CTU
    struct Piece
    {
#ifdef _DEBUG
#define DEBUG_PIECES
#endif
#ifdef DEBUG_PIECES
        residual_coding rc;
#endif
        uint16_t log2Size;
        int16_t i;

        void set(const residual_coding &rc, uint16_t log2Size, int16_t i)
        {
#ifdef DEBUG_PIECES
            this->rc = rc;
#endif
            this->log2Size = log2Size;
            this->i = i;
        }
    };

    Piece allocateBlock(int log2Size, int z)
    {
        Piece piece;
        piece.set(residual_coding(), uint16_t(log2Size), this->findAvailableBlock(log2Size, z));
        return piece;
    }

    void freeBlock(Piece piece)
    {
        auto i = piece.i;

        if (i < 0) return;

        switch (piece.log2Size - this->log2SizeMin)
        {
            default:
                ASSERT(0);
            case 0: // 1x1
            {
                uint8_t *p = reinterpret_cast<uint8_t *>(this->busyMask);
                ASSERT(p[i / 8] | (1 << i % 8));
                p[i / 8] &= ~(1 << i % 8);
            }
            break;
            case 1: // 2x2
            {
                ASSERT(i % 4 == 0);
                uint8_t *p = reinterpret_cast<uint8_t *>(this->busyMask);
                ASSERT((p[i / 8] & (0x0f << i % 8)) == (0x0f << i % 8));
                p[i / 8] &= ~(0x0f << i % 8);
            }
            break;
            case 2: // 4x4
            {
                ASSERT(i % 16 == 0);
                uint16_t *p = reinterpret_cast<uint16_t *>(this->busyMask);
                ASSERT(p[i / 16] == uint16_t(~0));
                p[i / 16] = 0;
            }
            break;
            case 3: // 8x8
            {
                ASSERT(i % 64 == 0);
                uint64_t *p = this->busyMask;
                ASSERT(p[i / 64] == ~uint64_t(0));
                p[i / 64] = 0;
            }
            break;
            case 4:
            {
                assert(i % 256 == 0);
                uint64_t *p = this->busyMask;
                ASSERT(p[i / 64 + 0] == ~uint64_t(0));
                p[i / 64 + 0] = 0;
                ASSERT(p[i / 64 + 1] == ~uint64_t(0));
                p[i / 64 + 1] = 0;
                ASSERT(p[i / 64 + 2] == ~uint64_t(0));
                p[i / 64 + 2] = 0;
                ASSERT(p[i / 64 + 3] == ~uint64_t(0));
                p[i / 64 + 3] = 0;
            }
            break;
        }
    }

    int step;

    int log2SizeMin;
    Sample *buffer;
    uint64_t *busyMask;
    uint64_t *busyMaskEnd;

    std::vector<Piece> pieces;
    Piece *piece;

    inline Raster<Sample> get(Piece piece)
    {
        return Raster<Sample>(&this->buffer[piece.i << 4], 1ll << piece.log2Size);
    }

    // recursively copy blocks of reconstructed samples from pieces' list to the reconstructed picture buffer.
    inline void commit(Raster<Sample> dst, int log2Size, Piece *&piece)
    {
        assert(piece->log2Size <= 6);
        assert(piece->log2Size >= 2);

        if (piece->log2Size == log2Size)
        {
            if (piece->i >= 0)
            {
                Raster<Sample> src = this->get(*piece);
                int const nCbS = 1 << log2Size;
                copyBlock<Sample>(dst, src, nCbS, nCbS);
            }
            ++piece;
        }
        else
        {
            assert(piece->log2Size < log2Size);
            --log2Size;
            int const nCbS = 1 << log2Size;
            this->commit(dst.offset(0, 0), log2Size, piece);
            this->commit(dst.offset(nCbS, 0), log2Size, piece);
            this->commit(dst.offset(0, nCbS), log2Size, piece);
            this->commit(dst.offset(nCbS, nCbS), log2Size, piece);
        }
    }
};


template <typename Sample>
struct StateReconstructionCache
{
    template <class H>
    StateReconstructionCache(H &h, int n)
    {
        this->components[0].allocateBuffer(h[MinTbLog2SizeY()], h[CtbLog2SizeY()], n);
        this->components[1].allocateBuffer(h[MinTbLog2SizeY()], h[CtbLog2SizeY()] - 1, n);
        this->components[2].allocateBuffer(h[MinTbLog2SizeY()], h[CtbLog2SizeY()] - 1, n);
    }

    ~StateReconstructionCache()
    {
        this->components[0].freeBuffer();
        this->components[1].freeBuffer();
        this->components[2].freeBuffer();
    }

    void reset()
    {
        this->components[0].reset();
        this->components[1].reset();
        this->components[2].reset();
    }

    ReconstructionCache<Sample> components[3];
};

#endif
