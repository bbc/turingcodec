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

#ifndef INCLUDED_CabacWriter_h
#define INCLUDED_CabacWriter_h

#include "Handlers.h"
#include "Global.h"

#include <vector>
#include <cstdint>


struct BitWriter
{
    std::vector<uint8_t> *data;
    std::vector<uint8_t> standby;
    int shift;

    BitWriter(std::vector<uint8_t> *buffer = 0)
        :
        shift(0),
        data(buffer ? buffer : &this->standby)
    {
    }
    inline void writeByte(int byte)
    {
        assert(this->shift == 0);
        this->data->push_back(byte);
    }
    inline void writeBit(int bit)
    {
        if (!this->shift)
        {
            this->data->push_back(0);
            this->shift = 8;
        }
        // For some reason, MSVC STL libraries don't like back() here.
        (*this->data)[this->data->size() - 1] |= (bit << --this->shift);
    }
    inline void writeBits(int n, std::uint64_t value)
    {
        std::uint64_t valueMask = 1;
        valueMask <<= n;
        valueMask >>= 1;
        for (; valueMask; valueMask >>= 1)
        {
            writeBit((value & valueMask) ? 1 : 0);
        }
    }
    inline size_t pos() const
    {
        return this->data->size();
    }
    static inline void insertEp3Bytes(std::vector<uint8_t> &buffer, size_t pos)
    {
        std::vector<uint8_t> data(buffer.begin() + pos, buffer.end());
        buffer.erase(buffer.begin() + pos, buffer.end());
        int zeroes = 0;
        for (auto i = data.begin(); i != data.end(); buffer.push_back(*i++))
        {
            if (zeroes == 2)
            {
                const bool ep3ByteNeeded = !(*i & ~0x03);
                if (ep3ByteNeeded)
                {
                    buffer.push_back(0x03);
                    zeroes = 0;
                }
            }
            zeroes = (*i == 0x00) ? (zeroes + 1) : 0;
        }
    }
};


struct NalWriter :
    BitWriter
{
};


struct CabacWriter :
    BitWriter
{
    unsigned ivlLow;
    unsigned ivlCurrRange;
    int bitsLeft;
    int numBufferedBytes;
    unsigned bufferedByte;

    CabacWriter(std::vector<uint8_t> &buffer)
        :
        BitWriter(&buffer)
    {
    }

    void restartCabac2()
    {
        this->ivlLow = 0;
        this->ivlCurrRange = 510;
        this->bitsLeft = 23;
        this->numBufferedBytes = 0;
        this->bufferedByte = 0xff;
    }
    void testAndWriteOut()
    {
        if (this->bitsLeft < 12) 
            this->writeOut();
    }
    void writeOut()
    {
        unsigned leadByte = this->ivlLow >> (24 - this->bitsLeft);
        this->bitsLeft += 8;
        this->ivlLow &= 0xffffffffu >> this->bitsLeft;

        if (leadByte == 0xff)
        {
            ++this->numBufferedBytes;
        }
        else
        {
            if (this->numBufferedBytes > 0)
            {
                unsigned carry = leadByte >> 8;
                unsigned byte = this->bufferedByte + carry;
                this->bufferedByte = leadByte & 0xff;
                this->writeByte(byte);
                byte = (0xff + carry) & 0xff;
                while (this->numBufferedBytes > 1)
                {
                    this->writeByte(byte);
                    --this->numBufferedBytes;
                }
            }
            else
            {
                this->numBufferedBytes = 1;
                this->bufferedByte = leadByte;
            }
        }
    }
    void finishCabac()
    {
        if (this->ivlLow >> (32 - this->bitsLeft))
        {
            this->writeByte(this->bufferedByte + 1);
            while (this->numBufferedBytes > 1)
            {
                this->writeByte(0x00);
                --this->numBufferedBytes;
            }
            this->ivlLow -= 1 << (32 - this->bitsLeft);
        }
        else
        {
            if (this->numBufferedBytes > 0)
            {
                this->writeByte(this->bufferedByte);
            }
            while (this->numBufferedBytes > 1)
            {
                this->writeByte(0xff);
                this->numBufferedBytes--;
            }
        }
        this->writeBits(24 - this->bitsLeft, this->ivlLow >> 8);
    }
};



#endif
