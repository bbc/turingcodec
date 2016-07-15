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

#ifndef INCLUDED_StreamReader_h
#define INCLUDED_StreamReader_h

#pragma once


#include "CodeCache.h"
#include "Global.h"
#include <iostream>


struct ExceptionOverrun { };


template <class Derived>
struct Bookmarkable
{
    struct Bookmark :
        Derived::State
        {
            Bookmark(Derived &s) :
            Derived::State(s.state),
            s(s)
            {
            }
            ~Bookmark()
            {
                this->s.state = *static_cast<typename Derived::State *>(this);
            }
            Derived &s;
        };
};


// Reads bitstream from a std::istream
struct StreamReader :
    Bookmarkable<StreamReader>
    {
        std::unique_ptr<CodeCache> codeCache;
        std::vector<CodeCache::Code>::iterator codeIt;

        void readRbspData(std::vector<uint8_t> &rbspData, size_t endPositionNalUnit, std::vector<size_t> &eb3pPositions)
        {
            size_t i = 0;
            this->adjustCodeIt();

            eb3pPositions.clear();

            while (this->state.pos < endPositionNalUnit)
            {
                size_t nBytes = size_t(endPositionNalUnit - this->state.pos);

                while (this->codeIt != this->codeCache->codes.end())
                {
                    if (this->codeIt->byte == 0x03)
                    {
                        const size_t nBytesToNextEp3b = size_t(this->codeIt->position + this->codeIt->numberOfZeroes) - this->state.pos;

                        if (nBytes > nBytesToNextEp3b)
                        {
                            nBytes = nBytesToNextEp3b;
                            this->codeIt++;
                        }

                        break;
                    }
                    this->codeIt++;
                }

                is->clear();
                is->seekg(this->state.pos);
                is->read(reinterpret_cast<char *>(&rbspData[i]), nBytes);

                i += static_cast<size_t>(nBytes);
                eb3pPositions.push_back(i);

                this->state.pos += nBytes;
                if (this->state.pos < endPositionNalUnit)++this->state.pos; // escaped byte
            }
            rbspData.resize(i);

            // preload the next byte to leave the class in expected state
            is->clear();
            is->seekg(this->state.pos);
            is->read(reinterpret_cast<char *>(&this->state.byte), 1);
        }

        // move codeIt to the correct position in codeCache according to current stream position
        // after this call, codeIt will indicate the next code following  stream position pos
        void adjustCodeIt(bool truncate = false)
        {
            const std::streamsize pos = static_cast<std::streamsize>(this->state.pos);

            const bool byteAligned = this->state.mask == 0x80;
            assert(byteAligned);

            if (truncate)
            {
                // if pos is mid-code-pattern, set codeIt to that code
                while (this->codeIt > this->codeCache->codes.begin() && ((this->codeIt - 1)->position + std::max<std::streamsize>(0, (this->codeIt - 1)->numberOfZeroes - 3)) >= pos) --this->codeIt;
                while (this->codeIt < this->codeCache->codes.end() && (this->codeIt->position + std::max<std::streamsize>(0, this->codeIt->numberOfZeroes - 3)) < pos) ++this->codeIt;
            }
            else
            {
                // if pos is mid-code-pattern, set codeIt to the next code
                while (this->codeIt > this->codeCache->codes.begin() && (this->codeIt - 1)->position >= pos)--this->codeIt;
                while (this->codeIt < this->codeCache->codes.end() && this->codeIt->position < pos) ++this->codeIt;
            }
        }

        size_t numBytesInNalUnit()
        {
            StreamReader &stream = *this;
            StreamReader::Bookmark mark(stream);

            try
            {
                // Fast and simple parsing of the first part of byte_stream_nal_unit()

                if (stream.nextBytes(4) == 0x00000001)
                {
                    stream.readBytes(1);
                }

                if (stream.readBytes(3) != 0x000001) return 0;
            }
            catch (ExceptionOverrun &)
            {
                return 0;
            }

            auto const nalUnitStartPosition = stream.state.pos;

            for (this->adjustCodeIt(); this->codeIt != this->codeCache->codes.end(); ++this->codeIt)
            {
                if (this->codeIt->byte == 0x01 || this->codeIt->byte == CodeCache::eos)
                {
                    ASSERT((size_t)this->codeIt->position >= nalUnitStartPosition);
                    return static_cast<int>(this->codeIt->position - nalUnitStartPosition);
                }
            }
            ASSERT(!"this should not occur due to end-of-stream codes entry with byte=0xff");
            ASSERT(stream.len > nalUnitStartPosition);
            return stream.len - nalUnitStartPosition;
        }

        void open(std::istream &is, size_t truncate, bool useCodeCache = true, size_t nCodes = 0, CodeCache::Code *codes = 0)
        {
            if (!is.good()) throw std::exception();
            this->is = &is;

            if (useCodeCache)
            {
                this->codeCache.reset(new CodeCache());

                if (!codes)
                {
                    this->len = this->codeCache->load(is, truncate);
                }
                else
                {
                    this->codeCache->codes.insert(this->codeCache->codes.begin(), codes, codes + nCodes);

                    this->is->seekg(0, std::ios::end);
                    this->len = static_cast<size_t>(is.tellg());
                }
                this->codeIt = this->codeCache->codes.begin();
            }
            else
            {
                this->is->seekg(0, std::ios::end);
                this->len = static_cast<size_t>(is.tellg());
            }

            this->state.pos = -1;
            this->state.mask = 0;
            this->bit();
        }

        std::istream *is;
        size_t len;

        struct State
        {
            unsigned char mask;
            unsigned char byte;
            size_t pos;
        };

        State state;

        int bit()
        {
            this->checkOverrun();
            char val = this->state.byte & this->state.mask;
            this->state.mask >>= 1;
            if (!this->state.mask)
            {
                this->state.mask = 0x80;
                is->clear();
                is->seekg(++this->state.pos);
                is->read(reinterpret_cast<char *>(&this->state.byte), 1);
            }
            return val ? 1 : 0;
        }
        int byte()
        {
            assert(this->state.mask == 0x80);
            this->checkOverrun();
            const int val = this->state.byte;
            is->clear();
            is->seekg(++this->state.pos);
            is->read(reinterpret_cast<char *>(&this->state.byte), 1);
            return val;
        }
        int readBytes(int n)
        {
            int v = 0;
            while (n--)
            {
                v <<= 8;
                v |= this->byte();
            }
            return v;
        }
        unsigned nextBytes(int n)
        {
            if (this->state.pos + n >= this->len) return ~0;

            Bookmark mark(*this);
            unsigned v = 0;
            while (n--)
            {
                v <<= 8;
                v |= this->byte();
            }
            return v;
        }
    private:
        void checkOverrun() const
        {
            if (this->state.pos == this->len)
            {
                throw ExceptionOverrun();
            }
        }
    };

template <class T, class S>
T bitPos(S &s);

template <class T>
static T bitPos(StreamReader& streamReader)
{
    T pos = 0;
    unsigned char m = 0x80;
    while (m != streamReader.state.mask)
    {
        m >>= 1;
        ++pos;
    }
    return pos + static_cast<T>(8 * streamReader.state.pos);
};

template <class T = double>
struct Bits;

template <class T = size_t>
struct Bytes;

template <class Units>
struct Position { };

// Type Stream is used to identify the current stream. It may be NAL stream or RBSP stream (after emulation prevention bytes removal)
struct Stream { };

template <class S>
struct Access<Stream, S, typename std::enable_if<std::is_base_of<StreamReader, S>::value>::type>
{
    typedef StreamReader SetType;
    typedef StreamReader &Type;
    static Type get(Stream, StreamReader &s)
    {
        return s;
    }
};

template <class S>
struct Access<Position<Bytes<size_t>>, S, typename std::enable_if<std::is_base_of<StreamReader, S>::value>::type>
{
    typedef size_t Type;
    static Type get(Position<Bytes<size_t>>, StreamReader &s)
    {
        return bitPos<size_t>(s) / 8;
    }
};

static size_t bitLen(const StreamReader& streamReader)
{
    return 8 * streamReader.len;
};

template <class S>
struct Access<more_data_in_byte_stream, S, typename std::enable_if<std::is_base_of<StreamReader, S>::value>::type>
{
    typedef bool Type;
    static Type get(more_data_in_byte_stream, StreamReader &streamReader)
    {
        return streamReader.state.mask != 0x80 || (bitPos<int>(streamReader) != bitLen(streamReader));
    }
};

struct BitReader :
    Bookmarkable<BitReader>
    {
        BitReader() { }
        template <class T>
        BitReader(const T *begin, const T *end)
        :
        begin(reinterpret_cast<const unsigned char*>(begin)),
        end(reinterpret_cast<const unsigned char*>(end))
        {
            this->state.p = this->begin;
            this->state.mask = 0x80;
            //	this->ivlCurrRange = 510;
            //	this->bitsNeeded = -8;
        }
        const unsigned char *begin;
        const unsigned char *end;
        struct State
        {
            const unsigned char *p;
            int mask;
        };
        State state;
        int bit()
        {
            this->checkOverrun();
            int val = *this->state.p & this->state.mask;
            this->state.mask >>= 1;
            if (!this->state.mask)
            {
                this->state.mask = 0x80;
                ++this->state.p;
            }
            return val ? 1 : 0;
        }
        unsigned byte()
        {
            assert(this->state.mask == 0x80);
            this->checkOverrun();
            return *this->state.p++;
        }
        int readBytes(int n)
        {
            int v = 0;
            while (n--)
            {
                if (this->state.p == this->end)
                {
                    throw ExceptionOverrun();
                }
                v <<= 8;
                v |= this->byte();
            }
            return v;
        }
        unsigned nextBytes(int n)
        {
            if (this->state.p + n >= this->end) return ~0;

            Bookmark mark(*this);
            unsigned v = 0;
            while (n--)
            {
                v <<= 8;
                v |= this->byte();
            }
            return v;
        }

        void rewind(int n)
        {
            assert(sizeof(this->state.mask) > 1); // current rewind implementation will not work with char-size mask
            while (n--)
            {
                this->state.mask <<= 1;
                if (this->state.mask == 0x100)
                {
                    this->state.mask = 0x1;
                    --this->state.p;
                    assert(this->state.p >= this->begin);
                }
            }
        }

    private:
        void checkOverrun()
        {
            if (this->state.p == this->end) throw ExceptionOverrun();
        }
    };

static size_t bitLen(const BitReader& bitReader)
{
    return 8 * (bitReader.end - bitReader.begin);
};


static void seekStream(BitReader& bitReader, size_t bitPos)
{
    bitReader.state.p = bitReader.begin + bitPos / 8;
    bitReader.state.mask = 0x80 >> (bitPos % 8);
}


static void seekStream(StreamReader& streamReader, size_t bitPos)
{
    streamReader.state.pos = bitPos / 8;
    streamReader.state.mask = 0x80 >> (bitPos % 8);
    streamReader.is->seekg(streamReader.state.pos);
    streamReader.is->read(reinterpret_cast<char *>(&streamReader.state.byte), 1);
}


template <class H>
void seek(H &h, size_t bitPos)
{
    seekStream(h[::Stream()], bitPos);
}


static bool endOfStream(StreamReader& streamReader)
{
    return streamReader.state.pos >= streamReader.len;
}


template <class S>
struct Access<byte_aligned, S, typename std::enable_if<std::is_base_of<BitReader, S>::value>::type>
{
    typedef bool Type;
    static Type get(byte_aligned, BitReader &bitReader)
    {
        return bitReader.state.mask == 0x80;
    }
};


template <class S>
struct Access<more_rbsp_data, S, typename std::enable_if<std::is_base_of<BitReader, S>::value>::type>
{
    typedef bool Type;
    static Type get(more_rbsp_data, BitReader &bitReader)
    {
        if (bitReader.state.p == bitReader.end)
        {
            return false;
        }
        else if (bitReader.state.p + 1 == bitReader.end)
        {
            const int lastByte = *bitReader.state.p;
            return !!(lastByte & (bitReader.state.mask - 1));
        }
        else
        {
            return true;
        }
    }
};


template <class S>
struct Access<more_rbsp_trailing_data, S, typename std::enable_if<std::is_base_of<BitReader, S>::value>::type>
{
    typedef bool Type;
    static Type get(more_rbsp_trailing_data, BitReader &bitReader)
    {
        return bitReader.state.p != bitReader.end;
    }
};


template <class S>
struct Access<more_data_in_payload, S, typename std::enable_if<std::is_base_of<BitReader, S>::value>::type>
{
    typedef bool Type;
    static Type get(more_data_in_payload, BitReader &bitReader)
    {
        return bitReader.state.mask != 0x80 || bitReader.state.p < bitReader.end;
    }
};


template <class S>
struct Access<payload_extension_present, S, typename std::enable_if<std::is_base_of<BitReader, S>::value>::type>
{
    typedef bool Type;
    static Type get(payload_extension_present, BitReader &bitReader)
    {
        if (bitReader.state.p + 1 == bitReader.end)
        {
            int pattern = (bitReader.state.mask << 1) - 1;
            return ((pattern & *bitReader.state.p) != bitReader.state.mask);
        }
        else
        {
            return true;
        }
    }
};

#endif
