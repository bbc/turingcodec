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
#include "../md5.h"

struct hash_type { };

enum HashType
{
    MD5,
    CRC,
    CHKSUM,
};

DEFINE_STRUCT_ARITY_2(picture_md5, cIdx, i);
template<> struct ValueType<picture_md5> { typedef std::uint8_t Type; };
template <> struct IsValueArray<picture_md5> : std::true_type {};
template <> struct ValueHolder<picture_md5>
{
    std::array<std::array<ValueType<picture_md5>::Type, 16>, 3> value;
    ValueType<picture_md5>::Type &get(picture_md5 v)
    {
        return value[v.cIdx][v.i];
    }
};

DEFINE_VALUE_ARRAY_1(picture_crc, cIdx, 3);
DEFINE_VALUE_ARRAY_1(picture_checksum, cIdx, 3);

template <>
struct Syntax<decoded_picture_hash>
{
    template <class H> static void go(decoded_picture_hash fun, H &hh);
};

template <class H>
void Syntax<decoded_picture_hash>::go(decoded_picture_hash fun, H &hh)
{
    StateSlice *stateSlice = hh;
    auto h = hh;//.extend(stateSlice, hh);

    if (!h[Active<Sps>()])
    {
        h(Violation("D.3.19", "no active SPS so chroma_format_idc unknown"));
        throw Abort();
    }

    h(hash_type(), u(8));
    for (int cIdx = 0; cIdx < (h[chroma_format_idc()] == 0 ? 1 : 3); cIdx++)
    {
        if (h[hash_type()] == 0)
            for (int i = 0; i < 16; i++)
                h(picture_md5(cIdx, i), b(8));
        else if (h[hash_type()] == 1)
            h(picture_crc(cIdx), u(16));
        else if (h[hash_type()] == 2)
            h(picture_checksum(cIdx), u(32));
    }
}

struct DecodedPictureHash :
    ValueHolder<hash_type>,
    ValueHolder<picture_md5>,
    ValueHolder<picture_crc>,
    ValueHolder<picture_checksum>
    {
        void computeMd5Sum(md5_byte_t digest[16], unsigned char* bufferPtr, intptr_t stride, int height, int width)
        {
            md5_state_t md5Picture;
            md5_init(&md5Picture);
            for (int r = 0; r < height; r++, bufferPtr += stride)
                md5_append(&md5Picture, bufferPtr, width);
            md5_finish(&md5Picture, digest);
        }

        void computeMd5Sum(md5_byte_t digest[16], unsigned short* bufferPtr, intptr_t stride, int height, int width)
        {
            md5_state_t md5Picture;
            md5_init(&md5Picture);
            unsigned char *rowOfPixels = new unsigned char[width << 1];
            for (int r = 0; r < height; r++, bufferPtr += stride)
            {
                int idx = 0;
                for(int c = 0; c < width; c++)
                {
                    rowOfPixels[idx++] = bufferPtr[c] & 0xff;
                    rowOfPixels[idx++] = bufferPtr[c] >> 8;
                }
                md5_append(&md5Picture, rowOfPixels, width << 1);
            }
            md5_finish(&md5Picture, digest);
            delete[] rowOfPixels;
        }

        template<class Ptr>
        int computeCrc(Ptr* bufferPtr, intptr_t stride, int height, int width, int bitDepth)
        {
            int crcMsb;
            int bitVal;
            int crcVal = 0xffff;
            int bitIdx;
            for (int y = 0; y < height; y++, bufferPtr += stride)
            {
                for (int x = 0; x < width; x++)
                {
                    // take CRC of first pictureData byte
                    for (bitIdx = 0; bitIdx < 8; bitIdx++)
                    {
                        crcMsb = (crcVal >> 15) & 1;
                        bitVal = (bufferPtr[x] >> (7 - bitIdx)) & 1;
                        crcVal = (((crcVal << 1) + bitVal) & 0xffff) ^ (crcMsb * 0x1021);
                    }
                    if(bitDepth > 8)
                    {
                      for(bitIdx=0; bitIdx<8; bitIdx++)
                      {
                        crcMsb = (crcVal >> 15) & 1;
                        bitVal = (bufferPtr[x] >> (15 - bitIdx)) & 1;
                        crcVal = (((crcVal << 1) + bitVal) & 0xffff) ^ (crcMsb * 0x1021);
                      }
                    }
                }
            }
            for (bitIdx = 0; bitIdx < 16; bitIdx++)
            {
                crcMsb = (crcVal >> 15) & 1;
                crcVal = ((crcVal << 1) & 0xffff) ^ (crcMsb * 0x1021);
            }
            return crcVal;
        }

        template<class Ptr>
        int computeCheckSum(Ptr* bufferPtr, intptr_t stride, int height, int width, int bitDepth)
        {
            int checksum = 0;
            unsigned char xor_mask;

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    xor_mask = (x & 0xff) ^ (y & 0xff) ^ (x >> 8) ^ (y >> 8);
                    checksum = (checksum + ((bufferPtr[y*stride + x] & 0xff) ^ xor_mask)) & 0xffffffff;
                    if( bitDepth > 8 )
                    {
                        checksum = (checksum + ((static_cast<unsigned>(bufferPtr[y*stride + x]) >> 8) ^ xor_mask )) & 0xffffffff;
                    }
                }
            }

            return checksum;
        }
    };


template <class H> void Read<decoded_picture_hash>::go(decoded_picture_hash f, H &h)
{
    DecodedPictureHash decodedPictureHash;
    auto h3 = h.extend(&decodedPictureHash);

    Syntax<decoded_picture_hash>::go(f, h3);
}


