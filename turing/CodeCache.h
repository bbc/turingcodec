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

#ifndef INCLUDED_CodeCache_H
#define INCLUDED_CodeCache_H

// Readr management/caching of start_code_prefix_one_3bytes and nal_unit_header() information.


#include <vector>
#include <stdint.h>
#include <iostream>


struct CodeCache
{
    struct Code
    {
        std::streamsize position;
        std::streamsize numberOfZeroes;
        uint8_t byte; // if in range 0 to 3, this indicates pattern 00 00 <byte>. 0xff indicates end of stream.
        uint8_t nextByte; // the value in pattern 00 00 <byte> <nextByte>

        bool endOfStream() const { return this->byte == 0xff; }
    };

    std::vector<Code> codes;

    // magic number - entry represents the end of the stream, not a startcode
    static const uint8_t eos = 0xff;

    size_t load(std::istream &is, size_t truncate)
    {
        std::vector<uint8_t> buffer(4096);
        size_t pos = 0;
        size_t consecutiveZeroes = 0;

        // loop over chunks of bitstream
        while (true)
        {
            // read bitstream from the istream
            is.read(reinterpret_cast<char *>(&buffer[0]), buffer.size());
            size_t n = static_cast<size_t>(is.gcount());

            if (truncate && pos + n > truncate) n = truncate - pos;
            const bool last = n != buffer.size();

            is.clear();

            // loop over all bytes in the chunk
            for (size_t i = 0; i < n; ++i)
            {
                const auto &byte = buffer[i];
                if (byte == 0x00)
                {
                    ++consecutiveZeroes;
                }
                else
                {
                    if ((consecutiveZeroes >= 2 && byte >= 0x01 && byte <= 0x03)
                            || (consecutiveZeroes >= 3))
                    {
                        // append a new entry to the list of codes
                        Code code;
                        code.numberOfZeroes = consecutiveZeroes;
                        code.position = pos + i - consecutiveZeroes;
                        code.byte = byte;
                        code.nextByte = i + 1 < buffer.size() ? buffer[i + 1] : 0;
                        this->codes.push_back(code);
                    }
                    consecutiveZeroes = 0;
                }
            }
            pos += n;
            if (last)
            {
                Code code;
                code.numberOfZeroes = consecutiveZeroes;
                code.position = pos - consecutiveZeroes;
                code.byte = eos;
                code.nextByte = 0;
                this->codes.push_back(code);
                break;
            }
        }


        return pos;
    }
};

#endif
