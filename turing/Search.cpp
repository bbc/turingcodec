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

// Explicit template instantiations for Search.hpp
#include "Search.hpp"

template void Search<coding_quadtree>::go<Handler<Search<Mode<0>>, Candidate<uint8_t>, StateEncodeSubstream<uint8_t>, StateEncodePicture2<uint8_t>, StateEncode> >(coding_quadtree const &, Handler<Search<Mode<0>>, Candidate<uint8_t>, StateEncodeSubstream<uint8_t>, StateEncodePicture2<uint8_t>, StateEncode> &);
template void Search<coding_quadtree>::go<Handler<Search<Mode<0>>, Candidate<uint16_t>, StateEncodeSubstream<uint16_t>, StateEncodePicture2<uint16_t>, StateEncode> >(coding_quadtree const &, Handler<Search<Mode<0>>, Candidate<uint16_t>, StateEncodeSubstream<uint16_t>, StateEncodePicture2<uint16_t>, StateEncode> &);
template void Search<sao>::go<Handler<Search<void>, Candidate<uint8_t>, StateEncodeSubstream<uint8_t>, StateEncodePicture2<uint8_t>, StateEncode> >(sao const&, Handler<Search<void>, Candidate<uint8_t>, StateEncodeSubstream<uint8_t>, StateEncodePicture2<uint8_t>, StateEncode>&);
template void Search<sao>::go<Handler<Search<void>, Candidate<uint16_t>, StateEncodeSubstream<uint16_t>, StateEncodePicture2<uint16_t>, StateEncode> >(sao const&, Handler<Search<void>, Candidate<uint16_t>, StateEncodeSubstream<uint16_t>, StateEncodePicture2<uint16_t>, StateEncode>&);
