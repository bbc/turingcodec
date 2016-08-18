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

#include "Global.h"
#include "Read.h"
#include "StateDecode.h"
#include "Decode.h"
#include "GlobalState.h"
#include "Syntax.h"
#include "SyntaxRbsp.hpp"
#include "Read.hpp"


template void Syntax<struct slice_segment_layer_rbsp>::go<struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> >(struct slice_segment_layer_rbsp const &, struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> &);
template void Syntax<struct access_unit_delimiter_rbsp>::go<struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> >(struct access_unit_delimiter_rbsp const &, struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> &);
template void Syntax<struct end_of_bitstream_rbsp>::go<struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> >(struct end_of_bitstream_rbsp const &, struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> &);
template void Syntax<struct filler_data_rbsp>::go<struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> >(struct filler_data_rbsp const &, struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> &);
template void Syntax<struct video_parameter_set_rbsp>::go<struct Handler<struct Decode<void>, struct HrdArray, struct ProfileTierLevel, struct Vps, struct RbspState, struct StatePicture, struct StateDecode> >(struct video_parameter_set_rbsp const &, struct Handler<struct Decode<void>, struct HrdArray, struct ProfileTierLevel, struct Vps, struct RbspState, struct StatePicture, struct StateDecode> &);
template void Syntax<struct seq_parameter_set_rbsp>::go<struct Handler<struct Decode<void>, struct HrdArray, struct ProfileTierLevel, struct ScalingListState, struct Sps, struct RbspState, struct StatePicture, struct StateDecode> >(struct seq_parameter_set_rbsp const &, struct Handler<struct Decode<void>, struct HrdArray, struct ProfileTierLevel, struct ScalingListState, struct Sps, struct RbspState, struct StatePicture, struct StateDecode> &);
template void Syntax<struct pic_parameter_set_rbsp>::go<struct Handler<struct Decode<void>, struct ScalingListState, struct Pps, struct RbspState, struct StatePicture, struct StateDecode> >(struct pic_parameter_set_rbsp const &, struct Handler<struct Decode<void>, struct ScalingListState, struct Pps, struct RbspState, struct StatePicture, struct StateDecode> &);
template void Syntax<struct end_of_seq_rbsp>::go<struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> >(struct end_of_seq_rbsp const &, struct Handler<struct Decode<void>, struct RbspState, struct StatePicture, struct StateDecode> &);
