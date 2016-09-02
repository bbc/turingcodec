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
#include "SyntaxCtu.hpp"
#include "Read.hpp"


template void Syntax<struct coding_tree_unit>::go<struct Handler<struct Decode<void>,struct StateSubstream,struct StateReconstructedPicture<unsigned short>,struct RbspState,struct StatePicture,struct StateDecode> >(struct coding_tree_unit const &,struct Handler<struct Decode<void>,struct StateSubstream,struct StateReconstructedPicture<unsigned short>,struct RbspState,struct StatePicture,struct StateDecode> &);
template void Syntax<struct coding_tree_unit>::go<struct Handler<struct Decode<void>,struct StateSubstream,struct StateReconstructedPicture<unsigned char>,struct RbspState,struct StatePicture,struct StateDecode> >(struct coding_tree_unit const &,struct Handler<struct Decode<void>,struct StateSubstream,struct StateReconstructedPicture<unsigned char>,struct RbspState,struct StatePicture,struct StateDecode> &);
