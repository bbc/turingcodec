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

#ifndef INCLUDED_StateFunctionTable_h
#define INCLUDED_StateFunctionTable_h

#pragma once

#include "havoc/pred_inter.h"
#include "havoc/pred_intra.h"
#include "havoc/quantize.h"
#include "havoc/transform.h"
#include "havoc/sad.h"
#include "havoc/ssd.h"
#include "havoc/hadamard.h"
#include "havoc/havoc.h"


 // State comprising tables of pointers for optimised functions.
struct StateFunctionTables :
    HavocTablePredUni<uint16_t>,
    HavocTablePredBi<uint16_t>,
    havoc::table_inverse_transform_add<uint16_t>,
    havoc::table_inverse_transform_add<uint8_t>,
    havoc::table_inverse_transform,
    HavocTablePredUni<uint8_t>,
    HavocTablePredBi<uint8_t>,
    havoc_table_quantize_inverse,
    havoc_table_quantize,
    havoc_table_quantize_reconstruct,
    havoc::intra::Table<uint16_t>,
    havoc::intra::Table<uint8_t>,
    havoc::table_transform<8>, // encode only
    havoc::table_transform<10>, // encode only
    havoc_table_ssd<uint8_t>, // encode only
    havoc_table_ssd<uint16_t>, // encode only
    havoc_table_sad<uint8_t>, // encode only
    havoc_table_sad<uint16_t>, // encode only
    havoc_table_sad_multiref<uint8_t>, // encode only
    havoc_table_sad_multiref<uint16_t>, // encode only
    havoc_table_hadamard_satd<uint8_t>, // encode only
    havoc_table_hadamard_satd<uint16_t>, // encode only
    havoc::TableSubtractBi<uint8_t>, // encode only
    havoc::TableSubtractBi<uint16_t> // encode only
{
    StateFunctionTables(bool encoder, havoc_instruction_set mask = havoc_instruction_set_support()) 
        :
        instruction_set_support(mask)
    {
        this->code = havoc_new_code(mask, 12000000);
        havocPopulatePredUni<uint8_t>(this, this->code);
        havocPopulatePredUni<uint16_t>(this, this->code);
        havocPopulatePredBi<uint8_t>(this, this->code);
        havocPopulatePredBi<uint16_t>(this, this->code);
        havoc_populate_quantize_inverse(this, this->code);
        havoc_populate_quantize(this, this->code);
        havoc_populate_quantize_reconstruct(this, this->code);
        havoc::populate_inverse_transform(this, this->code, encoder ? 1 : 0);
        havoc::populate_inverse_transform_add<uint8_t>(this, this->code, encoder ? 1 : 0);
        havoc::populate_inverse_transform_add<uint16_t>(this, this->code, encoder ? 1 : 0);
        this->havoc::intra::Table<uint16_t>::populate(this->code);
        this->havoc::intra::Table<uint8_t>::populate(this->code);
        havoc::populate_transform<8>(this, this->code); // encode
        havoc::populate_transform<10>(this, this->code); // encode
        havoc_populate_ssd<uint8_t>(this, this->code); // encode
        havoc_populate_ssd<uint16_t>(this, this->code); // encode
        havoc_populate_sad<uint8_t>(this, this->code); // encode
        havoc_populate_sad<uint16_t>(this, this->code); // encode
        havoc_populate_sad_multiref<uint8_t>(this, this->code); // encode
        havoc_populate_sad_multiref<uint16_t>(this, this->code); // encode
        havoc_populate_hadamard_satd<uint8_t>(this, this->code); // encode
        havoc_populate_hadamard_satd<uint16_t>(this, this->code); // encode
        havoc::populateSubtractBi<uint8_t>(this, this->code); // encode
        havoc::populateSubtractBi<uint16_t>(this, this->code); // encode
    }

    ~StateFunctionTables()
    {
        havoc_delete_code(this->code);
    }

    havoc_code code;

    havoc_instruction_set const instruction_set_support;
};

#endif
