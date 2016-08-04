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

#ifndef INCLUDED_BitField_h
#define INCLUDED_BitField_h

// BitField type inspired by http://blog.codef00.com/2014/12/06/portable-bitfields-using-c11/


// Bit field descriptor designed to be grouped together in a union. 
// Raw is the underlying type - an unsigned integer type.
// offset represents the shifted position of the bit field
// size is the number of bits in the field.
template <typename Raw, int offset, int size>
struct BitField
{
    BitField &operator=(int other)
    {
        this->value &= ~mask;
        this->value |= other << offset;
        return *this;
    }

    operator int() const { return (this->value & mask) >> offset;  }

    BitField const &operator |=(BitField other)
            {
        this->value |= other.value & mask;
        return *this;
            }

private:

    static const Raw mask = ((1 << size) - 1) << offset;

    Raw value;
};

template <typename Raw, int offset, int size, int n>
struct BitFieldArray
{
    struct Entry
    {
        Entry &operator=(int other)
        {
            this->value &= ~this->mask();
            this->value |= other << this->off;
            return *this;
        }

        operator int() const { return (this->value &  this->mask()) >> this->off; }

        Raw &value;
        const int off;

        Raw mask() const { return ((1 << size) - 1) << this->off; }
    };

    Entry operator[](int i)
    {
        return Entry{ value, offset + i * size };
    }

    Raw value;
};

#endif	
