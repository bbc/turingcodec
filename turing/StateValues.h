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

// Constructs to manage the storage and manipulation of state values. Such values
// may be syntax elements, derived values or other intermediate values used
// during processing. Values are identified using a C++ type, e.g. "end_of_slice_segment_flag".
// review: rename file - looks too much like it should contains a "struct StateValues"

#ifndef INCLUDED_StateValue_h
#define INCLUDED_StateValue_h

#pragma once

#include <type_traits>
#include <cstdint>
#include <array>
#include <cassert>


// Base type used via std::is_base_of<> to detirmine that value access is not possible.
struct NoAccess
{
};


// Default behaviour for Access is to derive a value from others
// Default derivation behaviour is NoAccess (i.e. cannot access value V via H)
template <class V, class H, class Enable = void>
struct Derive :
    NoAccess
{
};


// Provides read (and optionally write) access to a value V contained in structure S.
// Default behaviour of an Accessor is to derive the value from other values accessible via H.
// An example of a derived value in HEVC is PicSizeInSamplesY.
template <class V, class S, class Enable=void>
struct Access :
    Derive<V, S>
    {
    };

// Metafunction: does structure type S have value V, i.e. does it support Access<V, S> ?
// Returns true_type if handler H can handle requests for value V.
// If true_type, then Access<V, H> can be used to access that value.
template <class V, class S, class Enable = void>
struct Accessible :
    std::integral_constant<bool, !std::is_base_of<NoAccess, Access<V, S>>::value>
    {
    };

// Metafunction that returns storage type for variable type V.
template <class V>
struct ValueType
{
    typedef int Type;
};

// default implementation, different for arrays, etc.
template <class V>
struct ValueHolder
{
    typedef typename ValueType<V>::Type Type;
    ValueHolder(Type value=Type()) : value(value) { }
    Type value;
    Type &get(V v) { return this->value; }
    const Type &get(V v) const { return this->value; }
};

template <class V, int n>
struct Index
{
    static const int value = n;
};

template <class V>
struct IsValueArray : std::false_type { };

template <class V, int size>
struct ValueArray1
{
    typedef typename ValueType<V>::Type Type;
    std::array<Type, size> value;
    ValueArray1()
    {
        std::fill(std::begin(this->value), std::end(this->value), Type());
    }
    Type &get(V v)
    {
        auto const index = v.argValue(Index<V, 0>::value);
        return value[index];
    }
    void set(Type t);
};


template <class V, int size0, int size1>
struct ValueArray2
{
    typename ValueType<V>::Type value[size0][size1];
    typename ValueType<V>::Type &get(V v)
    {
        return value[v.argValue(Index<V, 0>::value)][v.argValue(Index<V, 1>::value)];
    }
    void set(typename ValueType<V>::Type t);
};


template <class V, int size0, int size1, int size2>
struct ValueArray3
{
    std::array<std::array<std::array<typename ValueType<V>::Type, size2>, size1>, size0> value;
    typename ValueType<V>::Type &get(V v)
    {
        return value[v.argValue(Index<V, 0>::value)][v.argValue(Index<V, 1>::value)][v.argValue(Index<V, 2>::value)];
    }
    void set(typename ValueType<V>::Type t);
};


// todo: potentially replace this preprocessing with C++ type magic
#define DEFINE_VALUE_ARRAY_1(name, index, size) \
        \
        DEFINE_STRUCT_ARITY_1(name, index); \
        template <> struct IsValueArray<name> : std::true_type { }; \
        template <> struct ValueHolder<name> : ValueArray1<name, size> { };

// todo: potentially replace this preprocessing with C++ type magic
#define DEFINE_VALUE_ARRAY_2(name, index, size, index1, size1) \
        \
        DEFINE_STRUCT_ARITY_2(name, index, index1); \
        template <> struct IsValueArray<name> : std::true_type { }; \
        template <> struct ValueHolder<name> : ValueArray2<name, size, size1> { };

// todo: potentially replace this preprocessing with C++ type magic
#define DEFINE_VALUE_ARRAY_3(name, index, size, index1, size1, index2, size2) \
        \
        DEFINE_STRUCT_ARITY_3(name, index, index1, index2); \
        template <> struct IsValueArray<name> : std::true_type { }; \
        template <> struct ValueHolder<name> : ValueArray3<name, size, size1, size2> { };


// Used by  to access value V according to Verb and Tuple.
template <class V, class Verb, class Tuple, class Enable = void>
struct HandleValue;


template <class V, class S>
struct Access<V, S, typename std::enable_if<!std::is_const<S>::value && std::is_base_of<ValueHolder<V>, S>::value>::type>
{
    typedef typename ValueType<V>::Type &Type;
    typedef typename ValueType<V>::Type SetType;
    static Type &get(V v, ValueHolder<V> &vh)
    {
        return vh.get(v);
    }
    static void set(V v, typename ValueType<V>::Type t, ValueHolder<V> &vh)
    {
        vh.get(v) = t;
    }
};

template <class V, class S>
struct Access<V, const S, typename std::enable_if<std::is_base_of<ValueHolder<V>, S>::value>::type>
{
    typedef const typename ValueType<V>::Type &Type;
    static Type get(V v, ValueHolder<V> const &vh)
    {
        return vh.get(v);
    }
};


// A ValueCache captures the current value of a parameter for onward processing.
// Useful to avoid repeated re-evaluation of a derived value, for example ScanIdx.
template <class V>
struct ValueCache
{
    template <class H>
    ValueCache(H &h)
    :
    value(h[V()])
    {
    }
    typename ValueType<V>::Type value;
};


template <class V, class S>
struct Access<V, S, typename std::enable_if<std::is_base_of<ValueCache<V>, S>::value>::type>
:
ValueType<V>
{
    static typename ValueType<V>::Type get(V, S &s)
    {
        return s.ValueCache<V>::value;
    }
};


template <class V>
struct ConstValue
{
    typename ValueType<V>::Type v;
};

template <class V>
struct ValueConstBase
{
};


// ValueConst fixes the specified value to a compile-time constant.
// Useful to statically disable unused parts of the specification
template <class V, typename ValueType<V>::Type v>
struct ValueConst :
    ValueConstBase<V>
    {
        constexpr operator ConstValue<V>() const { return ConstValue<V>{v}; }

        // Use this constructor to verify that parent H has same value for V as we are setting const
        template <class H>
        ValueConst(H &h)
        {
            assert(h[V()] == v || !"h has wrong value of V - hardcoded value will be v");
        }

        ValueConst() { }
    };

template <typename Tag, int v>
constexpr int getVal(ValueConst<Tag, v> *)
{
    return v;
}

template <typename Tag, unsigned N>
constexpr auto getConstVal(ValueConst<Tag, N>)
{
    return std::integral_constant<unsigned, N>{};
};


template <class V, class S>
struct Access<V, S, typename std::enable_if<std::is_base_of<ValueConstBase<V>, S>::value>::type> :
    ValueType<V>
    {
        static constexpr unsigned value = decltype(getConstVal<V>(std::declval<S>()))::value;

        constexpr static typename ValueType<V>::Type get(V, S const&)
        {
            return getVal<V>((S*)nullptr);
        }
    };

#endif
