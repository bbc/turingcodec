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

// Constructs to facilitate access to all current state via single "handler" variable.

#ifndef INCLUDED_Handlers_h
#define INCLUDED_Handlers_h

#pragma once

#include "StateValues.h"
#include <cassert>


// Value tag Struct<S> provides a reference to one of the consituent structures that make up a handler.
// Value tag Struct<S *> provides a reference to a pointer to S. This allows another instance of S to be swapped into the handler.
template <class S>
struct Struct
{
};


template <class T, class S>
struct Access<Struct<T>, S, typename std::enable_if<std::is_base_of<T, S>::value>::type>
{
    typedef T &Type;
    static Type get(Struct<T>, S &s)
    {
        return static_cast<T &>(s);
    }
};


// Value tag providing access to the concrete state structure that inherits from Base
template <class Base>
struct Concrete
{
};


template <class Base, class S>
struct Access<Concrete<Base>, S, typename std::enable_if<std::is_base_of<Base, S>::value>::type>
{
    typedef S &Type;
    typedef S ActualType;
    static Type get(Concrete<Base>, S &s)
    {
        return s;
    }
};


// A syntax element: a pair of Value type and a Mnenomic type.
template <class V, class M>
struct Element
{
    V v;
    M m;
};


// If the state conveyed by H includes a copy of the current syntax function, F, update that copy.
template <class F, class H, class Enable = void>
struct CopyValueToState
{
    CopyValueToState(H &, F const &)
    {
    }
};


// A base class for use in the CRTP pattern. Provides [] operators for easy access to the base class's data values.
template <class Derived>
struct AccessOperators
{
    template <class V>
    inline typename Access<V, Derived>::Type operator[](V v)
    {
        return Access<V, Derived>::get(v, *static_cast<Derived*>(this));
    }

    template <class V>
    inline typename Access<V, const Derived>::Type const operator[](V v) const
    {
        return Access<V, const Derived>::get(v, *static_cast<const Derived*>(this));
    }
};


// A heterogenous tuple of pointers to state structures.
template <class... States>
struct PointerTuple
{
};


// Combines a PointerTuple and a Verb which describes what action will be taken when the syntax tree is walked.
// Verb<void> is used as a tag type to represent the action to be invoked, e.g.
// Write<void>, Read<>.
template <class Verb, class ...States>
struct Handler :
    PointerTuple<>
    {
    };


// Specialisation for tuple size > 0.
template <class State, class... States>
struct PointerTuple<State, States...> :
PointerTuple<States...>
{
    State *state;

    // Access to a value whose state is represented somewhere within the tuple
    template <class V>
    typename Access<V, PointerTuple>::Type operator[](V v)
    {
        return Access<V, PointerTuple>::get(v, *this);
    }

    // Convenient access to a structure somewhere within the tuple
    template <class S>
    operator S *()
    {
        return &Access<Struct<S>, PointerTuple>::get(Struct<S>(), *this);
    }

    // Returns a new handler with a different verb.
    template <class NewVerb>
    Handler<NewVerb, State, States...> change()
    {
        Handler<NewVerb, State, States...> h;
        static_cast<PointerTuple<State, States...> &>(h) = *this;
        return h;
    }
};


// Specialisation for tuple size > 1.
template <template <typename> class Verb, class Specialisation, class State, class ...States>
struct Handler<Verb<Specialisation>, State, States...> :
    PointerTuple<State, States...>
{
    // Creates a new handler containing all of our state pointers plus one more.
    template <class StateExtra>
    Handler<Verb<Specialisation>, StateExtra, State, States...> extend(StateExtra *stateExtra)
    {
        Handler<Verb<Specialisation>, StateExtra, State, States...> h;
        static_cast<PointerTuple<State, States...> &>(h) = *this;
        h.state = stateExtra;
        return h;
    }

    // Invokes a syntax element having value V and mnenomic M
    template <class V, class M>
    inline void operator()(V v, M m)
    {
        Verb<Element<V, M>>::go(Element<V, M>{v, m}, *this);
    }

    // Invokes syntax construct F according to Verb
    template <class F>
    inline void operator()(F fun)
    {
        CopyValueToState<F, Handler> copy{ *this, fun };
        Verb<F>::go(fun, *this);
    }

    template <class V>
    inline typename HandleValue<V, Verb<Specialisation>, PointerTuple<State, States...>>::Type operator[](V v)
    {
        return HandleValue<V, Verb<Specialisation>, PointerTuple<State, States...>>::get(v, *this);
    }

    typedef Verb<void> Tag;
    typedef Specialisation Mode;
};



// Construct allowing Verb to force values to constants
template <class V, class Verb, class Tuple, class Enable>
struct HandleValue :
    Access<V, Tuple>
    {
    };


// Specialisation that allows Verb to forces values to constants
template <class V, class Verb, class Tuple>
struct HandleValue<V, Verb, Tuple, typename std::enable_if<std::is_base_of<ValueConstBase<V>, Verb>::value>::type>
{
    template <class H>
    static constexpr typename ValueType<V>::Type get(V v, H &h)
    {
        return ConstValue<V>(*((Verb *)0)).v;
    }

    typedef int Type;
};


// Access allowing a struct pointer somewhere in the tuple to be modified
// Review: opaque, consider removal
template <class S, class ...States>
struct Access<Struct<S*>, PointerTuple<S, States...>>
{
    typedef S *&Type;
    static Type get(Struct<S*>, PointerTuple<S, States...> &tuple)
    {
        return tuple.state;
    }
};


// Access to values in a handler is the same as to values in a tuple (may be overriden by ConstValue/HandleValue according to Verb).
template <class V, class Verb, class ...States>
struct Access<V, Handler<Verb, States...>> :
Access<V, PointerTuple<States...>>
{
};


// Helper macros to aid definition of defined values.
#define DEFINE_DERIVED_LONG(V) \
    template <class H> \
    struct Derive<V, H> \
        : \
        NoAccess \
    { \
        typedef ValueType<V>::Type Type; \
        static inline Type get(V v, H &h)


#define DEFINE_DERIVED(V, F) \
    struct V { }; \
    DEFINE_DERIVED_LONG(V) \
        { \
            return F; \
        } \
    }; \


// Convenience function
template <class V, class H>
static inline typename ValueType<V>::Type derive(V v, H &h)
{
    return Derive<V, H>::get(v, h);
}


// A metafunction to return the expected value of elements whose value
// is fixed by the HEVC standard.
template <class V> struct Fixed;


// Accessor for fixed-value elements
template <class V, class H>
struct Derive <V, H, typename std::enable_if<Fixed<V>::value == Fixed<V>::value>::type> :
    NoAccess
{
    typedef typename ValueType<V>::Type Type;

    static void set(V, Type x, H &h)
    {
        assert(x == Fixed<V>::value);
    }
};


// If value is found in State, access it there
template <class V, class State, class ...States>
struct Access<V, PointerTuple<State, States...>, typename std::enable_if<Accessible<V, State>::value>::type> :
    Access<V, State>
    {
        static typename Access<V, State>::Type get(V v, PointerTuple<State, States...> &pointers)
        {
            return Access<V, State>::get(v, *pointers.state);
        }
        static void set(V v, typename Access<V, State>::Type x, PointerTuple<State, States...> &pointers)
        {
            Access<V, State>::set(v, x, *pointers.state);
        }
    };


// If value is not found in State but is somewhere within States, access it via recursion
template <class V, class State, class ...States>
struct Access<V, PointerTuple<State, States...>, typename std::enable_if<!Accessible<V, State>::value && Accessible<V, PointerTuple<States...>>::value>::type> :
    Access<V, PointerTuple<States...>>
    {
        static typename Access<V, PointerTuple<States...>>::Type get(V v, PointerTuple<State, States...> &pointers)
        {
            return Access<V, PointerTuple<States...>>::get(v, static_cast<PointerTuple<States...> &>(pointers));
        }
        static void set(V v, typename Access<V, PointerTuple<States...>>::Type x, PointerTuple<State, States...> &pointers)
        {
            Access<V, PointerTuple<States...>>::set(v, x, static_cast<PointerTuple<States...> &>(pointers));
        }
    };


template <class F> struct IsRecursive : std::false_type { };
template <> struct IsRecursive<struct coding_quadtree> : std::true_type {};
template <> struct IsRecursive<struct transform_tree> : std::true_type {};


// Construct to store syntax function values to state so that, for example, we always know the current coding_unit()
template <class F, class H>
struct CopyValueToState<F, H, typename std::enable_if<Accessible<Struct<F>, H>::value && !IsRecursive<F>::value>::type>
{
    CopyValueToState(H &h, F const &f)
    {
        h[Struct<F>()] = f;
    }
};


// Specialisation for recursive syntax functions
template <class F, class H>
struct CopyValueToState<F, H, typename std::enable_if<Accessible<Struct<F>, H>::value && IsRecursive<F>::value>::type >
{
    CopyValueToState(H &h, F const &f) :
        copy(h[Struct<F>()])
    {
        this->previous = this->copy;
        h[Struct<F>()] = f;
    }
    ~CopyValueToState()
    {
        this->copy = this->previous;
    }
    F &copy;
    F previous;
};


// Do nothing
template <class F> struct Null
{
    template <class H> static void go(F, H&)
    {
    }
};

#endif
