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

#ifndef INCLUDED_Violation_h
#define INCLUDED_Violation_h

#pragma once

#include <string>
#include <map>

#ifdef HAVE_BOOST
    #ifdef __INTEL_COMPILER
// boost_1_55_0\boost/format/alt_sstream.hpp(65): warning #809: exception specification for virtual function "boost::io::basic_altstringbuf<Ch, Tr, Alloc>::~basic_altstringbuf [with Ch=char, Tr=std::char_traits<char>, Alloc=std::allocator<char>]" is incompatible with that of overridden function "std::basic_streambuf<_Elem, _Traits>::~basic_streambuf [with _Elem=char, _Traits=std::char_traits<char>]"
        #pragma warning disable 809
    #endif
    #include <boost/format.hpp>
#else
    #include <sstream>
    #include <string>
#endif


struct Violation
{
    Violation(std::string clause, const char *format, bool decodable = false, int profileIdc = -1) :
        clause(clause),
        format(format),
        message(format),
        decodable(decodable),
        profileIdc(profileIdc)
        {
        }
    const char *format;
#ifdef HAVE_BOOST_FORMAT
    boost::format message;
#else
    std::string message;
#endif
    std::string clause;
    bool decodable;
    int profileIdc; // -1 means non-profile-specific violation

    template <class T>
    Violation &operator%(const T &t)
    {
#ifdef HAVE_BOOST_FORMAT
        this->message = this->message % t;
#else
        std::ostringstream oss(this->message);
    oss << " " << t;
    this->message = oss.str();
#endif
    return *this;
    }
};


static bool trim2(std::string &s)
{
    size_t pos = s.find_first_of(" :<");
    if (pos == std::string::npos) return false;
    if (s[pos] == '<')
    {
        s.erase(pos);
        return false;
    }
    pos += (s[pos] == ':') ? 2 : 1;
    s.erase(0, pos);
    return true;
}

static std::string trim(std::string name)
{
    while (trim2(name));
    while (name.find_first_of("0123456789") == 0) name = name.substr(1);
    return name;
}

template <class T>
struct TypeName
{
    static const std::string value;
};

template <class T>
const std::string TypeName<T>::value(trim(typeid(T).name()));


#endif
