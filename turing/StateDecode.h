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

#ifndef INCLUDED_StateDecode_h
#define INCLUDED_StateDecode_h

#pragma once

#include "StreamReader.h"
#include "StatePicture.h"
#include "ProgressReporter.h"
#include "GlobalState.h"
#include "md5.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <string>


struct StateRbspData
{
    std::vector<uint8_t> rbspData;
    std::vector<size_t> eb3pPositions;
};


template <class S>
struct Access<rbsp_byte, S, typename std::enable_if<std::is_base_of<StateRbspData, S>::value>::type>
{
    typedef char Type;
    static void set(rbsp_byte e, char c, S &s)
    {
        assert(s.rbspData.size() == e.k);
        s.rbspData.push_back(c);
    }
    static Type get(rbsp_byte e, S &s)
    {
        return s.rbspData[e.k];
    }
};


struct StateDec :
    StateSequence,
    Neighbourhood,
    StreamReader,
    StateRbspData,
    ValueHolder<pic_type>,
    ValueHolder<Stop>
{
    StateDec()
    {
        StatePicturesBase *statePicturesBase = this;
        statePicturesBase->streamType = StatePicturesBase::streamTypeNal();
    }
};


struct StateDecode :
    StateDec,
    StateFunctionTables
{
    StateDecode(const boost::program_options::variables_map &vm, std::ostream &cout, std::ostream &cerr, size_t maxPictures = 0)
        :
#ifdef VALGRIND_FRIENDLY
        StateFunctionTables(false, havoc_instruction_set(HAVOC_C_OPT | HAVOC_C_REF)),
#else
        StateFunctionTables(false),
#endif
        vm(vm),
        ofs(getOptionParameter(vm, "output-file"), std::ios_base::binary),
        ifs(getOptionParameter(vm, "input-file"), std::ios_base::binary),
        maxPictures(maxPictures),
        n(0),
        progressReporter(vm.count("no-progress") ? 0 : &cerr, "decoded", "pictures", maxPictures)
    {
        if (!this->ifs)
        {
            throw std::runtime_error("could not open input file");
        }

        this->StreamReader::open(ifs, false);

        if (vm.count("output-file") != 0 && !this->ofs)
        {
            throw std::runtime_error("could not open output file");
        }
        if (vm.count("md5") != 0)
        {
            std::ifstream i(vm["md5"].as<std::string>());
            while (i)
            {
                std::string s;
                std::getline(i, s);
                std::transform(s.begin(), s.end(), s.begin(), ::tolower);
                if (s.length() >= 3 && s.substr(0, 3) == "md5")
                {
                    if (s.find("yuv") == std::string::npos)
                        continue;
                    s = s.substr(s.find(" = ") + 3, 32);
                }
                if (s.length() >= 16 && s[0] != '#')
                {
                    this->md5Expected = s.substr(0, 32);
                    break;
                }
            }
            if (this->md5Expected.length() != 32)
                throw std::runtime_error("could not read md5 digest file");
            md5_init(&this->md5Sum);
        }
    }

    struct Finished { };

    static std::string md5DigestToStr(md5_byte_t const digest[16])
    {
        std::ostringstream oss;
        for (int i = 0; i < 16; ++i)
            oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<unsigned>(digest[i]);
        return oss.str();
    }

    void finish()
    {
        this->progressReporter.finish();

        // review: remove md5 from decode() - can be done by caller
        if (vm.count("md5") != 0)
        {
            md5_byte_t md5Actual[16];
            md5_finish(&this->md5Sum, md5Actual);
            auto actual = md5DigestToStr(md5Actual);
            std::cout << actual << " actual\n";
            std::cout << this->md5Expected << " expected\n";

            if (actual != this->md5Expected)
            {
                throw std::runtime_error("MD5 mismatch");
            }
        }
    }

    static const std::string &getOptionParameter(const boost::program_options::variables_map &vm, const char *name)
    {
        static const std::string empty;
        return vm.count(name) ? vm[name].as<std::string>() : empty;
    }

    bool stop() const
    {
        return this->maxPictures != 0 && this->n == this->maxPictures;
    }

    std::ifstream ifs;
    std::ofstream ofs;
    const boost::program_options::variables_map &vm;
    std::string md5Expected;
    md5_state_t md5Sum;

    ProgressReporter progressReporter;
    size_t maxPictures;
    size_t n;
};

#endif
