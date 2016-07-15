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
#include "Read.h"
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
    StateSlice,
    ValueHolder<pic_type>,
    ValueHolder<Stop>
    {
        StateDec()
        {
            StatePicturesBase *statePicturesBase = this;
            statePicturesBase->streamType = StatePicturesBase::streamTypeNal();
            statePicturesBase->deblockWholeFrame = true;
        }
    };


    struct StateDecode :
        StateDec,
        StateFunctionTables
        {
            StateDecode(const boost::program_options::variables_map &vm, size_t maxPictures = 0)
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
            progressReporter(vm.count("no-progress") ? 0 : &std::cerr, "decoded", "pictures", maxPictures)
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
                    if (!i) throw std::runtime_error("could not read md5 digest file");
                    i >> this->md5Expected;
                    md5_init(&this->md5Sum);
                }
            }

            struct Finished { };

            template <class StatePicture> void deliver(StatePicture &dp, int chromaFormat, int bitDepthY, int bitDepthC)
            {
                if (!this->stop())
                {
                    if (this->ofs)
                    {
                        typedef typename StatePicture::Sample Sample;

                        // if requested and necessary, round picture to 8-bit precision
                        if (this->vm["8-bit"].as<bool>() && std::is_same<Sample, uint16_t>::value)
                        {
                            static std::ifstream ifs("hm.yuv", std::ios::binary);

                            std::shared_ptr<Picture<uint16_t>> pictureHm(new Picture<uint16_t>(
                                    dp.conformanceWindow[0].width,
                                    dp.conformanceWindow[0].height,
                                    chromaFormat, 0, 0, 32));

                            ifs >> *pictureHm;

                            auto *picture16 = reinterpret_cast<ThreePlanes<uint16_t> *>(&dp.conformanceWindow);
                            std::shared_ptr<Picture<uint8_t>> picture8(new Picture<uint8_t>(
                                    dp.conformanceWindow[0].width,
                                    dp.conformanceWindow[0].height,
                                    chromaFormat, 0, 0, 32));

                            for (int cIdx = 0; cIdx < 3; ++cIdx)
                            {
                                auto const bitDepth = cIdx ? bitDepthC : bitDepthY;
                                auto const add = 1 << bitDepth >> 8;
                                auto const shift = bitDepth - 8;
                                for (int y = 0; y < (*picture8)[cIdx].height; ++y)
                                {
                                    for (int x = 0; x < (*picture8)[cIdx].width; ++x)
                                    {
                                        int value = dp.conformanceWindow[cIdx](x, y);
                                        value += add;
                                        value >>= shift;
                                        if (value > 255) value = 255;
                                        (*picture8)[cIdx](x, y) = value;
                                    }
                                }
                            }

                            this->ofs << *picture8;
                        }
                        else
                        {
                            this->ofs << dp.conformanceWindow;
                        }
                    }

                    if (vm.count("md5"))
                    {
                        std::ostringstream oss;
                        oss << dp.conformanceWindow;
                        std::string const s = oss.str();
                        md5_append(&this->md5Sum, &reinterpret_cast<md5_byte_t const&>(s.front()), static_cast<int>(s.size()));
                    }

                    this->progressReporter.progress(++this->n);
                }

                if (this->stop())
                {
                    // Review: not really appropriate to use an exception for this purpose:
                    // may be neater to add logic into Bitstream and CodedVideoSequence processing to exit loops.
                    throw Finished();
                }
            }

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
