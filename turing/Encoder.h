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

#ifndef INCLUDED_Encoder_h
#define INCLUDED_Encoder_h

#pragma once

#include "Picture.h"
#include "AdaptiveQuantisation.h"
#include "turing.h"
#include "StateEncode.h"
#include "Levels.h"
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/tokenizer.hpp>
#include <cstddef>
#include <stdint.h>
#include <map>
#include <exception>
#include <vector>
#include <memory>


struct Encoder
{
    Encoder(boost::program_options::variables_map &vm);

    // prints summary information before encoding
    void printHeader(std::ostream &cout, std::string const &inputFile, std::string const &OutputFile);

    // prints summary information after encoding
    void printFooter(std::ostream &cout);

    // Writes bitstream headers (i.e. SPS, PPS, etc.)
    void headers(std::vector<uint8_t> &bitstream);

    struct PictureMetadata
    {
        int64_t pts;
        int64_t dts;
        bool keyframe;
    };

    // Synchronous, deterministic frame encode function.
    // Returns true if new encoded bitstream generated.
    // Call with pictureWrapper empty to flush at end of sequence.
    bool encodePicture(std::shared_ptr<PictureWrapper> pictureWrapper, std::vector<uint8_t> &bitstream, PictureMetadata &metadata);

    void parseInputRes();

    boost::program_options::variables_map &vm;
    StateEncode stateEncode;
    int frameWidth;
    int frameHeight;
    int pictureWidth;
    int pictureHeight;
    int confWinBottomOffset;
    int bitDepth;
    int externalBitDepth;
    boost::timer::cpu_timer cpuTimer;
    boost::timer::cpu_timer frameCpuTimer;
    size_t frameCount = 0;
    size_t byteCount = 0;
    double frameRate;
    Level bitstreamLevel;


    // Check whether one of the VUI parameters is set, so that VUI writing can be happen
    bool writeVui();

    template <class H>
    void setupPtl(H &h);

    template <class H>
    void setupVps(H &hhh, ProfileTierLevel *ptl);

    static int log2Size(int size, int min, int max, const char *errorMessage)
    {
        int shift = 0;
        while ((1 << shift) != size)
        {
            if ((1 << shift) > max) throw std::runtime_error(errorMessage);
            ++shift;
        }
        if ((1 << shift) < min) throw std::runtime_error(errorMessage);
        return shift;
    }

    template <class H>
    ProfileTierLevel *setupSps(H &hhh);

    template <class H>
    void setupPps(H &hh);

    template <class H>
    void setupVui(H &h);

    template <class H>
    void setupHrd(H &h);

    // returns parsed value of boolean settings of form "--no-thing" / "--thing"
    bool booleanSwitchSetting(std::string name, bool defaultValue);
};

#endif
