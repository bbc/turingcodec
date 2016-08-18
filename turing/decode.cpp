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
#include "Read.hpp"
#include "StateDecode.h"
#include "Decode.h"
#include "GlobalState.h"
#include "Syntax.h"
#include "SyntaxNal.hpp"
#include "Read.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>


namespace po = boost::program_options;


int parseDecodeOptions(po::variables_map &vm, int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr)
{
    po::options_description options("Options");
    options.add_options()
                ("output-file,o", po::value<std::string>(), "reconstructed YUV file name")
                ("8-bit,8", po::bool_switch(), "round output samples and write an 8-bit YUV file")
                ("frames", po::value<size_t>(), "number of frames to decode")
                ("no-progress", "suppress progress reporting to stderr")
                ("help,h", "display help message");

    po::options_description hidden("Hidden options");
    hidden.add_options()
                ("input-file", po::value<std::string>(), "input bitstream file name")
                ("md5", po::value<std::string>(), "file containing expected yuv output md5");

    po::options_description all;
    all.add(options).add(hidden);

    po::positional_options_description positional;
    positional.add("input-file", 1);

    try
    {
        po::store(po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);

        if (vm.count("help"))
        {
            cout << "decode: \n\n";
            cout << "usage: " << argv[0] << " [options] input-file\n\n";
            cout << options << "\n";
            return 1;
        }

        po::notify(vm);

        if (vm.count("input-file") != 1)
            throw std::runtime_error("no input file specified");
    }
    catch (std::exception & e)
    {
        cerr << argv[0] << ": unable to parse command line - " << e.what() << "\n";
        return 1;
    }

    return 0;
}



int decode(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr)
{
    po::variables_map vm;

    int rv = parseDecodeOptions(vm, argc, argv, cout, cerr);
    if (rv) return rv;

    const size_t nPictures = vm.count("frames") ? vm["frames"].as<size_t>() : 0;

    if (!boost::filesystem::exists(vm["input-file"].as<std::string>()))
    {
        cerr << argv[0] << ": cannot open input file " << vm["input-file"].as<std::string>() << "\n";
        return -1;
    }

    try
    {
        StateDecode stateDecode(vm, cout, cerr, nPictures);

        Handler<Decode<void>, StateDecode> h;
        h.state = &stateDecode;

        h(Bitstream(0));

        stateDecode.finish();
    }
    catch (Abort &)
    {
        cerr << argv[0] << ": problem decoding " << vm["input-file"].as<std::string>() << "\n";
        return 1;
    }
    catch (std::exception &e)
    {
        cerr << argv[0] << ": problem decoding " << vm["input-file"].as<std::string>() << " - " << e.what() << "\n";
        return 1;
    }
    catch (StateDecode::Finished &)
    {
    }

    return 0;
}


int dpbIndexPlus1(struct StatePicture *p, int refList, int refIdx)
{
    return p->dpbIndexPlus1[refList][refIdx];
}
