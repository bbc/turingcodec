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

// Self-test code relating to the "turing testdecode" command.

#include "ConformanceStreams.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/size.hpp>
#include <fstream>
#include <cassert>
#include <string>


int decode(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);


namespace po = boost::program_options;


int parseTestDecodeOptions(po::variables_map &vm, int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr)
{
    po::options_description options("Options");
    options.add_options()
        ("bitstreams", po::value<std::string>()->default_value("./bitstreams/"), "location of folder containing test clips from https://github.com/kupix/bitstreams")
        ("help,h", "display help message");

    po::options_description all;
    all.add(options);

    try
    {
        po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
    
        if (vm.count("help"))
        {
            cout << "Usage: " << argv[0] << " [options]\n\n";
            cout << options << "\n";
            std::exit(EXIT_SUCCESS);
        }

        po::notify(vm);
    }
    catch (std::exception & e)
    {
        cerr << argv[0] << ": unable to parse command line - " << e.what() << "\n";
        return 1;
    }
    return 0;
}


template <class Args>
void printArgs(std::ostream &os, Args &args)
{
    os << "\n";
    for (auto a : args) 
        os << a << " ";
    os << "\n";
}


void decodeStream(boost::filesystem::path bitstream, boost::filesystem::path md5, std::ostream &cout, std::ostream &cerr, int &total, int &passed)
{
    std::vector<const char *> argvDecode;

    argvDecode.push_back("decode");

    argvDecode.push_back("--md5");
    auto md5string = md5.string();
    argvDecode.push_back(md5string.c_str());
    
    auto bitstreamstring = bitstream.string();
    argvDecode.push_back(bitstreamstring.c_str());

    printArgs(cout, argvDecode);

    int rv = decode((int)argvDecode.size(), &argvDecode[0], cout, cerr);

    ++total;

    if (rv)
        cout << "testdecode failed when decoding " << bitstream << "\n";
    else
        ++passed;

    cout << passed << "/" << total << " streams passed\n";
}


int decodeAllStreams(const po::variables_map& vm, std::ostream &cout, std::ostream &cerr)
{
    cout << "\nTest decoding of conformance streams\n";

    std::vector<ConformanceStreams::Stream> streams;
    ConformanceStreams::find(streams, vm["bitstreams"].as<std::string>(), cout, cerr);

    int total = 0;
    int passed = 0;

    for (auto const&stream : streams)
        decodeStream(stream.bitstream, stream.md5, cout, cerr, total, passed);

    const int rv = static_cast<int>(total - passed);

    if (rv)
        cerr << "decodeConformanceStreams(): one or more bitstreams failed to decode correctly\n";

    return rv;
}


int testdecode(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr)
{
    po::variables_map vm;

    int rv = parseTestDecodeOptions(vm, argc, argv, cout, cerr);

    try
    {
        if (!rv) rv = decodeAllStreams(vm, cout, cerr);
    }
    catch (std::runtime_error& e)
    {
        cerr << "\nstd::runtime_error exception: " << e.what() << "\n";
        rv = -1;
    }

    if (rv)
        cout << "\n\n** testdecode failed **\n";
    else
        cout << "\n\n** testdecode completed successfully **\n";

    return rv;
}
