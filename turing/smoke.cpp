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

// Self-test code relating to the "turing-exe smoke" command.

#include "ConformanceStreams.h"
#include "havoc/havoc.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/size.hpp>
#include <fstream>
#include <cassert>
#include <string>

int decode(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);
int encode(int argc, const char* const argv[]);
int psnr(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);
int signature(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);


namespace po = boost::program_options;


int parseSmokeOptions(po::variables_map &vm, int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr)
{
    po::options_description options("Options");
    options.add_options()
		        ("neptune", po::value<std::string>()->default_value(neptuneDefault), "location of user's \"neptune\" folder containing test clips, etc.")
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
            return 0;
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

struct NeptuneFile
{
    NeptuneFile(po::variables_map vm, const char *name)
    {
        path = vm["neptune"].as<std::string>();
        path += name;
    }
    const char *name() const
    {
        return path.c_str();
    }
private:
    std::string path;
};

struct YuvFile
{
    YuvFile(const char *root, int bitDepth, int width, int height)
    {
        std::ostringstream oss;
        oss << root << bitDepth << "bit_" << width << "x" << height << ".yuv";
        this->path = oss.str();
    }
    const char *name() const
    {
        return path.c_str();
    }
private:
    std::string path;
};


template <class Args>
void printArgs(std::ostream &os, Args &args)
{
    os << "\n";
    for (auto a : args) os << a << " ";
    os << "\n";
}


int decodeEncodeDecode(const po::variables_map& vm, std::ostream &cout, std::ostream &cerr, int bitDepth)
{
    const char *filename = "conformance/SAO_A_MediaTek_4.bit";
    int width = 416;
    int height = 240;
    int frames = 60;

    if (bitDepth == 10)
    {
        filename = "conformance/WPP_C_ericsson_MAIN10_2.bit";
        frames = 48;
    }

    NeptuneFile bitstreamFile(vm, filename);

    YuvFile decodedYuv("decoded_", bitDepth, width, height);

    {
        boost::filesystem::remove(decodedYuv.name());

        const char *args[] =
        {
                "decode",
                "-o", decodedYuv.name(),
                bitstreamFile.name()
        };

        printArgs(cout, args);

        int rv = decode((int)boost::size(args), args, cout, cerr);

        if (rv)
        {
            cerr << "decodeEncodeDecode failed\n";
            return rv;
        }

        if (boost::filesystem::file_size(decodedYuv.name()) != width * height * frames * 3 / 2 * (bitDepth > 9 ? 2 : 1))
        {
            cerr << "decodeEncodeDecode failed to output " << frames << " frames\n";
            return -1;
        }
    }

    YuvFile encodedYuv("encoded_", bitDepth, width, height);
    boost::system::error_code ec;
    boost::filesystem::remove(encodedYuv.name(), ec);
    boost::filesystem::remove("encoded.bin", ec);

    std::ostringstream o;
    o << "--frames=" << frames;
    std::string framesOption = o.str();

    o.str("");
    o << "--input-res=" << width << "x" << height;
    std::string resOption = o.str();

    o.str("");
    o << "--bit-depth=" << bitDepth;
    std::string bitDepthOptions = o.str();

    //if (bitDepth == 9) bitDepthOptions = "--force-16";

    o.str("");
    o << "encoded" << bitDepth << ".bin";
    std::string bitstreamName = o.str();

    {
        const char *args[] =
        {
                "encode",
                "--ctu=16",
                "--wpp",
                //			"--asm=0",
                "--concurrent-frames=1",
                "--max-gop-n=24",
                "--speed=slow",
                "--amp",
                "--deblock",
                "--qp=36",
                framesOption.c_str(),
                resOption.c_str(),
                bitDepthOptions.c_str(),
                "--psnr",
                "--frame-rate=25",
                "--dump-frames", encodedYuv.name(),
                "-o",  bitstreamName.c_str(),
                decodedYuv.name()
        };

        printArgs(cout, args);

        int rv = encode((int)boost::size(args), args);
        if (rv)
        {
            cerr << "decodeEncodeDecode failed\n";
            return rv;
        }
    }

    cout << "encoded filesize is " << boost::filesystem::file_size(bitstreamName.c_str()) << "\n";

    if (bitDepth == 8 || bitDepth == 10)
    {
        std::ostringstream bitDepthOption;
        bitDepthOption << "--bit-depth=" << bitDepth;
        auto s = bitDepthOption.str();

        const char *args[] =
        {
                "psnr",
                s.c_str(),
                encodedYuv.name(),
                decodedYuv.name()
        };

        printArgs(cout, args);

        std::ostringstream oss;
        int rv = psnr((int)boost::size(args), args, oss, oss);
        cout << oss.str();
    }

    YuvFile encodedDecodedYuv("encoded_decoded_", bitDepth, width, height);
    boost::filesystem::remove(encodedDecodedYuv.name());

    {
        const char *args[] =
        {
                "decode",
                "-o", encodedDecodedYuv.name(),
                bitstreamName.c_str()
        };

        printArgs(cout, args);

        int rv = decode((int)boost::size(args), args, cout, cerr);

        if (rv)
        {
            cerr << "decodeEncodeDecode failed\n";
            return rv;
        }

        if (boost::filesystem::file_size(encodedDecodedYuv.name()) != width * height * frames * 3 / 2 * (bitDepth > 9 ? 2 : 1))
        {
            cerr << "decodeEncodeDecode failed to output " << frames << " frames\n";
            return -1;
        }
    }

    if (bitDepth == 8 || bitDepth == 10)
    {
        std::ostringstream bitDepthOption;
        bitDepthOption << "--bit-depth=" << bitDepth;
        auto s = bitDepthOption.str();

        const char *args[] =
        {
                "psnr",
                s.c_str(),
                encodedYuv.name(),
                encodedDecodedYuv.name()
        };

        printArgs(cout, args);

        std::ostringstream oss;
        int rv = psnr((int)boost::size(args), args, oss, oss);
        cout << oss.str();

        if (oss.str().find("identical") == std::string::npos)
        {
            cerr << "encoder and decoder yuv reconstructions do not match\n";
            return -1;
        }
    }

    return 0;
}

int decodeConformanceStreams(const po::variables_map& vm, std::ostream &cout, std::ostream &cerr)
{
    cout << "\nTest decoding of conformance streams\n";

    size_t total = 0, passed = 0;

    for (auto const entry : ConformanceStreams::streams)
    {
        //		if (entry.profile != 2 && entry.profile != 2 && entry.profile != 4) continue;

        ConformanceStreams::Locator locator(entry.name, vm["neptune"].as<std::string>().c_str());
        std::string md5 = locator.md5.string();
        std::string bitstream = locator.bitstream.string();

        std::vector<const char *> argvDecode;

        argvDecode.push_back("decode");
        argvDecode.push_back("--md5");
        argvDecode.push_back(md5.c_str());
        argvDecode.push_back(bitstream.c_str());

        printArgs(cout, argvDecode);

        int rv = decode((int)argvDecode.size(), &argvDecode[0], cout, cerr);

        ++total;

        if (rv)
        {
            cout << "decodeConformanceStreams smoke test failed decoding " << entry.name << "\n";
            continue;
        }

        ++passed;

        cout << passed << "/" << total << " streams passed\n";
    }

    const int rv = static_cast<int>(total - passed);

    if (rv)
    {
        cerr << "decodeConformanceStreams one or more bitstreams failed to decode correctly\n";
    }

    return rv;
}


int smoke(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr)
{
    po::variables_map vm;

    int rv = parseSmokeOptions(vm, argc, argv, cout, cerr);

    try
    {
        if (!rv) rv = decodeEncodeDecode(vm, cout, cerr, 8);
        if (!rv) rv = decodeEncodeDecode(vm, cout, cerr, 10);
        if (!rv) rv = decodeConformanceStreams(vm, cout, cerr);
        if (!rv) rv = havoc_main(0, 0);
    }
    catch (std::runtime_error& e)
    {
        cerr << "\nstd::runtime_error exception: " << e.what() << "\n";
        rv = -1;
    }

    if (rv)
        cout << "\n\n** smoke test failed **\n";
    else
        cout << "\n\n** smoke test completed successfully **\n";

    return rv;
}
