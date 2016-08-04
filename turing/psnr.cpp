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


#include "havoc/diff.h"
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdint.h>
#include <cassert>
#include <string>
#include <cmath>


namespace po = boost::program_options;


int psnr(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr)
{
    po::options_description options("Options");
    options.add_options()
                ("bit-depth", po::value<int>()->default_value(8), "video bit depth (of both input YUV and output stream)");
    ("help,h", "display help message");

    po::options_description hidden("Hidden options");
    hidden.add_options()
                ("input-file", po::value<std::vector<std::string>>());

    po::options_description all;
    all.add(options).add(hidden);

    po::positional_options_description positional;
    positional.add("input-file", -1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);

        if (vm.count("help"))
        {
            cout << "psnr: \n\n";
            cout << "usage: " << argv[0] << " [options] input-file1 input-file2\n\n";
            cout << options << "\n";
            std::exit(EXIT_SUCCESS);
        }

        po::notify(vm);

        if (vm["input-file"].as<std::vector<std::string>>().size() != 2)
        {
            throw std::runtime_error("two input files not specified");
        }
    }
    catch (std::exception & e)
    {
        cerr << argv[0] << ": unable to parse command line - " << e.what() << "\n";
        return 1;
    }

    std::vector<std::shared_ptr<std::ifstream>> files;

    for (auto &name : vm["input-file"].as<std::vector<std::string>>())
    {
        files.push_back(std::shared_ptr<std::ifstream>(new std::ifstream(name.c_str(), std::ios::binary | std::ios::in)));
        if (!*files.back())
        {
            cout << "failed to open input file: " << name << "\n";
            return -1;
        }
    }

    double sumOfSquareErrors = 0.0;
    double totalSamples = 0.0;

    auto const bitDepth = vm["bit-depth"].as<int>();

    havoc_code code = havoc_new_code(havoc_instruction_set(~0), 4000);

    while (true)
    {
        const size_t blockSize = 4096;
        HAVOC_ALIGN(32, uint8_t, buffer[2][blockSize]);

        int n = blockSize;
        for (int i=0; i<2; ++i)
        {
            buffer[i];
            files[i]->read(reinterpret_cast<char *>(&buffer[i][0]), blockSize);
            std::streamsize k = files[i]->gcount();
            if (int(k) < n) n = int(k);
        }

        if (bitDepth == 8)
        {
            havoc_ssd_linear *ssd = havoc_get_ssd_linear(n, code);
            sumOfSquareErrors += ssd(&buffer[0][0], &buffer[1][0], n);

            totalSamples += n;
        }
        else
        {
            auto const *p0 = reinterpret_cast<uint16_t *>(&buffer[0][0]);
            auto const *p1 = reinterpret_cast<uint16_t *>(&buffer[1][0]);
            for (int i = 0; i < n / 2; ++i)
            {
                int diff = p0[i] - p1[i];
                sumOfSquareErrors += diff * diff;
            }

            totalSamples += n / 2;
        }

        if (n != blockSize) break;
    }

    havoc_delete_code(code);

    if (totalSamples == 0.0)
    {
        cout << "no samples processed (one or both input files is empty)\n";
        return -1;
    }

    if (sumOfSquareErrors == 0.0)
    {
        cout << "files are identical\n";
    }
    else
    {
        const double max = bitDepth == 8 ? 255.0 : 1023.0;
        double mse = sumOfSquareErrors / totalSamples;

        double psnr = 10.0 * log10( max * max / mse );

        cout << "PSNR: " << psnr << "dB\n";
    }

    cout << "Total samples: " << totalSamples << "\n";

    return 0;
}
