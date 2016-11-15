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

// command-line entry point and dispatch to sub commands such as encode, decode, etc.

#ifdef CUSTOM
#include "custom.h"
#else
#define COMMAND_DECODE
#define COMMAND_ENCODE
#define COMMAND_PSNR
#define COMMAND_TESTDECODE
#define COMMAND_HAVOC
#define COMMAND_SIGNATURE
#define COMMAND_VERSION
#endif

#include "havoc/havoc.h"
#include "turing.h"

#include <string>
#include <vector>
#include <iostream>
#include <string>


int decode(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);
int psnr(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);
int testdecode(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);
int encode(int argc, const char* const argv[]);
int signature(int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr);



const char *gitDescribe();

int help(const char *programName,bool badArguments)
{
    std::ostream &o = badArguments ? std::cerr : std::cout;

    o <<
#ifdef CUSTOM
            CUSTOM
#else
            "Turing codec version: " << std::string(turing_version()) << "\n"
#endif
            "\n"
            "usage: " << programName << " <command> [<args>]\n"
            "\n"
            "Supported " << programName << " commands are:\n"
#ifdef COMMAND_ENCODE
            "   encode     Encode bitstream from input YUV\n"
#endif
#ifdef COMMAND_DECODE
            "   decode     Decode bitstream to output YUV\n"
#endif
#ifdef COMMAND_PSNR
            "   psnr       Report PSNR between two YUV files\n"
#endif
#ifdef COMMAND_TESTDECODE
            "   testdecode Test that decoding conformance streams produces correct output\n"
#endif
#ifdef COMMAND_SIGNATURE
            "   signature  Checks encoder integrity against expected outputs\n"
#endif
#ifdef COMMAND_HAVOC
            "   havoc      Runs havoc self test\n"
#endif
#ifdef COMMAND_VERSION
            "   version    Concise version description\n"
#endif
            "\n"
            "Use '" << programName << " help <command>' for more information on a specific command.\n\n";

    return badArguments ? 1 : 0;
}


int main(int argc, const char *argv[])
{
    std::vector<const char *> arguments(&argv[0], &argv[argc]);

    std::string programName = arguments[0];
    auto slashPos = programName.find_last_of("/\\");
    if (slashPos != std::string::npos) programName.erase(0, slashPos+1);

    arguments[0] = programName.c_str();

    if (arguments.size() <= 1) return help(arguments[0], true);

    std::string command = arguments[1];

    if (command == "help")
    {
        if (arguments.size() == 2) return help(arguments[0], false);

        // translate "turing help command" to "turing command --help"
        arguments.resize(3);
        command = arguments[1] = arguments[2];
        arguments[2] =  "--help";
    }

    std::string composite = arguments[0];
    composite += " ";
    composite += arguments[1];

    arguments[1] = composite.c_str();

    argc = int(arguments.size()) - 1;
    argv = &arguments[1];

#ifdef COMMAND_ENCODE
    if (command == "encode")
        return encode(argc, argv);
#endif
#ifdef COMMAND_DECODE
    if (command == "decode")
        return decode(argc, argv, std::cout, std::cerr);
#endif
#ifdef COMMAND_PSNR
    if (command == "psnr")
        return psnr(argc, argv, std::cout, std::cerr);
#endif
#ifdef COMMAND_TESTDECODE
    if (command == "testdecode")
        return testdecode(argc, argv, std::cout, std::cerr);
#endif
#ifdef COMMAND_HAVOC
    if (command == "havoc")
        return havoc_main(argc, argv);
#endif
#ifdef COMMAND_SIGNATURE
    if (command == "signature")
        return signature(argc, argv, std::cout, std::cerr);
#endif
#ifdef COMMAND_VERSION
    if (command == "version")
    {
        std::cout << std::string(turing_version());
        return 0;
    }
#endif

    return help(programName.c_str(), true);
}
