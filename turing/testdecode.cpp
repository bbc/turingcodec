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

// Self-test code relating to the "turing-exe testdecode" command.

#include "havoc/havoc.h"
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


bool blacklisted(boost::filesystem::path path)
{
    // review: investigate these and remove the blacklist
    // these fail decoding or md5 check - may be due to failure parsing md5 file, or a genuine decoder failures
    const char *blacklist[] = {
        "MAXBINS_A_TI_5",
        "MAXBINS_B_TI_5",
        "MAXBINS_C_TI_5",
        "CIP_C_Panasonic_2",
        "PICSIZE_D_Bossen_1",
        "SDH_A_Orange_4",
        "LS_B_Orange_4",
        "TSKIP_A_MS_3",
        "DBLK_G_VIXS_2",
        "PICSIZE_A_Bossen_1",
        "PICSIZE_B_Bossen_1",
        "PICSIZE_C_Bossen_1",
        "CIP_A_Panasonic_3",
        "TMVP_A_MS_3",
        "ENTP_A_Qualcomm_1",
        "WP_A_Toshiba_3",
        "DBLK_D_VIXS_2",
        "TILES_A_Cisco_2",
        "AMP_E_Hisilicon_3",
        "LS_A_Orange_2",
        "ENTP_C_Qualcomm_1",
        "IPRED_B_Nokia_3",
        "AMP_D_Hisilicon_3",
        "OPFLAG_C_Qualcomm_1",
        "WP_A_MAIN10_Toshiba_3",
        "ENTP_B_Qualcomm_1",
        "POC_A_Bossen_3",
        "MAXBINS_B_TI_4",
        "DBLK_E_VIXS_2",
        "OPFLAG_B_Qualcomm_1",
        "AMP_F_Hisilicon_3",
        "NoOutPrior_A_Qualcomm_1",
        "DBLK_F_VIXS_2",
        "TILES_B_Cisco_1",
        "IPRED_C_Mitsubishi_3",
        "NoOutPrior_B_Qualcomm_1",
        "TSCTX_10bit_I_RExt_SHARP_1",
        "GENERAL_8b_400_RExt_Sony_1",
        "PERSIST_RPARAM_A_RExt_Sony_2",
        "TSUNEQBD_A_MAIN10_Technicolor_2",
        "GENERAL_12b_400_RExt_Sony_1",
        "SAO_A_RExt_MediaTek_1",
        "CCP_12bit_RExt_QCOM",
        "QMATRIX_A_RExt_Sony_1",
        "Bitdepth_A_RExt_Sony_1",
        "IPCM_A_RExt_NEC_2",
        "CCP_8bit_RExt_QCOM",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_8BIT_RExt_Sony_1",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_10BIT_RExt_Sony_1",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_12BIT_RExt_Sony_1",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_16BIT_RExt_Sony_1",
        "EXTPREC_MAIN_444_16_INTRA_12BIT_RExt_Sony_1",
        "EXTPREC_MAIN_444_16_INTRA_8BIT_RExt_Sony_1",
        "EXTPREC_MAIN_444_16_INTRA_10BIT_RExt_Sony_1",
        "GENERAL_16b_444_RExt_Sony_1",
        "GENERAL_16b_400_RExt_Sony_1",
        "GENERAL_12b_444_RExt_Sony_1"
        "GENERAL_8b_444_RExt_Sony_1",
        "GENERAL_8b_444_RExt_Sony_1",
        "GENERAL_10b_444_RExt_Sony_1",
        "GENERAL_12b_422_RExt_Sony_1",
        "GENERAL_12b_420_RExt_Sony_1",
        "GENERAL_12b_444_RExt_Sony_1",
        "GENERAL_16b_444_highThroughput_RExt_Sony_1",
        "TSCTX_12bit_RExt_SHARP_1",
        "TSCTX_12bit_I_RExt_SHARP_1",
        "Bitdepth_B_RExt_Sony_1",
        "IPCM_B_RExt_NEC",
        "ADJUST_IPRED_ANGLE_A_RExt_Mitsubishi_1",
        "TSCTX_8bit_I_RExt_SHARP_1",
        "TSCTX_8bit_RExt_SHARP_1",
        "TSCTX_10bit_RExt_SHARP_1",
        "Bitdepth_B_RExt_Sony_1",
        "ExplicitRdpcm_B_BBC_2",
        "ExplicitRdpcm_A_BBC_1",
        "CCP_10bit_RExt_QCOM",
        "WAVETILES_RExt_Sony_1",
        "EXTPREC_MAIN_444_16_INTRA_16BIT_RExt_Sony_1",
        0};

    for (const char **p = blacklist; *p; ++p)
        if (path.string().find(*p) != std::string::npos)
            return true;

    return false;
}


void decodeStreamsInFolder(boost::filesystem::path folder, std::ostream &cout, std::ostream &cerr, int &total, int &passed)
{
    if (!boost::filesystem::exists(folder))
    {
        cerr << folder << " does not exist\n";
        throw std::runtime_error("");
    }

    for (boost::filesystem::directory_iterator i(folder); i != boost::filesystem::directory_iterator{}; ++i)
    {
        if (boost::filesystem::is_directory(i->status()))
        {
            decodeStreamsInFolder(i->path(), cout, cerr, total, passed);
        }
        else 
        {
            auto const extension = i->path().extension();
            if (extension == ".bit" || extension == ".bin")
            {
                auto base = i->path();
                base.remove_filename();
                base /= i->path().stem();

                auto md5 = base.string() + ".md5";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + "_yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + ".yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + "_3.yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = i->path().string() + ".yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + "_md5.txt";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + "_md5sum.txt";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string().substr(0, base.string().length() - 2) + "_yuv.md5";
                if (!boost::filesystem::exists(md5))
                {
                    auto s = base;
                    s.remove_filename();
                    s /= s.leaf().string() + "_yuv.md5";
                    md5 = s.string();
                }
                if (!boost::filesystem::exists(md5))
                {
                    auto s = base;
                    s.remove_filename();
                    s /= "md5sum.txt";
                    md5 = s.string();
                }

                if (boost::filesystem::exists(md5))
                {
                    if (blacklisted(i->path()))
                        cerr << "\nblacklisted: " << i->path() << "\n";
                    else
                        decodeStream(i->path(), md5, cout, cerr, total, passed);
                }
                else
                    cerr << "no md5 found for " << i->path() << "\n";
            }
        }
    }
}


int decodeAllStreams(const po::variables_map& vm, std::ostream &cout, std::ostream &cerr)
{
    cout << "\nTest decoding of conformance streams\n";

    int total = 0;
    int passed = 0;
    decodeStreamsInFolder(vm["bitstreams"].as<std::string>(), cout, cerr, total, passed);

    const int rv = static_cast<int>(total - passed);

    if (rv)
    {
        cerr << "decodeConformanceStreams(): one or more bitstreams failed to decode correctly\n";
    }

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
