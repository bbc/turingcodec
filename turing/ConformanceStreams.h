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

#ifndef INCLUDED_ConformanceStreams
#define INCLUDED_ConformanceStreams

#include <boost/filesystem.hpp>
#include <vector>


namespace ConformanceStreams {

bool blacklisted(boost::filesystem::path path)
{
    // review: decoder currently supports Main, Main10 and Main 4:2:2 10 profiles only - these streams do not fall within these profiles
    const char *blacklist[] = {
        "TSCTX_10bit_I_RExt_SHARP",
        "PERSIST_RPARAM_A_RExt_Sony",
        "SAO_A_RExt_MediaTek",
        "CCP_12bit_RExt_QCOM",
        "QMATRIX_A_RExt_Sony",
        "Bitdepth_A_RExt_Sony",
        "IPCM_A_RExt_NEC",
        "CCP_8bit_RExt_QCOM",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_8BIT_RExt_Sony",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_10BIT_RExt_Sony",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_12BIT_RExt_Sony",
        "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_16BIT_RExt_Sony",
        "EXTPREC_MAIN_444_16_INTRA_12BIT_RExt_Sony",
        "EXTPREC_MAIN_444_16_INTRA_8BIT_RExt_Sony",
        "EXTPREC_MAIN_444_16_INTRA_10BIT_RExt_Sony",
        "EXTPREC_MAIN_444_16_INTRA_16BIT_RExt_Sony",
        "GENERAL_16b_444_RExt_Sony",
        "GENERAL_16b_400_RExt_Sony",
        "GENERAL_12b_444_RExt_Sony",
        "GENERAL_8b_400_RExt_Sony",
        "GENERAL_8b_444_RExt_Sony",
        "GENERAL_10b_444_RExt_Sony",
        "GENERAL_12b_400_RExt_Sony",
        "GENERAL_12b_422_RExt_Sony",
        "GENERAL_12b_420_RExt_Sony",
        "GENERAL_12b_444_RExt_Sony",
        "GENERAL_16b_444_highThroughput_RExt_Sony",
        "TSCTX_12bit_RExt_SHARP",
        "TSCTX_12bit_I_RExt_SHARP",
        "Bitdepth_B_RExt_Sony",
        "IPCM_B_RExt_NEC",
        "ADJUST_IPRED_ANGLE_A_RExt_Mitsubishi",
        "TSCTX_8bit_I_RExt_SHARP",
        "TSCTX_8bit_RExt_SHARP",
        "TSCTX_10bit_RExt_SHARP",
        "Bitdepth_B_RExt_Sony",
        "ExplicitRdpcm_B_BBC",
        "ExplicitRdpcm_A_BBC",
        "CCP_10bit_RExt_QCOM",
        "WAVETILES_RExt_Sony",
        "SAODBLK_C_MainConcept_1",
        "3D-HEVC",
        "MV-HEVC",
        0};

    for (const char **p = blacklist; *p; ++p)
        if (path.string().find(*p) != std::string::npos)
            return true;

    return false;
}


struct Stream
{
    boost::filesystem::path bitstream;
    boost::filesystem::path md5;
};


void find(std::vector<Stream> &streams, boost::filesystem::path folder, std::ostream &cout, std::ostream &cerr)
{
    if (!boost::filesystem::exists(folder))
    {
        cerr << folder << " does not exist\n";
        throw std::runtime_error("");
    }

    for (boost::filesystem::directory_iterator i(folder); i != boost::filesystem::directory_iterator{}; ++i)
    {
        if (boost::filesystem::is_directory(i->status()))
            find(streams, i->path(), cout, cerr);
        else
        {
            auto const extension = i->path().extension();
            if (extension == ".bit" || extension == ".bin")
            {
                auto base = i->path();
                base.remove_filename();
                base /= i->path().stem();

                auto md5 = base.string() + ".yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + "_yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + "_3.yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = i->path().string() + ".yuv.md5";
                if (!boost::filesystem::exists(md5))
                    md5 = base.string() + ".md5";
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
                        streams.push_back({ i->path(), md5 });
                }
                else
                    cerr << "no md5 found for " << i->path() << "\n";
            }
        }
    }
}

}
#endif