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

#ifndef INCLUDED_ConformanceStream_h
#define INCLUDED_ConformanceStream_h

#pragma once


#ifdef WIN32
const char *neptuneDefault = "C:\\neptune\\";
#endif

#ifdef __GNUC__
const char *neptuneDefault = "./neptune/";
#endif

#include <boost/filesystem.hpp>

#include <iostream>


namespace ConformanceStreams {

    struct Piece
    {
        int profile;
        int level;
        int levelTenths;
        int fps;
        const char *name;
    };

    static const Piece streams[] =
    {
            { 1, 5, 1, 30, "AMP_A_Samsung_6" },
            { 1, 5, 1, 30, "AMP_B_Samsung_6" },
            { 1, 6, 2, 24, "AMP_D_Hisilicon_3" },
            { 1, 6, 2, 50, "AMP_E_Hisilicon_3" },
            { 1, 6, 2, 60, "AMP_F_Hisilicon_3" },
            { 1, 4, 0, 50, "AMVP_A_MTK_4" },
            { 1, 4, 0, 50, "AMVP_B_MTK_4" },
            { 1, 5, 1, 30, "AMVP_C_Samsung_6" },
            { 1, 5, 1, 30, "BUMPING_A_ericsson_1" },
            { 1, 3, 0, 50, "CAINIT_A_SHARP_4" },
            { 1, 3, 0, 50, "CAINIT_B_SHARP_4" },
            { 1, 3, 0, 50, "CAINIT_C_SHARP_3" },
            { 1, 3, 0, 50, "CAINIT_D_SHARP_3" },
            { 1, 3, 0, 50, "CAINIT_E_SHARP_3" },
            { 1, 3, 0, 50, "CAINIT_F_SHARP_3" },
            { 1, 3, 1, 50, "CAINIT_G_SHARP_3" },
            { 1, 3, 1, 50, "CAINIT_H_SHARP_3" },
            { 1, 4, 0, 30, "CIP_A_Panasonic_3" },
            { 1, 2, 0, 30, "cip_B_NEC_3" },
            { 1, 4, 0, 30, "CIP_C_Panasonic_2" },
            { 1, 4, 0, 30, "CONFWIN_A_Sony_1" },
            { 2, 4, 0, 30, "DBLK_A_MAIN10_VIXS_3" },
            { 1, 4, 0, 30, "DBLK_A_SONY_3" },
            { 1, 4, 0, 30, "DBLK_B_SONY_3" },
            { 1, 4, 0, 30, "DBLK_C_SONY_3" },
            { 1, 4, 1, 60, "DBLK_D_VIXS_2" },
            { 1, 4, 1, 60, "DBLK_E_VIXS_2" },
            { 1, 4, 1, 60, "DBLK_F_VIXS_2" },
            { 1, 4, 1, 60, "DBLK_G_VIXS_2" },
            { 1, 5, 0, 24, "DELTAQP_A_BRCM_4" },
            { 1, 4, 0, 30, "DELTAQP_B_SONY_3" },
            { 1, 4, 0, 30, "DELTAQP_C_SONY_3" },
            { 1, 4, 0, 24, "DSLICE_A_HHI_5" },
            { 1, 4, 0, 24, "DSLICE_B_HHI_5" },
            { 1, 4, 0, 24, "DSLICE_C_HHI_5" },
            { 1, 6, 2, 60, "ENTP_A_Qualcomm_1" },
            { 1, 6, 2, 60, "ENTP_B_Qualcomm_1" },
            { 1, 6, 2, 60, "ENTP_C_Qualcomm_1" },
            { 1, 3, 0, 30, "EXT_A_ericsson_4" },
            { 1, 6, 2, 50, "FILLER_A_Sony_1" },
            { 1, 6, 2, 50, "HRD_A_Fujitsu_3" },
            { 1, 2, 0, 30, "INITQP_A_Sony_1" },
            { 2, 2, 0, 30, "INITQP_B_Main10_Sony_1" },
            { 1, 2, 0, 30, "ipcm_A_NEC_3" },
            { 1, 2, 0, 30, "ipcm_B_NEC_3" },
            { 1, 2, 0, 30, "ipcm_C_NEC_3" },
            { 1, 2, 0, 30, "ipcm_D_NEC_3" },
            { 1, 2, 0, 30, "ipcm_E_NEC_2" },
            { 1, 5, 1, 30, "IPRED_A_docomo_2" },
            { 1, 4, 0, 30, "IPRED_B_Nokia_3" },
            { 1, 3, 0, 30, "IPRED_C_Mitsubishi_3" },
            { 1, 5, 0, 30, "LS_A_Orange_2" },
            { 1, 5, 0, 30, "LS_B_Orange_4" },
            { 1, 6, 2, 50, "LTRPSPS_A_Qualcomm_1" },
            { 1, 2, 0, 30, "MAXBINS_A_TI_5" },
            { 1, 2, 0, 30, "MAXBINS_B_TI_5" },
            { 1, 2, 0, 30, "MAXBINS_C_TI_5" },
            { 1, 2, 0, 30, "MERGE_A_TI_3" },
            { 1, 2, 0, 30, "MERGE_B_TI_3" },
            { 1, 2, 0, 30, "MERGE_C_TI_3" },
            { 1, 2, 0, 30, "MERGE_D_TI_3" },
            { 1, 2, 0, 30, "MERGE_E_TI_3" },
            { 1, 4, 0, 50, "MERGE_F_MTK_4" },
            { 1, 3, 1, 60, "MERGE_G_HHI_4" },
            { 1, 2, 0, 30, "MVCLIP_A_qualcomm_3" },
            { 1, 4, 0, 60, "MVDL1ZERO_A_docomo_4" },
            { 1, 2, 0, 30, "MVEDGE_A_qualcomm_3" },
            { 1, 3, 0, 30, "NoOutPrior_A_Qualcomm_1" },
            { 1, 3, 0, 30, "NoOutPrior_B_Qualcomm_1" },
            { 1, 3, 0, 30, "NUT_A_ericsson_5" },
            { 1, 6, 2, 50, "OPFLAG_A_Qualcomm_1" },
            { 1, 6, 2, 50, "OPFLAG_B_Qualcomm_1" },
            { 1, 6, 2, 50, "OPFLAG_C_Qualcomm_1" },
            { 1, 5, 1, 50, "PICSIZE_A_Bossen_1" },
            { 1, 5, 1, 50, "PICSIZE_B_Bossen_1" },
            { 1, 5, 1, 50, "PICSIZE_C_Bossen_1" },
            { 1, 5, 1, 50, "PICSIZE_D_Bossen_1" },
            { 1, 2, 0, 30, "PMERGE_A_TI_3" },
            { 1, 2, 0, 30, "PMERGE_B_TI_3" },
            { 1, 2, 0, 30, "PMERGE_C_TI_3" },
            { 1, 2, 0, 30, "PMERGE_D_TI_3" },
            { 1, 2, 0, 30, "PMERGE_E_TI_3" },
            { 1, 4, 0, 50, "POC_A_Bossen_3" },
            { 1, 6, 2, 30, "PPS_A_qualcomm_7" },
            { 1, 2, 1, 50, "PS_B_VIDYO_3" },
            { 1, 2, 0, 30, "RAP_A_docomo_6" },
            { 1, 6, 2, 60, "RAP_B_Bossen_2" },
            { 1, 2, 0, 30, "RPLM_A_qualcomm_4" },
            { 1, 2, 0, 30, "RPLM_B_qualcomm_4" },
            { 1, 2, 0, 30, "RPS_A_docomo_5" },
            { 1, 3, 0, 30, "RPS_B_qualcomm_5" },
            { 1, 3, 0, 30, "RPS_C_ericsson_5" },
            { 1, 3, 0, 30, "RPS_D_ericsson_6" },
            { 1, 3, 0, 30, "RPS_E_qualcomm_5" },
            { 1, 6, 2, 60, "RPS_F_docomo_2" },
            { 1, 3, 1, 60, "RQT_A_HHI_4" },
            { 1, 3, 1, 60, "RQT_B_HHI_4" },
            { 1, 3, 1, 60, "RQT_C_HHI_4" },
            { 1, 3, 1, 60, "RQT_D_HHI_4" },
            { 1, 3, 1, 60, "RQT_E_HHI_4" },
            { 1, 3, 1, 60, "RQT_F_HHI_4" },
            { 1, 3, 1, 60, "RQT_G_HHI_4" },
            { 1, 4, 0, 60, "SAO_A_MediaTek_4" },
            { 1, 4, 0, 60, "SAO_B_MediaTek_5" },
            { 1, 4, 1, 60, "SAO_C_Samsung_5" },
            { 1, 4, 1, 60, "SAO_D_Samsung_5" },
            { 1, 4, 0, 50, "SAO_E_Canon_4" },
            { 1, 4, 0, 50, "SAO_F_Canon_3" },
            { 1, 4, 0, 50, "SAO_G_Canon_3" },
            { 1, 4, 0, 50, "SAO_H_Parabola_1" },
            { 1, 4, 1, 50, "SDH_A_Orange_4" },
            { 1, 6, 2, 30, "SLICES_A_Rovi_3" },
            { 1, 4, 0, 50, "SLIST_A_Sony_5" },
            { 1, 4, 0, 50, "SLIST_B_Sony_9" },
            { 1, 4, 0, 50, "SLIST_C_Sony_4" },
            { 1, 4, 0, 50, "SLIST_D_Sony_9" },
            { 1, 6, 2, 30, "SLPPLP_A_VIDYO_2" },
            { 1, 5, 1, 50, "STRUCT_A_Samsung_6" },
            { 1, 5, 1, 50, "STRUCT_B_Samsung_6" },
            { 1, 4, 1, 60, "TILES_A_Cisco_2" },
            { 1, 4, 1, 60, "TILES_B_Cisco_1" },
            { 1, 2, 0, 30, "TMVP_A_MS_3" },
            { 1, 2, 1, 50, "TSCL_A_VIDYO_5" },
            { 1, 2, 1, 50, "TSCL_B_VIDYO_4" },
            { 1, 3, 1, 30, "TSKIP_A_MS_3" },
            { 2, 5, 1, 30, "TSUNEQBD_A_MAIN10_Technicolor_2" },
            { 1, 5, 0, 30, "TUSIZE_A_Samsung_1" },
            { 1, 6, 2, 30, "VPSID_A_VIDYO_2" },
            { 2, 6, 0, 60, "WP_A_MAIN10_Toshiba_3" },
            { 1, 6, 0, 60, "WP_A_Toshiba_3" },
            { 1, 6, 0, 60, "WP_B_Toshiba_3" },
            { 2, 6, 0, 60, "WP_MAIN10_B_Toshiba_3" },
            { 2, 2, 0, 50, "WPP_A_ericsson_MAIN10_2" },
            { 1, 2, 0, 50, "WPP_A_ericsson_MAIN_2" },
            { 2, 2, 0, 50, "WPP_B_ericsson_MAIN10_2" },
            { 1, 2, 0, 50, "WPP_B_ericsson_MAIN_2" },
            { 2, 2, 0, 50, "WPP_C_ericsson_MAIN10_2" },
            { 1, 2, 0, 50, "WPP_C_ericsson_MAIN_2" },
            { 2, 2, 0, 50, "WPP_D_ericsson_MAIN10_2" },
            { 1, 2, 0, 50, "WPP_D_ericsson_MAIN_2" },
            { 2, 2, 0, 50, "WPP_E_ericsson_MAIN10_2" },
            { 1, 2, 0, 50, "WPP_E_ericsson_MAIN_2" },
            { 2, 2, 0, 50, "WPP_F_ericsson_MAIN10_2" },
            { 1, 2, 0, 50, "WPP_F_ericsson_MAIN_2" },
            //{ 4, 2, 0, 24, "ADJUST_IPRED_ANGLE_A_RExt_Mitsubishi_1" },
            //{ 4, 4, 1, 60, "Bitdepth_A_RExt_Sony_1" },
            //{ 4, 4, 1, 60, "Bitdepth_B_RExt_Sony_1" },
            //{ 4, 2, 0, 24, "CCP_10bit_RExt_QCOM_1" },
            //{ 4, 2, 0, 30, "CCP_12bit_RExt_QCOM_1" },
            //{ 4, 2, 0, 30, "CCP_8bit_RExt_QCOM_1" },
            //{ 4, 6, 2, 50, "ExplicitRdpcm_A_BBC_1" },
            //{ 4, 6, 2, 30, "ExplicitRdpcm_B_BBC_2" },
            //{ 4, 0, 0, 0, "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_10BIT_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_12BIT_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_16BIT_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "EXTPREC_HIGHTHROUGHPUT_444_16_INTRA_8BIT_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "EXTPREC_MAIN_444_16_INTRA_10BIT_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "EXTPREC_MAIN_444_16_INTRA_12BIT_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "EXTPREC_MAIN_444_16_INTRA_16BIT_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "EXTPREC_MAIN_444_16_INTRA_8BIT_RExt_Sony_1" },
            { 4, 0, 0, 0, "GENERAL_10b_420_RExt_Sony_1" },
            { 4, 0, 0, 0, "GENERAL_10b_422_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_10b_444_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_12b_400_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_12b_420_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_12b_422_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_12b_444_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_16b_400_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_16b_444_highThroughput_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_16b_444_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_8b_400_RExt_Sony_1" },
            { 4, 0, 0, 0, "GENERAL_8b_420_RExt_Sony_1" },
            //{ 4, 0, 0, 0, "GENERAL_8b_444_RExt_Sony_1" },
            //{ 4, 2, 0, 30, "IPCM_A_RExt_NEC_1" },
            //{ 4, 2, 0, 30, "IPCM_B_RExt_NEC_1" },
            { 4, 4, 0, 24, "Main_422_10_A_RExt_Sony_2" },
            { 4, 5, 0, 30, "Main_422_10_B_RExt_Sony_2" },
            //{ 4, 2, 0, 30, "PERSIST_RPARAM_A_RExt_Sony_2" },
            //{ 4, 4, 0, 20, "QMATRIX_A_RExt_Sony_1" },
            //{ 4, 6, 2, 30, "SAO_A_RExt_MediaTek_1" },
            //{ 4, 2, 0, 30, "TSCTX_10bit_RExt_SHARP_1" },
            //{ 4, 2, 0, 30, "TSCTX_10bit_I_RExt_SHARP_1" },
            //{ 4, 2, 0, 30, "TSCTX_12bit_RExt_SHARP_1" },
            //{ 4, 2, 0, 30, "TSCTX_12bit_I_RExt_SHARP_1" },
            //{ 4, 2, 0, 30, "TSCTX_8bit_RExt_SHARP_1" },
            //{ 4, 2, 0, 30, "TSCTX_8bit_I_RExt_SHARP_1" },
            //{ 4, 0, 0, 0, "WAVETILES_RExt_Sony_1" },
    };

    struct Locator
    {
        Locator(const char *streamName, const char *neptunePath)
        {
            boost::filesystem::path conformanceStreamPath = neptunePath;
            conformanceStreamPath /= "conformance";
            conformanceStreamPath /= streamName;

            this->bitstream = conformanceStreamPath;
            this->bitstream += ".bin";
            if (!boost::filesystem::exists(this->bitstream))
            {
                this->bitstream = conformanceStreamPath;
                this->bitstream += ".bit";
            }

            if (!boost::filesystem::exists(this->bitstream))
            {
                std::cerr << "can't find " << streamName << "\n";
                throw std::runtime_error("conformance stream missing");
            }

            this->md5 = this->bitstream.parent_path() / "md5" / this->bitstream.filename();
            this->md5 += ".md5";

            if (!boost::filesystem::exists(this->md5))
            {
                std::cout << this->md5 << " missing\n";
            }
        }
        boost::filesystem::path bitstream;
        boost::filesystem::path md5;
    };

}

#endif
