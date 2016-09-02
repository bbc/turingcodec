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

#ifndef INCLUDED_HevcTypes_h
#define INCLUDED_HevcTypes_h

#pragma once

// Some HEVC-specific enumerations, e.g. NalUnitType and associated functions
// C-preprocessor "X-macros" are used here, more information here: http://stackoverflow.com/questions/6635851/real-world-use-of-x-macros

#include <cassert>
#include <iostream>


#define NAL_UNIT_TYPES \
        X(0, TRAIL_N, "Coded slice segment of a non-TSA, non-STSA trailing picture", slice_segment_layer_rbsp, VCL) \
        X(1, TRAIL_R, "Coded slice segment of a non-TSA, non-STSA trailing picture", slice_segment_layer_rbsp, VCL) \
        X(2, TSA_N, "Coded slice segment of a TSA picture", slice_segment_layer_rbsp, VCL) \
        X(3, TSA_R, "Coded slice segment of a TSA picture", slice_segment_layer_rbsp, VCL) \
        X(4, STSA_N, "Coded slice segment of a STSA picture", slice_segment_layer_rbsp, VCL) \
        X(5, STSA_R, "Coded slice segment of a STSA picture", slice_segment_layer_rbsp, VCL) \
        X(6, RADL_N, "Coded slice segment of a RADL picture", slice_segment_layer_rbsp, VCL) \
        X(7, RADL_R, "Coded slice segment of a RADL picture", slice_segment_layer_rbsp, VCL) \
        X(8, RASL_N, "Coded slice segment of a RASL picture", slice_segment_layer_rbsp, VCL) \
        X(9, RASL_R, "Coded slice segment of a RASL picture", slice_segment_layer_rbsp, VCL) \
        X(10, RSV_VCL_N10, "Reserved non-IRAP non-reference VCL NAL unit types", RbspReserved, VCL) \
        X(12, RSV_VCL_N12, "Reserved non-IRAP non-reference VCL NAL unit types", RbspReserved, VCL) \
        X(14, RSV_VCL_N14, "Reserved non-IRAP non-reference VCL NAL unit types", RbspReserved, VCL) \
        X(11, RSV_VCL_N11, "Reserved non-IRAP reference VCL NAL unit types", RbspReserved, VCL) \
        X(13, RSV_VCL_N13, "Reserved non-IRAP reference VCL NAL unit types", RbspReserved, VCL) \
        X(15, RSV_VCL_N15, "Reserved non-IRAP reference VCL NAL unit types", RbspReserved, VCL) \
        X(16, BLA_W_LP, "Coded slice segment of a BLA picture", slice_segment_layer_rbsp, VCL) \
        X(17, BLA_W_RADL, "Coded slice segment of a BLA picture", slice_segment_layer_rbsp, VCL) \
        X(18, BLA_N_LP, "Coded slice segment of a BLA picture", slice_segment_layer_rbsp, VCL) \
        X(19, IDR_W_RADL, "Coded slice segment of an IDR picture", slice_segment_layer_rbsp, VCL) \
        X(20, IDR_N_LP, "Coded slice segment of an IDR picture", slice_segment_layer_rbsp, VCL) \
        X(21, CRA_NUT, "Coded slice segment of a CRA picture", slice_segment_layer_rbsp, VCL) \
        X(22, RSV_IRAP_VCL22, "Reserved IRAP VCL NAL unit types", RbspReserved, VCL) \
        X(23, RSV_IRAP_VCL23, "Reserved IRAP VCL NAL unit types", RbspReserved, VCL) \
        X(24, RSV_VCL24, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(25, RSV_VCL25, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(26, RSV_VCL26, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(27, RSV_VCL27, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(28, RSV_VCL28, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(29, RSV_VCL29, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(30, RSV_VCL30, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(31, RSV_VCL31, "Reserved non-IRAP VCL NAL unit ", RbspReserved, VCL) \
        X(32, VPS_NUT, "Video parameter set", video_parameter_set_rbsp, non-VCL) \
        X(33, SPS_NUT, "Sequence parameter set", seq_parameter_set_rbsp, non-VCL) \
        X(34, PPS_NUT, "Picture<uint8_t> parameter set", pic_parameter_set_rbsp, non-VCL) \
        X(35, AUD_NUT, "Access unit delimiter", access_unit_delimiter_rbsp, non-VCL) \
        X(36, EOS_NUT, "End of sequence", end_of_seq_rbsp, non-VCL) \
        X(37, EOB_NUT, "End of bitsteam", end_of_bitstream_rbsp, non-VCL) \
        X(38, FD_NUT, "Filler data", filler_data_rbsp, non-VCL) \
        X(39, PREFIX_SEI_NUT, "Supplemental enhancement information (SEI)", sei_rbsp, non-VCL) \
        X(40, SUFFIX_SEI_NUT, "Supplemental enhancement information (SEI)", sei_rbsp, non-VCL) \
        X(41, RSV_NVCL41, "Reserved", RbspReserved, non-VCL) \
        X(42, RSV_NVCL42, "Reserved", RbspReserved, non-VCL) \
        X(43, RSV_NVCL43, "Reserved", RbspReserved, non-VCL) \
        X(44, RSV_NVCL44, "Reserved", RbspReserved, non-VCL) \
        X(45, RSV_NVCL45, "Reserved", RbspReserved, non-VCL) \
        X(46, RSV_NVCL46, "Reserved", RbspReserved, non-VCL) \
        X(47, RSV_NVCL47, "Reserved", RbspReserved, non-VCL) \
        X(48, UNSPEC48, "Unspecified", RbspUnspecified, non-VCL) \
        X(49, UNSPEC49, "Unspecified", RbspUnspecified, non-VCL) \
        X(50, UNSPEC50, "Unspecified", RbspUnspecified, non-VCL) \
        X(51, UNSPEC51, "Unspecified", RbspUnspecified, non-VCL) \
        X(52, UNSPEC52, "Unspecified", RbspUnspecified, non-VCL) \
        X(53, UNSPEC53, "Unspecified", RbspUnspecified, non-VCL) \
        X(54, UNSPEC54, "Unspecified", RbspUnspecified, non-VCL) \
        X(55, UNSPEC55, "Unspecified", RbspUnspecified, non-VCL) \
        X(56, UNSPEC56, "Unspecified", RbspUnspecified, non-VCL) \
        X(57, UNSPEC57, "Unspecified", RbspUnspecified, non-VCL) \
        X(58, UNSPEC58, "Unspecified", RbspUnspecified, non-VCL) \
        X(59, UNSPEC59, "Unspecified", RbspUnspecified, non-VCL) \
        X(60, UNSPEC60, "Unspecified", RbspUnspecified, non-VCL) \
        X(61, UNSPEC61, "Unspecified", RbspUnspecified, non-VCL) \
        X(62, UNSPEC62, "Unspecified", RbspUnspecified, non-VCL) \
        X(63, UNSPEC63, "Unspecified", RbspUnspecified, non-VCL)


enum NalUnitType
{
#define X(TYPE, NAME, DESCRIPTION, ELEMENT, VCLMSG) NAME = TYPE,
    NAL_UNIT_TYPES
#undef X
};


//review: antisocial function name
static const char *toString(int nut)
{
    switch (nut)
    {
#define X(TYPE, NAME, DESCRIPTION, ELEMENT, VCLMSG) case TYPE: return #NAME;
        NAL_UNIT_TYPES
#undef X
    }
    assert(false);
    return "";
#undef X
}

static const char *nalUnitTypeDescription(int nut)
{
    switch (nut)
    {
#define X(TYPE, NAME, DESCRIPTION, ELEMENT, VCLMSG) case TYPE: return DESCRIPTION;
        NAL_UNIT_TYPES
#undef X
    }
    assert(false);
    return "";
#undef X
}

static bool isVcl(int nut)
{
    const int VCL = 1;
    const int non = 1;
    switch (nut)
    {
#define X(TYPE, NAME, DESCRIPTION, ELEMENT, VCLMSG) case TYPE: return !!(VCLMSG);
        NAL_UNIT_TYPES
#undef X
    }
    assert(false);
    return false;
}

static bool isSliceSegment(int nut)
{
    return nut >= TRAIL_N && nut <= RSV_VCL31;
}

static bool isIdr(int nut)
{
    return nut == IDR_W_RADL || nut == IDR_N_LP;
}

static bool isIrap(int nut)
{
    return nut >= BLA_W_LP && nut <= RSV_IRAP_VCL23;
}

static bool isBla(int nut)
{
    return nut == BLA_W_LP || nut == BLA_W_RADL || nut == BLA_N_LP;
}

static bool isCra(int nut)
{
    return nut == CRA_NUT;
}


static bool isRasl(int nut)
{
    return nut == RASL_R || nut == RASL_N;
}

static bool isRadl(int nut)
{
    return nut == RADL_R || nut == RADL_N;
}

static bool isTsa(int nut)
{
    return nut == TSA_R || nut == TSA_N;
}

static bool isStsa(int nut)
{
    return nut == STSA_R || nut == STSA_N;
}

static bool isLeadingPicture(int nut)
{
    return isRasl(nut) || isRadl(nut);
}

static bool isSubLayerNonReferencePicture(int nut)
{
    switch (nut)
    {
        case TRAIL_N:
        case TSA_N:
        case STSA_N:
        case RADL_N:
        case RASL_N:
        case RSV_VCL_N10:
        case RSV_VCL_N12:
        case RSV_VCL_N14:
            return true;
        default:
            return false;
    }
}

static bool startsNewAccessUnit(int nut)
{
    switch (nut)
    {
        case AUD_NUT:
        case VPS_NUT:
        case SPS_NUT:
        case PPS_NUT:
        case PREFIX_SEI_NUT:
        case RSV_NVCL41:
        case RSV_NVCL42:
        case RSV_NVCL43:
        case RSV_NVCL44:
        case UNSPEC48:
        case UNSPEC49:
        case UNSPEC50:
        case UNSPEC51:
        case UNSPEC52:
        case UNSPEC53:
        case UNSPEC54:
        case UNSPEC55:
            return true;
        default:
            return false;
    }
}


#define REFERENCE_TYPES \
    X(UNUSED, "unused for reference") \
    X(SHORT_TERM, "used for short-term reference") \
    X(LONG_TERM, "used for long-term reference")

enum Reference
{
#define X(A, B) A,
    REFERENCE_TYPES
#undef X
};	

static char const* referenceName(Reference reference)
{
    switch (reference)
    {
    default: return 0;
#define X(A, B) case A: return B;
        REFERENCE_TYPES
#undef X
    }
}


#define SLICE_TYPES \
        X(I, 2) \
        X(P, 1) \
        X(B, 0)

#define X(A, B) static const int A = B;
SLICE_TYPES
#undef X

static char sliceTypeToChar(int sliceType)
{
#define X(A, B) \
        if (sliceType == B) return #A [0];

    SLICE_TYPES
#undef X
    assert(0); return 0;
}


#define CU_PRED_MODE_VALUES \
        X(MODE_INTER, 0) \
        X(MODE_INTRA, 1) \
        X(MODE_SKIP, 2)

#define X(A, B) \
        static const int A = B;

CU_PRED_MODE_VALUES
#undef X

#define PART_MODE_VALUES \
        X(PART_2Nx2N, 0) \
        X(PART_2NxN, 1) \
        X(PART_Nx2N, 2) \
        X(PART_NxN, 3) \
        X(PART_2NxnU, 4) \
        X(PART_2NxnD, 5) \
        X(PART_nLx2N, 6) \
        X(PART_nRx2N, 7)

enum PartModeType
{
#define X(A, B) A = B,
    PART_MODE_VALUES
#undef X
};

static bool isAmp(PartModeType partMode)
{
    return partMode >= PART_2NxnU;
}



// inter_pred_idc name association
static const int PRED_L0 = 0;
static const int PRED_L1 = 1;
static const int PRED_BI = 2;



#define INTRA_PRED_MODE_VALUES \
        X(INTRA_PLANAR, 0) \
        X(INTRA_DC, 1) \
        X(INTRA_ANGULAR2, 2) \
        X(INTRA_ANGULAR3, 3) \
        X(INTRA_ANGULAR4, 4) \
        X(INTRA_ANGULAR5, 5) \
        X(INTRA_ANGULAR6, 6) \
        X(INTRA_ANGULAR7, 7) \
        X(INTRA_ANGULAR8, 8) \
        X(INTRA_ANGULAR9, 9) \
        X(INTRA_ANGULAR10, 10) \
        X(INTRA_ANGULAR11, 11) \
        X(INTRA_ANGULAR12, 12) \
        X(INTRA_ANGULAR13, 13) \
        X(INTRA_ANGULAR14, 14) \
        X(INTRA_ANGULAR15, 15) \
        X(INTRA_ANGULAR16, 16) \
        X(INTRA_ANGULAR17, 17) \
        X(INTRA_ANGULAR18, 18) \
        X(INTRA_ANGULAR19, 19) \
        X(INTRA_ANGULAR20, 20) \
        X(INTRA_ANGULAR21, 21) \
        X(INTRA_ANGULAR22, 22) \
        X(INTRA_ANGULAR23, 23) \
        X(INTRA_ANGULAR24, 24) \
        X(INTRA_ANGULAR25, 25) \
        X(INTRA_ANGULAR26, 26) \
        X(INTRA_ANGULAR27, 27) \
        X(INTRA_ANGULAR28, 28) \
        X(INTRA_ANGULAR29, 29) \
        X(INTRA_ANGULAR30, 30) \
        X(INTRA_ANGULAR31, 31) \
        X(INTRA_ANGULAR32, 32) \
        X(INTRA_ANGULAR33, 33) \
        X(INTRA_ANGULAR34, 34)

#define X(A, B) \
        static const int A = B;

INTRA_PRED_MODE_VALUES
#undef X


static const char *partModeStr(int partMode)
{
    switch (partMode)
    {
#define X(A, B) case B: return #A;
        PART_MODE_VALUES
#undef X
        default:;
    }
    return 0;
};


static const char *cuPredModeStr(int CuPredMode)
{
    switch (CuPredMode)
    {
#define X(A, B) case B: return #A;
        CU_PRED_MODE_VALUES
#undef X
        default:;
    }
    return 0;
};


#define CHROMA_FORMATS \
        X(0, 0, "monochrome", 1, 1) \
        X(1, 0, "4:2:0", 2, 2) \
        X(2, 0, "4:2:2", 2, 1) \
        X(3, 0, "4:4:4", 1, 1) \
        X(3, 1, "4:4:4", 1, 1) \


static int subWidthC(int i)
{
#define X(chroma_format_idc, separate_color_plane_flag, name, SubWidthC, SubHeightC) \
        if (i == chroma_format_idc) return SubWidthC;

    CHROMA_FORMATS
#undef X
    assert(0);
    return 0;
}

static int subHeightC(int i)
{
#define X(chroma_format_idc, separate_color_plane_flag, name, SubWidthC, SubHeightC) \
        if (i == chroma_format_idc) return SubHeightC;

    CHROMA_FORMATS
#undef X
    assert(0);
    return 0;
}

// represents reference list L0
static const int L0 = 0;

// represents reference list L1
static const int L1 = 1;

#endif
