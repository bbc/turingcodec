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

#ifndef INCLUDED_StatePictures_h
#define INCLUDED_StatePictures_h

#pragma once

#include "StatePicture.h"
#include "StateCollocatedMotion.h"
#include "Handlers.h"
#include "Global.h"
#include "Padding.h"
#include "LoopFilter.h"
#include "Profiler.h"
#include <boost/iterator/filter_iterator.hpp>
#include <sstream>


struct Abort;

template <class> struct Table;
struct Sps;
struct Vps;
struct Pps;

// invoked as soon as a picture is decoded
struct JustDecoded
{
    JustDecoded(PicOrderCnt poc) : poc(poc) { }
    PicOrderCnt poc;
};

// invoked when a picture is to be displayed
struct OutputPicture
{
    OutputPicture(PicOrderCnt poc) : poc(poc) { }
    PicOrderCnt poc;
};

// invoked just before a picture is removed from the DPB
struct DeletePicture
{
    DeletePicture(PicOrderCnt poc) : poc(poc) { }
    PicOrderCnt poc;
};

// invoked when DPB has been emptied
struct DpbClear { };

// Request a newly allocated picture: this will be returned as a std::shared_ptr<S> where S is derived from StatePicture
struct NewPicture
{
};


static inline bool comparePicOrderCntVal(std::shared_ptr<StatePicture> a, std::shared_ptr<StatePicture> b)
{
    return (*a)[PicOrderCntVal()] < (*b)[PicOrderCntVal()];
}


//#define DEBUG_DPB

struct StatePicturesBase
{
    typedef StatePicture Pic;

    StatePicturesBase() :
        sliceHeaderValid(),
        codedVideoSequenceId(-1),
        picture0(true),
        justSeenEndOfSeq(false),
        skipSliceSegmentData(false),
        sequenceDecodeOrder(0),
        pictureWanted(true),
        deblockWholeFrame(false)
    {
    }

    typedef std::array<int, 16> PocList;
    PocList PocStCurrBefore, PocStCurrAfter, PocStFoll, PocLtCurr, PocLtFoll;
    int NumPocStCurrBefore, NumPocStCurrAfter, NumPocStFoll, NumPocLtCurr, NumPocLtFoll;
    std::vector<int> PocLsbLt;
    std::vector<int> UsedByCurrPicLt;
    std::vector<int> DeltaPocMSBCycleLt;
    std::vector<int> CurrDeltaPocMsbPresentFlag;
    std::vector<int> FollDeltaPocMsbPresentFlag;

    // review: remove these two:
    bool pictureWanted;
    bool deblockWholeFrame;

    uint64_t sequenceDecodeOrder;
    size_t posEndCvs;

    bool sliceHeaderValid;
    int sliceNalUnitType;
    int sliceType;
    int picWidth;
    int picHeight;
    int idVps;
    int idSps;
    int idPps;
    //int sliceCount[3];
    bool picture0; // should be able to use codedVideoSequenceId to avoid need for this extra variable
    bool justSeenEndOfSeq;

    int QPY;
    int qPY_PREV;
    int qPY_A;
    int qPY_B;
    int qPY_PRED;

    bool skipSliceSegmentData;

    struct PrevTid0Pic
    {
        PrevTid0Pic() : good(false) { }
        bool good;
        int pic_order_cnt_lsb;
        int PicOrderCntMsb;
        int PicOrderCntVal;
    };
    PrevTid0Pic prevTid0Pic;

    typedef int CodedVideoSequenceId;
    CodedVideoSequenceId codedVideoSequenceId;

    int NoRaslOutputFlag;

    int pocCollocated;

    const char *streamType;
    static const char *streamTypeNal() { static const char *n = "NAL"; return n; }
    static const char *streamTypeRbsp() { static const char *n = "RBSP"; return n; }
    static const char *streamTypeCabac() { static const char *n = "CABAC"; return n; }
    static const char *streamTypeSei() { static const char *n = "SEI"; return n; }
};


template <class Sample>
void fillRectangle(Raster<Sample> p, Sample value, int width, int height)
{
    for (int y = 0; y<height; ++y)
    {
        for (int x = 0; x<width; ++x)
        {
            p(x, y) = value;
        }
    }
}


template <class Sample>
static void padPicture(Picture<Sample> &picture)
{
    const int pad = 80; // todo: be more intelligent
    Padding::padImage<Sample>(&picture[0](0, 0), picture[0].width, picture[0].height, (int)picture[0].stride, pad);
    Padding::padImage<Sample>(&picture[1](0, 0), picture[1].width, picture[1].height, (int)picture[1].stride, pad / 2);
    Padding::padImage<Sample>(&picture[2](0, 0), picture[2].width, picture[2].height, (int)picture[2].stride, pad / 2);
}


template <class Picture, class Enable=void>
struct SetupReconstructedPicture
{
    template <class H> static void go(Picture &dp, H &h) { };
};

template <class Picture, class Enable = void>
struct SetupSaoPicture
{
    template <class H> static void go(Picture &dp, H &h) { };
};

template <class Picture, class Enable = void>
struct SetupDeblockPicture
{
    template <class H> static void go(Picture &dp, H &h) { };
};

template <class Picture>
struct SetupReconstructedPicture<Picture, typename std::enable_if<std::is_base_of<ReconstructedPictureBase, Picture>::value>::type>
{
    typedef typename Picture::Sample Sample;

    template <class H> static void go(StateReconstructedPicture<Sample> &dp, H &h)
    {
        const int pad = 96;// h[CtbSizeY()] + 16; // review: less padding will suffice
        ::Picture<Sample> *picture = new ::Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);
        dp.picture.reset(picture);
    }
};

template <class Picture>
struct SetupSaoPicture<Picture, typename std::enable_if<std::is_base_of<ReconstructedPictureBase, Picture>::value>::type>
{
    typedef typename Picture::Sample Sample;

    template <class H> static void go(StateReconstructedPicture<Sample> &dp, H &h)
    {
        const int pad = 96;// h[CtbSizeY()] + 16; // review: less padding will suffice
        ::Picture<Sample> *picture = new ::Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);
        dp.saoPicture.reset(picture);
    }
};

template <class Picture>
struct SetupDeblockPicture<Picture, typename std::enable_if<std::is_base_of<ReconstructedPictureBase, Picture>::value>::type>
{
    typedef typename Picture::Sample Sample;

    template <class H> static void go(StateReconstructedPicture<Sample> &dp, H &h)
    {
        const int pad = 96;// h[CtbSizeY()] + 16; // review: less padding will suffice
        ::Picture<Sample> *picture = new ::Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);
        dp.deblockPicture.reset(picture);
    }
};

template <class Picture, class H>
void setupReconstructedPicture(Picture &dp, H &h)
{
    SetupReconstructedPicture<Picture>::go(dp, h);
    if (h[sample_adaptive_offset_enabled_flag()])
    {
        SetupSaoPicture<Picture>::go(dp, h);
        SetupDeblockPicture<Picture>::go(dp, h);
    }
};

template <class Picture> struct GeneratePicture
{
    template <class H> static void go(Picture &dp, H &h)
    {
    }
};

template <typename Sample> struct GeneratePicture<StateReconstructedPicture<Sample>>
{
    template <class H> static void go(StateReconstructedPicture<Sample> &dp, H &h)
    {
        const int pad = 96;// h[CtbSizeY()] + 16; // review: less padding will suffice
        Picture<Sample> *picture = new Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);
        dp.reconstructedPicture.reset(picture);

        for (int cIdx = 0; cIdx<3; ++cIdx)
        {
            int bitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
            Raster<Sample> samples = (*dp.reconstructedPicture)[cIdx];
            fillRectangle<Sample>(samples, 1 << (bitDepth - 1), (*dp.reconstructedPicture)[cIdx].width, (*dp.reconstructedPicture)[cIdx].height);
        }
    }
};


struct StatePictures :
    StatePicturesBase
{
    // Derived classes may override this function to return a different picture type
    virtual StatePicture *newPicture()
    {
        return new StatePicture;
    }

    std::array<std::shared_ptr<StatePicture>, 16> RefPicSetStCurrBefore, RefPicSetStCurrAfter, RefPicSetLtCurr;

    template <class H>
    void generatePicture(H &h, int poc, Reference reference)
    {
        // std::cout << "generating POC=" << poc << "\n";
        this->dpb.push_back(std::shared_ptr<StatePicture>(this->newPicture()));
        Pic &dp = *this->dpb.back();
        this->setupDecodedPicture(dp, h);
        dp[PicOrderCntVal()] = poc;
        dp.reference = reference;
        dp.neededForOutput = false;
        dp.reconstructed = true;
        dp.generatedUnavailable = true;

        auto *p = &h[Concrete<StatePicture>()];
        GeneratePicture<decltype(*p)>::go(*p, h);
    }

    // Decoding process for picture order count
    template <class H>
    void decodePictureOrderCount(H &h)
    {
        int PicOrderCntMsb = 0;

        // suspect this is never set - problem?
        const bool firstPictureInStream = (this->codedVideoSequenceId < 0);


        if (
                !isIdr(h[nal_unit_type()]) &&
                !isBla(h[nal_unit_type()]) &&
                !(isIdr(h[nal_unit_type()]) && firstPictureInStream) &&
                this->prevTid0Pic.good)
        {
            assert(prevTid0Pic.pic_order_cnt_lsb >= 0);
            const int prevPicOrderCntLsb  = prevTid0Pic.pic_order_cnt_lsb;
            const int prevPicOrderCntMsb = prevTid0Pic.PicOrderCntMsb;

            if ( ( h[slice_pic_order_cnt_lsb()] <  prevPicOrderCntLsb )  &&
                    ( ( prevPicOrderCntLsb -  h[slice_pic_order_cnt_lsb()] )  >=  ( h[MaxPicOrderCntLsb()] / 2 ) ) )
            {
                PicOrderCntMsb = prevPicOrderCntMsb + h[MaxPicOrderCntLsb()];
            }
            else if( (h[slice_pic_order_cnt_lsb()]  >  prevPicOrderCntLsb )  &&
                    ( (h[slice_pic_order_cnt_lsb()] - prevPicOrderCntLsb )  >  ( h[MaxPicOrderCntLsb()] / 2 ) ) )
            {
                PicOrderCntMsb = prevPicOrderCntMsb - h[MaxPicOrderCntLsb()];
            }
            else
            {
                PicOrderCntMsb = prevPicOrderCntMsb;
            }
        }

        int64_t poc = int64_t(PicOrderCntMsb) + int64_t(h[slice_pic_order_cnt_lsb()]);

        if (poc >= (int64_t(1) << 31) || poc < -(int64_t(1) << 31))
        {
            h(Violation("8.3.1", "PicOrderCntVal{%1%} is not in the range of -2^31 to 2^31 - 1 , inclusive") % poc); // CondCheck 8.3.1-A
        }

        h[PicOrderCntVal()] = PicOrderCntMsb + h[slice_pic_order_cnt_lsb()];

        if (h[TemporalId()] == 0
                && !isRasl(h[nal_unit_type()])
                && !isRadl(h[nal_unit_type()])
                && !isSubLayerNonReferencePicture(h[nal_unit_type()]))
        {
            this->prevTid0Pic.pic_order_cnt_lsb = h[slice_pic_order_cnt_lsb()];
            this->prevTid0Pic.PicOrderCntMsb = PicOrderCntMsb;
            this->prevTid0Pic.PicOrderCntVal = h[PicOrderCntVal()];
            this->prevTid0Pic.good = true;
        }
    }

    // Decoding process for reference picture set
    template <class H>
    void decodeReferencePictureSet(H &h)
    {
        // When the current picture is an IRAP picture with NoRaslOutputFlag equal to 1,
        if (isIrap(h[nal_unit_type()]) && NoRaslOutputFlag == 1)
        {
            // all reference pictures currently in the DPB (if any) are marked as "unused for reference".
            for (auto &dp : this->dpb)
            {
                dp->reference = UNUSED;
            }
        }

        {
            if (isIdr(h[nal_unit_type()]) || isBla(h[nal_unit_type()]))
            {
                NumPocStCurrBefore = NumPocStCurrAfter = NumPocStFoll = NumPocLtCurr = NumPocLtFoll = 0;
            }
            else
            {
                int i, j, k;
                for (i = 0, j = 0, k = 0; i < h[NumNegativePics(h[CurrRpsIdx()])]; i++)
                {
                    if (h[UsedByCurrPicS0(h[CurrRpsIdx()], i)])
                    {
                        PocStCurrBefore[j++] = h[PicOrderCntVal()] + h[DeltaPocS0(h[CurrRpsIdx()], i)];
                    }
                    else
                    {
                        PocStFoll[k++] = h[PicOrderCntVal()] + h[DeltaPocS0(h[CurrRpsIdx()], i)];
                    }
                }
                NumPocStCurrBefore = j;
                for (i = 0, j = 0; i < h[NumPositivePics(h[CurrRpsIdx()])]; i++)
                {
                    if (h[UsedByCurrPicS1(h[CurrRpsIdx()], i)])
                    {
                        int x = PocStCurrAfter[j++] = h[PicOrderCntVal()] + h[DeltaPocS1(h[CurrRpsIdx()], i)];
                        if (x == 53)
                        {
                            x = 53;
                        }
                    }
                    else
                    {
                        PocStFoll[k++] = h[PicOrderCntVal()] + h[DeltaPocS1(h[CurrRpsIdx()], i)];
                    }
                }
                NumPocStCurrAfter = j;
                NumPocStFoll = k;

                PocLsbLt.resize(h[num_long_term_sps()] + h[num_long_term_pics()]);
                UsedByCurrPicLt.resize(h[num_long_term_sps()] + h[num_long_term_pics()]);
                DeltaPocMSBCycleLt.resize(h[num_long_term_sps()] + h[num_long_term_pics()]);
                CurrDeltaPocMsbPresentFlag.resize(h[num_long_term_sps()] + h[num_long_term_pics()]);
                FollDeltaPocMsbPresentFlag.resize(h[num_long_term_sps()] + h[num_long_term_pics()]);

                for (i = 0, j = 0, k = 0; i < h[num_long_term_sps()] + h[num_long_term_pics()]; i++)
                {
                    if (i < h[num_long_term_sps()])
                    {
                        PocLsbLt[i] = h[lt_ref_pic_poc_lsb_sps(h[lt_idx_sps(i)])];
                        UsedByCurrPicLt[i] = h[used_by_curr_pic_lt_sps_flag(h[lt_idx_sps(i)])];
                    }
                    else
                    {
                        PocLsbLt[i] = h[poc_lsb_lt(i)];
                        UsedByCurrPicLt[i] = h[used_by_curr_pic_lt_flag(i)];
                    }
                    if (i == 0 || (i == h[num_long_term_sps()] && !h[delta_poc_msb_present_flag(i)]))
                    {
                        DeltaPocMSBCycleLt[i] = h[delta_poc_msb_cycle_lt(i)];
                    }
                    else
                    {
                        DeltaPocMSBCycleLt[i] = h[delta_poc_msb_cycle_lt(i)] + DeltaPocMSBCycleLt[i - 1];
                    }

                    int pocLt = PocLsbLt[i];
                    if (h[delta_poc_msb_present_flag(i)])
                    {
                        pocLt += h[PicOrderCntVal()] - DeltaPocMSBCycleLt[i] * h[MaxPicOrderCntLsb()] - h[slice_pic_order_cnt_lsb()];
                    }
                    else
                    {
                        // Not in standard.
                        // pocLt is not a PicOrderCntVal, it is the LSB part.
                        // Look in DPB to find actual PicOrderCntVal
                        StatePicture *p = getPicByPoc(pocLt, h[MaxPicOrderCntLsb()]).get();
                        if (p)
                        {
                            pocLt = (*p)[PicOrderCntVal()];
                        }
                    }

                    if (UsedByCurrPicLt[i])
                    {
                        PocLtCurr[j] = pocLt;
                        CurrDeltaPocMsbPresentFlag[j++] = h[delta_poc_msb_present_flag(i)];
                    }
                    else
                    {
                        PocLtFoll[k] = pocLt;
                        FollDeltaPocMsbPresentFlag[k++] = h[delta_poc_msb_present_flag(i)];
                    }
                }
                NumPocLtCurr = j;
                NumPocLtFoll = k;
            }
        }

        RefPicSetStCurrBefore = std::array<std::shared_ptr<Pic>, 16>();
        for (int i = 0; i < NumPocStCurrBefore; ++i)
        {
            RefPicSetStCurrBefore[i] = this->getPicByPoc(PocStCurrBefore[i]);
            //			assert(RefPicSetStCurrBefore[i]);
        }

        RefPicSetStCurrAfter = std::array<std::shared_ptr<Pic>, 16>();
        for (int i = 0; i < NumPocStCurrAfter; ++i)
        {
            RefPicSetStCurrAfter[i] = this->getPicByPoc(PocStCurrAfter[i]);
            //			assert(RefPicSetStCurrAfter[i]);
        }

        RefPicSetLtCurr = std::array<std::shared_ptr<Pic>, 16>();
        for (int i = 0; i < NumPocLtCurr; ++i)
        {
            RefPicSetLtCurr[i] = this->getPicByPoc(PocLtCurr[i]);
            //			assert(RefPicSetLtCurr[i]);
        }

        //std::ostringstream o;
        // All reference pictures in the decoded picture buffer that are not included in RefPicSetLtCurr, RefPicSetLtFoll, RefPicSetStCurrBefore, RefPicSetStCurrAfter or RefPicSetStFoll are marked as "unused for reference".
        for (auto &p : this->dpb)
        {
            auto &dp = *p;

            //o << dp[PicOrderCntVal()] << ", ";
            dp.reference = UNUSED;
            bool usedByCurrentPicture = false;

            for (int i=0; i<NumPocStCurrBefore; ++i)
            {
                if (dp[PicOrderCntVal()] == PocStCurrBefore[i])
                {
                    dp.reference = SHORT_TERM;
                    usedByCurrentPicture = true;
                }
            }
            for (int i=0; i<NumPocStCurrAfter; ++i)
            {
                if (dp[PicOrderCntVal()] == PocStCurrAfter[i])
                {
                    dp.reference = SHORT_TERM;
                    usedByCurrentPicture = true;
                }
            }
            for (int i=0; i<NumPocStFoll; ++i)
            {
                if (dp[PicOrderCntVal()] == PocStFoll[i])
                {
                    dp.reference = SHORT_TERM;
                }
            }
            for (int i=0; i<NumPocLtCurr; ++i)
            {
                if (dp[PicOrderCntVal()] == PocLtCurr[i])
                {
                    dp.reference = LONG_TERM;
                    usedByCurrentPicture = true;
                }
            }
            for (int i=0; i<NumPocLtFoll; ++i)
            {
                if (dp[PicOrderCntVal()] == PocLtFoll[i])
                {
                    dp.reference = LONG_TERM;
                }
            }
        }
    }

    // Decoding process for generating unavailable reference pictures
    template <class H>
    void generateUnavailablePictures(H &h)
    {
        {
            //std::cout << "Generate: ";
            {
                for (int i=0; i<NumPocStFoll; ++i)
                {
                    if (!getPicByPoc(PocStFoll[i]))
                    {
                        this->generatePicture(h, PocStFoll[i], SHORT_TERM);
                        //std::cout << "G";
                    }
                    //std::cout << PocStFoll[i] << " ";
                }
                for (int i=0; i<NumPocLtFoll; ++i)
                {
                    if (!getPicByPoc(PocLtFoll[i]))
                    {
                        this->generatePicture(h, PocLtFoll[i], LONG_TERM);
                        //std::cout << "G";
                    }
                    //std::cout << PocLtFoll[i] << " ";
                }
            }
            //std::cout << "\n";
        }
    }

    template <class H>
    void constructReferencePictureLists(H &h)
    {
        StatePicture *statePicture2 = h;
        statePicture2->allBackwards = true;

        // Decoding process for reference picture lists construction
        // This process is invoked at the beginning of the decoding process for each P or B slice.
        // At the beginning of the decoding process for each slice, the reference picture list RefPicList0, and for B slices RefPicList1, are derived as follows.
        if (h[slice_type()] == P || h[slice_type()] == B)
        {
            std::array<std::shared_ptr<Pic>, 16> RefPicListTemp0;
            int NumRpsCurrTempList0 = std::max<>(h[num_ref_idx_l0_active_minus1()] + 1, h[NumPocTotalCurr()] );

            int rIdx = 0;
            while( rIdx < NumRpsCurrTempList0 )
            {
                // review - better understanding and handling of this condition:
                if (NumPocStCurrBefore + NumPocStCurrAfter + NumPocLtCurr == 0) throw Abort();

                for( int i = 0; i < NumPocStCurrBefore && rIdx < NumRpsCurrTempList0; rIdx++, i++ )
                {
                    RefPicListTemp0[ rIdx ] = RefPicSetStCurrBefore[ i ];
                }
                for( int i = 0;  i < NumPocStCurrAfter && rIdx < NumRpsCurrTempList0; rIdx++, i++ )
                {
                    RefPicListTemp0[ rIdx ] = RefPicSetStCurrAfter[ i ];
                }
                for( int i = 0; i < NumPocLtCurr && rIdx < NumRpsCurrTempList0; rIdx++, i++ )
                {
                    RefPicListTemp0[ rIdx ] = RefPicSetLtCurr[ i ];
                }
            }

            for( rIdx = 0; rIdx <= h[num_ref_idx_l0_active_minus1()]; rIdx++)
            {
                h[RefPicList(L0)][rIdx].dp = RefPicListTemp0[h[ref_pic_list_modification_flag_l0()] ? h[list_entry_l0(rIdx)] : rIdx];

                if (!h[RefPicList(L0)][rIdx].dp)
                {
                    // review: this is not a problem if the current picture will not be output or used for reference
                    h(Violation("8.3.2", "RefPicList0 contains a \"no reference picture\" entry"));
                    throw Abort();
                }

                h[RefPicList(L0)][rIdx].reference = h[RefPicList(L0)][rIdx].dp->reference;

                static_cast<StatePicture *>(h)->dpbIndexPlus1[L0][rIdx] = positionOfPictureInDpb(this->dpb, h[RefPicList(L0)][rIdx].dp.get()) + 1;

                const int poc = (*h[RefPicList(L0)][rIdx].dp)[PicOrderCntVal()];
                if (DiffPicOrderCnt(h, poc, h[PicOrderCntVal()]) > 0)
                {
                    statePicture2->allBackwards = false;
                }
            }
        }

        if (h[slice_type()] == B)
        {
            std::array<std::shared_ptr<Pic>, 16> RefPicListTemp1;
            int NumRpsCurrTempList1 = std::max<>(h[num_ref_idx_l1_active_minus1()] + 1, h[NumPocTotalCurr()] );

            int rIdx = 0;
            while( rIdx < NumRpsCurrTempList1 )
            {
                for( int i = 0;  i < NumPocStCurrAfter && rIdx < NumRpsCurrTempList1; rIdx++, i++ )
                {
                    RefPicListTemp1[ rIdx ] = RefPicSetStCurrAfter[ i ];
                }
                for( int i = 0; i < NumPocStCurrBefore && rIdx < NumRpsCurrTempList1; rIdx++, i++ )
                {
                    RefPicListTemp1[ rIdx ] = RefPicSetStCurrBefore[ i ];
                }
                for( int i = 0; i < NumPocLtCurr && rIdx < NumRpsCurrTempList1; rIdx++, i++ )
                {
                    RefPicListTemp1[ rIdx ] = RefPicSetLtCurr[ i ];
                }
            }

            for( rIdx = 0; rIdx <= h[num_ref_idx_l1_active_minus1()]; rIdx++)
            {
                h[RefPicList(L1)][rIdx].dp = RefPicListTemp1[h[ref_pic_list_modification_flag_l1()] ? h[list_entry_l1(rIdx)] : rIdx];

                if (!h[RefPicList(L1)][rIdx].dp)
                {
                    // review: this is not a problem if the current picture will not be output or used for reference
                    h(Violation("8.3.2", "RefPicList1 contains a \"no reference picture\" entry"));
                    throw Abort();
                }

                h[RefPicList(L1)][rIdx].reference = h[RefPicList(L1)][rIdx].dp->reference;

                static_cast<StatePicture *>(h)->dpbIndexPlus1[L1][rIdx] = positionOfPictureInDpb(this->dpb, h[RefPicList(L1)][rIdx].dp.get()) + 1;

                const int poc = h[RefPicList(L1)][rIdx].dp ? (*h[RefPicList(L1)][rIdx].dp)[PicOrderCntVal()] : 0;
                if (DiffPicOrderCnt(h, poc, h[PicOrderCntVal()]) > 0)
                {
                    statePicture2->allBackwards = false;
                }
            }
        }

        // Move to verifier
        if ( (h[slice_type()] == P || h[slice_type()] == B) && h[slice_temporal_mvp_enabled_flag()] )
        {
            if (h[collocated_ref_idx()] >= 0 || h[collocated_ref_idx()] < 16)
            {
                auto dp = h[RefPicList(h[collocated_from_l0_flag()] ? L0 : L1)][h[collocated_ref_idx()]].dp;
                if (dp)
                {
                    const int poc = (*dp)[PicOrderCntVal()];
                    StatePictures& state = *this;

                    if (h[first_slice_segment_in_pic_flag()])
                    {
                        state.pocCollocated = poc;
                    }
                    else
                    {
                        if (state.pocCollocated != poc)
                        {
                            h(Violation("7.4.7", "picture referred to by collocated_ref_idx (POC=%1%) is different to that in a previous slice (POC=%2%)")
                              % poc
                              % state.pocCollocated); // CondCheck 7.4.7-AG
                        }
                    }
                }
            }
        }

        if (h[slice_temporal_mvp_enabled_flag()] == 0 && h[TemporalId()] == 0)
        {
            // When both slice_temporal_mvp_enabled_flag and TemporalId are equal to 0, the syntax elements for all coded pictures that follow the current picture in decoding order shall be constrained such that no temporal motion vector from any picture that precedes the current picture in decoding order is used in decoding of any coded picture that follows the current picture in decoding order .
            for (auto dp : this->dpb)
            {
                dp->motion.reset();
            }
        }
    }

    template <class H>
    void outputAndRemoveDpbPictures(H &h)
    {
        {
            // C.5.2  Operation of the output order DPB
            // C.5.2.1  General
            // C.5.2.2  Output and removal of pictures from the DPB
            struct NotNeededForOutputAndUnusedForReference
            {
                static bool f(std::shared_ptr<StatePicture> p)
                {
                    return !p->neededForOutput && p->reference == UNUSED;
                }
            };

            struct PicLatencyCountHigh
            {
                PicLatencyCountHigh(int SpsMaxLatencyPictures) : SpsMaxLatencyPictures(SpsMaxLatencyPictures) { }
                int SpsMaxLatencyPictures;
                bool operator()(std::shared_ptr<StatePicture> dp)
                {
                    return dp->neededForOutput && dp->PicLatencyCount >= this->SpsMaxLatencyPictures;
                }
            };

            if (isIrap(h[nal_unit_type()]) && NoRaslOutputFlag == 1 && !picture0)
            {
                int NoOutputOfPriorPicsFlag = h[no_output_of_prior_pics_flag()];

                if (isCra(h[nal_unit_type()]))
                {
                    NoOutputOfPriorPicsFlag = 1;
                }

                if (NoOutputOfPriorPicsFlag)
                {
                    // If NoOutputOfPriorPicsFlag is equal to 1, all picture storage buffers in the DPB are emptied without output of the pictures they contain, and the DPB fullness is set equal to 0.
#ifdef DEBUG_DPB
                    std::cout << "clearing DPB\n";
#endif
                    while (!this->dpb.empty())
                    {
                        auto i = std::min_element(this->dpb.begin(), this->dpb.end(), comparePicOrderCntVal);
                        std::shared_ptr<Pic> dp = *i;
                        const auto poc = (*dp)[PicOrderCntVal()];
                        h(DeletePicture(poc));
                        this->dpb.erase(i);
                    }
                }
                else
                {
#ifdef DEBUG_DPB
                    std::cout << "remove NotNeededForOutputAndUnusedForReference\n";
#endif

                    // Otherwise (NoOutputOfPriorPicsFlag is equal to 0), all picture storage buffers containing a picture that is marked as "not needed for output" and "unused for reference" are emptied (without output),
                    for (auto dp : this->dpb)
                    {
                        if (NotNeededForOutputAndUnusedForReference::f(dp))
                        {
                            auto const poc = (*dp)[PicOrderCntVal()];
                            h(DeletePicture(poc));
                        }
                    }

                    this->dpb.erase(std::remove_if(this->dpb.begin(), this->dpb.end(), &NotNeededForOutputAndUnusedForReference::f), this->dpb.end());

                    //  and all non-empty picture storage buffers in the DPB are emptied by repeatedly invoking the "bumping" process specified in subclause C.5.2.4, and the DPB fullness is set equal to 0.
                    while (!this->dpb.empty()) this->bumpingProcess(h, true);
                }
            }
            else
            {
#ifdef DEBUG_DPB
                std::cout << "remove NotNeededForOutputAndUnusedForReference\n";
#endif
                // Otherwise (the current picture is not an IRAP picture with NoRaslOutputFlag equal to 1), all picture storage buffers containing a picture which are marked as "not needed for output" and "unused for reference" are emptied (without output).
                for (auto dp : this->dpb)
                {
                    if (NotNeededForOutputAndUnusedForReference::f(dp))
                    {
                        auto const poc = (*dp)[PicOrderCntVal()];
                        h(DeletePicture(poc));
                    }
                }

                this->dpb.erase(std::remove_if(this->dpb.begin(), this->dpb.end(), &NotNeededForOutputAndUnusedForReference::f), this->dpb.end());

                while (true)
                {
                    const int HighestTid= h[sps_max_sub_layers_minus1()];
                    // When one or more of the following conditions are true, the "bumping" process specified in subclause C.5.2.4 is invoked repeatedly while further decrementing the DPB fullness by one for each additional picture storage buffer that is emptied, until none of the following conditions are true
                    const int numberOfPicturesMarkedForOutput = static_cast<int>(std::count_if(dpb.begin(), dpb.end(), IsMarkedForOutput()));

                    //std::cout << "HighestTid=" << HighestTid << " numberOfPicturesMarkedForOutput=" << numberOfPicturesMarkedForOutput << "\n";

                    const bool c1 = numberOfPicturesMarkedForOutput > h[sps_max_num_reorder_pics( HighestTid ) ];
                    if (!c1)
                    {
                        const bool c2 = h[sps_max_latency_increase_plus1( HighestTid ) ] != 0 && std::count_if(dpb.begin(), dpb.end(), PicLatencyCountHigh( h[SpsMaxLatencyPictures(HighestTid)] ) );
                        //std::cout << "c2=" << c2 << "\n";
                        if (!c2)
                        {
                            //const int numberOfPicturesWithLowerOrEqualTemporalId = static_cast<int>(std::count_if(dpb.begin(), dpb.end(), TemporalIdLowerThanOrEqualTo(HighestTid)));
                            //const bool c3 = numberOfPicturesWithLowerOrEqualTemporalId > h[sps_max_dec_pic_buffering_minus1(HighestTid)] + 1;
                            //const int numberOfPicturesWithLowerOrEqualTemporalId = static_cast<int>(std::count_if(dpb.begin(), dpb.end(), TemporalIdLowerThanOrEqualTo(HighestTid)));
                            const bool c3 = dpb.size() >= size_t(h[sps_max_dec_pic_buffering_minus1(HighestTid)] + 1);
                            //std::cout << "c3=" << c3 << "\n";
                            if (!c3) break;
                        }
                    }

                    this->bumpingProcess(h);
                }
            }
        }

        if (this->dpb.empty())
        {
            h(DpbClear());
        }
    }

    template <class H>
    void sliceHeaderDone(H &h)
    {
        auto *decodedPicture = &h[Concrete<StatePicture>()];

        this->sliceHeaderValid = true;

        sliceHeaderDone2(h);

        if (h[first_slice_segment_in_pic_flag()])
        {
            if (isRasl(h[nal_unit_type()]))
            {
                //std::cout << "RASL ";
            }
            if (isIrap(h[nal_unit_type()]))
            {
                //std::cout << "IRAP ";
                this->NoRaslOutputFlag = (isIdr(h[nal_unit_type()]) ||  isBla(h[nal_unit_type()]) || this->picture0 || this->justSeenEndOfSeq) ? 1 : 0;
            }

            // Variables and functions relating to picture order count are derived. This needs to be invoked only for the first slice segment of a picture.
            // 8.3.1
            this->decodePictureOrderCount(h);

            int poc = h[PicOrderCntVal()];
#ifdef DEBUG_DPB
            std::cout << "\nPOC is " << h[PicOrderCntVal()] << " NoRaslOutputFlag=" << NoRaslOutputFlag << "\n";
#endif
            // The decoding process for RPS is invoked, wherein reference pictures may be marked as "unused for reference" or "used for long-term reference".
            this->decodeReferencePictureSet(h); // 8.3.2

            this->outputAndRemoveDpbPictures(h); // Annex C

            int n = 0;
            for (auto &dp : this->dpb)
            {
                dp->n = n++;
            }

            // When the current picture is a BLA picture or is a CRA picture with NoRaslOutputFlag equal to 1,
            if (isBla(h[nal_unit_type()]) || (isCra(h[nal_unit_type()]) && this->NoRaslOutputFlag == 1))
            {
                //the decoding process  for generating unavailable reference pictures specified in subclause 8.3.3 is invoked.
                this->generateUnavailablePictures(h); // 8.3.3
            }

            this->setupDecodedPicture(*decodedPicture, h);

            // PicOutputFlag is set as follows.
            decodedPicture->PicOutputFlag = (isRasl(h[nal_unit_type()]) && this->NoRaslOutputFlag) ? 0 : h[pic_output_flag()];

            static_cast<StatePicture *>(h)->loopFilterPicture.reset(new LoopFilter::Picture(h));
        }

        // At the beginning of the decoding process for each P or B slice, the decoding process for reference picture lists construction
        // is invoked for derivation of reference picture list 0 (RefPicList0) and, when decoding a B slice, reference picture list 1 (RefPicList1).
        if (!h[dependent_slice_segment_flag()])
        {
            this->constructReferencePictureLists(h);
        }

        this->picture0 = false;
        this->justSeenEndOfSeq = false;

        if (!h[dependent_slice_segment_flag()])
        {
            this->qPY_PREV = h[SliceQpY()];
        }
    }

    std::shared_ptr<Pic> getPicByPoc(int poc, int modulo = 0)
    {
        const int mask = modulo - 1;

        assert((poc & mask) == poc);

        for (auto &dp : this->dpb)
        {
            if (((*dp)[PicOrderCntVal()] & mask) == poc)
            {
                return dp;
            }
        }
        return 0;
    }

    template <class H>
    void accessUnitDone(H &h)
    {
#ifdef DEBUG_DPB
        std::cout << "access unit done\n";
#endif

        auto ppsEntry = h[Table<Pps>()].find(h[slice_pic_parameter_set_id()]);
        if (ppsEntry == h[Table<Pps>()].end()) return;

        auto spsEntry = h[Table<Sps>()].find(h[pps_seq_parameter_set_id()]);
        if (spsEntry == h[Table<Sps>()].end()) return;

        const int HighestTid= h[sps_max_sub_layers_minus1()];

        // The processes specified in this subclause happen instantaneously when the last decoding unit of access unit n containing the current picture is removed from the CPB.

        // For each picture in the DPB that is marked as "needed for output", the associated variable PicLatencyCount is set equal to PicLatencyCount + 1.
        for(auto &dp : this->dpb)
        {
            if (dp->neededForOutput) ++dp->PicLatencyCount;
        }

#ifdef DEBUG_DPB
        std::cout << "just decoded " << h[PicOrderCntVal()] << "\n";
#endif

        // The current picture is considered as decoded after the last decoding unit of the picture is
        // decoded. The current decoded picture is stored in an empty picture storage buffer in the
        // DPB, and the following applies:
        // - If the current decoded picture has PicOutputFlag equal to 1, it is marked as "needed for
        // output" and its associated variable PicLatencyCount is set equal to 0.
        // - Otherwise (the current decoded picture has PicOutputFlag equal to 0), it is marked as
        // "not needed for output".
        // The current decoded picture is marked as "used for short-term reference".

        StatePicture *statePicture = h;
        this->dpb.push_back(statePicture->shared_from_this());//this->currPic);

        {
            auto hPic = h.extend(this);
            hPic(JustDecoded((*statePicture)[PicOrderCntVal()]));
        }

        if (this->dpb.back()->PicOutputFlag)
        {
            this->dpb.back()->neededForOutput = true;
            this->dpb.back()->PicLatencyCount  = 0;
        }
        else
        {
            this->dpb.back()->neededForOutput = false;
        }

#ifdef DEBUG_DPB
        std::cout << "dpb[" << h[sps_max_dec_pic_buffering_minus1(0)] + 1 << "] contains: ";
        for (auto &i : this->dpb)
        {
            std::cout << " " << i->PicOrderCntVal;
            if (i->neededForOutput) std::cout << "o";
            if (i->reference == SHORT_TERM) std::cout << "s";
            if (i->reference == LONG_TERM) std::cout << "l";
        }
        std::cout << "\n";
#endif

        {
            while (true)
            {
                // When one or more of the following conditions are true, the "bumping" process specified in subclause C.5.2.4 is invoked repeatedly until none of the following conditions are true:
                const auto numberOfPicturesMarkedForOutput = std::count_if(dpb.begin(), dpb.end(), IsMarkedForOutput());
                const bool c1 = static_cast<int>(numberOfPicturesMarkedForOutput) > h[sps_max_num_reorder_pics(HighestTid)];

                if (!c1)
                {
                    const bool c2 = h[sps_max_latency_increase_plus1(HighestTid) ] != 0 && std::count_if(dpb.begin(), dpb.end(), IsMarkedForOutputAndPicLatencyCountGreaterThanOrEqualTo(h[SpsMaxLatencyPictures(HighestTid)]));
                    if (!c2)
                    {
                        break;
                    }
                }
                this->bumpingProcess(h);
            }
        }
    }

    template <class H>
    bool bumpingProcess(H &h, bool flush = false)
    {
#ifdef DEBUG_DPB
        std::cout << "bumpingProcess(" << flush << ")\n";
#endif
        typename DecodedPictureBuffer::iterator i = this->bump();
        if (i != this->dpb.end())
        {
            std::shared_ptr<StatePicture> dp = *i;
            assert(dp->neededForOutput);
            const auto poc = (*dp)[PicOrderCntVal()];
            h(OutputPicture(poc));

#ifdef DEBUG_DPB
            std::cout << "  outputting " << (*dp)[PicOrderCntVal()] << "\n";
#endif
            dp->neededForOutput = false;
            if (dp->reference == UNUSED || flush)
            {
#ifdef DEBUG_DPB
                std::cout << (flush ? "  flushing " : "  removing ") << (*dp)[PicOrderCntVal()] << "\n";
#endif
                h(DeletePicture(poc));
                this->dpb.erase(i);
            }
            return true;
        }
        return false;
    }

    template <class H, class Picture>
    void setupDecodedPicture(Picture &dp, H &h)
    {
        dp.neededForOutput = true;
        dp.reference = SHORT_TERM;
        dp.TemporalId = h[TemporalId()];
        dp.codedVideoSequenceId = this->codedVideoSequenceId;
        dp.PicLatencyCount = 0;
        dp.generatedUnavailable = false;
        dp.sequenceDecodeOrder = this->sequenceDecodeOrder++;
        dp.nal_unit_type = h[nal_unit_type()];
        dp.reconstructed = false;

        auto &picture = h[Concrete<StatePicture>()];
        setupReconstructedPicture(picture, h); // defunct

        //h(InstantiatePicture());

#ifdef VALGRIND_FRIENDLY
        for (int cIdx=0; cIdx<3; ++cIdx)
        {
            int bitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
            Raster<std::uint8_t> samples = (*dp.reconstructedPicture)[cIdx].offset(-1, -1);
            fillRectangle<std::uint8_t>(samples, 1 << (bitDepth - 1), (*dp.reconstructedPicture)[cIdx].width+2, (*dp.reconstructedPicture)[cIdx].height+2);
        }
#endif

        dp.motion.reset(new StateCollocatedMotion(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], this->dpb, h[PicOrderCntVal()]));

        {
            int bitmask = 0;
            for (const auto i : this->dpb)
            {
                bitmask |= 1 << i->notionalPositionInDpb;
            }

            dp.notionalPositionInDpb = 0;
            while ((1 << dp.notionalPositionInDpb) & bitmask) ++dp.notionalPositionInDpb;

            assert(dp.notionalPositionInDpb < 32); // actually should be <16 in HEVC
        }
    }

    DecodedPictureBuffer dpb;

    struct IsMarkedForOutput
    {
        bool operator()(const std::shared_ptr<StatePicture> dp) const
        {
            return dp->neededForOutput;
        }
    };

    struct IsMarkedForOutputAndPicLatencyCountGreaterThanOrEqualTo
    {
        IsMarkedForOutputAndPicLatencyCountGreaterThanOrEqualTo(int n) : n(n) { }
        int n;

        bool operator()(const std::shared_ptr<Pic> dp) const
        {
            return dp->neededForOutput && dp->PicLatencyCount >= this->n;
        }
    };

    struct TemporalIdLowerThanOrEqualTo
    {
        TemporalIdLowerThanOrEqualTo(int n) : n(n) { }
        int n;
        bool operator()(const std::shared_ptr<Pic> dp) const
        {
            return dp->TemporalId <= this->n;
        }
    };

    DecodedPictureBuffer::iterator bump()
    {
        // review: use of STL and filter_iterator may be elegant but makes this code opaque
        typedef boost::filter_iterator<IsMarkedForOutput, typename DecodedPictureBuffer::iterator> IteratorMarkedForOutput;
        IteratorMarkedForOutput itBegin(IsMarkedForOutput(), this->dpb.begin(), this->dpb.end());
        IteratorMarkedForOutput itEnd(IsMarkedForOutput(), this->dpb.end(), this->dpb.end());

        return std::min_element(itBegin, itEnd, comparePicOrderCntVal).base();
    }
};

template <typename Sample, class H> static void finishPicture(H h)
{
    StatePicture *statePicture2 = h;

    StateReconstructedPicture<Sample> *reconstructedPicture = h;
    auto &picture = *reconstructedPicture->picture;

    Profiler::Timers *timers = h;

    timers->processed(picture[0].width * picture[0].height * 3 / 2);

    StatePictures *statePictures = h;
    if (statePictures->pictureWanted  && statePictures->deblockWholeFrame)
    {
        if (!h[pps_deblocking_filter_disabled_flag()])
        {
            auto &recPictureL = picture[0];
            auto &recPictureCb = picture[1];
            auto &recPictureCr = picture[2];

            timers->deblock.start();
            statePicture2->loopFilterPicture->template deblock<EDGE_VER>(h, recPictureL, recPictureCb, recPictureCr, 0, 0, h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()]);
            statePicture2->loopFilterPicture->template deblock<EDGE_HOR>(h, recPictureL, recPictureCb, recPictureCr, 0, 0, h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()]);
            timers->deblock.stop();
        }

        if (h[sample_adaptive_offset_enabled_flag()])
        {
            const int pad = 96;// h[CtbSizeY()] + 16; // less padding will suffice
            Picture<Sample> *pictureSao = new Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);

            statePicture2->loopFilterPicture->applySao2<Sample>(h, *pictureSao, picture);

            reconstructedPicture->picture.reset(pictureSao);
        }

        if (statePictures->deblockWholeFrame)
        {
            auto &picture = *reconstructedPicture->picture;
            timers->pad.start();
            padPicture<Sample>(picture);
            timers->pad.stop();
        }
    }

    // review: remove clunky destruction, reconstruction
    reconstructedPicture->conformanceWindow.~ThreePlanes<Sample>();
    new (&reconstructedPicture->conformanceWindow) ThreePlanes<Sample>(
            *reconstructedPicture->picture,
            h[SubWidthC()] * h[conf_win_left_offset()],
            h[SubHeightC()] * h[conf_win_top_offset()],
            h[SubWidthC()] * h[conf_win_right_offset()],
            h[SubHeightC()] * h[conf_win_bottom_offset()]);
}

#endif
