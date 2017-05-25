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

#include "Picture.h"
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

// invoked after headers parsed but before beginning to process the current picture
struct PictureBegin { };


// invoked after the current picture has been completely decoded
struct PictureDone { };


// invoked when a picture is output from the DPB for display
struct PictureOutput { std::shared_ptr<StatePicture> statePicture; };

template <> struct Syntax<PictureOutput> : Null<PictureOutput> { };


// invoked when a picture is removed from the DPB
struct PictureDelete { std::shared_ptr<StatePicture> statePicture; };

template <> struct Syntax<PictureDelete> : Null<PictureDelete> { };


// invoked when DPB has been emptied
struct DpbClear {};

template <> struct Syntax<DpbClear> : Null<DpbClear> { };


// Request a newly allocated picture: this will be returned as a std::shared_ptr<S> where S is derived from StatePicture
struct NewPicture {};


static inline bool comparePicOrderCntVal(std::shared_ptr<StatePicture> a, std::shared_ptr<StatePicture> b)
{
    return (*a)[PicOrderCntVal()] < (*b)[PicOrderCntVal()];
}


struct StatePicturesBase
{
    StatePicturesBase() :
        sliceHeaderValid(),
        codedVideoSequenceId(-1),
        picture0(true),
        justSeenEndOfSeq(false),
        sequenceDecodeOrder(0)
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


template <class Picture, class Enable = void>
struct SetupStateReconstructedPicture
{
    template <class H> static void go(Picture &dp, H &h, bool saoslow) 
    { };
};

template <class Picture>
struct SetupStateReconstructedPicture<Picture, typename std::enable_if<std::is_base_of<StateReconstructedPictureBase, Picture>::value>::type>
{
    typedef typename Picture::Sample Sample;

    template <class H> static void go(StateReconstructedPicture<Sample> &dp, H &h, bool saoslow)
    {
        const int pad = 96;// h[CtbSizeY()] + 16; // review: less padding will suffice
        ::Picture<Sample> *picture = new ::Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);
        dp.picture.reset(picture);
        if(h[sample_adaptive_offset_enabled_flag()])
        {
            ::Picture<Sample> *saoPicture = new ::Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);
            dp.saoPicture.reset(saoPicture);
            if (saoslow)
            {
                ::Picture<Sample> *deblockPicture = new ::Picture<Sample>(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], h[chroma_format_idc()], pad, pad, 32);
                dp.deblockPicture.reset(deblockPicture);
            }
        }
    }
};

template <class Picture, class H>
void setupStateReconstructedPicture(Picture &dp, H &h, bool saoslow = false)
{
    auto &picture = h[Concrete<StatePicture>()];
    SetupStateReconstructedPicture<Picture>::go(dp, h, saoslow);
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

        for (int cIdx = 0; cIdx < 3; ++cIdx)
        {
            int bitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
            Raster<Sample> samples = (*dp.reconstructedPicture)[cIdx];
            fillRectangle<Sample>(samples, 1 << (bitDepth - 1), (*dp.reconstructedPicture)[cIdx].width, (*dp.reconstructedPicture)[cIdx].height);
        }
    }
};

struct PictureNew {};


struct PictureGenerate
{
    int poc;
    Reference reference;
};


template <> struct Syntax<PictureGenerate>
{
    template <class H> static void go(PictureGenerate pg, H &h);
};


struct StatePictures :
    StatePicturesBase
{
    std::array<std::shared_ptr<StatePicture>, 16> RefPicSetStCurrBefore, RefPicSetStCurrAfter, RefPicSetLtCurr;

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
            const int prevPicOrderCntLsb = prevTid0Pic.pic_order_cnt_lsb;
            const int prevPicOrderCntMsb = prevTid0Pic.PicOrderCntMsb;

            if ((h[slice_pic_order_cnt_lsb()] < prevPicOrderCntLsb) &&
                ((prevPicOrderCntLsb - h[slice_pic_order_cnt_lsb()]) >= (h[MaxPicOrderCntLsb()] / 2)))
            {
                PicOrderCntMsb = prevPicOrderCntMsb + h[MaxPicOrderCntLsb()];
            }
            else if ((h[slice_pic_order_cnt_lsb()] > prevPicOrderCntLsb) &&
                ((h[slice_pic_order_cnt_lsb()] - prevPicOrderCntLsb) > (h[MaxPicOrderCntLsb()] / 2)))
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
                        PocStCurrAfter[j++] = h[PicOrderCntVal()] + h[DeltaPocS1(h[CurrRpsIdx()], i)];
                    else
                        PocStFoll[k++] = h[PicOrderCntVal()] + h[DeltaPocS1(h[CurrRpsIdx()], i)];
                }
                NumPocStCurrAfter = j;
                NumPocStFoll = k;

                this->PocLsbLt.resize(h[num_long_term_sps()] + h[num_long_term_pics()]);
                this->UsedByCurrPicLt.resize(h[num_long_term_sps()] + h[num_long_term_pics()]);
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
                        DeltaPocMSBCycleLt[i] = h[delta_poc_msb_cycle_lt(i)];
                    else
                        DeltaPocMSBCycleLt[i] = h[delta_poc_msb_cycle_lt(i)] + DeltaPocMSBCycleLt[i - 1];

                    int pocLt = PocLsbLt[i];
                    if (h[delta_poc_msb_present_flag(i)])
                        pocLt += h[PicOrderCntVal()] - DeltaPocMSBCycleLt[i] * h[MaxPicOrderCntLsb()] - h[slice_pic_order_cnt_lsb()];
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

        RefPicSetStCurrBefore = std::array<std::shared_ptr<StatePicture>, 16>();
        for (int i = 0; i < NumPocStCurrBefore; ++i)
            RefPicSetStCurrBefore[i] = this->getPicByPoc(PocStCurrBefore[i]);

        RefPicSetStCurrAfter = std::array<std::shared_ptr<StatePicture>, 16>();
        for (int i = 0; i < NumPocStCurrAfter; ++i)
            RefPicSetStCurrAfter[i] = this->getPicByPoc(PocStCurrAfter[i]);

        RefPicSetLtCurr = std::array<std::shared_ptr<StatePicture>, 16>();
        for (int i = 0; i < NumPocLtCurr; ++i)
            RefPicSetLtCurr[i] = this->getPicByPoc(PocLtCurr[i]);

        // All reference pictures in the decoded picture buffer that are not included in RefPicSetLtCurr, RefPicSetLtFoll, RefPicSetStCurrBefore, RefPicSetStCurrAfter or RefPicSetStFoll are marked as "unused for reference".
        for (auto &p : this->dpb)
        {
            auto &dp = *p;

            dp.reference = UNUSED;

            for (int i = 0; i < NumPocStCurrBefore; ++i)
                if (dp[PicOrderCntVal()] == PocStCurrBefore[i])
                    dp.reference = SHORT_TERM;

            for (int i = 0; i < NumPocStCurrAfter; ++i)
                if (dp[PicOrderCntVal()] == PocStCurrAfter[i])
                    dp.reference = SHORT_TERM;

            for (int i = 0; i < NumPocStFoll; ++i)
                if (dp[PicOrderCntVal()] == PocStFoll[i])
                    dp.reference = SHORT_TERM;

            for (int i = 0; i < NumPocLtCurr; ++i)
                if (dp[PicOrderCntVal()] == PocLtCurr[i])
                    dp.reference = LONG_TERM;

            for (int i = 0; i < NumPocLtFoll; ++i)
                if (dp[PicOrderCntVal()] == PocLtFoll[i])
                    dp.reference = LONG_TERM;
        }
    }

    // Decoding process for generating unavailable reference pictures
    template <class H>
    void generateUnavailablePictures(H &h)
    {
        for (int i = 0; i < NumPocStFoll; ++i)
            if (!getPicByPoc(PocStFoll[i]))
                h(PictureGenerate{ PocStFoll[i], SHORT_TERM });

        for (int i = 0; i < NumPocLtFoll; ++i)
            if (!getPicByPoc(PocLtFoll[i]))
                h(PictureGenerate{ PocLtFoll[i], LONG_TERM });
    }

    template <class H>
    void constructReferencePictureLists(H &h)
    {
        StatePicture *statePicture = h;
        statePicture->allBackwards = true;

        // Decoding process for reference picture lists construction
        // This process is invoked at the beginning of the decoding process for each P or B slice.
        // At the beginning of the decoding process for each slice, the reference picture list RefPicList0, and for B slices RefPicList1, are derived as follows.
        if (h[slice_type()] == P || h[slice_type()] == B)
        {
            std::array<std::shared_ptr<StatePicture>, 16> RefPicListTemp0;
            int NumRpsCurrTempList0 = std::max<>(h[num_ref_idx_l0_active_minus1()] + 1, h[NumPocTotalCurr()]);

            int rIdx = 0;
            while (rIdx < NumRpsCurrTempList0)
            {
                // review - better understanding and handling of this condition:
                if (NumPocStCurrBefore + NumPocStCurrAfter + NumPocLtCurr == 0) 
                    throw Abort();

                for (int i = 0; i < NumPocStCurrBefore && rIdx < NumRpsCurrTempList0; rIdx++, i++)
                    RefPicListTemp0[rIdx] = RefPicSetStCurrBefore[i];

                for (int i = 0; i < NumPocStCurrAfter && rIdx < NumRpsCurrTempList0; rIdx++, i++)
                    RefPicListTemp0[rIdx] = RefPicSetStCurrAfter[i];

                for (int i = 0; i < NumPocLtCurr && rIdx < NumRpsCurrTempList0; rIdx++, i++)
                    RefPicListTemp0[rIdx] = RefPicSetLtCurr[i];
            }

            for (rIdx = 0; rIdx <= h[num_ref_idx_l0_active_minus1()]; rIdx++)
            {
                h[RefPicList(L0)][rIdx].dp = RefPicListTemp0[h[ref_pic_list_modification_flag_l0()] ? h[list_entry_l0(rIdx)] : rIdx];

                if (!h[RefPicList(L0)][rIdx].dp)
                {
                    // review: this is not a problem if the current picture will not be output or used for reference
                    h(Violation("8.3.2", "RefPicList0 contains a \"no reference picture\" entry"));
                    throw Abort();
                }

                h[RefPicList(L0)][rIdx].reference = h[RefPicList(L0)][rIdx].dp->reference;

                statePicture->dpbIndexPlus1[L0][rIdx] = positionOfPictureInDpb(this->dpb, h[RefPicList(L0)][rIdx].dp.get()) + 1;

                const int poc = (*h[RefPicList(L0)][rIdx].dp)[PicOrderCntVal()];

                if (DiffPicOrderCnt(h, poc, h[PicOrderCntVal()]) > 0)
                    statePicture->allBackwards = false;
            }
        }

        if (h[slice_type()] == B)
        {
            std::array<std::shared_ptr<StatePicture>, 16> RefPicListTemp1;
            int NumRpsCurrTempList1 = std::max<>(h[num_ref_idx_l1_active_minus1()] + 1, h[NumPocTotalCurr()]);

            int rIdx = 0;
            while (rIdx < NumRpsCurrTempList1)
            {
                for (int i = 0; i < NumPocStCurrAfter && rIdx < NumRpsCurrTempList1; rIdx++, i++)
                    RefPicListTemp1[rIdx] = RefPicSetStCurrAfter[i];

                for (int i = 0; i < NumPocStCurrBefore && rIdx < NumRpsCurrTempList1; rIdx++, i++)
                    RefPicListTemp1[rIdx] = RefPicSetStCurrBefore[i];

                for (int i = 0; i < NumPocLtCurr && rIdx < NumRpsCurrTempList1; rIdx++, i++)
                    RefPicListTemp1[rIdx] = RefPicSetLtCurr[i];
            }

            for (rIdx = 0; rIdx <= h[num_ref_idx_l1_active_minus1()]; rIdx++)
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
                    statePicture->allBackwards = false;
            }
        }

        // Move to verifier
        if ((h[slice_type()] == P || h[slice_type()] == B) && h[slice_temporal_mvp_enabled_flag()])
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
                    NoOutputOfPriorPicsFlag = 1;

                if (NoOutputOfPriorPicsFlag)
                {
                    // If NoOutputOfPriorPicsFlag is equal to 1, all picture storage buffers in the DPB are emptied without output of the pictures they contain, and the DPB fullness is set equal to 0.
                    // all picture storage buffers in the DPB are emptied without output of the pictures they contain, and the DPB fullness is set equal to 0.
                    while (!this->dpb.empty())
                    {
                        auto i = this->dpb.begin();
                        h(PictureDelete{ *i });
                        this->dpb.erase(i);
                    }
                }
                else
                {
                    // all picture storage buffers containing a picture that is marked as "not needed for output" and "unused for reference" are emptied (without output),
                    for (auto i = this->dpb.begin(); i != this->dpb.end(); )
                    {
                        auto &dp = **i;
                        if (!dp.neededForOutput && dp.reference == UNUSED)
                        {
                            h(PictureDelete{ *i });
                            i = this->dpb.erase(i);
                        }
                        else
                            ++i;
                    }

                    //  and all non-empty picture storage buffers in the DPB are emptied by repeatedly invoking the "bumping" process specified in subclause C.5.2.4, and the DPB fullness is set equal to 0.
                    while (!this->dpb.empty()) 
                        this->bumpingProcess(h, true);
                }
            }
            else
            {
                // Otherwise (the current picture is not an IRAP picture with NoRaslOutputFlag equal to 1), all picture storage buffers containing a picture which are marked as "not needed for output" and "unused for reference" are emptied (without output).
                for (auto i = this->dpb.begin(); i != this->dpb.end(); )
                {
                    auto &dp = **i;
                    if (!dp.neededForOutput && dp.reference == UNUSED)
                    {
                        h(PictureDelete{ *i });
                        i = this->dpb.erase(i);
                    }
                    else
                        ++i;
                }

                while (true)
                {
                    const int HighestTid = h[sps_max_sub_layers_minus1()];
                    // When one or more of the following conditions are true, the "bumping" process specified in subclause C.5.2.4 is invoked repeatedly while further decrementing the DPB fullness by one for each additional picture storage buffer that is emptied, until none of the following conditions are true
                    const int numberOfPicturesMarkedForOutput = static_cast<int>(std::count_if(dpb.begin(), dpb.end(), IsMarkedForOutput()));

                    //std::cout << "HighestTid=" << HighestTid << " numberOfPicturesMarkedForOutput=" << numberOfPicturesMarkedForOutput << "\n";

                    const bool c1 = numberOfPicturesMarkedForOutput > h[sps_max_num_reorder_pics(HighestTid)];
                    if (!c1)
                    {
                        const bool c2 = h[sps_max_latency_increase_plus1(HighestTid)] != 0 && std::count_if(dpb.begin(), dpb.end(), PicLatencyCountHigh(h[SpsMaxLatencyPictures(HighestTid)]));
                        if (!c2)
                        {
                            const bool c3 = dpb.size() >= size_t(h[sps_max_dec_pic_buffering_minus1(HighestTid)] + 1);

                            if (!c3)
                                break;
                        }
                    }

                    this->bumpingProcess(h);
                }
            }
        }

        if (this->dpb.empty())
            h(DpbClear());
    }

    template <class H>
    void sliceHeaderDone(H &h)
    {
        auto *decodedPicture = &h[Concrete<StatePicture>()];
        this->sliceHeaderValid = true;

        if (h[first_slice_segment_in_pic_flag()])
            h(PictureBegin());

        if (!h[dependent_slice_segment_flag()])
        {
            h[SliceAddrRs()] = h[slice_segment_address()];

            // At the beginning of the decoding process for each P or B slice, the decoding process for reference picture lists construction
            // is invoked for derivation of reference picture list 0 (RefPicList0) and, when decoding a B slice, reference picture list 1 (RefPicList1).
            this->constructReferencePictureLists(h);
        }

        this->picture0 = false;
        this->justSeenEndOfSeq = false;
        this->sliceHeaderValid = true;
    }

    std::shared_ptr<StatePicture> getPicByPoc(int poc, int modulo = 0)
    {
        const int mask = modulo - 1;

        assert((poc & mask) == poc);

        for (auto &dp : this->dpb)
            if (((*dp)[PicOrderCntVal()] & mask) == poc)
                return dp;

        return 0;
    }

    template <class H>
    bool bumpingProcess(H &h, bool flush = false)
    {
        typename DecodedPictureBuffer::iterator i = this->bump();
        if (i != this->dpb.end())
        {
            std::shared_ptr<StatePicture> statePicture = *i;

            assert(statePicture->neededForOutput);
            h(PictureOutput{ statePicture });

            statePicture->neededForOutput = false;
            statePicture->hasBeenOutput = true;

            if (statePicture->reference == UNUSED || flush)
            {
                h(PictureDelete{ statePicture });

                this->dpb.erase(i);
            }

            return true;
        }
        return false;
    }

    template <class H>
    void setupDecodedPicture(StatePicture *statePicture, H &h)
    {
        StatePictures *statePictures = h;
        statePicture->neededForOutput = true;
        statePicture->reference = SHORT_TERM;
        statePicture->TemporalId = h[TemporalId()];
        statePicture->codedVideoSequenceId = this->codedVideoSequenceId;
        statePicture->PicLatencyCount = 0;
        statePicture->generatedUnavailable = false;
        statePicture->sequenceDecodeOrder = this->sequenceDecodeOrder++;
        statePicture->nal_unit_type = h[nal_unit_type()];
        statePicture->reconstructed = false;

        statePicture->motion.reset(new StateCollocatedMotion(h[pic_width_in_luma_samples()], h[pic_height_in_luma_samples()], statePictures->dpb, h[PicOrderCntVal()]));

        {
            int bitmask = 0;
            for (const auto i : statePictures->dpb)
                bitmask |= 1 << i->notionalPositionInDpb;

            statePicture->notionalPositionInDpb = 0;
            while ((1 << statePicture->notionalPositionInDpb) & bitmask)
                ++statePicture->notionalPositionInDpb;

            assert(statePicture->notionalPositionInDpb < 32); // actually should be <16 in HEVC
        }

        statePictures->dpb.push_back(statePicture->shared_from_this());
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

        bool operator()(const std::shared_ptr<StatePicture> p) const
        {
            return p->neededForOutput && p->PicLatencyCount >= this->n;
        }
    };

    struct TemporalIdLowerThanOrEqualTo
    {
        TemporalIdLowerThanOrEqualTo(int n) : n(n) { }
        int n;
        bool operator()(const std::shared_ptr<StatePicture> p) const
        {
            return p->TemporalId <= this->n;
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

template <> struct Syntax<PictureBegin>
{
    template <class H> static void go(PictureBegin, H &h)
    {
        StatePictures *statePictures = h;

        std::vector<int> colWidth(h[num_tile_columns_minus1()] + 1);
        if (h[uniform_spacing_flag()])
        {
            for (int i = 0; i <= h[num_tile_columns_minus1()]; i++)
                colWidth[i] = ((i + 1) * h[PicWidthInCtbsY()]) / (h[num_tile_columns_minus1()] + 1) -
                    (i * h[PicWidthInCtbsY()]) / (h[num_tile_columns_minus1()] + 1);
        }
        else
        {
            colWidth[h[num_tile_columns_minus1()]] = h[PicWidthInCtbsY()];
            for (int i = 0; i < h[num_tile_columns_minus1()]; i++) 
            {
                colWidth[i] = h[column_width_minus1(i)] + 1;
                colWidth[h[num_tile_columns_minus1()]] -= colWidth[i];
            }
        }

        std::vector<int> rowHeight(h[num_tile_rows_minus1()] + 1);
        if (h[uniform_spacing_flag()])
        {
            for (int i = 0; i <= h[num_tile_rows_minus1()]; i++)
                rowHeight[i] = ((i + 1) * h[PicHeightInCtbsY()]) / (h[num_tile_rows_minus1()] + 1) -
                    (i * h[PicHeightInCtbsY()]) / (h[num_tile_rows_minus1()] + 1);
        }
        else
        {
            rowHeight[h[num_tile_rows_minus1()]] = h[PicHeightInCtbsY()];
            for (int i = 0; i < h[num_tile_rows_minus1()]; i++) 
            {
                rowHeight[i] = h[row_height_minus1(i)] + 1;
                rowHeight[h[num_tile_rows_minus1()]] -= rowHeight[i];
            }
        }

        h[ContainerOf<colBd>()].resize(h[num_tile_columns_minus1()] + 2);
        for (int i = 0; i <= h[num_tile_columns_minus1()]; i++)
        {
            h[colBd(i + 1)] = h[colBd(i)] + colWidth[i];
        }

        h[ContainerOf<rowBd>()].resize(h[num_tile_rows_minus1()] + 2);
        for (int i = 0; i <= h[num_tile_rows_minus1()]; i++)
        {
            h[rowBd(i + 1)] = h[rowBd(i)] + rowHeight[i];
        }

        h[ContainerOf<CtbAddrRsToTs>()].resize(h[PicSizeInCtbsY()] + 2);
        for (int ctbAddrRS = 0; ctbAddrRS < h[PicSizeInCtbsY()]; ctbAddrRS++)
        {
            const int tbX = ctbAddrRS %  h[PicWidthInCtbsY()];
            const int tbY = ctbAddrRS / h[PicWidthInCtbsY()];
            int tileX;
            for (int i = 0; i <= h[num_tile_columns_minus1()]; i++)
                if (tbX >= h[colBd(i)])
                    tileX = i;
            int tileY;
            for (int j = 0; j <= h[num_tile_rows_minus1()]; j++)
                if (tbY >= h[rowBd(j)])
                    tileY = j;
            h[CtbAddrRsToTs(ctbAddrRS)] = 0;
            for (int i = 0; i < tileX; i++)
                h[CtbAddrRsToTs(ctbAddrRS)] += rowHeight[tileY] * colWidth[i];
            for (int j = 0; j < tileY; j++)
                h[CtbAddrRsToTs(ctbAddrRS)] += h[PicWidthInCtbsY()] * rowHeight[j];
            h[CtbAddrRsToTs(ctbAddrRS)] += (tbY - h[rowBd(tileY)]) * colWidth[tileX] + tbX - h[colBd(tileX)];
        }

        h[ContainerOf<CtbAddrTsToRs>()].resize(h[PicSizeInCtbsY()] + 1);
        for (int ctbAddrRS = 0; ctbAddrRS < h[PicSizeInCtbsY()]; ctbAddrRS++)
            h[CtbAddrTsToRs(h[CtbAddrRsToTs(ctbAddrRS)])] = ctbAddrRS;
        h[CtbAddrTsToRs(h[PicSizeInCtbsY()])] = h[PicSizeInCtbsY()];

        h[ContainerOf<TileId>()].resize(h[PicSizeInCtbsY()] + 2);
        h[TileId(0)] = -2;
        for (int j = 0, tIdx = 0; j <= h[num_tile_rows_minus1()]; j++)
            for (int i = 0; i <= h[num_tile_columns_minus1()]; i++, tIdx++)
                for (int y = h[rowBd(j)]; y < h[rowBd(j + 1)]; y++)
                    for (int x = h[colBd(i)]; x < h[colBd(i + 1)]; x++)
                        h[TileId(h[CtbAddrRsToTs(y* h[PicWidthInCtbsY()] + x)])] = tIdx;
        h[TileId(h[PicSizeInCtbsY()])] = -1;

        auto const nut = h[nal_unit_type()];

        if (isIrap(nut))
            statePictures->NoRaslOutputFlag = (isIdr(nut) || isBla(nut) || statePictures->picture0 || statePictures->justSeenEndOfSeq) ? 1 : 0;

        // Variables and functions relating to picture order count are derived. This needs to be invoked only for the first slice segment of a picture. 
        statePictures->decodePictureOrderCount(h); // 8.3.1

        // The decoding process for RPS is invoked, wherein reference pictures may be marked as "unused for reference" or "used for long-term reference".
        statePictures->decodeReferencePictureSet(h); // 8.3.2

        statePictures->outputAndRemoveDpbPictures(h); // Annex C

        // When the current picture is a BLA picture or is a CRA picture with NoRaslOutputFlag equal to 1,
        if (isBla(nut) || (isCra(nut) && statePictures->NoRaslOutputFlag == 1))
        {
            //the decoding process for generating unavailable reference pictures specified in subclause 8.3.3 is invoked.
            statePictures->generateUnavailablePictures(h); // 8.3.3
        }

        // The processes specified in this subclause happen instantaneously when the last decoding unit of access unit n containing the current picture is removed from the CPB.
        // For each picture in the DPB that is marked as "needed for output", the associated variable PicLatencyCount is set equal to PicLatencyCount + 1.
        int n = 0;
        for (auto &dp : statePictures->dpb)
        {
            dp->n = n++;
            if (dp->neededForOutput)
                ++dp->PicLatencyCount;
        }

        StatePicture *statePicture = h;
        statePictures->setupDecodedPicture(statePicture, h);

        statePicture->PicOutputFlag = (isRasl(nut) && statePictures->NoRaslOutputFlag) ? 0 : h[pic_output_flag()];
        statePicture->neededForOutput = !!statePictures->dpb.back()->PicOutputFlag;
        statePicture->PicLatencyCount = 0;

        statePicture->loopFilterPicture.reset(new LoopFilter::Picture(h));
    }
};


template <class Verb, class Tuple>
struct HandleValue<NewPicture, Verb, Tuple>
{
    typedef std::shared_ptr<StatePicture> Type;
    template <class H> static Type get(NewPicture, H &)
    {
        return Type{ new StatePicture };
    }
};


template <class H> void Syntax<PictureGenerate>::go(PictureGenerate pg, H &h)
{
    StatePictures *statePictures = h;

    std::shared_ptr<StatePicture> statePicture = h[NewPicture()];

    (*statePicture)[PicOrderCntVal()] = pg.poc;

    statePictures->setupDecodedPicture(statePicture.get(), h);

    statePicture->reference = pg.reference;
    statePicture->neededForOutput = false;
    statePicture->reconstructed = true;
    statePicture->generatedUnavailable = true;
}


template <> struct Syntax<PictureDone>
{
    template <class H> static void go(PictureDone, H &h)
    {
        StatePictures *statePictures = h;
        auto &dpb = statePictures->dpb;

        const int HighestTid = h[sps_max_sub_layers_minus1()];

        while (true)
        {
            // When one or more of the following conditions are true, the "bumping" process specified in subclause C.5.2.4 is invoked repeatedly until none of the following conditions are true:
            const auto numberOfPicturesMarkedForOutput = std::count_if(dpb.begin(), dpb.end(), StatePictures::IsMarkedForOutput());
            const bool c1 = static_cast<int>(numberOfPicturesMarkedForOutput) > h[sps_max_num_reorder_pics(HighestTid)];

            if (!c1)
            {
                const bool c2 = h[sps_max_latency_increase_plus1(HighestTid)] != 0 && std::count_if(dpb.begin(), dpb.end(), StatePictures::IsMarkedForOutputAndPicLatencyCountGreaterThanOrEqualTo(h[SpsMaxLatencyPictures(HighestTid)]));
                if (!c2)
                    break;
            }
            statePictures->bumpingProcess(h);
        }
    }
};


static const int tablesDs = 0;
static const int tablesWpp = 1;


template <>
struct Syntax<ContextsInitialize>
{
    template <class H> static void go(ContextsInitialize, H &h)
    {
        Contexts& contexts = *static_cast<Contexts *>(h);
        contexts.initialize(h[SliceQpY()], h[initType()]);
        h[StatCoeff(0)] = 0;
        h[StatCoeff(1)] = 0;
        h[StatCoeff(2)] = 0;
        h[StatCoeff(3)] = 0;
    }
};

template <>
struct Syntax<ContextsSave>
{
    template <class H> static void go(ContextsSave e, H &h)
    {
        Contexts *contexts = h;
        StateSlice *stateSlice = h;

        stateSlice->savedContexts[e.i] = *contexts;
    }
};

template <>
struct Syntax<ContextsRestore>
{
    template <class H> static void go(ContextsRestore e, H &h)
    {
        Contexts *contexts = h;
        StateSlice *stateSlice = h;

        *contexts = stateSlice->savedContexts[e.i];
    }
};


template <class H>
void preCtu(H &h)
{
    Neighbourhood *neighbourhood = h;
    StateSpatial *stateSpatial = h;
    AvailabilityCtu *availabilityCtu = h;

    availabilityCtu->init(h);

    auto const rx = (h[CtbAddrInRs()] % h[PicWidthInCtbsY()]);
    auto const ry = (h[CtbAddrInRs()] / h[PicWidthInCtbsY()]);

    bool const codingTreeUnitIsFirstInSlice = h[SliceAddrRs()] == h[CtbAddrInRs()];
    bool const codingTreeUnitIsFirstInTile = h[CtbAddrInTs()] == 0 || h[TileId(h[CtbAddrInTs()])] != h[TileId(h[CtbAddrInTs()] - 1)];
    bool const codingTreeUnitIsFirstInRow = rx == 0;

    bool const resetQp =
        codingTreeUnitIsFirstInSlice ||
        codingTreeUnitIsFirstInTile ||
        (codingTreeUnitIsFirstInRow && h[entropy_coding_sync_enabled_flag()] == 1);

    if (resetQp)
        h[QpY()] = h[SliceQpY()];

    const bool initialiseCabac =
        h[CtbAddrInRs()] == h[slice_segment_address()] ||
        codingTreeUnitIsFirstInTile ||
        (codingTreeUnitIsFirstInRow && h[entropy_coding_sync_enabled_flag()] == 1);

    if (initialiseCabac)
    {
        h(CabacRestart());

        if (codingTreeUnitIsFirstInTile)
        {
            h(ContextsInitialize());
        }
        else if (h[entropy_coding_sync_enabled_flag()] == 1 && codingTreeUnitIsFirstInRow)
        {
            // start of WPP row
            auto const y0 = ry << h[CtbLog2SizeY()];
            AvailabilityCtu *availabilityCtu = h;
            auto const availableFlagT = availabilityCtu->available(0, y0, h[CtbSizeY()], y0 - h[CtbSizeY()], h[CtbLog2SizeY()]);
            if (availableFlagT)
                h(ContextsRestore(tablesWpp));
            else
                h(ContextsInitialize());
        }
        else if (h[CtbAddrInRs()] == h[slice_segment_address()] && h[dependent_slice_segment_flag()])
        {
            // start of depdendent slice
            h(ContextsRestore(tablesDs));
        }
        else
        {
            // start of slice
            h(ContextsInitialize());
        }
    }

    if (codingTreeUnitIsFirstInTile || codingTreeUnitIsFirstInSlice)
    {
        Turing::Rectangle rectangle;
        rectangle.x0 = 0;
        rectangle.y0 = 0;
        rectangle.width = h[PicWidthInCtbsY()];
        rectangle.height = h[PicHeightInCtbsY()];

        stateSpatial->snakeSaoCtuData.foreach(rectangle, 0, 0, 0, [](SaoCtuData &data) { new (&data) SaoCtuData; });

        rectangle.width <<= h[CtbLog2SizeY()];
        rectangle.height <<= h[CtbLog2SizeY()];

        neighbourhood->snakeMerge.foreach(rectangle, h[MinCbLog2SizeY()] - 1, 1, 1, [](BlockData &data) { data.reset(); });
        neighbourhood->snake.foreach(rectangle, h[MinCbLog2SizeY()] - 1, 1, 1, [](BlockData &data) { data.reset(); });
    }

    SaoCtuData saoCtuData = SaoCtuData();
    stateSpatial->snakeSaoCtuData.commitRectangle(Turing::Rectangle{ rx, ry, 1, 1 }, saoCtuData, 0);
}


template <class H>
void postCtu(H &h)
{
    int const rx = (h[CtbAddrInRs()] % h[PicWidthInCtbsY()]);

    if (h[entropy_coding_sync_enabled_flag()] == 1 && rx == 1)
        h(ContextsSave(tablesWpp));
}

#endif
