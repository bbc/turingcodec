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

#ifndef INCLUDED_StatePicture_h
#define INCLUDED_StatePicture_h

#pragma once

#include "StateSpatial.h"
#include "LoopFilter.h"
#include "StateCollocatedMotion.h"
#include "StateParameterSets.h"
#include "StateValues.h"
#include "ScalingMatrices.h"
#include "Cabac.h"
#include <memory>
#include <array>


struct StateReconstructedPictureBase
{
    virtual ~StateReconstructedPictureBase() { };
};


template <typename T>
struct StateReconstructedPicture :
    StateReconstructedPictureBase
{
    typedef T Sample;
    std::shared_ptr<Picture<T>> picture;
    std::shared_ptr<Picture<T>> saoPicture;
    std::shared_ptr<Picture<T>> deblockPicture;
    ThreePlanes<T> conformanceWindow;
};

// metafunction to retrieve sample type from handler type
template <class H>
struct SampleType
{
    typedef typename Access<Concrete<StateReconstructedPictureBase>, H>::ActualType::Sample Type;
    static_assert(std::is_same<Type, uint8_t>::value || std::is_same<Type, uint16_t>::value, "");
};


// State that relates to a NAL unit
struct StateNalUnit :
    ValueHolder<nal_unit_type>,
    ValueHolder<nuh_layer_id>,
    ValueHolder<nuh_temporal_id_plus1>
{
    StateNalUnit()
    {
        this->ValueHolder<nuh_temporal_id_plus1>::value = 1;
    }
};


// Decoder or encoder state that persists while working on a particular slice of video
struct StateSlice :
    StateNalUnit,
    ValueHolder<SliceAddrRs>,
    ValueHolder<CtbAddrInTs>,
    ValueHolder<QpY>,
    SliceSegmentHeader,
    ScalingMatrices,
    ActiveParameterSets,
    short_term_ref_pic_set
{
    Contexts savedContexts[2];
    int savedStatCoeff[2][4];
};


// Decoder or encoder state that persists while working on a particular video frame
struct StatePicture :
    std::enable_shared_from_this<StatePicture>,
    AccessOperators<StatePicture>,
    ValueHolder<PicOrderCntVal>,
    ValueHolder<colBd>,
    ValueHolder<rowBd>,
    ValueHolder<CtbAddrRsToTs>,
    ValueHolder<CtbAddrTsToRs>,
    ValueHolder<TileId>,
    StateSlice,
    StateSpatial
    {
        using AccessOperators<StatePicture>::operator[];

        virtual ~StatePicture() { }

        std::shared_ptr<StateReconstructedPictureBase> reconstructedPicture;

        std::shared_ptr<LoopFilter::Picture> loopFilterPicture;

        Reference reference;

        bool neededForOutput;
        bool hasBeenOutput = false;
        int nal_unit_type;
        int TemporalId;
        int codedVideoSequenceId;
        int PicLatencyCount;
        int PicOutputFlag;
        bool generatedUnavailable;
        bool reconstructed;
        std::uint64_t sequenceDecodeOrder;
        int notionalPositionInDpb;
        int n;

        std::shared_ptr<struct StateCollocatedMotion> motion;

        struct RplEntry
        {
            std::shared_ptr<struct StatePicture> dp;

            operator PicOrderCnt() { return (*this->dp)[PicOrderCntVal()]; }
            // This may not be the same as dp->reference - it is constant for the life of the RPL whereas dp->reference may change (in another thread).
            Reference reference;
        };

        typedef std::array<RplEntry, 16> Rpl;

        Rpl refPicList[2];

        bool allBackwards; // true if DiffPicOrderCnt( aPic, currPic ) is less than or equal to 0 for every picture aPic in every reference picture list of the current slice

        std::uint8_t dpbIndexPlus1[2][16];

        void clearRpls()
        {
            this->refPicList[0] = Rpl();
            this->refPicList[1] = Rpl();
        }
    };


typedef StatePicture CurrPic;

static bool LongTermRefPic(StatePicture::RplEntry picX)
{
    return picX.reference == LONG_TERM;
}

typedef std::vector<std::shared_ptr<StatePicture>> DecodedPictureBuffer;

static int positionOfPictureInDpb(DecodedPictureBuffer const &dpb, StatePicture const *dp)
{
    assert(dp == dpb[dp->n].get());
    return dp->n;
}

static StatePicture const *getDpbPictureByIndex(DecodedPictureBuffer const &dpb, int dpbIndex)
{
    return dpb[dpbIndex].get();
}

template <class S>
struct Access<RefPicList, S, typename std::enable_if<std::is_base_of<StatePicture, S>::value>::type>
{
    typedef StatePicture::Rpl &Type;
    static Type get(RefPicList rpl, StatePicture &s)
    {
        return s.refPicList[rpl.x];
    }
};


template <class S>
struct Access<ReconstructedSamples, S, typename std::enable_if<std::is_base_of<StateReconstructedPictureBase, S>::value>::type>
{
    typedef typename S::Sample Sample;
    typedef Raster<Sample> Type;
    static Type get(ReconstructedSamples e, StateReconstructedPicture<Sample> &s)
    {
        return (*s.picture)(e.x, e.y, e.cIdx);
    }
};

template <class H>
void commitCu(H &h, const coding_quadtree &cqt)
{
    Neighbourhood *neighbourhood = h;
    Snake<BlockData>::Cursor *cursor = h;

    BlockData &blockData = cursor->current(0, 0, neighbourhood->MinCbLog2SizeYMinus1);
    const bool isIntra = !blockData.predFlag(L0) && !blockData.predFlag(L1);

    neighbourhood->recordMerge(h, cqt, isIntra);

    if (isIntra)
    {
        StatePicture *statePicture = h;
        statePicture->motion->fillRectangleIntra(cqt);
    }
}

template <class H>
void commitPu(H &h, prediction_unit pu, const PuData &puData)
{
    StatePicture *decodedPicture = h;
    decodedPicture->motion->fillRectangle(pu, puData);

    Neighbourhood *neighbourhood = h;
    neighbourhood->recordMerge(h, pu);
}

#endif
