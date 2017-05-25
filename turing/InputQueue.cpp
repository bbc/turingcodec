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

#include "InputQueue.h"
#include "HevcTypes.h"
#include <cassert>
#include <deque>

struct LookaheadAnalysisResults
{
};


namespace {

struct Piece
{
    Piece(std::shared_ptr<PictureWrapper> picture)
        :
        picture(picture)
    {
    }

    std::shared_ptr<InputQueue::Docket> docket;
    std::shared_ptr<PictureWrapper> picture;
    std::shared_ptr<AdaptiveQuantisation> aqInfo;
    std::shared_ptr<LookaheadAnalysisResults> lookaheadAnalysisResults;

    bool done() const
    {
        return !this->picture;
    }
};

void addReference(InputQueue::Docket &docket, int ref)
{
    if (ref < 0) docket.references.negative.insert(ref);
    if (ref > 0) docket.references.positive.insert(ref);
}

}


struct InputQueue::State
{
    State(int maxGopN, int maxGopM, bool fieldCoding, bool shotChange, int segmentLength, int baseQP)
        :
        maxGopN(maxGopN),
        maxGopM(maxGopM),
        fieldCoding(fieldCoding),
        shotChange(shotChange),
        segmentLength(segmentLength),
        baseQP(baseQP),
        scd(new ShotChangeDetection())
    {
    }

    const int maxGopN;
    const int maxGopM;
    const bool fieldCoding;
    const bool shotChange;
    const int segmentLength;
    const int baseQP;

    // list of pictures in preanalysis stage
    std::deque<Piece> entriesPreanalysis;

    // Shot change detection engine
    std::unique_ptr<ShotChangeDetection> scd;

    //list of pictures use during SOP planning - ~8 pictures
    std::deque<Piece> entries;

    struct Timestamp
    {
        int64_t timestamp;
        size_t sequence;
    };

    std::deque<Timestamp> timestamps;

    size_t pictureInputCount = 0;
    size_t picturePreprocessCount = 0;
    int segmentFront = 0;
    int sequenceFront = 0;
    int sequenceIFrame = 0;
    int sequenceIdr = 0;
    int sequenceSeg = 0;
    int hierarchyLevelFront = 0;
    bool finish = false;
    int gopSize;
    int currentIntraPoc = 0;
    int previousIntraPoc = 0;
    int currentSegmentPoc = 0;
    int previousSegmentPoc = 0;
    void process();

    void createDocket(int size, int i, int nut, int qpOffset, double qpFactor, int sopLevel, int hierarchyLevel, int numSameHierarchyLevel, int ref1 = 0, int ref2 = 0, int ref3 = 0, int ref4 = 0);

    bool isValidReference(int i, int delta) const
    {
        assert(i > 0);

        if (delta < 0) 
            return true;
        if (delta == 0) 
            return false;
        auto const  offset =(sequenceIFrame < segmentFront)? this->entries.size() : (sequenceIFrame - segmentFront + 1);
        auto const  deltaLimit = (static_cast<int>(this->entries.size() > offset)) ? offset : static_cast<int>(this->entries.size());
        if (i - 1 + delta >= static_cast<int>(deltaLimit)) 
            return false;
        return !!this->entries.at(i - 1 + delta).docket;
    }

    int computeNextIdr(int segmentFront, int sequenceFront, int currDefaultIdr, int intraPeriod, bool& isShotChange)
    {
        int nextIdr = currDefaultIdr + intraPeriod;
        for (int i = segmentFront; i < (currDefaultIdr + intraPeriod); i++)
        {
            int index = i + sequenceFront - segmentFront;
            index = (index < 0) ? 0 : index;
            if (index >= this->scd->m_shotChangeList.size())
            {
                nextIdr = currDefaultIdr;
                break;
            }
                
            if (this->scd->m_shotChangeList[index])
            {
                nextIdr = i;
                isShotChange = true;
                break;
            }
        }
        return nextIdr;
    }
};


InputQueue::InputQueue(int maxGopN, int maxGopM, bool fieldCoding, bool shotChange, int segmentLength, int baseQP)
    :
    state(new State(maxGopN, maxGopM, fieldCoding, shotChange, segmentLength, baseQP))
{
}


InputQueue::~InputQueue() = default;


void InputQueue::State::createDocket(int max, int i, int nut, int qpOffset, double qpFactor, int sopLevel, int hierarchyLevel, int numSameHierarchyLevel, int ref1, int ref2, int ref3, int ref4)
{
    assert(i > 0);
    if (i <= max)
    {
        Docket *docket = new Docket();
        docket->poc = this->segmentFront + i - 1;
        docket->absolutePoc = this->sequenceFront + i - 1;
        docket->nut = nut;
        docket->qpOffset = qpOffset;
        docket->qpFactor = qpFactor;
        docket->currentGopSize = gopSize;
        docket->sopLevel = sopLevel;
        docket->pocInSop = i;
        docket->sopId = this->segmentFront;
        docket->hierarchyLevel = this->hierarchyLevelFront + hierarchyLevel;
        docket->numSameHierarchyLevel = numSameHierarchyLevel;
        docket->intraFramePoc  = (this->currentIntraPoc + (this->sequenceFront - this->segmentFront)) <= docket->absolutePoc? this->currentIntraPoc : this->previousIntraPoc;

        if ((this->sequenceIdr == this->sequenceIFrame && docket->poc == this->sequenceIdr) || this->segmentFront == 0)
        {
            this->previousSegmentPoc = this->currentSegmentPoc;
            this->currentSegmentPoc = docket->absolutePoc;
        }

        docket->segmentPoc = this->currentSegmentPoc <= docket->absolutePoc ? this->currentSegmentPoc : this->previousSegmentPoc;

        if (this->isValidReference(i, ref1)) addReference(*docket, ref1);
        if (this->isValidReference(i, ref2)) addReference(*docket, ref2);
        if (this->isValidReference(i, ref3)) addReference(*docket, ref3);
        if (this->isValidReference(i, ref4)) addReference(*docket, ref4);

        docket->sliceType = isIrap(nut) ? I : B;
        this->entries.at(i - 1).docket.reset(docket);
    }
}


void InputQueue::State::process()
{
    if (this->entries.empty() ||  this->entries.front().docket) return;
    bool insertIdr = false; 
    if (this->shotChange)
    {
        if ((this->sequenceIFrame - this->segmentFront) < 0)
        {
            this->sequenceIFrame = computeNextIdr(this->segmentFront, this->sequenceFront, this->sequenceIFrame, this->maxGopN, insertIdr);
            if (insertIdr)
            {
                this->sequenceIdr = this->sequenceIFrame;

            }

        }
            
    }
    else
    {
        if (this->sequenceIFrame - this->segmentFront < 0)
            this->sequenceIFrame += this->maxGopN;
    }

    if (this->segmentLength != -1) 
    {
        if (this->sequenceSeg - this->segmentFront < 0)
            this->sequenceSeg += this->segmentLength;
    }


    gopSize = this->maxGopM;

    if (gopSize != 1 && gopSize != 8) 
        throw std::runtime_error("max-gop-m must be either 1 or 8"); // review: still the case?

    char lastPicture = 'P';

    if (this->finish && int(entries.size()) < gopSize)
    {
        gopSize = int(entries.size());
    }
    
    const int gopSizeIdr = this->sequenceIFrame - this->segmentFront + 1;
    if ((this->sequenceIFrame - this->segmentFront)>=0 && gopSizeIdr <= gopSize)
    {
        gopSize = gopSizeIdr;
        if (this->sequenceIFrame == this->sequenceIdr)
        {
            if (gopSize == 1)
            {
                lastPicture = 'R';
            }
            else
            {
                lastPicture = 'P';
                gopSize = gopSize - 1;
            }
        }
        else
        {
            lastPicture = 'I';
        }
    }

    if (this->segmentLength != -1)
    {
        const int gopSizeSeg = this->sequenceSeg - this->segmentFront + 1;
        if (gopSizeSeg <= gopSize) 
        {
            gopSize = gopSizeSeg;
            if (gopSizeSeg == 1) 
            {
                lastPicture = 'R';
            }
            else 
            {
                lastPicture = 'P';
                gopSize = gopSize - 1;
            }
        }
    }

  
    if (int(entries.size()) < gopSize) 
        return;

    int max = gopSize;

    auto nutR = TRAIL_R;
    auto nutN = TRAIL_N;

    if (lastPicture == 'I' || lastPicture == 'R')
    {
        this->previousIntraPoc = this->currentIntraPoc;
        this->currentIntraPoc = this->segmentFront + gopSize - 1;
    }

    if (lastPicture == 'I')
    {
        int nut = this->segmentFront ? CRA_NUT : IDR_N_LP;
        this->createDocket(gopSize, gopSize, nut, 0, 0.4420, 1, 0, 1,-gopSize);
        max = gopSize - 1;
        nutR = RASL_R;
        nutN = RASL_N;
    }
    else if (lastPicture == 'R') 
    {
        this->segmentFront = 0;
        this->sequenceIFrame = 0;
        this->hierarchyLevelFront = 0;
        this->currentIntraPoc = 0;
        this->previousIntraPoc = 0;
        int nut = IDR_N_LP;
        this->createDocket(gopSize, gopSize, nut, 0, 0.4420, 1, 0, 1, -gopSize);
        max = gopSize - 1;
        nutR = RASL_R;

        nutN = RASL_N;
    }
    else if (lastPicture == 'P')
    {
        this->createDocket(gopSize, gopSize, TRAIL_R,( (this->baseQP + 1) > 51 ? (51 - this->baseQP) : 1), 0.4420, 1, 0, 1, -gopSize, -gopSize);
        max = gopSize - 1;
    }

    if (gopSize == 2)
    {
        this->createDocket(max, 1, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.6800, 2, 1, 1, -1, 1);
    }
    else if (gopSize == 3)
    {
        this->createDocket(max, 2, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.3536, 2, 1, 1, -2, 1);
        this->createDocket(max, 1, nutN, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.6800, 3, 2, 1, -1, 2, 1);
    }
    else if (gopSize == 4)
    {
        this->createDocket(max, 2, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.3536, 2, 1, 1, -2, 2);
        this->createDocket(max, 1, nutN, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.6800, 3, 2, 1,-1, 3, 1);
        this->createDocket(max, 3, nutN, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.6800, 3, 3, 1,-1, 1);
    }
    else if (gopSize == 5)
    {
        this->createDocket(max, 3, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.3536, 2, 1, 1, -3, 2);
        this->createDocket(max, 1, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.3536, 2, 2, 1, -1, 4, 2);
        this->createDocket(max, 2, nutN, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.6800, 3, 3, 1, -2, 3, -1, 1);
        this->createDocket(max, 4, nutN, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.6800, 3, 4, 1, -4, 1, -1);
    }
    else if (gopSize == 6)
    {
        this->createDocket(max, 3, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.3536, 2, 1, 1, -3, 3);
        this->createDocket(max, 1, nutR, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.3536, 3, 2, 1, -1, 5, 2);
        this->createDocket(max, 2, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 3, 2, -2, 4, 1, -1);
        this->createDocket(max, 5, nutR, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.3536, 3, 3, 2, -5, 1, -2);
        this->createDocket(max, 4, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 4, 1, -4, 2, -1, 1);
    }
    else if (gopSize == 7)
    {
        this->createDocket(max, 4, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.3536, 2, 1, 1, -4, 3);
        this->createDocket(max, 2, nutR, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.3536, 3, 2, 1, -2, 5, 2);
        this->createDocket(max, 1, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 3, 1, -1, 6, 3, 1);
        this->createDocket(max, 3, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 4, 2, -3, 4, 1, -1);
        this->createDocket(max, 6, nutR, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.3536, 3, 4, 2, -2, 1);
        this->createDocket(max, 5, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 5, 1, -1, 2, 1);
    }
    else if (gopSize == 8)
    {
        this->createDocket(max, 4, nutR, ((this->baseQP + 2) > 51 ? (51 - this->baseQP) : 2), 0.3536, 2, 1, 1, -4, 4);
        this->createDocket(max, 2, nutR, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.3536, 3, 2, 1, -2, 2, 6);
        this->createDocket(max, 1, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 3, 1, -1, 1, 3, 7);
        this->createDocket(max, 3, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 4, 2, -1, 1, -3, 5);
        this->createDocket(max, 6, nutR, ((this->baseQP + 3) > 51 ? (51 - this->baseQP) : 3), 0.3536, 3, 4, 2, -2, 2, -6);
        this->createDocket(max, 5, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 5, 2, -1, 1, 3, -5);
        this->createDocket(max, 7, nutN, ((this->baseQP + 4) > 51 ? (51 - this->baseQP) : 4), 0.6800, 4, 5, 2, -1, 1, -7);
    }

    int sopSizeToMaxHierarchyLevel[8] = { 1, 2, 3, 4, 5, 5, 5, 6 };
    this->hierarchyLevelFront += sopSizeToMaxHierarchyLevel[gopSize - 1];
}


void InputQueue::append(std::shared_ptr<PictureWrapper> picture, std::shared_ptr<AdaptiveQuantisation> aqInfo)
{
    assert(!this->state->finish);
    this->state->entriesPreanalysis.push_back(Piece(picture));
    this->state->entriesPreanalysis.back().aqInfo = aqInfo;

    auto &timestamps = this->state->timestamps;

    size_t const reorderDelay = 3; // review: could reduce - this relates to *_max_num_reorder_pics

    if (this->state->picturePreprocessCount == 1)
    {
        // upon the second picture, we know the PTS period
        auto const period = picture->pts - timestamps.front().timestamp;
        for (size_t i=0; i<reorderDelay; ++i)
            timestamps.push_front({ timestamps.front().timestamp - period, i });
    }
    timestamps.push_back({ picture->pts, this->state->picturePreprocessCount + reorderDelay });
    ++this->state->picturePreprocessCount;
}


void InputQueue::endOfInput()
{
    this->state->finish = true;
}

void InputQueue::preanalyse()
{
    auto n = this->state->entriesPreanalysis.size();
    if (n == 0)
        return;
    auto const windowSize = this->state->shotChange ? 40 : this->state->maxGopM;

    if (n >= windowSize ||  this->state->finish)
    {
        if (n > windowSize)
            n = windowSize;

        if (this->state->shotChange)
            this->state->scd->processSeq(this->state->entriesPreanalysis, static_cast<int>(n));

        for (int i = 0; i < n; ++i)
        {
            this->state->entries.push_back(this->state->entriesPreanalysis.front());
            this->state->entriesPreanalysis.pop_front();
            ++this->state->pictureInputCount;
        }
    }
}

bool InputQueue::eos() const
{
    return this->state->finish;
}

std::shared_ptr<InputQueue::Docket> InputQueue::getDocket()
{
    if (this->state->pictureInputCount == 1 && !this->eos())
        return nullptr;
    
    this->state->process();

    int isIntraFrame = -1;
    bool encodeIntraFrame = false;
    int size = int(this->state->entries.size() > 8) ? 8 : int(this->state->entries.size());

    std::shared_ptr<InputQueue::Docket> docket;
    for (int i = 0; i < size; ++i)
    {
        docket = this->state->entries.at(i).docket;
        if (!docket) 
            break;

        bool canEncode = true;

        for (int deltaPoc : docket->references.positive)
            if (!this->state->entries.at(i + deltaPoc).done())
                canEncode = false;

        if (canEncode)
        {
            docket->picture = this->state->entries.at(i).picture;
            docket->aqInfo = this->state->entries.at(i).aqInfo;
            this->state->entries.at(i).picture = 0;
            break;
        }
    }

    while (!this->state->entries.empty() && this->state->entries.front().done())
    {
        this->state->entries.pop_front();
        ++this->state->segmentFront;
        ++this->state->sequenceFront;
    }

    if (docket)
    {
        docket->dts = this->state->timestamps.front().timestamp;
        this->state->timestamps.pop_front();
    }

    return docket;
}
