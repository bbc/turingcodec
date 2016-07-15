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


namespace {

    struct Piece
    {
        Piece(std::shared_ptr<PictureWrapper> picture) : picture(picture) { }

        std::shared_ptr<InputQueue::Docket> docket;
        std::shared_ptr<PictureWrapper> picture;
        std::shared_ptr<AdaptiveQuantisation> aqInfo;

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
    State(const boost::program_options::variables_map &vm) :
        vm(vm),
        sequenceFront(0),
        sequenceIdr(0),
        finish(false)
        {
        }

    const boost::program_options::variables_map &vm;
    std::deque<Piece> entries;
    int sequenceFront;
    int sequenceIdr;
    bool finish;
    int gopSize;

    void process();

    void createDocket(int size, int i, int nut, int qpOffset, double qpFactor, int ref1 = 0, int ref2 = 0, int ref3 = 0, int ref4 = 0);

    bool isValidReference(int i, int delta) const
    {
        assert(i > 0);
        if (delta < 0) return true;
#ifdef FORCE_P_SLICES
        return false;
#endif
        if (delta == 0) return false;
        int offset = sequenceIdr - sequenceFront + 1;
        int deltaLimit = (static_cast<int>(this->entries.size() > offset)) ? offset : static_cast<int>(this->entries.size());
        if (i - 1 + delta >= deltaLimit) return false;
        return !!this->entries.at(i - 1 + delta).docket;
    }

    std::vector<int>        m_shotChangeList;
    void setShotChangeList(std::vector<int>& shotChangeList) { if (shotChangeList.size()) m_shotChangeList.swap(shotChangeList); };
    int computeNextIdr(int sequenceFront, int nextDefaultIdr, bool fieldCoding = false)
    {
        int nextIdr = nextDefaultIdr;
        int scale = fieldCoding ? 1 : 0;
        for (int i = sequenceFront; i < nextDefaultIdr; i++)
        {
            int index = (i >> scale);
            index = (index < 0) ? 0 : index;
            if (index >= m_shotChangeList.size())
                break;

            if (i % (fieldCoding + 1) == 0 && m_shotChangeList[index])
            {
                nextIdr = i;
                break;
            }
        }
        return nextIdr;
    }
};


InputQueue::InputQueue(const boost::program_options::variables_map &vm)
:
	        state(new State(vm))
{
}


void InputQueue::State::createDocket(int max, int i, int nut, int qpOffset, double qpFactor, int ref1, int ref2, int ref3, int ref4)
{
    assert(i > 0);
    if (i <= max)
    {
        Docket *docket = new Docket();
        docket->poc = this->sequenceFront + i - 1;
        docket->nut = nut;
        docket->qpOffset = qpOffset;
        docket->qpFactor = qpFactor;
        docket->currentGopSize = gopSize;

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

    if (this->entries.empty() || this->entries.front().docket) return;
    if (this->vm.at("shot-change").as<bool>())
    {
        if ((this->sequenceIdr - this->sequenceFront) < 0)
            this->sequenceIdr = computeNextIdr(this->sequenceFront, (this->sequenceIdr + this->vm.at("max-gop-n").as<int>()), this->vm.at("field-coding").as<bool>());
    }
    else
        if (this->sequenceIdr - this->sequenceFront < 0) this->sequenceIdr += this->vm.at("max-gop-n").as<int>();

    gopSize = this->vm.at("max-gop-m").as<int>();

    if (gopSize != 1 && gopSize != 8) throw std::runtime_error("max-gop-m must be either 1 or 8");
    char lastPicture = 'P';

    if (this->finish && int(entries.size()) < gopSize)
    {
        gopSize = int(entries.size());
        lastPicture = 0;
    }

    const int gopSizeIdr = this->sequenceIdr - this->sequenceFront + 1;
    if (gopSizeIdr <= gopSize)
    {
        gopSize = gopSizeIdr;
        lastPicture = 'I';
    }
    if (int(entries.size()) < gopSize) return;

    int max = gopSize;

    auto nutR = TRAIL_R;
    auto nutN = TRAIL_N;

    if (lastPicture == 'I')
    {
        int nut = this->sequenceFront ? CRA_NUT : IDR_N_LP;
        this->createDocket(gopSize, gopSize, nut, 0, 0.4420, -gopSize);
        max = gopSize - 1;
        nutR = RASL_R;
        nutN = RASL_N;
    }
    else if (lastPicture == 'P')
    {
        this->createDocket(gopSize, gopSize, TRAIL_R, 1, 0.4420, -gopSize, -gopSize);
        max = gopSize - 1;
    }

    if (!this->finish && gopSize != 8)
    {
        if (gopSize == 2)
        {
            this->createDocket(max, 1, nutR, 2, 0.6800, -1, 1);
        }
        else if (gopSize == 3)
        {
            this->createDocket(max, 2, nutR, 2, 0.3536, -2, 1);
            this->createDocket(max, 1, nutN, 3, 0.6800, -1, 2, 1);
        }
        else if (gopSize == 4)
        {
            this->createDocket(max, 2, nutR, 2, 0.3536, -2, 2);
            this->createDocket(max, 1, nutN, 3, 0.6800, -1, 3, 1);
            this->createDocket(max, 3, nutN, 3, 0.6800, -1, 1);
        }
        else if (gopSize == 5)
        {
            this->createDocket(max, 3, nutR, 2, 0.3536, -3, 2);
            this->createDocket(max, 1, nutR, 2, 0.3536, -1, 4, 2);
            this->createDocket(max, 2, nutN, 3, 0.6800, -2, 3, -1, 1);
            this->createDocket(max, 4, nutN, 3, 0.6800, -4, 1, -1);
        }
        else if (gopSize == 6)
        {
            this->createDocket(max, 3, nutR, 2, 0.3536, -3, 3);
            this->createDocket(max, 1, nutR, 3, 0.3536, -1, 5, 2);
            this->createDocket(max, 2, nutN, 4, 0.6800, -2, 4, 1, -1);
            this->createDocket(max, 5, nutR, 3, 0.3536, -5, 1, -2);
            this->createDocket(max, 4, nutN, 4, 0.6800, -4, 2, -1, 1);
        }
        else if (gopSize == 7)
        {
            this->createDocket(max, 4, nutR, 2, 0.3536, -4, 3);
            this->createDocket(max, 2, nutR, 3, 0.3536, -2, 5, 2);
            this->createDocket(max, 1, nutN, 4, 0.6800, -1, 6, 3, 1);
            this->createDocket(max, 3, nutN, 4, 0.6800, -3, 4, 1, -1);
            this->createDocket(max, 6, nutR, 3, 0.3536, -2, 1);
            this->createDocket(max, 5, nutN, 4, 0.6800, -1, 2, 1);
        }
    }
    else
    {
        this->createDocket(max, 4, nutR, 2, 0.3536, -4, 4);
        this->createDocket(max, 2, nutR, 3, 0.3536, -2, 2, 6);
        this->createDocket(max, 1, nutN, 4, 0.6800, -1, 1, 3, 7);
        this->createDocket(max, 3, nutN, 4, 0.6800, -1, 1, -3, 5);
        this->createDocket(max, 6, nutR, 3, 0.3536, -2, 2, -6);
        this->createDocket(max, 5, nutN, 4, 0.6800, -1, 1, 3, -5);
        this->createDocket(max, 7, nutN, 4, 0.6800, -1, 1, -7);
    }
}


void InputQueue::append(std::shared_ptr<PictureWrapper> picture, std::shared_ptr<AdaptiveQuantisation> aqInfo)
{
    assert(!this->state->finish);
    this->state->entries.push_back(Piece(picture));
    this->state->entries.back().aqInfo = aqInfo;
}


void InputQueue::endOfInput()
{
    this->state->finish = true;
}


bool InputQueue::eos() const
{
    return this->state->finish;
}


std::shared_ptr<InputQueue::Docket> InputQueue::getDocket()
{

    int currentGopSize;
    this->state->setShotChangeList(m_shotChangeList);
    this->state->process();

    std::shared_ptr<InputQueue::Docket> docket;
    int isIntraFrame = -1;
    bool encodeIntraFrame = false;
    int size = int(this->state->entries.size() > 8) ? 8 : int(this->state->entries.size());
    for (int i = 0; i < size; ++i)
    {
        //docket = this->state->entries.at(i).docket;
        if (!this->state->entries.at(i).docket || this->state->entries.at(i).done()) break;
        if (this->state->entries.at(i).docket->sliceType == I)
        {
            isIntraFrame = i;
            encodeIntraFrame = true;
        }
    }

    for (int i = 0; i < int(this->state->entries.size()); ++i)
    {
        docket = this->state->entries.at(i).docket;
        if (!docket) break;

        bool canEncode = true;
        for (int deltaPoc : docket->references.positive)
        {
            if (!this->state->entries.at(i + deltaPoc).done()) canEncode = false;
        }
        if (encodeIntraFrame && i != isIntraFrame)
        {
            canEncode = false;
        }
        if (encodeIntraFrame && i == isIntraFrame)
        {
            canEncode = true;
        }
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
        ++this->state->sequenceFront;
    }

    return docket;
}

