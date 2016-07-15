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

// Deblocking top-level task.

#include "TaskDeblock.h"
#include "Global.h"
#include "Handlers.h"
#include "Write.h"
#include "Padding.h"


template <class H>
TaskDeblock<H>::TaskDeblock(H &h, std::shared_ptr<StateEncodePicture> stateEncodePicture, int begin, int end) :
    h(h),
    stateEncodePicture(stateEncodePicture),
    begin(begin),
    end(end),
    position(begin)
    {
        StateWavefront *stateWavefront = h;

        this->syncIn = &stateWavefront->encoded;
        this->syncOut = &stateWavefront->deblocked;
    }


template <class H>
bool TaskDeblock<H>::blocked()
{
    int rx = position % h[PicWidthInCtbsY()];
    int ry = position / h[PicWidthInCtbsY()];

    // Blocked by previous wavefront deblock task? (above CTU not yet processed)
    if (!(*this->syncOut)(rx, ry - 1)) return true;

    if (rx < h[PicWidthInCtbsY()] - 1) ++rx;
    if (ry < h[PicHeightInCtbsY()] - 1) ++ry;

    // Blocked by reconstruction task? (bottom-right CTU not ready)
    if (!(*this->syncIn)(rx, ry)) return true;

    return false;
}


template <class H>
bool TaskDeblock<H>::run()
{
    typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
    static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");

    Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->postprocess);

    ThreadPool *threadPool = h;

    //TODO: Review invocation order of this statement in SyntaxRbsb.hpp
    h[slice_deblocking_filter_disabled_flag()] = h[pps_deblocking_filter_disabled_flag()];

    Picture<Sample> *picture = &h[ReconstructedPicture()];

    while (true)
    {
        const int rx = this->position % h[PicWidthInCtbsY()];
        const int ry = this->position / h[PicWidthInCtbsY()];

        if (this->end != h[PicSizeInCtbsY()] && rx == 0)
        {
            assert(h[entropy_coding_sync_enabled_flag()]);

            // Create the next deblock wavefront task
            threadPool->add(newTaskDeblock(this->h, this->stateEncodePicture, this->end, this->end + h[PicWidthInCtbsY()]));
        }

        threadPool->lock();
        assert(!this->blocked());
        threadPool->unlock();

        auto &recPictureL = (*picture)[0];
        auto &recPictureCb = (*picture)[1];
        auto &recPictureCr = (*picture)[2];

        // deblocking
        if (!h[slice_deblocking_filter_disabled_flag()])
        {
            Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->deblock);

            {
                int xBegin = rx << h[CtbLog2SizeY()];
                int yBegin = ry << h[CtbLog2SizeY()];
                if (rx) xBegin += 8;
                if (ry) yBegin += 8;
                int xEnd = std::min(((rx + 1) << h[CtbLog2SizeY()]) + 8, h[pic_width_in_luma_samples()]);
                int yEnd = std::min(((ry + 1) << h[CtbLog2SizeY()]) + 8, h[pic_height_in_luma_samples()]);

                static_cast<StatePicture *>(h)->loopFilterPicture->deblock<EDGE_VER>(h, recPictureL, recPictureCb, recPictureCr, xBegin, yBegin, xEnd, yEnd);
            }

            {
                int xBegin = rx << h[CtbLog2SizeY()];
                int yBegin = ry << h[CtbLog2SizeY()];
                if (ry) yBegin += 8;
                int xEnd = std::min((rx + 1) << h[CtbLog2SizeY()], h[pic_width_in_luma_samples()]);
                int yEnd = std::min(((ry + 1) << h[CtbLog2SizeY()]) + 8, h[pic_height_in_luma_samples()]);

                static_cast<StatePicture *>(h)->loopFilterPicture->deblock<EDGE_HOR>(h, recPictureL, recPictureCb, recPictureCr, xBegin, yBegin, xEnd, yEnd);
            }
        }

        // todo: SAO here (and then may want to rename the task appropriately)

        // padding
        if (isSubLayerNonReferencePicture(h[nal_unit_type()]))
        {
            // No need to pad non-reference pictures (note: in fact they don't need even
            // to be allocated with a padding boundary).
        }
        else
        {
            Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->pad);

            const bool left = rx == 0;
            const bool top = ry == 0;
            const bool right = rx == h[PicWidthInCtbsY()] - 1;
            const bool bottom = ry == h[PicHeightInCtbsY()] - 1;

            int x0 = rx << h[CtbLog2SizeY()];
            int y0 = ry << h[CtbLog2SizeY()];

            const int pad = 80; // todo: be more intelligent

            int width = std::min(h[CtbSizeY()], h[pic_width_in_luma_samples()] - x0);
            int height = std::min(h[CtbSizeY()], h[pic_height_in_luma_samples()] - y0);
            Padding::padBlock(&recPictureL(x0, y0), width, height, recPictureL.stride, pad, top, bottom, left, right);

            x0 >>= 1;
            y0 >>= 1;
            width >>= 1;
            height >>= 1;

            Padding::padBlock(&recPictureCb(x0, y0), width, height, recPictureCb.stride, pad >> 1, top, bottom, left, right);
            Padding::padBlock(&recPictureCr(x0, y0), width, height, recPictureCr.stride, pad >> 1, top, bottom, left, right);
        }

        threadPool->lock();
        this->syncOut->set(rx, ry);
        const bool done = (++this->position == this->end);
        const bool blocked = done ? false : this->blocked();
        threadPool->nudge();
        threadPool->unlock();

        if (done) break;
        if (blocked) return true;
    }

    if (this->end == h[PicSizeInCtbsY()])
    {
        StateEncode *stateEncode = h;
        StatePicture *statePicture = h;
        ReconstructedPicture2<Sample> *currPic = h;

        if (stateEncode->psnrAnalysis)
        {
            StateEncodePicture *stateEncodePicture = h;
            PictureWrapper &pictureWrapper = *stateEncodePicture->docket->picture;
            typedef typename Access<Concrete<ReconstructedPictureBase>, H>::ActualType::Sample Sample;
            static_assert(std::is_same<Sample, uint8_t>::value || std::is_same<Sample, uint16_t>::value, "");
            auto &picture = dynamic_cast<Picture<Sample> &>(pictureWrapper);
            stateEncode->psnrAnalysis->analyse(*currPic->picture, picture);
        }

        if (stateEncode->decodedHashSei)
        {
            stateEncode->hashElement.clear();
            auto &picture = *currPic->picture;

            if (stateEncode->hashType == MD5)
            {
                for (int cIdx = 0; cIdx < 3; cIdx++)
                {
                    md5_byte_t digest[16];
                    int currentBitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
                    if(currentBitDepth == 8)
                    {
                        stateEncode->decodedPictureHash.computeMd5Sum(
                                digest,
                                (unsigned char*)picture[cIdx].p,
                                picture[cIdx].stride,
                                picture[cIdx].height,
                                picture[cIdx].width);
                    }
                    else
                    {
                        stateEncode->decodedPictureHash.computeMd5Sum(
                                digest,
                                (unsigned short*)picture[cIdx].p,
                                picture[cIdx].stride,
                                picture[cIdx].height,
                                picture[cIdx].width);
                    }

                    for (int i = 0; i < 16; i++)
                    {
                        h[picture_md5(cIdx, i)] = digest[i];
                        stateEncode->hashElement.push_back(digest[i]);
                    }
                }
            }
            else if (stateEncode->hashType == CRC)
            {
                for (int cIdx = 0; cIdx < 3; cIdx++)
                {
                    int currentBitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
                    int crc;
                    if(currentBitDepth == 8)
                    {
                        crc = stateEncode->decodedPictureHash.computeCrc((unsigned char*)picture[cIdx].p,
                                                                             picture[cIdx].stride,
                                                                             picture[cIdx].height,
                                                                             picture[cIdx].width,
                                                                             currentBitDepth);
                    }
                    else
                    {
                        crc = stateEncode->decodedPictureHash.computeCrc((unsigned short*)picture[cIdx].p,
                                                                             picture[cIdx].stride,
                                                                             picture[cIdx].height,
                                                                             picture[cIdx].width,
                                                                             currentBitDepth);
                    }
                    h[picture_crc(cIdx)] = crc;
                    stateEncode->hashElement.push_back((crc >> 8) & 0xff);
                    stateEncode->hashElement.push_back(crc & 0xff);
                }
            }
            else if (stateEncode->hashType == CHKSUM)
            {
                for (int cIdx = 0; cIdx < 3; cIdx++)
                {
                    int currentBitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
                    int checkSum;
                    if(currentBitDepth == 8)
                    {
                        checkSum = stateEncode->decodedPictureHash.computeCheckSum((unsigned char*)picture[cIdx].p,
                                                                                       picture[cIdx].stride,
                                                                                       picture[cIdx].height,
                                                                                       picture[cIdx].width,
                                                                                       currentBitDepth
                                                                                       );
                    }
                    else
                    {
                        checkSum = stateEncode->decodedPictureHash.computeCheckSum((unsigned short*)picture[cIdx].p,
                                                                                       picture[cIdx].stride,
                                                                                       picture[cIdx].height,
                                                                                       picture[cIdx].width,
                                                                                       currentBitDepth
                                                                                       );
                    }
                    h[picture_checksum(cIdx)] = checkSum;
                    stateEncode->hashElement.push_back((checkSum >> 24) & 0xff);
                    stateEncode->hashElement.push_back((checkSum >> 16) & 0xff);
                    stateEncode->hashElement.push_back((checkSum >> 8) & 0xff);
                    stateEncode->hashElement.push_back(checkSum & 0xff);
                }
            }
        }

        // Signal that picture is now completely reconstructed and postprocessed.
        // Note: could signal this before padding.
        {

            std::unique_lock<std::mutex> lock(threadPool->mutex());

            statePicture->reconstructed = true;
            stateEncode->responsesAvailable.notify_all();
        }

    }

    delete this;
    return false;
}


// explict template instantiations
template struct TaskDeblock<Handler<Encode<void>, StateEncodePicture2<uint8_t>, StateEncode>>;
template struct TaskDeblock<Handler<Encode<void>, StateEncodePicture2<uint16_t>, StateEncode>>;
