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

#include "TaskSao.h"
#include "Global.h"
#include "Handlers.h"
#include "Write.h"
#include "Padding.h"

template <class H>
TaskSao<H>::TaskSao(H &h, std::shared_ptr<StateEncodePicture> stateEncodePicture, int begin, int end)
    :
    h(h),
    stateEncodePicture(stateEncodePicture),
    begin(begin),
    end(end),
    position(begin)
{
    StateWavefront *stateWavefront = h;

    this->syncIn = &stateWavefront->deblocked;
    this->syncOut = &stateWavefront->saoed;
}


template <class H>
bool TaskSao<H>::blocked()
{
    int rx = position % h[PicWidthInCtbsY()];
    int ry = position / h[PicWidthInCtbsY()];

    // Blocked by previous wavefront sao task? (above CTU not yet processed)
    if (!(*this->syncOut)(rx, ry - 1)) return true;

    if (rx < h[PicWidthInCtbsY()] - 1) ++rx;
    if (ry < h[PicHeightInCtbsY()] - 1) ++ry;

    // Blocked by deblock task? (bottom-right CTU not ready)
    if (!(*this->syncIn)(rx, ry)) 
        return true;

    return false;
}


template <class H>
bool TaskSao<H>::run()
{
    using Sample = typename SampleType<H>::Type;

    Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->postprocess);

    ThreadPool *threadPool = h;

    StateReconstructedPicture<Sample> *stateReconstructedPicture = h;
    Picture<Sample> *picture = stateReconstructedPicture->picture.get();
    Picture<Sample> *saoPicture = stateReconstructedPicture->saoPicture.get();
    while (true)
    {
        const int rx = this->position % h[PicWidthInCtbsY()];
        const int ry = this->position / h[PicWidthInCtbsY()];

        if (this->end != h[PicSizeInCtbsY()] && rx == 0)
        {
            assert(h[entropy_coding_sync_enabled_flag()]);
            // Create the next deblock sao task
            threadPool->add(newTaskSao(this->h, this->stateEncodePicture, this->end, this->end + h[PicWidthInCtbsY()]));
        }

        threadPool->lock();
        assert(!this->blocked());
        threadPool->unlock();

        auto &recPictureL = (*picture)[0];
        auto &recPictureCb = (*picture)[1];
        auto &recPictureCr = (*picture)[2];

        if (h[sample_adaptive_offset_enabled_flag()])
        {
            for (int cIdx = 0; cIdx < 3; cIdx++)
            {
                int xBegin = rx << h[CtbLog2SizeY()];
                int yBegin = ry << h[CtbLog2SizeY()];
                int xEnd = std::min(((rx + 2) << h[CtbLog2SizeY()]), h[pic_width_in_luma_samples()]);
                int yEnd = std::min(((ry + 2) << h[CtbLog2SizeY()]), h[pic_height_in_luma_samples()]);
                if (cIdx != 0)
                {
                    xBegin >>= 1;
                    yBegin >>= 1;
                    xEnd >>= 1;
                    yEnd >>= 1;
                }
                for (int y = yBegin; y < yEnd; ++y)
                {
                    for (int x = xBegin; x < xEnd; ++x)
                    {
                        (*saoPicture)[cIdx](x, y) = (*picture)[cIdx](x, y);
                    }
                }
            }

            static_cast<StatePicture *>(h)->loopFilterPicture->applySaoCTU<Sample>(h, rx, ry);
        }
    
        // padding
        if (isSubLayerNonReferencePicture(h[nal_unit_type()]))
        {
            // No need to pad non-reference pictures (note: in fact they don't need even
            // to be allocated with a padding boundary).
        }
        else if (h[sample_adaptive_offset_enabled_flag()])
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
        StateReconstructedPicture<Sample> *currPic = h;
        const int currentPoc = stateEncodePicture->docket->absolutePoc;
        if (stateEncode->psnrAnalysis)
        {
            StateEncodePicture *stateEncodePicture = h;
            auto &pictureInput = static_cast<PictureWrap<Sample> &>(*static_cast<StateEncodePicture *>(h)->docket->picture);
            stateEncode->psnrAnalysis->analyse(currentPoc, *currPic->picture, pictureInput);
        }

        if (stateEncode->decodedHashSei)
        {
            StateEncode::FrameHash currentFrameHash;
            auto &picture = *currPic->picture;

            if (stateEncode->hashType == MD5)
            {
                for (int cIdx = 0; cIdx < 3; cIdx++)
                {
                    md5_byte_t digest[16];
                    int currentBitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
                    if (currentBitDepth == 8)
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
                        currentFrameHash.hash.push_back(digest[i]);
                    }
                }
                stateEncode->addFrameHash(currentPoc, currentFrameHash);
            }
            else if (stateEncode->hashType == CRC)
            {
                for (int cIdx = 0; cIdx < 3; cIdx++)
                {
                    int currentBitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
                    int crc;
                    if (currentBitDepth == 8)
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
                    currentFrameHash.hash.push_back((crc >> 8) & 0xff);
                    currentFrameHash.hash.push_back(crc & 0xff);
                }
                stateEncode->addFrameHash(currentPoc, currentFrameHash);
            }
            else if (stateEncode->hashType == CHKSUM)
            {
                for (int cIdx = 0; cIdx < 3; cIdx++)
                {
                    int currentBitDepth = cIdx ? h[BitDepthC()] : h[BitDepthY()];
                    int checkSum;
                    if (currentBitDepth == 8)
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
                    currentFrameHash.hash.push_back((checkSum >> 24) & 0xff);
                    currentFrameHash.hash.push_back((checkSum >> 16) & 0xff);
                    currentFrameHash.hash.push_back((checkSum >> 8) & 0xff);
                    currentFrameHash.hash.push_back(checkSum & 0xff);
                }
                stateEncode->addFrameHash(currentPoc, currentFrameHash);
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

// explicit template instantiations
template struct TaskSao<Handler<Encode<void>, StateEncodePicture2<uint8_t>, StateEncode>>;
template struct TaskSao<Handler<Encode<void>, StateEncodePicture2<uint16_t>, StateEncode>>;
