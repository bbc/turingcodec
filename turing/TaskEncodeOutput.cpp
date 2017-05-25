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

// Encoder output task responsible for writing all bitstream.


#include "TaskEncodeOutput.h"
#include "Write.h"
#include "StateEncode.h"
#include "SyntaxRbsp.hpp"
#include "SyntaxNal.hpp"
#include "sei/all.h"
#include "RateControl.h"
#include <string>


template <class H>
TaskEncodeOutput<H>::TaskEncodeOutput(H &h)
:
h(h)
{
}


template <class H>
bool TaskEncodeOutput<H>::blocked()
{
    StateEncode *stateEncode = this->h;

    for (auto &response : stateEncode->responses)
    {
        if (response.hungry) continue;

        if (response.eos)
        {
            return false;
        }

        StateWavefront *stateWavefront = response.picture.get();

        if (stateWavefront->nSubstreamsUnfinished != 0)
        {
            // Task is blocked: its next output picture is still being encoded by TaskEncodeSubstream.
            return true;
        }
        else if (stateEncode->decodedHashSei && !response.picture->reconstructed)
        {
            // Task is blocked.
            return true;
        }
        else if (!response.done)
        {
            return false;
        }
    }

    return true;
}


// emits an access unit, returns true if it's a keyframe
template <class H>
bool writeOut(H &h)
{
    StateEncode *stateEncode = h;
    StateEncodePicture *stateEncodePicture = h;

    {
        SliceSegmentHeaderIndependent *header = h;
        header->populateEntryPoints(stateEncodePicture->substreams);
    }

    auto const nut = h[nal_unit_type()];
    auto const isKeyframe = isIdr(nut) || isBla(nut) || isCra(nut);
    const bool writeParameterSets = isKeyframe && (stateEncodePicture->sequenceDecodeOrder == 0 || stateEncode->repeatHeaders);
    
    if (writeParameterSets)
    {
        writeHeaders(h);
    }
    else
    {
        // zero_byte is required on first NALU in AU
        h(zero_byte()  /* equal to 0x00 */, f(8));
    }

    // Alternative transfer characteristics
    bool writeAtcSei = (stateEncode->preferredTransferCharacteristics > -1) && (isIrap(nut) || h[PicOrderCntVal()] == 0);
    if (writeAtcSei)
    {
        h[nal_unit_type()] = PREFIX_SEI_NUT;
        h[last_payload_type_byte()] = PayloadTypeOf<alternative_transfer_characteristics>::value;
        h[preferred_transfer_characteristics()] = stateEncode->preferredTransferCharacteristics;
        h(byte_stream_nal_unit(0));
    }

    if(!stateEncode->userDataUnregSeiWritten)
    {
        // Review: add command line options
        auto const message = "Turing codec version " + std::string(turing_version());

        StateWriteUserDataUnregistered *stateWriteUserDataUnregistered = h;
        stateWriteUserDataUnregistered->p = message.c_str();
        stateWriteUserDataUnregistered->end = stateWriteUserDataUnregistered->p + message.length();

        h[uuid_iso_iec_11578()] = boost::uuids::string_generator()("ac9e584e2d484bdd8ccbf3ff9a878e69");

        h[nal_unit_type()] = PREFIX_SEI_NUT;
        h[last_payload_type_byte()] = PayloadTypeOf<user_data_unregistered>::value;
        h(byte_stream_nal_unit(0));

        stateEncode->userDataUnregSeiWritten = true;
    }

    const bool writeMasterDisplayInfoSei = stateEncode->masteringDisplayInfoPresent && (isIrap(nut) || h[PicOrderCntVal()] == 0);
    if(writeMasterDisplayInfoSei)
    {
        for(int c = 0; c < 3; c++)
        {
            h[display_primaries_x(c)] = stateEncode->masterDisplayInfo.displayPrimariesX[c];
            h[display_primaries_y(c)] = stateEncode->masterDisplayInfo.displayPrimariesY[c];
        }
        h[white_point_x()] = stateEncode->masterDisplayInfo.whitePointX;
        h[white_point_y()] = stateEncode->masterDisplayInfo.whitePointY;
        h[max_display_mastering_luminance()] = stateEncode->masterDisplayInfo.maxDisplayMasteringLuma;
        h[min_display_mastering_luminance()] = stateEncode->masterDisplayInfo.minDisplayMasteringLuma;

        h[nal_unit_type()] = PREFIX_SEI_NUT;
        h[last_payload_type_byte()] = PayloadTypeOf<mastering_display_colour_volume>::value;
        h(byte_stream_nal_unit(0));
    }

    if (stateEncode->fieldcoding)
    {
        PictureWrapper picturewrapper = *static_cast<StateEncodePicture *>(h)->docket->picture;

        // active_parameter_sets
        h[nal_unit_type()] = PREFIX_SEI_NUT;
        h[last_payload_type_byte()] = PayloadTypeOf<active_parameter_sets>::value;
        h[active_video_parameter_set_id()] = 0;
        h[self_contained_cvs_flag()] = 1;
        h(byte_stream_nal_unit(0));

        // pic_timing
        h[nal_unit_type()] = PREFIX_SEI_NUT;
        h[last_payload_type_byte()] = PayloadTypeOf<pic_timing>::value;
        h[frame_field_info_present_flag()] = 1;
        h[pic_struct()] = picturewrapper.fieldTB;
        h[source_scan_type()] = 2;
        h[duplicate_flag()] = 0;
        h(byte_stream_nal_unit(0));
    }
    else if (stateEncode->framedoubling)
    {
        // active_parameter_sets
        h[nal_unit_type()] = PREFIX_SEI_NUT;
        h[last_payload_type_byte()] = PayloadTypeOf<active_parameter_sets>::value;
        h[active_video_parameter_set_id()] = 0;
        h[self_contained_cvs_flag()] = 1;
        h(byte_stream_nal_unit(0));

        // pic_timing
        h[nal_unit_type()] = PREFIX_SEI_NUT;
        h[last_payload_type_byte()] = PayloadTypeOf<pic_timing>::value;
        h[frame_field_info_present_flag()] = 1;
        h[pic_struct()] = 7;
        h(byte_stream_nal_unit(0));
    }

    // Restore slice segment nal_unit_type
    h[nal_unit_type()] = nut;

    // Write the slice segment NALU (slice_segment_data() and slice_segment_trailing_bits() suppressed here)
    NalWriter *nalWriter = h;
    h(byte_stream_nal_unit(0));

    {
        // Write out substream data (already with emulation prevention applied)
        for (const auto &substream : stateEncodePicture->substreams)
        {
            nalWriter->data->insert(nalWriter->data->end(), substream.begin(), substream.end());
        }
    }

    if (stateEncode->decodedHashSei)
    {
        // Write the decoded_picture_hash
        h[nal_unit_type()] = SUFFIX_SEI_NUT;
        h[hash_type()] = stateEncode->hashType;
        h[last_payload_type_byte()] = PayloadTypeOf<decoded_picture_hash>::value;
        h(byte_stream_nal_unit(0));
    }

    if (stateEncode->useRateControl)
    {
        StateEncodePicture *stateEncodePicture = h;
        int currentPictureLevel = stateEncodePicture->docket->sopLevel;
        size_t codingBits = (nalWriter->data->size()) << 3;
        int segmentPoc = stateEncodePicture->docket->segmentPoc;
        stateEncode->rateControlParams->takeTokenOnRCLevel();
        stateEncode->rateControlMap.find(segmentPoc)->second->updateAfterEncoding(stateEncodePicture->docket, static_cast<int>(codingBits));
        stateEncode->rateControlMap.find(segmentPoc)->second->updateSequenceControllerFinishedFrames(stateEncodePicture->docket);
        stateEncode->rateControlParams->releaseTokenOnRCLevel();
    }

    return isKeyframe;
}


template <class H>
bool TaskEncodeOutput<H>::run()
{
    Profiler::Scope scope(static_cast<Profiler::Timers*>(h)->output);

    StateEncode *stateEncode = this->h;
    ThreadPool *threadPool = this->h;

    {
        std::unique_lock<std::mutex> lock(threadPool->mutex());
        for (auto &response : stateEncode->responses)
        {
            if (response.hungry)
            {
                continue;
            }

            if (response.eos)
            {
                response.done = true;
                stateEncode->responsesAvailable.notify_all();
                delete this;
                return false;
            }

            StateWavefront *stateWavefront = response.picture.get();

            if (stateWavefront->nSubstreamsUnfinished)
            {
                // blocked
                break;
            }

            if (stateEncode->decodedHashSei && !response.picture->reconstructed)
            {
                break;
            }

            if (!response.done)
            {
                auto h =  this->h.extend(&*response.picture);

                response.keyframe = writeOut(h);

                response.done = true;
                stateEncode->responsesAvailable.notify_all();
            }
        }

        return true;
    }
}


// Explict template instantiation
template struct TaskEncodeOutput<Handler<Encode<void>, StateEncode> >;
