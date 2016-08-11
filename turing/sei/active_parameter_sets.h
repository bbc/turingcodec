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

#pragma once

#include "../SyntaxSei.h"


//

struct active_video_parameter_set_id { };
struct self_contained_cvs_flag { };
struct no_parameter_set_update_flag { };
struct num_sps_ids_minus1 { };
DEFINE_VALUE_ARRAY_1(active_seq_parameter_set_id, i, 16);


// review: struct name
struct ActiveParameterSets2 :
    ValueHolder<active_video_parameter_set_id>,
    ValueHolder<self_contained_cvs_flag>,
    ValueHolder<no_parameter_set_update_flag>,
    ValueHolder<num_sps_ids_minus1>,
    ValueHolder<active_seq_parameter_set_id>
    {
    };

template <>
struct Syntax<active_parameter_sets>
{
    template <class H> static void go(active_parameter_sets fun, H &h)
    {
        h(active_video_parameter_set_id(), u(4));
        h(self_contained_cvs_flag(), u(1));
        h(no_parameter_set_update_flag(), u(1));
        h(num_sps_ids_minus1(), ue(v));
        for (int i = 0; i <= h[num_sps_ids_minus1()]; i++)
        {
            h(active_seq_parameter_set_id(i), ue(v));
        }
    }
};


template <class H> void Read<active_parameter_sets>::go(active_parameter_sets f, H &h)
{
    ActiveParameterSets2 activeParameterSets;
    auto h3 = h.extend(&activeParameterSets);

    Syntax<active_parameter_sets>::go(f, h3);

    h[Active<Vps>()] = h[Table<Vps>()][h3[active_video_parameter_set_id()]];
    h[Active<Sps>()] = h[Table<Sps>()][h3[active_seq_parameter_set_id(0)]];
}


#ifdef REVIEW_LATER

                template <class H> static void Read<buffering_period>::go(buffering_period f, H &h)
                {
                    BufferingPeriod bufferingPeriod;
                    auto h3 = h.extend(&bufferingPeriod);

                    auto h2 = h3.extend(getHrd(h));

                    Syntax<buffering_period>::go(f, h2);
                }


                template <class H> void Read<pic_timing>::go(pic_timing f, H &h)
        {
                    Hrd *hrd = getHrd(h);

                    if (hrd)
                    {
                        auto h2 = extendHandler(*hrd, h);
                        PicTiming picTiming;
                        auto h3 = extendHandler(picTiming, h2);

                        Syntax<pic_timing>::go(f, h3);
                    }
                    else
                    {
                        // Seek to end of SEI payload so as not to trigger a further error
                        seek(h, bitLen(h[::Stream()]));
                    }
        }


                template <class H> void Read<decoded_picture_hash>::go(decoded_picture_hash f, H &h)
        {
                    DecodedPictureHash decodedPictureHash;
                    auto h3 = h.extend(&decodedPictureHash);

                    Syntax<decoded_picture_hash>::go(f, h3);
        }


                template <class H> void Read<user_data_registered_itu_t_t35>::go(user_data_registered_itu_t_t35 f, H &h)
        {
                    UserDataRegistered userDataRegistered;
                    auto h3 = h.extend(&userDataRegistered);

                    Syntax<user_data_registered_itu_t_t35>::go(f, h3);
        }

                template <class H> void Read<user_data_unregistered>:: go(user_data_unregistered f, H &h)
        {
                    UserDataUnregistered userDataRegistered;
                    auto h3 = h.extend(&userDataRegistered);

                    Syntax<user_data_unregistered>::go(f, h3);
        }


                template <class H> void Read<recovery_point>::go(recovery_point f, H &h)
        {
                    RecoveryPoint recoveryPoint;
                    auto h3 = h.extend(&recoveryPoint);

                    Syntax<recovery_point>::go(f, h3);
        }


                template <class H> void Read<scene_info>::go(scene_info f, H &h)
        {
                    SceneInfo sceneInfo;
                    auto h3 = h.extend(&sceneInfo);

                    Syntax<scene_info>::go(f, h3);
        }


                template <class H> void Read<picture_snapshot>::go(picture_snapshot f, H &h)
        {
                    PictureSnapshot pictureSnapshot;
                    auto h3 = h.extend(&pictureSnapshot);

                    Syntax<picture_snapshot>::go(f, h3);
        }


                template <class H> void Read<progressive_refinement_segment_start>::go(progressive_refinement_segment_start f, H &h)
        {
                    ProgressiveRefinementSegmentStart progressiveRefinementSegmentStart;
                    auto h3 = h.extend(&progressiveRefinementSegmentStart);

                    Syntax<progressive_refinement_segment_start>::go(f, h3);
        }


                template <class H> void Read<progressive_refinement_segment_end>::go(progressive_refinement_segment_end f, H &h)
        {
                    ProgressiveRefinementSegmentEnd progressiveRefinementSegmentEnd;
                    auto h3 = h.extend(&progressiveRefinementSegmentEnd);

                    Syntax<progressive_refinement_segment_end>::go(f, h3);
        }


                template <class H> void Read<film_grain_characteristics>::go(film_grain_characteristics f, H &h)
        {
                    FilmGrainCharacteristics filmGrainCharacteristics;
                    auto h3 = h.extend(&filmGrainCharacteristics);

                    Syntax<film_grain_characteristics>::go(f, h3);
        }


                template <class H> void Read<post_filter_hint>::go(post_filter_hint f, H &h)
        {
                    PostFilterHint postFilterHint;
                    auto h3 = extendHandler(postFilterHint, h);

                    Syntax<post_filter_hint>::go(f, h3);
        }


                template <class H> void Read<tone_mapping_info>::go(tone_mapping_info f, H &h)
        {
                    ToneMappingInfo toneMappingInfo;
                    auto h3 = h.extend(&toneMappingInfo);

                    Syntax<tone_mapping_info>::go(f, h3);
        }




                template <class H> void Read<pan_scan_rect>:: go(pan_scan_rect f, H &h)
        {
                    PanScanRect panScanRect;
                    auto h3 = h.extend(&panScanRect);

                    Syntax<pan_scan_rect>::go(f, h3);
        }


                template <class H> void Read<frame_packing_arrangement>::go(frame_packing_arrangement f, H &h)
        {
                    FramePackingArrangement framePackingArrangement;
                    auto h3 = h.extend(&framePackingArrangement);

                    Syntax<frame_packing_arrangement>::go(f, h3);
        }


                template <class H> void Read<display_orientation>::go(display_orientation f, H &h)
        {
                    DisplayOrientation displayOrientation;
                    auto h3 = h.extend(&displayOrientation);

                    Syntax<display_orientation>::go(f, h3);
        }


                template <class H> void Read<structure_of_pictures_info>::go(structure_of_pictures_info f, H &h)
        {
                    StructureOfPicturesInfo structureOfPicturesInfo;
                    auto h3 = h.extend(&structureOfPicturesInfo);

                    Syntax<structure_of_pictures_info>::go(f, h3);
        }


                template <class H> void Read<active_parameter_sets>::go(active_parameter_sets f, H &h)
        {
                    ActiveParameterSets2 activeParameterSets;
                    auto h3 = h.extend(&activeParameterSets);

                    Syntax<active_parameter_sets>::go(f, h3);

                    h[Active<Vps>()] = h[Table<Vps>()][h3[active_video_parameter_set_id()]];
                    h[Active<Sps>()] = h[Table<Sps>()][h3[active_seq_parameter_set_id(0)]];
        }


                template <class H> void Read<decoding_unit_info>::go(decoding_unit_info f, H &h)
        {
                    Hrd *hrd = getHrd(h);

                    if (hrd)
                    {
                        auto h2 = h.extend(hrd);

                        DecodingUnitInfo decodingUnitInfo;
                        auto h3 = h2.extend(&decodingUnitInfo);

                        Syntax<decoding_unit_info>::go(f, h3);
                    }
                    else
                    {
                        // Seek to end of SEI payload so as not to trigger a further error
                        seek(h, bitLen(h[::Stream()]));
                    }
        }


                template <class H> static void Read<temporal_sub_layer_zero_index>:go(temporal_sub_layer_zero_index f, H &h)
        {
                    TemporalSubLayerZeroIndex temporalSubLayerZeroIndex;
                    auto h3 = h.extend(&temporalSubLayerZeroIndex);

                    Syntax<temporal_sub_layer_zero_index>::go(f, h3);
        }

                template <class H> static void Read<scalable_nesting>::go(scalable_nesting f, H &h)
        {
                    ScalableNesting *scalableNesting = h;

                    if (scalableNesting->nested)
                    {
                        h(Violation("D.3.1", "scalable_nesting() within scalable_nesting()")); // CondCheck D.3.1-M
                    }
                    else
                    {
                        scalableNesting->nested = true;
                        Syntax<scalable_nesting>::go(f, h);
                        scalableNesting->nested = false;
                    }
        }


#endif
