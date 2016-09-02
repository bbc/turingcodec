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

struct scene_info_present_flag { };
struct prev_scene_id_valid_flag { };
struct scene_id { };
struct scene_transition_type { };
struct second_scene_id { };

template <>
struct Syntax<scene_info>
{
    template <class H> static void go(scene_info fun, H &h);
};

template <class H>
void Syntax<scene_info>::go(scene_info fun, H &h)
{
    h(scene_info_present_flag(), u(1));
    if (h[scene_info_present_flag()])
    {
        h(prev_scene_id_valid_flag(), u(1));
        h(scene_id(), ue(v));
        h(scene_transition_type(), ue(v));
        if (h[scene_transition_type()] > 3)
            h(second_scene_id(), ue(v));
    }
}


struct SceneInfo :
    ValueHolder<scene_info_present_flag>,
    ValueHolder<prev_scene_id_valid_flag>,
    ValueHolder<scene_id>,
    ValueHolder<scene_transition_type>,
    ValueHolder<second_scene_id>
    {
    };


template <class H> void Read<scene_info>::go(scene_info f, H &h)
{
    SceneInfo sceneInfo;
    auto h3 = h.extend(&sceneInfo);

    Syntax<scene_info>::go(f, h3);
}
