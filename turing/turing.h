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

#ifndef INCLUDED_turing_h
#define INCLUDED_turing_h

#define TURING_API_VERSION 2

/* This file represents the public API of the Turing encoder */

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct turing_encoder_settings
{
    int argc;
    char const* const* argv;
} turing_encoder_settings;


typedef struct turing_raster
{
    uint8_t *p;
    size_t stride; // in bytes
} turing_raster;


typedef struct turing_picture
{
    turing_raster image[3];
    int64_t pts;
} turing_picture;


typedef struct turing_bitstream
{
    uint8_t *p;
    int size;
} turing_bitstream;

typedef struct turing_encoder_output
{
    turing_bitstream bitstream;
    int64_t pts;
    int64_t dts;
    int keyframe;
} turing_encoder_output;

typedef struct turing_encoder turing_encoder;


const char *turing_version(void);

turing_encoder *turing_create_encoder(turing_encoder_settings settings);

// returned pointer is good until the next encoder call
turing_bitstream const* turing_encode_headers(turing_encoder *encoder);

// Encode one picture
// Set 'headers' non-zero to allow encoder to write headers on keyframes
// returned pointer is good until the next encoder call
turing_encoder_output const* turing_encode_picture(turing_encoder *encoder, turing_picture *picture);

void turing_destroy_encoder(turing_encoder *encoder);

int turing_check_binary_option(const char *option);


#ifdef __cplusplus
}
#endif

#endif
