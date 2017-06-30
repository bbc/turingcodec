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

// Encoder API implementation and command-line encoder

#include "Encoder.h"
#include "ProgressReporter.h"
#include "Speed.h"
#include "Picture.h"
#include "SCDetection.h"
#include "turing.h"
#include "git-describe.h"
//#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <cstring>
#include <map>


#pragma optimize ("", off)

namespace po = boost::program_options;
using std::string;


// custom parsing for "--speed" command-line option
std::istream& operator >> (std::istream& is, Speed::Type &speed)
{
    std::string value;
    is >> value;
    bool valid = false;
#define X(name) if (value == #name) { speed = Speed::name; valid = true; }
    ENCODER_SPEED_PRESETS_XMACRO
#undef X
    if (!valid) throw po::invalid_option_value(value);
    return is;
}


int parseEncodeOptions(po::variables_map &vm, int argc, const char* const argv[], std::ostream &cout, std::ostream &cerr, std::string &inputFilename, std::string &outputFilename)
{
    po::options_description options;

    po::options_description optionsInput("Input options");
    optionsInput.add_options()
                ("input-res", po::value<std::string>()->required(), "video frame resolution <width>x<height>")
                ("seek", po::value<size_t>(0), "number of initial frames to be omitted before encoding starts")
                ("frames", po::value<size_t>()->required(), "number of frames to encode") // <review - would be better if this was not required so encoder can encode whole sequence
                ("frame-rate", po::value<double>()->required(), "sequence frame rate") // <review - this needs to be a rational number for, e.g. 30000/1001 framerates
                ("bit-depth", po::value<int>()->default_value(8), "video bit depth (of both input YUV and output stream)");
    options.add(optionsInput);

    po::options_description optionsOutput("Output options");
    optionsOutput.add_options()
                ("output-file,o", po::value<std::string>(&outputFilename), "output file name")
                ("dump-pictures", po::value<std::string>(), "reconstructed YUV file name (separate fields if field-coding is enabled)")
                ("dump-frames", po::value<std::string>(), "reconstructed YUV file name (interleaved frames if field-coding is enabled)")
                ("hash", po::value<int>(), "Decoded picture hash: 0 = MD5 sum, 1 = CRC, 2 = Check sum")
                ("atc-sei", po::value<int>()->default_value(-1), "Alternative transfer characteristics SEI message: --atc-sei ptc (preferred transfer characteristics)")
                ("mastering-display-info", po::value<std::string>(), "Mastering display colour volume SEI message: --mastering-display-info <string> "
                        "where string has the following format:\n"
                        "\"G(x,y)B(x,y)R(x,y)WP(x,y)L(max,min)\"\n"
                        "- with G(x,y) denoting the normalised x and y chromaticity coordinates for the green colour (and similarly for blue and red),\n"
                        "- WP(x,y) denoting the normalised x and y chromaticity coordinates, respectively, of the white point of the mastering display and\n"
                        "- L(max,min) denoting the nominal maximum and minimum display luminance, respectively, of the mastering display\n"
                        "See Sections D.2.27 and D.3.27 of the HEVC/H.265 spec. for more information");
    options.add(optionsOutput);

    po::options_description optionsRate("Rate control options");
    optionsRate.add_options()
                ("qp", po::value<int>()->default_value(26), "quantization parameter")
                ("aq", po::bool_switch()->default_value(false), "TM5-like adaptive quantisation based on psycho-visual model")
                // review: consider replacing aq-depth with aq-size - this is more consistent with other options
                ("aq-depth", po::value<int>()->default_value(3), "maximum depth at which adaptive quantisation can be performed")
                ("aq-range", po::value<int>()->default_value(6), "maximum range at which the qp can vary during adaptive quantisation")
                // review: consider replacing dqp-depth with dqp-size - this is more consistent with other options
                ("dqp-depth", po::value<int>()->default_value(-1), "cu depth where QP can be varied")
                ("bitrate", po::value<int>(), "cbr rate control based on lambda-rate model");
    options.add(optionsRate);

    po::options_description optionsStructure("Structure options");
    optionsStructure.add_options()
        ("shot-change", po::bool_switch()->default_value(false), "enable shot change detection")
        ("field-coding", po::bool_switch()->default_value(false), "enable field coding")
        ("frame-doubling", po::bool_switch()->default_value(false), "enable frame doubling")
        ("max-gop-n", po::value<int>()->default_value(250), "maximum intra picture interval")
        ("max-gop-m", po::value<int>()->default_value(8), "maximum anchor picture interval (1 or 8)")
        ("segment", po::value<int>()->default_value(-1), "Enable IDR segmentation (-1 for Disabled)")
        ("wpp", po::bool_switch()->default_value(true), "enable wave-front parallel processing (default enabled)")
        ("ctu", po::value<int>()->default_value(64), "CTU size in luma pixels")
        ("min-cu", po::value<int>()->default_value(8), "minimum CU size in luma pixels")
        ("repeat-headers", po::bool_switch(), "emit VPS/SPS/PPS on every keyframe (default first keyframe only)");
    options.add(optionsStructure);

    po::options_description optionsTools("Coding tool options");
    optionsTools.add_options()
                ("deblock", po::bool_switch()->default_value(true), "enable deblocking filter")
                ("sao", po::bool_switch(), "enable sample adaptive offset filter")
                ("strong-intra-smoothing", po::bool_switch(), "enable strong intra smoothing")
                ("rqt", po::bool_switch(), "enable one level of rqt (inter coding only)")
                ("amp", po::bool_switch(), "enable asymmetric motion partitions")
                ("smp", po::bool_switch()->default_value(false), "enable symmetric motion partitions")
                ("rdoq", po::bool_switch(), "enable rate distortion optimised quantisation")
                ("max-num-merge-cand", po::value<int>()->default_value(5), "maximum number of merge candidates tested")
                ("sdh", po::bool_switch(), "enable sign data hiding, only with --rdoq")
                ("tskip", po::bool_switch(), "enable transform skip");
    options.add(optionsTools);

    po::options_description optionsPerformance("Performance options");
    optionsPerformance.add_options()
                ("speed", po::value<Speed::Type>()->default_value(Speed::slow), "speed / efficiency tradeoff (one of "
#define X(name) " " #name
                        ENCODER_SPEED_PRESETS_XMACRO
#undef X
                        ")")
                        ("fdm", po::bool_switch(), "fast decision for merge")
                        ("fdam", po::bool_switch(), "fast decision for all modes")
                        ("ecu", po::bool_switch(), "enable early cu termination")
                        ("esd", po::bool_switch(), "enable early skip detection")
                        ("cfm", po::bool_switch(), "enable coding flag mode")
                        ("met", po::bool_switch(), "enable multiple early termination for motion estimation")
                        ("aps", po::bool_switch(), "enable adaptive partition selection")
                        ("rcudepth", po::bool_switch(), "enable rcu-depth algorithm")
                        ("sao-slow-mode", po::bool_switch()->default_value(false), "enable slow sao mode (more accurate)");
    options.add(optionsPerformance);

    po::options_description optionsOptimisation("Optimisation options");
    optionsOptimisation.add_options()
                ("threads", po::value<int>()->default_value(0), "size of thread pool (0=auto detect)")
                ("concurrent-frames", po::value<int>()->default_value(4), "maximum number of pictures that may encode in parallel")
                ("no-parallel-processing", po::bool_switch()->default_value(false), "disable parallelisation (wpp, multithreading and encoding of concurrent frames)")
                ("asm", po::value<int>()->default_value(1), "enable assembly language optimisations");
    options.add(optionsOptimisation);

    po::options_description optionsReporting("Reporting");
    optionsReporting.add_options()
                ("psnr", "measure PSNR and report")
                ("profiler", "profile CPU usage and report")
                ("verbosity", po::value<int>()->default_value(1), "output during encoding (0, 1 or 2)")
                ("help,h", "display help message");
    options.add(optionsReporting);

    po::options_description optionsVui("Video Usability Information (VUI)");
    optionsVui.add_options()
                ("sar", po::value<std::string>(), "sample aspect ratio <w:h>")
                ("display-window", po::value<std::string>(), "specify display window <left,right,top,bottom>")
                ("overscan", po::value<std::string>(), "indicate whether is appropriate to display the overscan area <show|crop>")
                ("video-format", po::value<std::string>(), "indicate the video format <integer|string>")
                ("range", po::value<std::string>(), "indicate the video full range <full|limited>")
                ("colourprim", po::value<std::string>(), "indicate the colour primaries <interger|string>")
                ("transfer-characteristics", po::value<std::string>(), "indicate the transfer function associated with the source material <integer>")
                ("colour-matrix", po::value<std::string>(), "indicate the colour matrix used to derive luma and chroma")
                ("chroma-loc", po::value<std::string>(), "indicate the location for chroma samples <0..5>");
    options.add(optionsVui);

    po::options_description hidden("Hidden options");
    hidden.add_options()
                ("no-rqt", po::bool_switch(), "disable one level of rqt")
                ("no-strong-intra-smoothing", po::bool_switch(), "disable strong intra smoothing")
                ("no-wpp", po::bool_switch(), "disable wpp")
                ("no-deblock", po::bool_switch(), "disable deblocking")
                ("no-sao", po::bool_switch(), "disable sao")
                ("no-rect", po::bool_switch(), "disable rectangular partitions for inter coding")
                ("no-amp", po::bool_switch(), "disable asymmetric motion partition")
                ("no-smp", po::bool_switch()->default_value(false), "disable symmetric motion partition")
                ("no-fdm", po::bool_switch(), "disabled fast decision for merge")
                ("no-fdam", po::bool_switch(), "disabled fast decision for all modes")
                ("no-ecu", po::bool_switch(), "disabled early cu termination")
                ("no-esd", po::bool_switch(), "disable early skip detection")
                ("no-cfm", po::bool_switch(), "disable coding flag mode")
                ("no-met", po::bool_switch(), "disable multiple early termination")
                ("no-sao-slow-mode", po::bool_switch(), "disable slow sao mode")
                ("no-rdoq", po::bool_switch(), "disable rate distortion optimised quantisation")
                ("no-rcudepth", po::bool_switch(), "disable rcu-depth algorithm")
                ("no-sdh", po::bool_switch(), "disable sign data hiding")
                ("no-tskip", po::bool_switch(), "disable transform skip")
                ("no-aps", po::bool_switch(), "disable adaptive partition selection")
                ("force-16", po::bool_switch(), "force usage of uint16_t sample types internally for 8-bit coding")
                ("internal-bit-depth", po::value<int>()->default_value(8), "internal bit depth")
                ("input-file", po::value<std::string>(&inputFilename), "input file name");

    po::options_description all;
    all.add(options).add(hidden);

    po::positional_options_description positional;
    positional.add("input-file", 1);

    try
    {
        po::store(po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);

        if (vm.count("help"))
        {
            cout << "Usage: " << argv[0] << " [options] input-file\n";
            cout << options << "\n";
            std::exit(EXIT_SUCCESS);
        }

        po::notify(vm);

        if (vm.count("input-file") != 1)
        {
            throw std::runtime_error("no input file specified");
        }
    }
    catch (std::exception & e)
    {
        cerr << argv[0] << ": unable to parse command line - " << e.what() << "\n";
        return 1;
    }

    return 0;
}


#define xstr(s) str(s)
#define str(s) #s


const char *gitDescribe()
{
    const char *s = xstr(GIT_DESCRIBE);
    return (s[0] >= '0' && s[0] <= '9') ? s : "<unknown>";
}

int turing_check_binary_option(const char *option)
{
    map<string, bool> supportedBinaryOptions; // option name and default value

    // review: duplication here, would be better to enumerate these options where they are originally defined
    // review: default values are not used
    supportedBinaryOptions["aq"] = false;
    supportedBinaryOptions["shot-change"] = false;
    supportedBinaryOptions["field-coding"] = false;
    supportedBinaryOptions["frame-doubling"] = false;
    supportedBinaryOptions["wpp"] = true;
    supportedBinaryOptions["repeat-headers"] = true;
    supportedBinaryOptions["deblock"] = true;
    supportedBinaryOptions["sao"] = false;
    supportedBinaryOptions["strong-intra-smoothing"] = true;
    supportedBinaryOptions["rqt"] = true;
    supportedBinaryOptions["amp"] = true;
    supportedBinaryOptions["smp"] = false;
    supportedBinaryOptions["rdoq"] = true;
    supportedBinaryOptions["sdh"] = true;
    supportedBinaryOptions["tskip"] = true;
    supportedBinaryOptions["fdm"] = true;
    supportedBinaryOptions["fdam"] = true;
    supportedBinaryOptions["ecu"] = true;
    supportedBinaryOptions["esd"] = true;
    supportedBinaryOptions["cfm"] = true;
    supportedBinaryOptions["met"] = true;
    supportedBinaryOptions["aps"] = true;
    supportedBinaryOptions["rcudepth"] = true;
    supportedBinaryOptions["sao-slow-mode"] = false;
    supportedBinaryOptions["no-parallel-processing"] = true;

    string currentOption(option);
    int isBinaryOption = 0;

    if(currentOption == "no-parallel-processing")
    {
        isBinaryOption = 1;
    }
    else
    {
        // Strip out any no- prefix for disabling-like options
        auto const idx = currentOption.find("no-", 0);
        if(idx != string::npos)
        {
            auto const length = currentOption.length() - 3;
            currentOption = currentOption.substr(idx+3, length);
        }

        auto optionPresent = supportedBinaryOptions.find(currentOption);
        if(optionPresent != supportedBinaryOptions.end())
        {
            isBinaryOption = 1;
        }
    }

    return isBinaryOption;
}


// Review - no need for Encoder and turing_encoder - merge these two?
struct turing_encoder
{
    turing_encoder(turing_encoder_settings const &settings)
    {
        int rv = parseEncodeOptions(this->vm, settings.argc, settings.argv, std::cout, std::cerr, this->inputFilename, this->outputFilename);

        this->verbosity = this->vm["verbosity"].as<int>();

        if (rv) throw std::runtime_error("parseEncodeOptions failed"); // review, error handling

        this->encoder.reset(new Encoder(this->vm));
    }

    size_t bytesPerSample() const
    {
        return 1 + (this->vm.at("bit-depth").as<int>() > 8 || this->vm.at("internal-bit-depth").as<int>() > 8);
    }

    size_t bytesPerInputSample() const
    {
        return 1 + (this->vm.at("bit-depth").as<int>() > 8);
    }

    size_t bytesPerInputPicture() const
    {
        return this->bytesPerInputSample() * this->encoder->pictureWidth * this->encoder->pictureHeight * 3 / 2;
    }

    size_t bytesPerInputFrame() const
    {
        return this->bytesPerInputSample() * this->encoder->frameWidth * this->encoder->frameHeight * 3 / 2;
    }

    turing_bitstream* headers()
    {
        this->bitstream.clear();

        this->encoder->headers(this->bitstream);

        this->output = {
                &this->bitstream.front() ,
                static_cast<int>(this->bitstream.size()) };

        return &this->output.bitstream;
    }

    static int line(int height, int y, bool fieldCoding, bool bottomField)
    {
        if (fieldCoding)
            y = (y << 1) + bottomField;
        if (y >= height)
            y = height - 1;
        return y;
    }

    turing_encoder_output* encode(turing_picture *picture)
    {
        this->bitstream.clear();
        this->output = { 0 };
        if (picture)
        {
            // note that encoder API is picture based but for interlaced coding it is being used to pass frames
            // review: move frame-to-fields conversion above the API?
            auto frame = picture->image;

            auto const fieldCoding = this->vm["field-coding"].as<bool>();

            bool b = false;

            Encoder::PictureMetadata metadata;

            for (int field = 0; field < (fieldCoding ? 2 : 1); ++field)
            {
                std::shared_ptr<PictureWrapper> pictureWrapper;

                if (this->vm.at("bit-depth").as<int>() > 8 || this->vm.at("internal-bit-depth").as<int>() > 8)
                {
                    auto pictureWrap = std::make_shared<PictureWrap<uint16_t>>(this->encoder->pictureWidth, this->encoder->pictureHeight, 1, 0, 0, 32);
                    pictureWrap->sampleSize = 16;
                    pictureWrap->fieldTB = (fieldCoding) ? (field+1) : 0;

                    if (this->vm.at("bit-depth").as<int>() == 8)
                    {
                        for (int cIdx = 0; cIdx < 3; ++cIdx)
                        {
                            for (int y = 0; y < ((this->encoder->frameHeight >> (fieldCoding ? 1 : 0)) >> (cIdx ? 1 : 0)); ++y)
                            {
                                auto p = frame[cIdx].p + line(this->encoder->frameHeight >> (cIdx ? 1 : 0), y, fieldCoding, !field) * frame[cIdx].stride;
                                for (int x = 0; x < (*pictureWrap)[cIdx].width; ++x)
                                    (*pictureWrap)[cIdx](x, y) = p[x] << 2;
                            }
                            if (fieldCoding)
                            {
                                for (int y = ((this->encoder->frameHeight >> (fieldCoding ? 1 : 0)) >> (cIdx ? 1 : 0)); y < (*pictureWrap)[cIdx].height; ++y)
                                {
                                    memset(&(*pictureWrap)[cIdx](0, y), 0, (*pictureWrap)[cIdx].width * sizeof(uint16_t));
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int cIdx = 0; cIdx < 3; ++cIdx)
                        {
                            for (int y = 0; y < ((this->encoder->frameHeight >> (fieldCoding ? 1 : 0)) >> (cIdx ? 1 : 0)); ++y)
                            {
                                auto p = frame[cIdx].p + line(this->encoder->frameHeight >> (cIdx ? 1 : 0), y, fieldCoding, !field) * frame[cIdx].stride;
                                memcpy(&(*pictureWrap)[cIdx](0, y), p, (*pictureWrap)[cIdx].width * sizeof(uint16_t));
                            }
                            if (fieldCoding)
                            {
                                for (int y = ((this->encoder->frameHeight >> (fieldCoding ? 1 : 0)) >> (cIdx ? 1 : 0)); y < (*pictureWrap)[cIdx].height; ++y)
                                {
                                    memset(&(*pictureWrap)[cIdx](0, y), 0, (*pictureWrap)[cIdx].width * sizeof(uint16_t));
                                }
                            }
                        }
                    }
                    pictureWrapper = pictureWrap;
                }
                else
                {
                    auto pictureWrap = std::make_shared<PictureWrap<uint8_t>>(this->encoder->pictureWidth, this->encoder->pictureHeight, 1, 0, 0, 32);
                    pictureWrap->sampleSize = 8;
                    pictureWrap->fieldTB = (fieldCoding) ? (field + 1) : 0;

                    for (int cIdx = 0; cIdx < 3; ++cIdx)
                    {
                        for (int y = 0; y < ((this->encoder->frameHeight>>(fieldCoding ? 1 : 0)) >> (cIdx ? 1 : 0)); ++y)
                        {
                            auto p = frame[cIdx].p + line(this->encoder->frameHeight >> (cIdx ? 1 : 0), y, fieldCoding, !field) * frame[cIdx].stride;
                            memcpy(&(*pictureWrap)[cIdx](0, y), p, (*pictureWrap)[cIdx].width * sizeof(uint8_t));
                        }
                        if(fieldCoding)
                        {
                            for (int y = ((this->encoder->frameHeight >> (fieldCoding ? 1 : 0)) >> (cIdx ? 1 : 0)); y < (*pictureWrap)[cIdx].height; ++y)
                            {
                                memset(&(*pictureWrap)[cIdx](0, y), 0, (*pictureWrap)[cIdx].width * sizeof(uint8_t));
                            }
                        }
                    }

                    pictureWrapper = pictureWrap;
                }

                pictureWrapper->pts = picture->pts;

                std::vector<uint8_t> temp;
                b |= this->encoder->encodePicture(pictureWrapper, temp, metadata);

                this->bitstream.insert(this->bitstream.end(), temp.begin(), temp.end());
            }

            if (b)
            {
                this->output.bitstream.p = &this->bitstream.front();
                this->output.bitstream.size = static_cast<int>(this->bitstream.size());
                this->output.pts = metadata.pts;
                this->output.dts = metadata.dts;
                this->output.keyframe = metadata.keyframe;
            }
        }
        else
        {
            while (true)
            {
                Encoder::PictureMetadata metadata;
                bool b = this->encoder->encodePicture(nullptr, this->bitstream, metadata);
                if (b)
                {
                    if (this->bitstream.empty())
                        continue;

                    this->output.bitstream.p = &this->bitstream.front();
                    this->output.bitstream.size = static_cast<int>(this->bitstream.size());
                    this->output.pts = metadata.pts;
                    this->output.dts = metadata.dts;
                    this->output.keyframe = metadata.keyframe;
                    break;
                }
                else
                    break;
            }
        }

        return &this->output;
    }

    std::unique_ptr<Encoder> encoder;
    po::variables_map vm;
    std::string inputFilename, outputFilename;
    std::vector<uint8_t> bitstream;
    turing_encoder_output output;
    int verbosity;
};


const char *turing_version()
{
    return "1.1";
}


turing_encoder *turing_create_encoder(turing_encoder_settings settings)
{
    turing_encoder *encoder = 0;

    try
    {
        encoder = new turing_encoder(settings);
        if (encoder->verbosity)
        {
            std::cout << "[turing] Project Turing HEVC Encoder Version (" << turing_version() << ") - hash commit: " << gitDescribe() << "\n";
            std::cout << "[turing] patents pending / copyright 2016 Project Turing contributors\n";
        }
    }
    catch (std::exception &e)
    {
        std::cerr << "[turing] " << e.what() << "\n";
    }

    return encoder;
}


turing_bitstream const* turing_encode_headers(turing_encoder *encoder)
{
    try
    {
        return encoder->headers();
    }
    catch (std::exception &)
    {
    }

    static turing_bitstream const error{ 0, -1 };
    return &error;
}


turing_encoder_output const* turing_encode_picture(turing_encoder *encoder, turing_picture *picture)
{
    try
    {
        return encoder->encode(picture);
    }
    catch (std::exception &)
    {
    }

    static turing_encoder_output const error{ 0, -1 };
    return &error;
}


void turing_destroy_encoder(turing_encoder *encoder)
{
    if (encoder)
        delete encoder;
}


std::ostream &operator<<(std::ostream &o, turing_bitstream const &bitstream)
{
    o.write(reinterpret_cast<char *>(bitstream.p), bitstream.size);
    return o;
}


int encode(int argc, const char* const argv[])
{
    auto* encoder = turing_create_encoder(turing_encoder_settings{ argc, argv });

    if (!encoder)
    {
        std::cerr << "[turing] failed to create encoder\n";
        return -1;
    }

    ProgressReporter progressReporter(&std::cout, "encoded ");

    try
    {
        std::ifstream ifs(encoder->inputFilename.c_str(), std::ios_base::binary);
        if (!ifs)
        {
            string errorMessage = "could not open input video file: ";
            errorMessage += encoder->inputFilename.c_str();
            throw std::runtime_error(errorMessage.c_str());
        }
        ifs.seekg(0, std::ios_base::end);

        auto const nInputFrames = static_cast<size_t>(ifs.tellg()) / encoder->bytesPerInputFrame();
        auto const firstFrame = encoder->vm.count("seek") ? encoder->vm["seek"].as<size_t>() : 0;
        auto const nFrames = encoder->vm.count("frames") ? encoder->vm["frames"].as<size_t>() : nInputFrames - firstFrame;

        if (firstFrame >= nInputFrames)
            throw std::runtime_error("cannot seek beyond end of input video file");
        if (firstFrame + nFrames > nInputFrames)
            throw std::runtime_error("encode would read beyond end of input video file");

        ifs.seekg(firstFrame * encoder->bytesPerInputFrame(), std::ios_base::beg);

        std::ofstream ofs(encoder->outputFilename.c_str(), std::ios_base::binary);
        if (!ofs)
            throw std::runtime_error("could not open output bitstream file");

        const bool progress = encoder->vm["verbosity"].as<int>() == 1;

        int nPicturesOutput = 0;

        {
            if (encoder->vm["verbosity"].as<int>())
                encoder->encoder->printHeader(std::cout, encoder->inputFilename, encoder->outputFilename);

            for (int i = 0; i < nFrames; ++i)
            {
                std::vector<uint8_t> buffer;
                buffer.resize(encoder->bytesPerInputFrame());
                ifs.read(reinterpret_cast<char *>(&buffer.front()), encoder->bytesPerInputFrame());

                turing_picture picture;
                picture.image[0].stride = encoder->encoder->frameWidth * encoder->bytesPerInputSample();
                picture.image[1].stride = encoder->encoder->frameWidth * encoder->bytesPerInputSample() / 2;
                picture.image[2].stride = encoder->encoder->frameWidth * encoder->bytesPerInputSample() / 2;
                picture.image[0].p = &buffer.front();
                picture.image[1].p = picture.image[0].p + encoder->bytesPerInputFrame() / 6 * 4;
                picture.image[2].p = picture.image[1].p + encoder->bytesPerInputFrame() / 6;
                picture.pts = i;

                auto *output = turing_encode_picture(encoder, &picture);
                if (output->bitstream.size)
                {
                    ofs << output->bitstream;
                    if(progress)
                    {
                        progressReporter.progress(++nPicturesOutput);
                    }
                }
            }

            while (true)
            {
                auto *output = turing_encode_picture(encoder, 0);
                if (output->bitstream.size)
                {
                    ofs << output->bitstream;
                    if(progress)
                    {
                        progressReporter.progress(++nPicturesOutput);
                    }
                }
                else
                    break;
            }

            if (encoder->vm["verbosity"].as<int>())
                encoder->encoder->printFooter(std::cout);
        }
    }
    catch (std::exception &e)
    {
        std::cerr << argv[0] << ": failed - " << e.what() << std::endl;
        return 1;
    }

    turing_destroy_encoder(encoder);

    return 0;
}
