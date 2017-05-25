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

#include "TaskEncodeInput.h"
#include "TaskEncodeOutput.h"
#include "TaskEncodeSubstream.h"
#include "TaskDeblock.h"
#include "TaskSao.h"
#include "Encoder.h"
#include "Write.h"
#include "Syntax.h"
#include "HevcTypes.h"
#include "InputQueue.h"
#include "ThreadPool.h"
#include "Levels.h"
#include "SyntaxNal.hpp"
#include "SyntaxRbsp.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/chrono.hpp>
#include <mutex>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include <stdio.h>


void Encoder::parseInputRes()
{
    auto resolution = this->vm.at("input-res").as<std::string>();
    auto x = resolution.find('x');
    if (x == std::string::npos)
        throw std::runtime_error("malformed --input-res option");
    resolution[x] = ' ';
    std::istringstream ss{ resolution };
    ss >> this->frameWidth >> this->frameHeight;
}


Encoder::Encoder(boost::program_options::variables_map &vm) :
    frameRate(this->vm["frame-rate"].as<double>()),
    bitDepth((this->vm["internal-bit-depth"].as<int>() > this->vm["bit-depth"].as<int>()) ? this->vm["internal-bit-depth"].as<int>() : this->vm["bit-depth"].as<int>()),
    externalBitDepth(this->vm["bit-depth"].as<int>()),
    vm(vm),
    stateEncode(vm)
{
    this->parseInputRes();

    Speed &speed = this->stateEncode;
    speed = this->vm["speed"].as<Speed::Type>();

    this->stateEncode.internalbitdepth = this->bitDepth;
    this->stateEncode.externalbitdepth = this->externalBitDepth;

    this->stateEncode.scd = this->vm["shot-change"].as<bool>();
    this->stateEncode.fieldcoding = this->vm["field-coding"].as<bool>();
    this->stateEncode.framedoubling = this->vm["frame-doubling"].as<bool>();
    this->stateEncode.amp = this->booleanSwitchSetting("amp", speed.useAmp());
    this->stateEncode.sao = this->booleanSwitchSetting("sao", speed.useSao());
    this->stateEncode.smp = this->vm["smp"].as<bool>();    //this->booleanSwitchSetting("smp", speed.useSmp());
    this->stateEncode.nosmp = this->vm["no-smp"].as<bool>();
    this->stateEncode.fdm = this->booleanSwitchSetting("fdm", speed.useFdm());
    this->stateEncode.fdam = this->booleanSwitchSetting("fdam", speed.useFdam());
    this->stateEncode.ecu = this->booleanSwitchSetting("ecu", speed.useEcu());
    this->stateEncode.esd = this->booleanSwitchSetting("esd", speed.useEsd());
    this->stateEncode.cfm = this->booleanSwitchSetting("cfm", speed.useCfm());
    this->stateEncode.met = this->booleanSwitchSetting("met", speed.useMet());
    this->stateEncode.rcudepth = this->booleanSwitchSetting("rcudepth", speed.useRcuDepth());
    this->stateEncode.rqt = this->booleanSwitchSetting("rqt", speed.useRqt());
    this->stateEncode.rdoq = this->booleanSwitchSetting("rdoq", speed.useRdoq());
    this->stateEncode.sdh = this->booleanSwitchSetting("sdh", speed.useSdh()) && this->stateEncode.rdoq;
    this->stateEncode.tskip = this->booleanSwitchSetting("tskip", speed.useTSkip());
    this->stateEncode.aps = this->booleanSwitchSetting("aps", speed.useAps());
    this->stateEncode.saoslow = this->booleanSwitchSetting("sao-slow-mode", speed.useSaoSlow());
    this->stateEncode.verbosity = this->vm["verbosity"].as<int>();
    this->stateEncode.gopM = std::min(this->vm["max-gop-m"].as<int>(), this->vm["max-gop-n"].as<int>());
    this->stateEncode.baseQp = this->vm["qp"].as<int>();
    this->stateEncode.maxnummergecand = (this->vm["max-num-merge-cand"].defaulted()) ? (speed.setMaxNumMergeCand()) : this->vm["max-num-merge-cand"].as<int>();
    this->stateEncode.preferredTransferCharacteristics = this->vm["atc-sei"].as<int>();

    this->stateEncode.wpp = this->booleanSwitchSetting("wpp", true);
    this->stateEncode.concurrentFrames = this->vm["concurrent-frames"].as<int>();
    if (this->vm["no-parallel-processing"].as<bool>())
    {
        this->stateEncode.wpp = false;
        this->stateEncode.concurrentFrames = 1;
    }
    else
    {
        if (this->stateEncode.wpp && vm["threads"].as<int>() == 1 && this->stateEncode.concurrentFrames == 1)
        {
            std::cout << "Warning: --wpp brings no advantage in single-threaded mode.\n" << std::endl;
        }
    }

    if (this->vm.count("dump-frames"))
    {
        this->stateEncode.fileOutYuvFrames.open(vm["dump-frames"].as<std::string>(), std::ios::binary);
        if (!this->stateEncode.fileOutYuvFrames)
        {
            throw std::runtime_error("failed to open reconstructed YUV file for writing");
        }
    }

    if (this->vm.count("dump-pictures"))
    {
        this->stateEncode.fileOutYuvPictures.open(vm["dump-pictures"].as<std::string>(), std::ios::binary);
        if (!this->stateEncode.fileOutYuvPictures)
        {
            throw std::runtime_error("failed to open reconstructed YUV file for writing");
        }
    }

    if (this->vm.count("psnr") || (this->stateEncode.verbosity == 2))
        this->stateEncode.psnrAnalysis.reset(new PsnrAnalysis((this->vm["internal-bit-depth"].as<int>() > this->vm["bit-depth"].as<int>()) ? this->vm["internal-bit-depth"].as<int>() : this->vm["bit-depth"].as<int>()));

    if (this->vm.count("hash"))
    {
        this->stateEncode.decodedHashSei = true;
        if (!(0 <= this->vm["hash"].as<int>() && this->vm["hash"].as<int>() <= 2))
        {
            throw std::runtime_error("Error: decoded hash type should be in the range [0, 2] inclusive\n");
        }
        this->stateEncode.hashType = HashType(this->vm["hash"].as<int>());
    }
    else
    {
        this->stateEncode.decodedHashSei = false;
    }

    int deltaQpDepth = this->vm.at("dqp-depth").as<int>();
    if (deltaQpDepth > -1)
    {
        int log2DiffMaxMin = static_cast<int>(log2Size(this->vm["ctu"].as<int>(), 8, 64, "invalid ctu value") - log2Size(this->vm["min-cu"].as<int>(), 8, 64, "invalid min-cu value"));
        std::string errMsg("Error the depth at which QP can be varied should be in the range [0, ");
        std::ostringstream numberToString;
        numberToString << log2DiffMaxMin;
        errMsg += numberToString.str();
        errMsg += "] inclusive\n";
        if (!(0 <= deltaQpDepth && deltaQpDepth <= log2DiffMaxMin))
            throw std::runtime_error(errMsg.c_str());
    }

    this->stateEncode.useAq = this->vm.at("aq").as<bool>();
    if (this->stateEncode.useAq)
    {
        int aqDepth = this->vm["aq-depth"].as<int>();
        int aqRange = this->vm["aq-range"].as<int>();
        int maxCuSize = this->vm["ctu"].as<int>();
        int minCuSize = this->vm["min-cu"].as<int>();
        const auto log2DiffMaxMin = static_cast<int>(log2Size(maxCuSize, 8, 64, "invalid ctu value") - log2Size(minCuSize, 8, 64, "invalid min-cu value"));

        std::string errMsg("Error the depth for adaptive quantisation should be in the range [0, ");
        std::ostringstream numberToString;
        numberToString << log2DiffMaxMin;
        errMsg += numberToString.str();
        errMsg += "] inclusive\n";

        if (!(0 <= aqDepth && aqDepth <= log2DiffMaxMin))
            throw std::runtime_error(errMsg.c_str());
    }

    this->stateEncode.enableProfiler = !!this->vm.count("profiler");

    Handler<Encode<void>, StateEncode> h;
    h.state = &this->stateEncode;

    this->pictureHeight = this->frameHeight;
    this->pictureWidth = this->frameWidth;

    if (this->vm["field-coding"].as<bool>())
    {
        this->pictureHeight /= 2;
        const auto minCuSize = this->vm["min-cu"].as<int>();
        this->pictureHeight += minCuSize - 1;
        this->pictureHeight /= minCuSize;
        this->pictureHeight *= minCuSize;
        // review - allow non-CU-multiple dimensions also for frame coding?
        // review - is it necessary to set the conformance window accordingly?
        if (this->pictureHeight*2 != this->frameHeight)
            std::cout << "Warning: picture height will be padded to " << (this->pictureHeight * 2) << " samples to allow field-coding with minimum cu size of " << minCuSize << "\n" << std::endl;

        this->confWinBottomOffset = (this->pictureHeight * 2 - this->frameHeight) / 2;
    }
    else
        this->confWinBottomOffset = 0;

    this->stateEncode.useRateControl = !!this->vm.count("bitrate");
    if (this->stateEncode.useRateControl)
    {
        int frames = static_cast<int>(this->vm["frames"].as<std::size_t>());
        if(this->vm["field-coding"].as<bool>())
            frames *= 2;

        this->stateEncode.rateControlParams.reset(new RateControlParameters(
                (double)this->vm["bitrate"].as<int>(),

                this->vm["frame-rate"].as<double>(),
                this->vm["max-gop-n"].as<int>(),
                this->vm["max-gop-m"].as<int>(),
                this->pictureHeight,
                this->pictureWidth,
                this->vm["ctu"].as<int>(),
                this->vm["qp"].as<int>(),
                this->stateEncode.wpp,
                static_cast<int>(this->stateEncode.concurrentFrames)));
    }

    this->stateEncode.repeatHeaders = this->vm["repeat-headers"].as<bool>();

    if(this->vm.count("mastering-display-info"))
    {
        this->stateEncode.masteringDisplayInfoPresent = true;
        string masteringDisplayInfo = this->vm["mastering-display-info"].as<string>();

        sscanf(masteringDisplayInfo.c_str(), "G(%hu,%hu)B(%hu,%hu)R(%hu,%hu)WP(%hu,%hu)L(%u,%u)",
                   &this->stateEncode.masterDisplayInfo.displayPrimariesX[0],
                   &this->stateEncode.masterDisplayInfo.displayPrimariesY[0],
                   &this->stateEncode.masterDisplayInfo.displayPrimariesX[1],
                   &this->stateEncode.masterDisplayInfo.displayPrimariesY[1],
                   &this->stateEncode.masterDisplayInfo.displayPrimariesX[2],
                   &this->stateEncode.masterDisplayInfo.displayPrimariesY[2],
                   &this->stateEncode.masterDisplayInfo.whitePointX,
                   &this->stateEncode.masterDisplayInfo.whitePointY,
                   &this->stateEncode.masterDisplayInfo.maxDisplayMasteringLuma,
                   &this->stateEncode.masterDisplayInfo.minDisplayMasteringLuma);
    }
    else
    {
        this->stateEncode.masteringDisplayInfoPresent = false;
    }

    this->setupPps(h);
    ProfileTierLevel *ptl = this->setupSps(h);
    this->setupVps(h, ptl);

    {
        // Kick off  task. This will block until the first output frame is encoded.
        ThreadPool *threadPool = &this->stateEncode;
        threadPool->add(newTaskEncodeInput(h));
        threadPool->add(newTaskEncodeOutput(h));
    }
}


// returns parsed value of boolean settings of form "--no-thing" / "--thing"
bool Encoder::booleanSwitchSetting(std::string name, bool defaultValue)
{
    if (this->vm["no-" + name].as<bool>()) return false;
    if (this->vm[name].as<bool>()) return true;
    return defaultValue;
}


void Encoder::printHeader(std::ostream &cout, std::string const &inputFile, std::string const &outputFile)
{
    bool const wpp = this->stateEncode.wpp;

    std::string reconFile =  this->vm.count("dump-frames") ? this->vm["dump-frames"].as<std::string>() :
            (this->vm.count("dump-pictures") ? this->vm["dump-pictures"].as<std::string>() : "none");
    auto const startFrame = this->vm.count("seek") ? this->vm.at("seek").as<size_t>() : 0;
    auto const nFrames = (this->stateEncode.fieldcoding ? (this->vm.at("frames").as<size_t>() << 1) : this->vm.at("frames").as<size_t>());
    auto const stopFrame = startFrame + nFrames - 1;
    cout << "Input          file             : " << inputFile << "\n";
    cout << "Bitstream      file             : " << outputFile << "\n";
    cout << "Reconstruction file             : " << reconFile << "\n";
    cout << "Input          resolution       : "
            << this->frameWidth << "x" << this->frameHeight << " "
            << this->vm.at("frame-rate").as<double>() << "Hz "
            << this->vm.at("bit-depth").as<int>() << "-bit";
    if (this->vm.at("internal-bit-depth").as<int>() == 10)
        cout << " (16-bit datapath)";
    cout << "\n";
    cout << "Frame          range            : " << startFrame << " - " << stopFrame << " ( " << stopFrame - startFrame + 1 << (this->stateEncode.fieldcoding ? " fields )\n" : " frames )\n");
    cout << "Frame/Field    coding           : " << (this->stateEncode.fieldcoding ? "Field based coding\n" : "Frame based coding\n");
    cout << "Coding unit size (min/max)      : " << this->vm.at("ctu").as<int>() << " / " << this->vm.at("min-cu").as<int>() << "\n";
    cout << "Intra period                    : " << this->vm.at("max-gop-n").as<int>() << "\n";
    if (this->vm.at("segment").as<int>() != -1)
        cout << "IDR Segment length              : " << this->vm.at("segment").as<int>() << "\n";
    if (this->stateEncode.scd)
        cout << "Shot change detection enabled. Warning: up to 48 frames will be kept in the buffer.\n";
    cout << "Structure Of Picture (SOP) size : " << this->vm.at("max-gop-m").as<int>() << "\n";
    cout << "Quantisation Parameter (QP)     : " << this->vm.at("qp").as<int>() << "\n";
    if (this->stateEncode.useRateControl)
        cout << "Target rate (CBR)               : " << this->vm.at("bitrate").as<int>() << " [kbit/s]\n";
    if (this->vm.at("dqp-depth").as<int>() >= 0)
        cout << "Maximum depth to select QP      : " << this->vm.at("dqp-depth").as<int>() << "\n";
    if (this->vm.at("aq").as<bool>())
        cout << "AQ enabled (depth/range)        : " << this->vm.at("aq-depth").as<int>() << " / " << this->vm.at("aq-range").as<int>() << "\n";
    cout << "Speed                           : " << Speed::Type(this->stateEncode) << "\n";
    cout << "Threads                         : " << this->stateEncode.ThreadPool::size() << "\n";
    cout << "Concurrent frames               : " << this->stateEncode.concurrentFrames << "\n";
    cout << "Optimisations                   :";
#define X(A, B, C) if (this->stateEncode.StateFunctionTables::instruction_set_support & (1 << A)) cout << " "<< #B;
    HAVOC_INSTRUCTION_SET_XMACRO
#undef X
    cout << "\n";
    cout << "Coding tools                    : ecu=" << this->stateEncode.ecu;
    cout << " fdm=" << this->stateEncode.fdm;
    cout << " fdam=" << this->stateEncode.fdam;
    cout << " cfm=" << this->stateEncode.cfm;
    cout << " esd=" << this->stateEncode.esd;
    cout << " met=" << this->stateEncode.met;
    cout << " rqt=" << this->stateEncode.rqt;
    cout << " rdoq=" << this->stateEncode.rdoq;
    cout << " amp=" << this->stateEncode.amp;

    cout << " smp=" << (!!this->stateEncode.smp ? "1" : (!!this->stateEncode.nosmp ? "0" : (!!static_cast<Speed &>(this->stateEncode).useSmp(6) ? "1" : "restricted")));
    cout << " deblock=" << !!(this->vm.at("deblock").as<bool>() && !this->vm.at("no-deblock").as<bool>());
    cout << " sao=" << this->stateEncode.sao;
    cout << " sao-slow-mode=" << this->stateEncode.saoslow;
    cout << " rcudepth=" << this->stateEncode.rcudepth;
    cout << " wpp=" << wpp;
    cout << " sdh=" << this->stateEncode.sdh;
    cout << " tskip=" << this->stateEncode.tskip;
    cout << " aps=" << this->stateEncode.aps;
    cout << "\n\n";
    cout.flush();

    this->cpuTimer.start();
}


void Encoder::printFooter(std::ostream &cout)
{
    cout << "\n";
    double frameRate = vm.at("frame-rate").as<double>();
    double bitrate = (double)(this->byteCount * 8) / (double)this->frameCount * frameRate / 1000.0;
    cout << "Encoded " << this->frameCount << (this->frameCount == 1 ? " picture" : " pictures") << " to " << bitrate << " (kbps)\n";
    cout << "Bytes written to file: " << this->byteCount << "\n";

    if (this->stateEncode.psnrAnalysis)
    {
        cout << "\n";
        this->stateEncode.psnrAnalysis->report(cout);
    }

    cout << "\n";
    this->cpuTimer.stop();
    cout << this->cpuTimer.format() << "\n";

    if (this->stateEncode.enableProfiler)
    {
#ifdef WIN32 // rdtsc seems unreliable (sensitive to context switches?) under Linux whereas under Windows (Hypervisor Win8) consistent results are obtained
        cout << "\n";

        Profiler::Timers *timers = &this->stateEncode;
        timers->report(cout);
#endif
    }
}


struct StateHeaders :
    NalWriter,
    StateSlice,
    StateEncodePictureSei,
    Strps
{
};


void Encoder::headers(std::vector<uint8_t> &bitstream)
{
    assert(bitstream.empty());

    Handler<Write<void>, StateEncode> h1;
    h1.state = &this->stateEncode;

    StateHeaders stateHeaders;

    auto h = h1.extend(&stateHeaders);

    // Activate parameter sets
    h[Active<Vps>()] = h[Table<Vps>()][0];
    h[Active<Sps>()] = h[Table<Sps>()][0];
    h[Active<Pps>()] = h[Table<Pps>()][0];

    writeHeaders(h);

    bitstream = *stateHeaders.data;
}


static std::string hashToString(std::vector<int> hashElement, int numberOfChar)
{
    const char* hex = "0123456789abcdef";
    std::string result;

    for (int pos = 0; pos<int(hashElement.size()); pos++)
    {
        if ((pos % numberOfChar) == 0 && pos != 0)
        {
            result += ',';
        }
        result += hex[hashElement[pos] >> 4];
        result += hex[hashElement[pos] & 0xf];
    }

    return result;
}


bool Encoder::encodePicture(std::shared_ptr<PictureWrapper> picture, std::vector<uint8_t> &bitstream, Encoder::PictureMetadata &metadata)
{
    //if (picture)
    //	std::cout << "  Encoder::encodePicture " << picture->pts << "\n";

    assert(bitstream.empty());

    ThreadPool *threadPool = &this->stateEncode;
    InputQueue *inputQueue = &this->stateEncode;

    if (this->stateEncode.verbosity == 2) this->frameCpuTimer.start();

    const bool eos = !picture;

    {
        threadPool->lock();
        if (eos)
        {
            inputQueue->endOfInput();
        }
        else
        {
            std::shared_ptr<AdaptiveQuantisation> aqInfo;
            if (this->stateEncode.useAq)
            {
                aqInfo.reset(new AdaptiveQuantisation(this->vm["aq-depth"].as<int>(), this->vm["aq-range"].as<int>(), this->pictureHeight, this->pictureWidth, this->vm["ctu"].as<int>()));
                if (picture->sampleSize == 8)
                    aqInfo->preAnalysis<uint8_t>(picture);
                else
                    aqInfo->preAnalysis<uint16_t>(picture);
            }
            inputQueue->append(picture, aqInfo);
        }
        threadPool->nudge();
        threadPool->unlock();
    }

    {
        StateEncode::Response response;

        {
            std::unique_lock<std::mutex> lock(threadPool->mutex());
            while (stateEncode.responses.empty() || !stateEncode.responses.front().done)
            {
                stateEncode.responsesAvailable.wait(lock);
            }
            response = stateEncode.responses.front();
            stateEncode.responses.pop_front();

            threadPool->nudge();
        }

        assert(response.done);

        if (response.hungry && !eos)
        {
            //std::cout << "  } Encoder::encodePicture false (hungry && !eos)\n";
            return false;
        }

        if (response.eos)
        {
            stateEncode.decodedPicturesFlush((this->bitDepth > 8) ? 16 : 8, (this->externalBitDepth > 8) ? 16 : 8);
            return false;
        }

        if (response.picture)
        {
            NalWriter *nalWriter = &*response.picture;

            auto const bytes = nalWriter->data->size();

            if (!nalWriter->data->empty())
            {
                bitstream = *nalWriter->data;
                nalWriter->data->clear();
            }
            else
                assert(false);

            metadata.pts = response.picture->docket->picture->pts;
            metadata.dts = response.picture->docket->dts;
            metadata.keyframe = response.keyframe;

            //std::cout << "response keyframe/DTS/PTS " << metadata.keyframe << " " << metadata.dts << " " << metadata.pts << "\n";
            std::unique_lock<std::mutex> lock(threadPool->mutex());
            while (!response.picture->reconstructed)
            {
                stateEncode.responsesAvailable.wait_for(lock, std::chrono::milliseconds(10));
            }
            stateEncode.decodedPicturesFlush((this->bitDepth > 8) ? 16 : 8, (this->externalBitDepth > 8) ? 16 : 8);

            if (this->stateEncode.verbosity == 2)
            {
                this->frameCpuTimer.stop();
                auto &h = *response.picture;
                int qp = h[init_qp_minus26()] + 26 + h[slice_qp_delta()];
                boost::chrono::duration<double> frameEncoderTimeSec = boost::chrono::nanoseconds(this->frameCpuTimer.elapsed().user);
                int sliceType = h[slice_type()];

                std::ostringstream oss;
                
                const int currentPoc = h[PicOrderCntVal()];
                StateEncodePicture* stateEncodePicture = &h;
                const int absolutePoc = stateEncodePicture->docket->absolutePoc;
                double psnrY, psnrU, psnrV;
                stateEncode.psnrAnalysis->getPsnrFrameData(absolutePoc, psnrY, psnrU, psnrV);
                oss << "POC" << std::setw(5) << currentPoc << " ( " << sliceTypeToChar(sliceType) << "-SLICE, QP " << std::setw(4) << qp << " ) ";
                oss << std::setw(12) << 8 * bytes << " bits ";
                oss << std::setprecision(4);
                oss << std::fixed;
                oss << "[ Y" << std::setw(8) << psnrY << " dB  ";
                oss << " U" << std::setw(8) << psnrU << " dB  ";
                oss << " V" << std::setw(8) << psnrV << " dB ]";
                oss << std::setprecision(1);
                oss << " [ ET " << frameEncoderTimeSec.count() << " ] ";
                stateEncode.psnrAnalysis->removePsnrFrameData(absolutePoc);

                if (this->stateEncode.decodedHashSei)
                {
                    StateEncode::FrameHash currentHash = stateEncode.getFrameHash(absolutePoc);
                    std::string digestHashString;
                    switch (this->stateEncode.hashType)
                    {
                        case MD5:
                            digestHashString = "[MD5:" + hashToString(currentHash.hash, 16) + "]";
                            break;
                        case CRC:
                            digestHashString = "[CRC:" + hashToString(currentHash.hash, 2) + "]";
                            break;
                        case CHKSUM:
                            digestHashString = "[CKS:" + hashToString(currentHash.hash, 4) + "]";
                            break;
                    }
                    oss << digestHashString << " ";
                    stateEncode.removeFrameHash(absolutePoc);
                }

                // Implemented print out of reference frames used
                std::cout << oss.str() << "\n";
                std::cout.flush();
            }

            ++this->frameCount;
            this->byteCount += bytes;

            assert(bytes);
        }
    }
    //std::cout << "  } Encoder::encodePicture true\n";
    return true;
}


template <class H>
void Encoder::setupPtl(H &h)
{
    h[general_profile_compatibility_flag(2)] = 1;

    if (this->bitDepth == 8)
    {
        h[general_profile_idc()] = 1;
        h[general_profile_compatibility_flag(1)] = 1;
    }
    else
    {
        h[general_profile_idc()] = 2;
    }

    int a = h[PicSizeInSamplesY()];
    int b = (int)ceil(h[PicSizeInSamplesY()] * this->frameRate);

    h[general_level_idc()] = 0;

    for (auto level : levels)
    {
        if (level.parameters[Level::MaxLumaPs] >= h[PicSizeInSamplesY()] && level.parameters[Level::MaxLumaSr] >= h[PicSizeInSamplesY()] * this->frameRate)
        {
            h[general_level_idc()] = level.level_idc();
            this->bitstreamLevel = level;
            if(this->stateEncode.useRateControl)
                this->stateEncode.rateControlParams->initCpbInfo(static_cast<int>(level.parameters[Level::MaxCPB]) * 1000);
            break;
        }
    }
}


template <class H>
void Encoder::setupVps(H &hhh, ProfileTierLevel *ptl)
{
    std::shared_ptr<Vps> vps(new Vps());
    auto hh = hhh.extend(&*vps);
    auto h = hh.extend(&vps->ptl);
    h[Table<Vps>()][0] = vps;

    vps->ptl = *ptl;

    h[vps_reserved_0xffff_16bits()] = 0xffff;
    h[vps_reserved_three_2bits()] = 3;


    h[vps_sub_layer_ordering_info_present_flag()] = 1;
    h[vps_max_dec_pic_buffering_minus1(0)] = 4;
    h[vps_max_num_reorder_pics(0)] = 3;
    h[vps_temporal_id_nesting_flag()] = 1;
}


template <class H>
ProfileTierLevel *Encoder::setupSps(H &hhh)
{
    std::shared_ptr<Sps> sps(new Sps());
    auto hh = hhh.extend(&*sps);
    auto h = hh.extend(&sps->ptl);
    h[Table<Sps>()][0] = sps;

    {
        // spatial configuration
        h[chroma_format_idc()] = 1;
        h[pic_width_in_luma_samples()] = this->pictureWidth;
        h[pic_height_in_luma_samples()] = this->pictureHeight;

        if (this->confWinBottomOffset)
        {
            h[conformance_window_flag()] = 1;
            h[conf_win_left_offset()] = 0;
            h[conf_win_right_offset()] = 0;
            h[conf_win_top_offset()] = 0;
            h[conf_win_bottom_offset()] = this->confWinBottomOffset / (h[SubHeightC()]);
        }

        h[bit_depth_luma_minus8()] = this->bitDepth - 8;
        h[bit_depth_chroma_minus8()] = this->bitDepth - 8;
        h[sps_temporal_id_nesting_flag()] = 1;
        h[log2_max_pic_order_cnt_lsb_minus4()] = 2;
        const int log2MinCuSize = log2Size(this->vm["min-cu"].as<int>(), 8, 64, "invalid min-cu value");
        h[log2_min_luma_coding_block_size_minus3()] = log2MinCuSize - 3;
        const int log2CtuSize = log2Size(this->vm["ctu"].as<int>(), 16, 64, "invalid ctu value");
        h[log2_diff_max_min_luma_coding_block_size()] = log2CtuSize - log2MinCuSize;
        h[log2_min_luma_transform_block_size_minus2()] = 0;
        int log2MaxTransformBlockSize = std::min(5, h[CtbLog2SizeY()]);
        h[log2_diff_max_min_luma_transform_block_size()] = log2MaxTransformBlockSize - h[MinTbLog2SizeY()];
        h[max_transform_hierarchy_depth_intra()] = 1;
        h[max_transform_hierarchy_depth_inter()] = (this->stateEncode.rqt) ? 1 : 0;
        if (h[pic_width_in_luma_samples()] % (1 << h[MinCbLog2SizeY()])) throw std::runtime_error("picture width is not a multiple of minimum CU size");
        if (h[pic_height_in_luma_samples()] % (1 << h[MinCbLog2SizeY()])) throw std::runtime_error("picture height is not a multiple of minimum CU size");

        // temporal configuration
        h[sps_max_dec_pic_buffering_minus1()] = 4;
        h[sps_max_num_reorder_pics()] = 3; // review - could be 2?
        h[sps_temporal_id_nesting_flag()] = 1;

        // coding configuration
#ifdef FORCE_PCM_PICTURES
        h[pcm_enabled_flag()] = 1;
#endif
        if (h[pcm_enabled_flag()])
        {
            h[log2_min_pcm_luma_coding_block_size_minus3()] = 0;
            auto const log2MaxPcmBlockSize = std::min(5, h[CtbLog2SizeY()]);
            h[log2_diff_max_min_pcm_luma_coding_block_size()] = log2MaxPcmBlockSize - (h[log2_min_pcm_luma_coding_block_size_minus3()] + 3);
            h[pcm_sample_bit_depth_luma_minus1()] = 7;
            h[pcm_sample_bit_depth_chroma_minus1()] = 7;
            h[pcm_loop_filter_disabled_flag()] = 1;
        }
        h[strong_intra_smoothing_enabled_flag()] = this->booleanSwitchSetting("strong-intra-smoothing", true);
        h[amp_enabled_flag()] = this->stateEncode.amp;
        h[sample_adaptive_offset_enabled_flag()] = this->stateEncode.sao;
        h[sps_temporal_mvp_enabled_flag()] = 1;
    }

    setupPtl(h);

    if (writeVui())
        setupVui(h);

    return &sps->ptl;
}


template <class H>
void Encoder::setupPps(H &hh)
{
    std::shared_ptr<Pps> pps(new Pps());
    auto h = hh.extend(&*pps);
    h[Table<Pps>()][0] = pps;

    int const qp = this->vm["qp"].as<int>();
    h[init_qp_minus26()] = qp - 26;

    h[entropy_coding_sync_enabled_flag()] = this->stateEncode.wpp;

    int const deblock = this->booleanSwitchSetting("deblock", true);

    h[deblocking_filter_control_present_flag()] = !deblock;
    h[pps_deblocking_filter_disabled_flag()] = !deblock;
    h[sign_data_hiding_enabled_flag()] = this->stateEncode.sdh;
    h[transform_skip_enabled_flag()] = this->stateEncode.tskip;

    if (this->vm.count("bitrate"))
    {
        h[cu_qp_delta_enabled_flag()] = 1;
        h[diff_cu_qp_delta_depth()] = 0;
    }

    if (this->vm.at("dqp-depth").as<int>() >= 0)
    {
        h[cu_qp_delta_enabled_flag()] = 1;
        h[diff_cu_qp_delta_depth()] = this->vm.at("dqp-depth").as<int>();
    }

    if (this->vm.at("aq").as<bool>())
    {
        h[cu_qp_delta_enabled_flag()] = 1;
        h[diff_cu_qp_delta_depth()] = this->vm.at("aq-depth").as<int>();
    }
}

bool Encoder::writeVui()
{
    // Check whether one of the VUI parameters is set, so that VUI writing can happen
    if (this->vm.count("sar"))                      return true;
    if (this->vm.count("display-window"))           return true;
    if (this->vm.count("overscan"))                 return true;
    if (this->vm.count("video-format"))             return true;
    if (this->vm.count("range"))                    return true;
    if (this->vm.count("colourprim"))               return true;
    if (this->vm.count("transfer-characteristics")) return true;
    if (this->vm.count("colour-matrix"))            return true;
    if (this->vm.count("chroma-loc"))               return true;
    if (this->vm["field-coding"].as<bool>())        return true;
    if (this->vm["frame-doubling"].as<bool>())      return true;
    return false;
}


template <class H>
void Encoder::setupVui(H &h)
{
    using namespace std;
    using namespace boost;

    h[vui_parameters_present_flag()] = 1;
    if (this->vm.count("sar"))
    {
        h[aspect_ratio_info_present_flag()] = 1;

        string sampleRatio = this->vm["sar"].as<string>();

        if (sampleRatio.find(":", 0) == string::npos)
        {
            // sar specified as idc
            int sarIdc = atoi(sampleRatio.c_str());
            h[aspect_ratio_idc()] = sarIdc;
        }
        else
        {
            h[aspect_ratio_idc()] = 255;
            // Tokenise command line to get w and h
            char_separator<char> separator(":");
            tokenizer<char_separator<char>> tokens(sampleRatio, separator);
            tokenizer<char_separator<char>>::iterator it = tokens.begin();
            string widthStr = *it;
            int width = atoi(widthStr.c_str());
            string heightStr = *(++it);
            int height = atoi(heightStr.c_str());
            h[sar_width()] = width;
            h[sar_height()] = height;
        }
    }
    else
        h[aspect_ratio_info_present_flag()] = 0;

    if (this->vm.count("overscan"))
    {
        h[overscan_info_present_flag()] = 1;
        if (this->vm["overscan"].as<string>() == "show")
            h[overscan_appropriate_flag()] = 1;
        else
            h[overscan_appropriate_flag()] = 0;
    }
    else
        h[overscan_info_present_flag()] = 0;

    bool videoSignalTypeInfoFlag = this->vm.count("video-format") ||
            this->vm.count("range") ||
            this->vm.count("colourprim") ||
            this->vm.count("transfer-characteristics") ||
            this->vm.count("colour-matrix");
    if (videoSignalTypeInfoFlag)
    {
        h[video_signal_type_present_flag()] = 1;
        if (this->vm.count("video-format"))
        {
            // Read video format from command line
            int videoFormat = atoi(this->vm["video-format"].as<string>().c_str());
            h[video_format()] = videoFormat;
        }
        else
            h[video_format()] = 5; // Unspecified

        if (this->vm.count("range"))
        {
            if (this->vm["range"].as<string>() == "full")
                h[video_full_range_flag()] = 1;
            else
                h[video_full_range_flag()] = 0;
        }
        else
            h[video_full_range_flag()] = 0;

        bool colourDescriptionPresentFlag = this->vm.count("colourprim") ||
                this->vm.count("transfer-characteristics") ||
                this->vm.count("colour-matrix");
        if (colourDescriptionPresentFlag)
        {
            h[colour_description_present_flag()] = 1;
            if (this->vm.count("colourprim"))
            {
                int colourPrimaries = atoi(this->vm["colourprim"].as<string>().c_str());
                h[colour_primaries()] = colourPrimaries;
            }
            else
                h[colour_primaries()] = 2; // Unspecified

            if (this->vm.count("transfer-characteristics"))
            {
                int transferCharacteristics = atoi(this->vm["transfer-characteristics"].as<string>().c_str());
                h[transfer_characteristics()] = transferCharacteristics;
            }
            else
                h[transfer_characteristics()] = 2; // Unspecified

            if (this->vm.count("colour-matrix"))
            {
                int colourMatrix = atoi(this->vm["colour-matrix"].as<string>().c_str());
                h[matrix_coeffs()] = colourMatrix;
            }
            else
                h[matrix_coeffs()] = 2; // Unspecified
        }
    }
    else
        h[video_signal_type_present_flag()] = 0;

    if (this->vm.count("chroma-loc"))
    {
        h[chroma_loc_info_present_flag()] = 1;
        int chromaLocation = atoi(this->vm["chroma-loc"].as<string>().c_str());
        h[chroma_sample_loc_type_top_field()] = chromaLocation;
        h[chroma_sample_loc_type_bottom_field()] = chromaLocation;
    }
    else
        h[chroma_loc_info_present_flag()] = 0;

    h[neutral_chroma_indication_flag()] = 0;

    if ((vm["field-coding"].as<bool>()))
    {
        h[vui_hrd_parameters_present_flag()] = 0;
        h[field_seq_flag()] = 0;
        h[frame_field_info_present_flag()] = 1;
    }
    else if ((vm["frame-doubling"].as<bool>()))
    {
        h[vui_hrd_parameters_present_flag()] = 0;
        h[field_seq_flag()] = 0;
        h[frame_field_info_present_flag()] = 1;
    }
    else
    {
        h[field_seq_flag()] = 0;
        h[frame_field_info_present_flag()] = 0;
    }

    if (this->vm.count("display-window"))
    {
        h[default_display_window_flag()] = 1;
        char_separator<char> separator(", ");
        tokenizer<char_separator<char>> tokens(this->vm["display-window"].as<string>(), separator);
        tokenizer<char_separator<char>>::iterator it = tokens.begin();
        string leftStr = *it++;
        string rightStr = *it++;
        string topStr = *it++;
        string bottomStr = *it++;
        int left = atoi(leftStr.c_str());
        int right = atoi(rightStr.c_str());
        int top = atoi(topStr.c_str());
        int bottom = atoi(bottomStr.c_str());
        h[def_disp_win_left_offset()] = left;
        h[def_disp_win_right_offset()] = right;
        h[def_disp_win_top_offset()] = top;
        h[def_disp_win_bottom_offset()] = bottom;
    }
    else
        h[default_display_window_flag()] = 0;

    h[vui_timing_info_present_flag()] = 1;
    h[vui_num_units_in_tick()] = 1000;
    h[vui_time_scale()] = static_cast<int>(this->vm["frame-rate"].as<double>() * 1000 + 0.5);
    h[bitstream_restriction_flag()] = 0;
    if(this->stateEncode.useRateControl)
    {
        h[vui_hrd_parameters_present_flag()] = 1;
        Hrd hrdDefault{};
        Hrd::SubLayer sublayerDefault{};
        auto currentSps = h[Table<Sps>()][0];
        auto h2 = h.extend(&hrdDefault);
        auto h3 = h2.extend(&sublayerDefault);
        setupHrd(h3);
        hrdDefault.sublayers.push_back(sublayerDefault);
        currentSps->hrdArray.hrd.push_back(hrdDefault);
    }
}

int computeScale(int value)
{
    uint32_t mask = 0xffffffff;
    int scaleValue = 32;

    while ((value & mask) != 0)
    {
      scaleValue--;
      mask = (mask >> 1);
    }

    return scaleValue;
}

template <class H>
void Encoder::setupHrd(H &h)
{
    int targetRate = this->vm["bitrate"].as<int>()*1000;
    int bitrateScale = Clip3(0, 15, computeScale(targetRate) - 6);
    int cpbSize = static_cast<int>(this->bitstreamLevel.parameters[Level::MaxCPB]*1000);
    int cpbScale = Clip3(0, 15, computeScale(cpbSize) - 4);
    h[nal_hrd_parameters_present_flag()] = 1;
    h[bit_rate_scale()] = bitrateScale;
    h[cpb_size_scale()] = cpbScale;
    h[initial_cpb_removal_delay_length_minus1()] = 15;
    h[au_cpb_removal_delay_length_minus1()] = 5;
    h[dpb_output_delay_du_length_minus1()] = 5;

    h[fixed_pic_rate_general_flag(0)] = 1;
    h[fixed_pic_rate_within_cvs_flag(0)] = 1;
    h[elemental_duration_in_tc_minus1(0)] = 0;

    h[bit_rate_value_minus1(0)] = targetRate >> (bitrateScale + 6);
    h[cpb_size_value_minus1(0)] = cpbSize >> (cpbScale + 4);
    h[cbr_flag(0)] = 1;

}
