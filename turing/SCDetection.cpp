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

#include "SCDetection.h"
#include <deque>
#include <string>
#include <fstream>
#include <algorithm>
#include <memory>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <numeric>
#include <cmath>


std::ifstream i_seq;
const int SHIFT_DOWN = 2;
const int THRESHOLD_CONST = 45;
const int LOWER_THRESHOLD = 7;
const int DELAY = 64;
const int WINDOW_SIZE = 12;
int last_shot_change = 0;
const double LIKE_THRESHOLD =1.1;


double calc_likelihood (double avg1, double var1, double avg2, double var2)
{
    double tmp = (avg2 - avg1) / 2.0;
    tmp = tmp*tmp;
    double tmpv = (var1+var2) / 2.0;
    tmp = (tmp+tmpv) * (tmp+tmpv);
    return tmp/(var1*var2);
}

double ShotChangeDetection::getLikelihood (FramePtr prev, FramePtr cur)
{
    double ret_l = 0.0;
    int block_w = width_ >> 3;
    int block_h = height_ >> 3;
    std::vector<double> prev_avg;
    std::vector<double> prev_var;
    //calculate avg and var for all blocks on prev frame
    const unsigned char* p_frame = &(*prev)[0];
    for (auto j=block_h; j<height_-block_h; j+=block_h)
    {
        for (auto i=block_w; i<width_-block_w; i+=block_w)
        {
            const unsigned char* p_block = p_frame+j*width_+i;
            int64_t avg = 0;
            for (auto h=0; h<block_h; ++h)
            {
                for (auto w=0; w<block_w; ++w)
                {
                    avg+=*(p_block+h*height_+w);
                }
            }
            prev_avg.push_back (static_cast<double> (avg)/(block_h*block_w));
        }
    }
    int k=0;
    for (auto j=block_h; j<height_-block_h; j+=block_h)
    {
        for (auto i=block_w; i<width_-block_w; i+=block_w)
        {
            const unsigned char* p_block = p_frame+j*width_+i;
            double var = 0.0;
            for (auto h=0; h<block_h; ++h)
            {
                for (auto w=0; w<block_w; ++w)
                {
                    double elem = static_cast<double> (*(p_block+h*height_+w));
                    var += (elem-prev_avg[k])*(elem-prev_avg[k]);
                }
            }
            var = var/(block_h * block_w);
            prev_var.push_back (var);
            ++k;
        }
    }
    // calculate avg and variance for some blocks of current frame
    std::vector<double> cur_avg;
    std::vector<double> cur_var;
    p_frame = &(*cur)[0];
    for (auto j=2*block_h; j<height_-(2*block_h); j+=block_h)
    {
        for (auto i=2*block_w; i<width_-(2*block_w); i+=block_w)
        {
            const unsigned char* p_block = p_frame+j*width_+i;
            int64_t avg = 0;
            for (auto h=0; h<block_h; ++h)
            {
                for (auto w=0; w<block_w; ++w)
                {
                    avg+=*(p_block+h*height_+w);
                }
            }
            cur_avg.push_back (static_cast<double> (avg)/(block_h*block_w));
        }
    }
    k=0;
    for (auto j=2*block_h; j<height_-(2*block_h); j+=block_h)
    {
        for (auto i=2*block_w; i<width_-(2*block_w); i+=block_w)
        {
            const unsigned char* p_block = p_frame+j*width_+i;
            double var = 0.0;
            for (auto h=0; h<block_h; ++h)
            {
                for (auto w=0; w<block_w; ++w)
                {
                    double elem = static_cast<double> (*(p_block+h*height_+w));
                    var += (elem-cur_avg[k])*(elem-cur_avg[k]);
                }
            }
            var = var/(block_h * block_w);
            cur_var.push_back (var);
            ++k;
        }
    }
    // calculate likelihood factor;
    int prev_w = 6;
    int cur_h = 4;
    int cur_w = 4;
    for (auto j=0; j<cur_h; ++j)
    {
        for (auto i=0; i<cur_w; ++i)
        {
            int prev_j = j+1;
            int prev_i = i+1;
            double l = 10000000.0;
            for (auto s=prev_j-1; s<prev_j+2; ++s)
            {
                for (auto k=prev_i-1; k<prev_i+2; ++k)
                {
                    int cur_idx = i+j*cur_w;
                    int prev_idx = k+s*prev_w;
                    double tmp_l = calc_likelihood (prev_avg[prev_idx], prev_var[prev_idx]
                                                                                 , cur_avg[cur_idx], cur_var[cur_idx]);
                    if (tmp_l < l)
                        l = tmp_l;
                }
            }
            ret_l += l;
        }
    }
    return ret_l/16;
}


void ShotChangeDetection::processSeq (std::vector<int>& shotChangeList)
{
    shotChangeList.resize (frameNum_);
    memset (&shotChangeList[0], 0, shotChangeList.size() * sizeof (shotChangeList[0]));
    int tot_window = WINDOW_SIZE * 2 + 1;
    if (frameNum_ < tot_window) return;
    try
    {
        i_seq.open (filename_, std::ios::binary);
        if (i_seq.fail())
        {
            std::cerr << "Failed to open input sequence " << filename_ << "\n";
            exit (0);
        }
    }
    catch (std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
    int lumaSize = width_*height_;// *((fileBitDepth_ > 8) ? 2 : 1);
    int frameOffset = frameSize_ - (lumaSize* ((fileBitDepth_ > 8) ? 2 : 1));
    i_seq.seekg (frameSkip_ * frameSize_, std::ios::cur);

    using FrameQueue = std::deque<std::shared_ptr<std::vector<unsigned char>>>;
    FrameQueue frame_queue;
    std::vector<int> cur (1<<(8-SHIFT_DOWN), 0), next (1<<(8-SHIFT_DOWN), 0);
    std::deque<int> dhist_queue;
    auto p_frame = std::make_shared<std::vector<unsigned char>> (lumaSize);
    // read in first frame
    if (fileBitDepth_ == 8)
    {
        i_seq.read(reinterpret_cast<char*> (&(*(p_frame->begin())))
                   , p_frame->size() * sizeof((*p_frame)[0]));
    }
    else
    {
        //std::vector<unsigned short> temp(lumaSize);
        auto temp = std::make_shared<std::vector<unsigned short>>(lumaSize);
        i_seq.read(reinterpret_cast<char*> (&(*(temp->begin()))), temp->size() * sizeof((*temp)[0]));
        for (int y = 0; y < lumaSize; ++y)
        {
            (*p_frame)[y] = ((*temp)[y] >>2);
        }
    }
    if (static_cast<long> (p_frame->size()*sizeof ((*p_frame)[0]) *  ((fileBitDepth_ > 8) ? 2 : 1)) != i_seq.gcount()
            || !i_seq)
    {
        std::cout << "error reading input file\n";
        exit (EXIT_FAILURE);
    }
    i_seq.seekg (frameOffset, std::ios::cur);
    for (auto& v:*p_frame)
    {
        unsigned char s = v>>SHIFT_DOWN;
        ++(cur)[s];
    }
    frame_queue.push_back (p_frame);
    for (auto i =1; i<tot_window; ++i)
    {
        auto p_frame = std::make_shared<std::vector<unsigned char>> (lumaSize);
        if (fileBitDepth_ == 8)
        {
            i_seq.read(reinterpret_cast<char*> (&(*(p_frame->begin())))
                       , p_frame->size() * sizeof((*p_frame)[0]));
        }
        else
        {
            //std::vector<unsigned short> temp(lumaSize);
            auto temp = std::make_shared<std::vector<unsigned short>>(lumaSize);
            i_seq.read(reinterpret_cast<char*> (&(*(temp->begin()))), temp->size() * sizeof((*temp)[0]));
            for (int y = 0; y < lumaSize; ++y)
            {
                (*p_frame)[y] = ((*temp)[y] >> 2);
            }
        }

        //i_seq.read (reinterpret_cast<char*> (&(*(p_frame->begin())))
        //    , p_frame->size() * sizeof ((*p_frame)[0]));
        if (static_cast<long> (p_frame->size()*sizeof ((*p_frame)[0])*  ((fileBitDepth_ > 8) ? 2 : 1)) != i_seq.gcount()
                || !i_seq)
        {
            std::cout << "error reading input file\n";
            exit (EXIT_FAILURE);
        }
        i_seq.seekg (frameOffset, std::ios::cur);
        frame_queue.push_back (p_frame);
        memset (&next[0], 0, next.size()*sizeof (next[0]));
        for (auto& v:*p_frame)
        {
            unsigned char s = v>>SHIFT_DOWN;
            ++(next)[s];
        }
        int dhist=0;
        for (unsigned int k=0; k<cur.size(); ++k)
        {
            dhist += abs (cur[k] - next[k]);
        }
        dhist_queue.push_back (dhist);
        next.swap (cur);
    }
    for (auto i= tot_window; i <frameNum_; ++i)
    {
        double left_avg = std::accumulate (dhist_queue.begin(), dhist_queue.begin()+WINDOW_SIZE, 0.0) / (double)WINDOW_SIZE;
        double right_avg = std::accumulate (dhist_queue.begin()+WINDOW_SIZE +1, dhist_queue.end(), 0.0) / (double)WINDOW_SIZE;
        double temp = std::accumulate (dhist_queue.begin(), dhist_queue.begin() + WINDOW_SIZE, 0.0
                                       , [left_avg] (double current_sum, double elem)
                                       { return current_sum += (elem-left_avg)*(elem-left_avg);});
        double left_stddev = std::sqrt (temp/ (double)WINDOW_SIZE);
        temp = std::accumulate (dhist_queue.begin() + WINDOW_SIZE + 1, dhist_queue.end(), 0.0
                                , [right_avg] (double current_sum, double elem)
                                { return current_sum += (elem-right_avg)*(elem-right_avg);});
        double right_stddev = std::sqrt (temp / (double)WINDOW_SIZE);
        int max = *(std::max_element (dhist_queue.begin(), dhist_queue.end()));
        double th_max = std::max (left_avg + THRESHOLD_CONST * left_stddev
                                  , right_avg + THRESHOLD_CONST * right_stddev);
        double th_min = std::max (left_avg + LOWER_THRESHOLD * left_stddev
                                  , right_avg + LOWER_THRESHOLD * right_stddev);
        if (!(dhist_queue[WINDOW_SIZE] < max))
        {
            if (dhist_queue[WINDOW_SIZE] > th_max && (i-(WINDOW_SIZE+1)-last_shot_change>DELAY))
            {
                last_shot_change = i- (WINDOW_SIZE);
                shotChangeList[i- (WINDOW_SIZE)] = 1;
            }
            else if (dhist_queue[WINDOW_SIZE] > th_min && (i- (WINDOW_SIZE + 1) -last_shot_change>DELAY))
            {
                //std::cout << i-11 << "   - possible shot change";
                double likelihood = getLikelihood (frame_queue[WINDOW_SIZE], frame_queue[(WINDOW_SIZE + 1)]);
                if (likelihood < LIKE_THRESHOLD)
                {
                    last_shot_change = i- (WINDOW_SIZE);
                    shotChangeList[i- (WINDOW_SIZE)] = 1;
                }
            }
        }
        std::shared_ptr<std::vector<unsigned char>> p_f = frame_queue.front();
        frame_queue.pop_front();
        if (fileBitDepth_ == 8)
        {
            i_seq.read(reinterpret_cast<char*> (&(*(p_f->begin())))
                       , p_f->size() * sizeof((*p_f)[0]));
        }
        else
        {
            //std::vector<unsigned short> temp(lumaSize);
            auto temp = std::make_shared<std::vector<unsigned short>>(lumaSize);
            i_seq.read(reinterpret_cast<char*> (&(*(temp->begin()))), temp->size() * sizeof((*temp)[0]));
            for (int y = 0; y < lumaSize; ++y)
            {
                (*p_f)[y] = ((*temp)[y] >> 2);
            }
        }
        //i_seq.read (reinterpret_cast<char*> (&(*(p_f->begin())))
        //   , p_f->size() * sizeof ((*p_f)[0]));
        if (static_cast<long> (p_f->size()*sizeof ((*p_f)[0])*  ((fileBitDepth_ > 8) ? 2 : 1)) != i_seq.gcount()
                || !i_seq)
        {
            std::cout << "error reading input file\n";
            exit (EXIT_FAILURE);
        }
        i_seq.seekg (frameOffset, std::ios::cur);
        frame_queue.push_back (p_f);
        memset (&next[0], 0, next.size()*sizeof (next[0]));
        for (auto& v:*p_f)
        {
            unsigned char s = v>>SHIFT_DOWN;
            ++(next)[s];
        }
        int dhist=0;
        for (unsigned int k=0; k<cur.size(); ++k)
        {
            dhist += abs (cur[k] - next[k]);
        }
        dhist_queue.pop_front();
        dhist_queue.push_back (dhist);
        next.swap (cur);
    }
    i_seq.close();
}
