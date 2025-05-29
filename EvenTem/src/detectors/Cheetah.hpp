/* Copyright (C) 2025 Thomas Friedrich, Chu-Ping Yu, Arno Annys
 * University of Antwerp - All Rights Reserved. 
 * You may use, distribute and modify
 * this code under the terms of the GPL3 license.
 * You should have received a copy of the GPL3 license with
 * this file. If not, please visit: 
 * https://www.gnu.org/licenses/gpl-3.0.en.html
 * 
 * Authors: 
 *   Thomas Friedrich <>
 *   Chu-Ping Yu <>
 *   Arno Annys <arno.annys@uantwerpen.be>
 */

#ifndef CHEETAH_H
#define CHEETAH_H

#ifdef __GNUC__
#define PACK(__Declaration__) __Declaration__ __attribute__((__packed__))
#endif

#ifdef _MSC_VER
#define PACK(__Declaration__) __pragma(pack(push, 1)) __Declaration__ __pragma(pack(pop))
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <atomic>
#include <vector>
#include <array>
#include <thread>
#include <chrono>
#include <set>

#include "FileConnector.h"
#include "Timepix.hpp"

namespace CHEETAH_ADDITIONAL
{
    // const size_t BUFFER_SIZE = 1024;
    // const size_t N_BUFFER = 512;
    const size_t BUFFER_SIZE = 4096;
    const size_t N_BUFFER = 2048;
    using EVENT = uint64_t;
}; 

#pragma pack(push, 1)
struct state_after_buffer
{
    uint64_t tdc_offset;
    uint64_t toa_offset;
    int line_count[4];
};
#pragma pack(pop)

template <typename event, int buffer_size, int n_buffer>
class CHEETAH : public TIMEPIX<event, buffer_size, n_buffer>
{
private:

    std::vector<state_after_buffer> state_after_buffer_list;

    // header
    int chip_id;
    uint64_t tpx_header = 861425748; //(b'TPX3', 'little')

    // TDC
    uint64_t rise_t[4];
    uint64_t fall_t[4];
    bool rise_fall[4] = {false, false, false, false};
    int line_count[4] = {0, 0, 0, 0};
    int most_advanced_line = 0;
    uint64_t line_interval;
    uint64_t dt;

    int type;

    // event
    uint64_t toa = 0;
    uint64_t pack_44;
    int address_multiplier[4] = {1,-1,-1,1};
    int address_bias_x[4] = {256, 511, 255, 0};
    int address_bias_y[4] = {0, 511, 511, 0};

    // overflow correction
    uint64_t prev_toa;
    std::atomic<uint64_t> toa_offset = 0;
    uint64_t toa_overflow_drop = 4294967296;
    int last_offset_line = 0;

    uint64_t prev_tdc = 0;
    uint64_t tdc_offset = 0;
    uint64_t tdc_overflow_drop = 17179869184; //half of 34359738368 = tdc range
    int last_offset_line_tdc = 0;
    std::thread check_overflow_thread;


    void parse_event(event *packet)
    {
        toa = ((((*packet & 0xFFFF) << 14) + ((*packet >> 30) & 0x3FFF)) << 4) + toa_offset;
        uint64_t _probe_position = ( toa - (rise_t[chip_id] * 2)) / dt;
        if (_probe_position < this->nx)
        {
            pack_44 = (*packet >> 44);
            _probe_position += (line_count[chip_id] % this->ny) * this->nx;
            uint16_t _kx = (address_multiplier[chip_id] * (((pack_44 & 0x0FE00) >> 8) + ((pack_44 & 0x00007) >> 2)) + address_bias_x[chip_id]);
            uint16_t _ky = (address_multiplier[chip_id] * (((pack_44 & 0x001F8) >> 1) + (pack_44 & 0x00003)) + address_bias_y[chip_id]);

            switch(this->functionType)
            {
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::vstem:
                    this->vstem(_probe_position,_kx,_ky, this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::multi_vstem:
                    this->multi_vstem(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::mask_vstem:
                    this->mask_vstem(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::com:
                    this->com(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_8: 
                    this->count_chunked_8(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_16:
                    this->count_chunked_16(_probe_position,_kx,_ky,this->id_image);
                    break;  
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_32:
                    this->count_chunked_32(_probe_position,_kx,_ky,this->id_image);
                    break;  
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::pacbed:
                    this->pacbed(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::var:
                    this->var(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi:
                    this->roi(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_mask:
                    this->roi_mask(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_4D:
                    this->roi_4D(_probe_position,_kx,_ky,this->id_image);
                    break;  
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_electron:
                    this->write_electron(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_declusterer_buffer:
                    this->write_declusterer_buffer(_probe_position,_kx,_ky,this->id_image,toa,0);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::information:
                    this->information(_probe_position,_kx,_ky,this->id_image);
                    break;
                #ifdef GPRI_OPTION_ENABLED
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::GPRI:
                    this->GPRI(_probe_position,_kx,_ky,this->id_image);
                    break;
                #endif
            }
            ++this->n_events_processed;
        }
    };

    void parse_event_w_tot(event *packet)
    {
        toa = ((((*packet & 0xFFFF) << 14) + ((*packet >> 30) & 0x3FFF)) << 4) - ((*packet >> 16) & 0xF) + toa_offset;
        this->tot = ((*packet) >> (16 + 4)) & 0x3ff;
        uint64_t _probe_position = ( toa - (rise_t[chip_id] * 2)) / dt;
        if (_probe_position < this->nx)
        {
            pack_44 = (*packet >> 44);
            _probe_position += (line_count[chip_id] % this->ny) * this->nx;
            uint16_t _kx = (address_multiplier[chip_id] * (((pack_44 & 0x0FE00) >> 8) + ((pack_44 & 0x00007) >> 2)) + address_bias_x[chip_id]);
            uint16_t _ky = (address_multiplier[chip_id] * (((pack_44 & 0x001F8) >> 1) + (pack_44 & 0x00003)) + address_bias_y[chip_id]);

            switch(this->functionType)
            {
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::vstem:
                    this->vstem(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::multi_vstem:
                    this->multi_vstem(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::mask_vstem:
                    this->mask_vstem(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::com:
                    this->com(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_8: 
                    this->count_chunked_8(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_16:
                    this->count_chunked_16(_probe_position,_kx,_ky,this->id_image);
                    break;  
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_32:
                    this->count_chunked_32(_probe_position,_kx,_ky,this->id_image);
                    break;  
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::pacbed:
                    this->pacbed(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::var:
                    this->var(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi:
                    this->roi_ToT(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_mask:
                    this->roi_mask(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_4D:
                    this->roi_4D(_probe_position,_kx,_ky,this->id_image);
                    break;  
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_electron:
                    this->write_electron(_probe_position,_kx,_ky,this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_declusterer_buffer:
                    this->write_declusterer_buffer(_probe_position,_kx,_ky,this->id_image,toa*25./16.,this->tot);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::information:
                    this->information(_probe_position,_kx,_ky,this->id_image);
                    break;
                #ifdef GPRI_OPTION_ENABLED
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::GPRI:
                    this->GPRI(_probe_position,_kx,_ky,this->id_image);
                    break;
                #endif
            }
            ++this->n_events_processed;
        }
    };

    void check_toa_overflow()
    {
        if ((prev_toa > toa + toa_overflow_drop) && (this->current_line > 1) && (last_offset_line != this->current_line)) // toa drop bigger than half of toa range --> toa must have overflowed
        {
            toa_offset += 17179869184; 
            last_offset_line = this->current_line;
            std::cout << "toa overflow at line " << this->current_line << std::endl;
        }
        prev_toa = toa;
    };

    void continous_check_toa_overflow(){
        while ((*this->p_processor_line)!=-1)
        {
            check_toa_overflow();
            std::this_thread::sleep_for(std::chrono::microseconds(1));
        }
    }

    inline void schedule_buffer()
    {
        int buffer_id;

        while ((*this->p_processor_line)!=-1)
        {
            if (this->n_buffer_processed < this->n_buffer_filled)
            {
                buffer_id = this->n_buffer_processed % this->n_buf;

                if (!this->repetitions_reached)
                { 
                    process_buffer(&(this->buffer[buffer_id]));
                    state_after_buffer_list.push_back({tdc_offset, toa_offset, line_count[0], line_count[1], line_count[2], line_count[3]});

                    if (this->decluster) this->declusterer.set_buffer_read();

                    ++this->n_buffer_processed;
                }
                *this->p_preprocessor_line = (int)this->current_line;
                check_toa_overflow();
            }
            else
            {
                this->process_wait++;
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    };

    inline void process_buffer(std::array<event, buffer_size> *p_buffer)
    {
        switch (this->b_tot)
        {
           case true:
           for (int j = 0; j < buffer_size; j++)
            {
                type = which_type(&(*p_buffer)[j]);
                if ((type == 2) && rise_fall[chip_id] && (!this->repetitions_reached)) 
                {
                    parse_event_w_tot(&(*p_buffer)[j]);
                }
            }
            break;
            case false:
            for (int j = 0; j < buffer_size; j++)
            {
                type = which_type(&(*p_buffer)[j]);
                if ((type == 2) && rise_fall[chip_id] && (!this->repetitions_reached)) 
                {
                    parse_event(&(*p_buffer)[j]);
                }
            }
        }
        if (this->repetitions_reached)  
        {
            this->probe_position_total = this->nxy*this->repetitions+1;
            this->id_image = this->repetitions;
        }
    };

    inline int which_type(event *packet)
    {
        if ((*packet & 0xFFFFFFFF) == tpx_header) // header
        {
            chip_id = (*packet >> 32) & 0xff;
            return 0;
        } 
        else if (*packet >> 60 == 0x6) // TDC
        {
            process_tdc(packet);
            return 1;
        } 
        else if (*packet >> 60 == 0xb) // event
        {
            return 2;
        }
        else if (*packet >> 60 == 0x4)
        {
            std::cout << "global time" << std::endl;
            return 3;
        }
        else // unknown
        {
            std::cout << "unknown packet type" << std::endl;
            return 3;
        }
    };

    inline void process_tdc(event *packet)
    {
        if (((*packet >> 56) & 0x0F) == 15) // TDC1 rise
        {
            rise_fall[chip_id] = true;
            rise_t[chip_id] = ((*packet >> 9) & 0x7FFFFFFFF) + tdc_offset;

            if ((prev_tdc > rise_t[chip_id] + tdc_overflow_drop) && (this->current_line > 1) && (last_offset_line_tdc != this->current_line)) // tdc drop bigger than half of tdc range --> tdc must have overflowed
            {
            tdc_offset += 34359738368; 
            last_offset_line_tdc = this->current_line;
            std::cout << "tdc overflow at line " << this->current_line << std::endl;
            }

            prev_tdc = rise_t[chip_id];
        }
        else if (((*packet >> 56) & 0x0F) == 10) // TDC1 fall
        {
            rise_fall[chip_id] = false;
            fall_t[chip_id] = ((*packet >> 9) & 0x7FFFFFFFF) + tdc_offset;

            if ((prev_tdc > fall_t[chip_id] + tdc_overflow_drop) && (this->current_line > 1) && (last_offset_line_tdc != this->current_line)) // tdc drop bigger than half of tdc range --> tdc must have overflowed
            {
            tdc_offset += 34359738368; 
            last_offset_line_tdc = this->current_line;
            }

            prev_tdc = fall_t[chip_id];

            ++line_count[chip_id];

            if ((line_count[chip_id] <= line_count[0]) & (line_count[chip_id] <= line_count[1]) & (line_count[chip_id] <= line_count[2]) & (line_count[chip_id] <= line_count[3]))
            { 
                this->current_line = line_count[chip_id];
            }
            else if (line_count[chip_id] >= most_advanced_line)
            {
                most_advanced_line = line_count[chip_id];
                if (most_advanced_line%this->ny == 0)
                {
                    this->id_image = most_advanced_line / this->ny ;
                    this->flush_image(this->id_image);
                }
            }

            line_interval = (fall_t[chip_id] - rise_t[chip_id]) * 2; //factor 2 for difference in time unit of tdc and toa
            // std::cout << "line interval: " << line_interval << std::endl;
            dt = line_interval / this->nx; //unit 1.5625 ns
            // std::cout << dt << std::endl;
        }
        if (this->current_line == this->ny * this->repetitions) this->repetitions_reached = true;
    };

    void skip_to_buffer(int buffer_id)
    {
        this->n_buffer_processed = buffer_id;
        int min_line = std::min({state_after_buffer_list[buffer_id].line_count[0], state_after_buffer_list[buffer_id].line_count[1], state_after_buffer_list[buffer_id].line_count[2], state_after_buffer_list[buffer_id].line_count[3]});
        this->current_line = min_line;
        toa_offset = state_after_buffer_list[buffer_id].toa_offset;
        tdc_offset = state_after_buffer_list[buffer_id].tdc_offset;
        line_count[0] = state_after_buffer_list[buffer_id].line_count[0];
        line_count[1] = state_after_buffer_list[buffer_id].line_count[1];
        line_count[2] = state_after_buffer_list[buffer_id].line_count[2];
        line_count[3] = state_after_buffer_list[buffer_id].line_count[3];
    }
    
    void reset()
    {
        TIMEPIX<event, buffer_size, n_buffer>::reset();
        for (int i = 0; i < 4; i++)
        {
            rise_fall[i] = false;
            line_count[i] = 0;
        }
        most_advanced_line = 0;
        this->id_image = 0;
    };

public:
    void run()
    {
        reset();
        switch (this->mode)
        {
            case 0:
            {
                this->file.path = this->file_path;
                this->file.open_file();
                this->read_thread = std::thread(&CHEETAH<event, buffer_size, n_buffer>::read_file, this);
                break;
            }
            case 1:
            {
                //socket connection handled through seperate funtions in python binding
                this->read_thread = std::thread(&CHEETAH<event, buffer_size, n_buffer>::read_socket, this);
                break;
            }
        }
        this->proc_thread = std::thread(&CHEETAH<event, buffer_size, n_buffer>::schedule_buffer, this);
        this->check_overflow_thread = std::thread(&CHEETAH<event, buffer_size, n_buffer>::continous_check_toa_overflow, this);
        this->starttime  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    };

    void terminate()
    {
        TIMEPIX<event, buffer_size, n_buffer>::terminate();
        while ((!check_overflow_thread.joinable()) ) { std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
        check_overflow_thread.join();
    };

    CHEETAH(
        int &nx,
        int &ny,
        int &dt, // unit: ns,
        bool *b_cumulative,
        int repetitions,
        int *p_processor_line,
        int *p_preprocessor_line,
        int &mode,
        std::string &file_path,  
        SocketConnector socket
    ) : TIMEPIX<event, buffer_size, n_buffer>(
            nx,
            ny,
            b_cumulative,
            repetitions,
            p_processor_line,
            p_preprocessor_line,
            mode,
            file_path,
            socket
        ), dt(dt*16/25)
        {
            this->n_cam = 512;
            if (dt == 0){
                std::cout << "Dwell time not provided! This means sacrificing the first line" << std::endl;
                dt = 1000; 
            }
        };
};
#endif // CHEETAH_H
