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

#ifndef CHEETAH_PIXELTRIG_H
#define CHEETAH_PIXELTRIG_H

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

#include "FileConnector.h"
#include "Timepix.hpp"

template <typename event, int buffer_size, int n_buffer>
class CHEETAH_pixeltrig : public TIMEPIX<event, buffer_size, n_buffer>
{
private:
    // header
    int chip_id;
    uint64_t tpx_header = 861425748; //(b'TPX3', 'little')

    std::string &pattern_file;
    std::vector<uint64_t> pattern;
    uint64_t probe_count;

    // TDC
    uint64_t rise_t[4];
    uint64_t fall_t[4];
    bool rise_fall[4] = {false, false, false, false};
    int line_count[4] = {0, 0, 0, 0};
    int probe_count_chip[4] = {0, 0, 0, 0};

    int most_advanced_line = 0;
    int most_advanced_probe_count = 0;

    //buffer
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

    void parse_event(event *packet)
    {
        toa = ((((*packet & 0xFFFF) << 14) + ((*packet >> 30) & 0x3FFF)) << 4) + toa_offset;
        if (this->probe_count_chip[chip_id]%this->nxy < pattern.size()-1)
        {
            pack_44 = (*packet >> 44);
            uint64_t _probe_position = pattern[this->probe_count_chip[chip_id]%this->nxy];
            uint16_t _kx = (address_multiplier[chip_id] * (((pack_44 & 0x0FE00) >> 8) + ((pack_44 & 0x00007) >> 2)) + address_bias_x[chip_id]);
            uint16_t _ky = (address_multiplier[chip_id] * (((pack_44 & 0x001F8) >> 1) + (pack_44 & 0x00003)) + address_bias_y[chip_id]);
            ++this->n_events_processed;

            switch(this->functionType){
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::vstem:
                    this->vstem(_probe_position,_kx,_ky, this->id_image);
                    break;
                case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::multi_vstem:
                    this->multi_vstem(_probe_position,_kx,_ky,this->id_image);
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

    inline void schedule_buffer()
    {
        int buffer_id;

        while ((*this->p_processor_line)!=-1)
        {
            if (this->n_buffer_processed < this->n_buffer_filled)
            {
                buffer_id = this->n_buffer_processed % this->n_buf;

                process_buffer(&(this->buffer[buffer_id]));

                if (this->decluster) this->declusterer.set_buffer_read();
 
                ++this->n_buffer_processed;
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
        for (int j = 0; j < buffer_size; j++)
        {
            type = which_type(&(*p_buffer)[j]);
            if ((type == 2) & rise_fall[chip_id]) //currently always lose first line because no dwelltime known yet
            {
                parse_event(&(*p_buffer)[j]);
            }
        }
    };

    inline int which_type(event *packet)
    {
        if ((*packet & 0xFFFFFFFF) == tpx_header)
        {
            chip_id = (*packet >> 32) & 0xff;
            return 0;
        } // header
        else if (*packet >> 60 == 0x6)
        {
            process_tdc(packet);
            return 1;
        } // TDC
        else if (*packet >> 60 == 0xb)
        {
            return 2;
        } // event
        else
        {
            return 3;
        } // unknown
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
            }

            prev_tdc = rise_t[chip_id];
        }
        else if (((*packet >> 56) & 0x0F) == 10)  // TDC1 fall
        //OUDS or QD scan engine TDC line setting  6 - end of pixel clock
        {
            rise_fall[chip_id] = false;
            fall_t[chip_id] = ((*packet >> 9) & 0x7FFFFFFFF) + tdc_offset;

            if ((prev_tdc > fall_t[chip_id] + tdc_overflow_drop) && (this->current_line > 1) && (last_offset_line_tdc != this->current_line)) // tdc drop bigger than half of tdc range --> tdc must have overflowed
            {
            tdc_offset += 34359738368; 
            last_offset_line_tdc = this->current_line;
            }

            prev_tdc = fall_t[chip_id];

            ++probe_count_chip[chip_id];

            if ((probe_count_chip[chip_id] <= probe_count_chip[0]) & (probe_count_chip[chip_id] <= probe_count_chip[1]) & (probe_count_chip[chip_id] <= probe_count_chip[2]) & (probe_count_chip[chip_id] <= probe_count_chip[3]))
            { 
                this->probe_count = probe_count_chip[chip_id];
                this->current_line = this->probe_count/this->nx;
            }
            else if (probe_count_chip[chip_id] >= most_advanced_probe_count)
            {
                most_advanced_probe_count = probe_count_chip[chip_id];
                most_advanced_line = most_advanced_probe_count/this->nx;
                
                if (most_advanced_probe_count%this->nxy == 0)
                {
                    this->id_image = most_advanced_line / this->ny;
                    this->flush_image(this->id_image);
                }
            }

        }
    };
    
    void reset()
    {
        TIMEPIX<event, buffer_size, n_buffer>::reset();
        for (int i = 0; i < 4; i++)
        {
            rise_fall[i] = false;
            line_count[i] = 0;
        }
        most_advanced_probe_count = 0;
        most_advanced_line = 0;
        this->id_image = 0;
    };

    void read_patten_file()
    {
        std::ifstream file(pattern_file);
        if (file.is_open())
        {
            double element;
            while (file >> element)
            {
                pattern.push_back(element);
            }
            file.close();
            std::cout << "Pattern file read: " << pattern.size() << " positions" << std::endl;
        }
        else std::cout << "Unable to open pattern file" << std::endl;   
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
                this->read_thread = std::thread(&CHEETAH_pixeltrig<event, buffer_size, n_buffer>::read_file, this);
                break;
            }
            case 1:
            {
                //socket connection handled through seperate funtions in python binding
                this->read_thread = std::thread(&CHEETAH_pixeltrig<event, buffer_size, n_buffer>::read_socket, this);
                break;
            }
        }
        this->proc_thread = std::thread(&CHEETAH_pixeltrig<event, buffer_size, n_buffer>::schedule_buffer, this);
        this->starttime  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    };

    void terminate()
    {
        TIMEPIX<event, buffer_size, n_buffer>::terminate();
    };

    CHEETAH_pixeltrig(
        int &nx,
        int &ny,
        bool *b_cumulative,
        int repetitions,
        int *p_processor_line,
        int *p_preprocessor_line,
        int &mode,
        std::string &file_path,  
        SocketConnector socket,
        std::string &pattern_file
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
        ), pattern_file(pattern_file)
        {
            this->n_cam = 512;
            read_patten_file();
        }
};
#endif // CHEETAH_H
