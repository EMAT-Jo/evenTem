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

#ifndef SIMULATED_H
#define SIMULATED_H

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
#include "BoundedThreadPool.hpp"

namespace SIMULATED_ADDITIONAL
{
    // const size_t BUFFER_SIZE = 6400; //match L3 cache
    // const size_t N_BUFFER = 512;
    // const size_t BUFFER_SIZE = 204800; //match L3 ans L2 cache
    // const size_t N_BUFFER = 64;
    const size_t BUFFER_SIZE = 115200;
    const size_t N_BUFFER = 8;
    PACK(struct EVENT
    {
        uint16_t kx;
        uint16_t ky;
        uint16_t rx;
        uint16_t ry;
        uint16_t id_image;
    });
};


template <typename event, int buffer_size, int n_buffer>
class SIMULATED : public TIMEPIX<event, buffer_size, n_buffer>
{ 
private:
    BoundedThreadPool *event_parsing_pool = new BoundedThreadPool;

    inline void schedule_buffer()
    {
        int buffer_id;

        while ((*this->p_processor_line)!=-1)
        {
            if (this->n_buffer_processed < this->n_buffer_filled)
            {
                buffer_id = this->n_buffer_processed % this->n_buf;

                if (!this->repetitions_reached){ 
                    process_buffer(&(this->buffer[buffer_id]));
                    // event_parsing_pool->push_task([=]{process_buffer(&(this->buffer[buffer_id]));});
                    // if (this->decluster) this->declusterer.set_buffer_read();
                    ++this->n_buffer_processed;
                }
                *this->p_preprocessor_line = (int)this->current_line;
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
            if (!this->repetitions_reached) 
            {
                process_event(&(*p_buffer)[j]);
            }
        }

        if (!this->repetitions_reached) { 
            this->current_line = (*p_buffer).back().ry + this->ny * (*p_buffer).back().id_image;
            this->probe_position_total = (*p_buffer).back().ry * this->nx + (*p_buffer).back().rx;
            this->id_image = (*p_buffer).back().id_image;
        }
        else {
            this->current_line =  this->ny * this->repetitions;
            this->probe_position_total = this->nxy*this->repetitions+1;
            this->id_image = this->repetitions;
        }

    }
    
    inline void process_event(event *packet)
    {
        if (packet->id_image >= this->repetitions) {
            this->repetitions_reached = true; 
            return;
        }
        // this->process[0]((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image); //way slower
        // (this->*this->ProcessFun)((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image); //quite slower
        switch(this->functionType){
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::vstem:
                this->vstem((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::atomic_vstem:
                this->atomic_vstem((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::multi_vstem:
                this->multi_vstem((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::mask_vstem:
                this->mask_vstem((uint64_t)(packet->ry * this->ny + packet->rx),packet->kx,packet->ky,this->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::com:
                this->com((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_8: 
                this->count_chunked_8((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_16:
                this->count_chunked_16((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_32:
                this->count_chunked_32((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::pacbed:
                this->pacbed((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::var:
                this->var((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi:
                this->roi((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_4D:
                this->roi_4D((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_electron:
                this->write_electron((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::information:
                this->information((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            #ifdef GPRI_OPTION_ENABLED
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::GPRI:
                this->GPRI((uint64_t)(packet->ry * this->ny + packet->rx), packet->kx, packet->ky, packet->id_image);
                break;
            #endif
        }
        ++this->n_events_processed;
    };
    
public:

    void run()
    {
        this->reset();
        switch (this->mode)
        {
            case 0:
            {
                this->file.path = this->file_path;
                this->file.open_file();
                this->read_thread = std::thread(&SIMULATED<event, buffer_size, n_buffer>::read_file, this);
                break;
            }
            case 1:
            {
                throw std::invalid_argument("simulated must be file mode");
            }
        }
        this->proc_thread = std::thread(&SIMULATED<event, buffer_size, n_buffer>::schedule_buffer, this);
        this->starttime  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    };

    SIMULATED(
        int &nx,
        int &ny,
        int _n_cam,
        bool *b_cumulative,
        int repetitions,
        int *p_processor_line,
        int *p_preprocessor_line,
        int &mode,
        std::string &file_path,  
        SocketConnector &socket
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
    ) 
    {
        this->n_cam = _n_cam;
        // event_parsing_pool->init(4, 16);
    }
};

#endif // SIMULATED_H