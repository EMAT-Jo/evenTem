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

#ifndef ADVAPIX_H
#define ADVAPIX_H

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


namespace ADVAPIX_ADDITIONAL
{   
    const size_t BUFFER_SIZE = 14400;
    const size_t N_BUFFER = 1024;

    PACK(struct EVENT
    {
        uint32_t index;
        uint64_t toa;
        uint8_t overflow;
        uint8_t ftoa;
        uint16_t tot;
    });
};

#ifdef PIXET_ENABLED
#include "pxcapi.h"
#define SINGLE_CHIP_PIXSIZE      65536
#define ERRMSG_BUFF_SIZE         512
#define PAR_DDBLOCKSIZE         "DDBlockSize"
#define PAR_DDBUFFSIZE          "DDBuffSize"
#define PAR_DUMMYSPEED          "DDDummyDataSpeed"
#define PAR_BLOCKCOUNT          "BlockCount"
#define PAR_PROCESSDATA         "ProcessData"
#define PAR_TRG_STG             "TrgStg"
#define BLOCKSIZE 50
#define BUFFSIZE  5000
struct CallbackData_t
{
    int deviceIndex;
    int *p_n_buffer_filled;
    int *p_n_buffer_processed;
    std::shared_ptr<Tpx3Pixel[]> (*p_ragged_buffer)[ADVAPIX_ADDITIONAL::N_BUFFER];
    int (*p_ragged_buffer_sizes)[ADVAPIX_ADDITIONAL::N_BUFFER];
};

extern void OnTpx3Data(intptr_t eventData, intptr_t userData);
#else
struct Tpx3Pixel
{
    uint32_t index;
    uint64_t toa;
    uint8_t overflow;
    uint8_t ftoa;
    uint16_t tot;
};
#endif

template <typename event, int buffer_size, int n_buffer>
class ADVAPIX : public TIMEPIX<event, buffer_size, n_buffer>
{ 
private:
    int dt;

    BoundedThreadPool *event_parsing_pool = new BoundedThreadPool;

    #ifdef PIXET_ENABLED
    unsigned deviceIndex = 0 ;
    CallbackData_t *CallbackData;
    intptr_t CallbackDataAsIntPtr;
    bool trig = true;
    int meausure_time = 5;

    void connect()
    {
        int rc = pxcInitialize();
        if (rc) printf("Could not initialize Pixet");

        int connectedDevicesCount = pxcGetDevicesCount();
        printf("Connected devices: %d\n", connectedDevicesCount);

        if (connectedDevicesCount == 0) int pxcExit();

        char deviceName[256];
        memset(deviceName, 0, 256);
        pxcGetDeviceName(deviceIndex, deviceName, 256);

        char chipID[256];
        memset(chipID, 0, 256);
        pxcGetDeviceChipID(0, 0, chipID, 256);

        printf("Device %d Name %s, (ChipID: %s)\n", deviceIndex, deviceName, chipID);
        
        double th;
        rc = pxcGetThreshold(deviceIndex,0,&th);
        std::cout<< "threshold: "<< th << std::endl;

        double bias;
        rc = pxcGetBias(deviceIndex,&bias);
        std::cout<< "bias voltage: "<< bias << std::endl;

        pxcSetDeviceParameter(deviceIndex, PAR_PROCESSDATA, 1);

        pxcSetTimepix3Mode(deviceIndex, PXC_TPX3_OPM_TOA);
        pxcSetDeviceParameter(deviceIndex, PAR_DDBLOCKSIZE, BLOCKSIZE); //in MegaBytes unlike what docs say
        pxcSetDeviceParameter(deviceIndex, PAR_DDBUFFSIZE, BUFFSIZE); // in MegaBytes 1000

        CallbackData = new CallbackData_t{(int)deviceIndex, &this->n_buffer_filled, &this->n_buffer_processed};
        p_ragged_buffer = new std::shared_ptr<Tpx3Pixel[]>[1][n_buffer];
        for (int i = 0; i < n_buffer; ++i){
            (*p_ragged_buffer)[i] = std::shared_ptr<Tpx3Pixel[]>(new Tpx3Pixel[1], std::default_delete<Tpx3Pixel[]>());
            (*CallbackData->p_ragged_buffer)[i] = (*p_ragged_buffer)[i];
        }
        p_ragged_buffer_sizes = new int[1][n_buffer];
        CallbackData->p_ragged_buffer_sizes = p_ragged_buffer_sizes;
        CallbackDataAsIntPtr = reinterpret_cast<intptr_t>(CallbackData);

        if (trig == true) rc = pxcMeasureTpx3DataDrivenMode(deviceIndex, meausure_time, "", PXC_TRG_HWSTART, OnTpx3Data,CallbackDataAsIntPtr);
        else rc = pxcMeasureTpx3DataDrivenMode(deviceIndex, meausure_time, "", PXC_TRG_NO, OnTpx3Data,CallbackDataAsIntPtr);
    };
    #endif
 
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
                    switch (this->mode)
                    {
                        case 0:
                            process_buffer(&(this->buffer[buffer_id]));
                            // event_parsing_pool->push_task([=]{process_buffer(&(this->buffer[buffer_id]));});
                            break;
                        case 1:
                            std::cout << (*p_ragged_buffer_sizes)[buffer_id] << std::endl;
                            process_ragged_buffer((*p_ragged_buffer)[buffer_id],  (*p_ragged_buffer_sizes)[buffer_id]);
                            break;
                    }
                    if (this->decluster) {
                        this->declusterer.set_buffer_read();
                        #ifdef DBG_LOG 
                            Logger::getInstance().log("filled buffer " + std::to_string(n_buffer_filled));
                        #endif
                    }
                    ++this->n_buffer_processed;
                    *this->p_preprocessor_line = (int)this->current_line;
                }
            }
            else
            {
                this->process_wait++;
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    };

    inline void process_ragged_buffer(std::shared_ptr<Tpx3Pixel[]> p_buffer, size_t size)
    {
        for (int j = 0; j < size; j++)
        {
            if (!this->repetitions_reached) process_event(&p_buffer[j]);
        }

        if (!this->repetitions_reached) {
            this->probe_position_total = p_buffer[size-1].toa * 25 / this->dt;
            this->current_line = floor(this->probe_position_total / this->nx);
            this->id_image = this->probe_position_total / this->nxy;
        }
        else 
        {
            this->current_line =  this->ny * this->repetitions;
            this->probe_position_total = this->nxy+1;
            this->id_image = this->repetitions;
        }
    };



    inline void process_buffer(std::array<event, buffer_size> *p_buffer)
    {
        for (int j = 0; j < buffer_size; j++)
        {
            if (!this->repetitions_reached) 
            {
                process_event(&(*p_buffer)[j]);
                // for (int i_proc=0; i_proc<this->n_proc; i_proc++) { this->process[i_proc](); }
            }
        }

        if (!this->repetitions_reached) 
        {
            this->probe_position_total = (*p_buffer).back().toa * 25 / this->dt;
            this->current_line = floor(this->probe_position_total / this->nx);
            this->id_image = this->probe_position_total / this->nxy;
        }
        else 
        {
            this->current_line =  this->ny * this->repetitions;
            this->probe_position_total = this->nxy*this->repetitions+1;
            this->id_image = this->repetitions;
        }
    };

    uint64_t latest_pp = 0;
    
    inline void process_event(event *packet)
    {
        uint64_t _probe_position_total = packet->toa * 25 / this->dt;
        uint16_t _kx = packet->index % this->n_cam;
        uint16_t _ky = packet->index / this->n_cam;
        uint16_t _id_image = _probe_position_total / this->nxy;

        // if (_probe_position_total >= this->nxy*this->repetitions) 
        // {
        //     this->repetitions_reached = true;
        //     return;
        // }

        switch(this->functionType)
        {
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::vstem:
                this->vstem(_probe_position_total%this->nxy,_kx,_ky, _id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::multi_vstem:
                this->multi_vstem(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::mask_vstem:
                this->mask_vstem(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::atomic_vstem:
                this->atomic_vstem(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::com:
                this->com(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_8: 
                this->count_chunked_8(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_16:
                this->count_chunked_16(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_32:
                this->count_chunked_32(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::pacbed:
                this->pacbed(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::var:
                this->var(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi:
                this->roi(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_mask:
                this->roi_mask(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_4D:
                this->roi_4D(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_electron:
                this->write_electron(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_declusterer_buffer:
                this->write_declusterer_buffer(_probe_position_total%this->nxy,_kx,_ky,_id_image,packet->toa*25,packet->tot);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::information:
                this->information(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            #ifdef GPRI_OPTION_ENABLED
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::GPRI:
                this->GPRI(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            #endif
        }
        ++this->n_events_processed;
    };

    inline void process_event(Tpx3Pixel *packet)
    {
        uint64_t _probe_position_total = packet->toa * 25 / this->dt;
        uint16_t _kx = packet->index % this->n_cam;
        uint16_t _ky = packet->index / this->n_cam;
        uint16_t _id_image = _probe_position_total / this->nxy;

        if (_probe_position_total >= this->nxy*this->repetitions) 
        {
            this->repetitions_reached = true;
            return;
        }

        switch(this->functionType)
        {
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::vstem:
                this->vstem(_probe_position_total%this->nxy,_kx,_ky, _id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::multi_vstem:
                this->multi_vstem(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::mask_vstem:
                this->mask_vstem(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::com:
                this->com(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_8: 
                this->count_chunked_8(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_16:
                this->count_chunked_16(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::count_chunked_32:
                this->count_chunked_32(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::pacbed:
                this->pacbed(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::var:
                this->var(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi:
                this->roi(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_mask:
                this->roi_mask(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::roi_4D:
                this->roi_4D(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;  
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_electron:
                this->write_electron(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::write_declusterer_buffer:
                this->write_declusterer_buffer(_probe_position_total%this->nxy,_kx,_ky,_id_image,packet->toa*25,packet->tot);
                break;
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::information:
                this->information(_probe_position_total%this->nxy,_kx,_ky,_id_image);
                break;
            #ifdef GPRI_OPTION_ENABLED
            case TIMEPIX<event, buffer_size, n_buffer>::FunctionType::GPRI:
                this->GPRI(_probe_position_total%this->nxy,_kx,_ky,_id_image);
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
                this->read_thread = std::thread(&ADVAPIX<event, buffer_size, n_buffer>::read_file, this);
                break;
            }
            case 1:
            {   
                #ifdef PIXET_ENABLED
                this->read_thread = std::thread(&ADVAPIX<event, buffer_size, n_buffer>::connect, this);
                #endif
                break;
            }
        }
        this->proc_thread = std::thread(&ADVAPIX<event, buffer_size, n_buffer>::schedule_buffer, this);
        this->starttime  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    };

    std::shared_ptr<Tpx3Pixel[]> (*p_ragged_buffer)[n_buffer];
    int (*p_ragged_buffer_sizes)[n_buffer];


    ADVAPIX(
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
        ) ,
        dt(dt)
        {this->n_cam = 256;
        if (dt == 0) std::cout << "Dwell time not provided!" << std::endl;
        // event_parsing_pool->init(1, 16);
        }

    #ifdef PIXET_ENABLED
    ~ADVAPIX(){
        int pxcExit();
    };
    #endif
};

#endif // ADVAPIX_H