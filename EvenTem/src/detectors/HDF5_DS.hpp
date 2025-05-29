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

#ifndef HDF5_DS_H
#define HDF5_DS_H

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
#include <functional>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <type_traits>
#include <cstdint> 

#include "SocketConnector.h"
#include "FileConnector.h"
#include "BoundedThreadPool.hpp"
#include "FrameBased.hpp"

#include "H5Cpp.h"

namespace HDF5_ADDITIONAL{
    const int N_CAM = 64;
    const int BUFFER_SIZE = 128;
    const int HEAD_SIZE = 0;
    const int N_BUFFER = 32;
    using PIXEL = uint8_t;
}; 

template <int n_cam,int buffer_size, int HEAD_SIZE, int n_buffer, typename pixel>
class HDF5 : public FRAMEBASED<n_cam,buffer_size,HEAD_SIZE,n_buffer,pixel>
{
protected:

    inline void buffer_reading()
    {
        int _buffer_id;
        int _frame_id;
        while (*this->p_processor_line!=-1)
        {
            if ((this->n_frame_filled < (buffer_size*n_buffer + this->n_frame_processed))) 
            {
                _frame_id = this->n_frame_filled % buffer_size; 
                _buffer_id = (this->n_frame_filled/buffer_size)%n_buffer;
                read_frame(this->frame_buffer[_buffer_id][_frame_id]);

                ++this->n_frame_filled;

                if ((this->n_frame_filled%buffer_size == 0) && (this->n_buffer_filled < (n_buffer + this->n_buffer_processed))) 
                {
                    ++this->n_buffer_filled;
                }
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                this->read_wait++;
            }
        }
        // this->file.close_file();
    };
    

    inline void read_frame(std::array<pixel,n_cam*n_cam> &data)
    {
        int data_size = static_cast<int>(this->framesize * sizeof(pixel));
        int rx = this->frame_index % this->nx;
        int ry = this->frame_index / this->nx;
        if (this->frame_index < this->nxy)
        {
            hsize_t offset[4] = {static_cast<hsize_t>(rx), static_cast<hsize_t>(ry), 0, 0};
            hsize_t count[4] = {1, 1, dims[2], dims[3]};             
            dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
            H5::DataSpace memspace(4, count);
            dataset.read(data.data(), H5::PredType::NATIVE_UINT8, memspace, dataspace);
            frame_index++;
        }
    };

    int pre_run()
    {
        framesize = n_cam*n_cam ;
        H5::H5File file(this->file_path, H5F_ACC_RDONLY);
        dataset = file.openDataSet("4D");
        dataspace = dataset.getSpace();
        dataspace.getSimpleExtentDims(dims, NULL);
        return 8;
    }

public:

    void run()
    {
        this->reset();

        switch (this->mode)
        {
            case 0:
            {
                this->data_depth = pre_run();
                this->read_thread = std::thread(&HDF5::buffer_reading, this);
                break;
            }
            case 1:
            {
                break;
            }
        }
        this->proc_thread = std::thread(&HDF5::schedule_buffer, this);
        this->starttime  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    };

//-------------------------------------------------------------------------------------------------
// Variables
//-------------------------------------------------------------------------------------------------

    int framesize;
    H5::DataSet dataset;
    H5::DataSpace dataspace;
    hsize_t dims[4];
    int frame_index = 0;

    HDF5(
        int &nx,
        int &ny,
        bool *b_cumulative,
        int repetitions,
        int *p_processor_line,
        int *p_preprocessor_line,
        int &mode,
        std::string &file_path,  
        SocketConnector socket
    ) : FRAMEBASED<n_cam,buffer_size,HEAD_SIZE,n_buffer,pixel>(
        nx,
        ny,
        b_cumulative,
        repetitions,
        p_processor_line,
        p_preprocessor_line,
        mode,
        file_path, 
        socket)
    {
    }

};
#endif // HDF5_DS_H