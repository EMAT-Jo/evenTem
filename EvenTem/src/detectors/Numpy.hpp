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

#ifndef NUMPY_H
#define NUMPY_H

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

#include "SocketConnector.h"
#include "FileConnector.h"
#include "BoundedThreadPool.hpp"
#include "FrameBased.hpp"

template <int n_cam,int buffer_size, int HEAD_SIZE, int n_buffer, typename pixel>
class NUMPY : public FRAMEBASED<n_cam,buffer_size,HEAD_SIZE,n_buffer,pixel>
{
protected:

    void parse_npy_header()
    {
        // Read magic string
        char magic[6];
        this->file.read_data(magic, 6);
        if (std::memcmp(magic, "\x93NUMPY", 6) != 0) {
            throw std::runtime_error("Not a valid .npy file");
        }

        // Read version number
        uint8_t major_version, minor_version;
        this->file.read_data(reinterpret_cast<char*>(&major_version), 1);
        this->file.read_data(reinterpret_cast<char*>(&minor_version), 1);

        int16_t header_len;
        if (major_version == 1){
            this->file.read_data(reinterpret_cast<char*>(&header_len), 2);
            header_len = static_cast<int>(header_len);        
        }
        else{
            throw std::runtime_error("Unsupported version number of npy file, should be version 1 for simple structures");
        }

        // Read header
        std::string header(header_len, ' ');
        this->file.read_data(&header[0], header_len);

        // Extract shape and dtype
        auto find_shape = header.find("'shape': (");
        auto find_descr = header.find("'descr': '");
        auto find_order = header.find("'fortran_order': ");

        if (find_shape == std::string::npos || find_descr == std::string::npos) {
            throw std::runtime_error("Header parsing failed");
        }

        // Parse shape
        size_t start = find_shape + 10;
        size_t end = header.find(")", start);
        std::string shape_str = header.substr(start, end - start);
        
        this->shape.clear();
        size_t pos = 0;
        while ((pos = shape_str.find(',')) != std::string::npos) {
            this->shape.push_back(std::stoi(shape_str.substr(0, pos)));
            shape_str.erase(0, pos + 1);
        }
        if (!shape_str.empty()) {
            this->shape.push_back(std::stoi(shape_str));  // Add the last dimension
        }

        // Parse dtype
        this->dtype = header.substr(find_descr + 10, header.find("'", find_descr + 10) - (find_descr + 10));

        // Data starts after the header
        this->data_offset = 10 + header_len; // 10 bytes for magic, version, and header length
        this->file.seek_to(this->data_offset);
    }

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
                    // if (this->GPRI_enabled || this->frame_torch_enabled) this->Binned_Tensor_buffer[_buffer_id] = ((torch::from_blob(this->frame_buffer[_buffer_id], {buffer_size,n_cam, n_cam}, torch::TensorOptions().dtype(torch::kUInt8)).to(this->device)).view({buffer_size,n_cam/this->GPRI_detector_bin, this->GPRI_detector_bin, n_cam/this->GPRI_detector_bin, this->GPRI_detector_bin}).sum({2, 4}, /*keepdim=*/false)).to(this->Tensortype);
                    ++this->n_buffer_filled;
                }
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                this->read_wait++;
            }

        }
        this->file.close_file();
    };

    
    inline void read_data_file(char *buffer, int data_size)
    {
        this->file.read_data(buffer, data_size);
    };

    inline void read_frame(std::array<pixel,n_cam*n_cam> &data)
    {
        int data_size = static_cast<int>(this->framesize * sizeof(pixel));
        char *buffer = reinterpret_cast<char *>(&data[0]);
        read_data_file(buffer, data_size);
    };

    int pre_run()
    {
        this->parse_npy_header();
        if (this->shape.size() != 4) throw std::runtime_error("Numpy dataset is not 4D");
        if (this->shape[2] != this->shape[3]) throw std::runtime_error("Numpy dataset does not have square detector images");

        std::cout << "shape: scan: " << this->shape[0] << "x" << this->shape[1] << ", detector: " << this->shape[2] << "x" << this->shape[3] << std::endl;

        this->framesize = this->shape[2] * this->shape[3];

        if (this->dtype == "|u1") {
            std::cout << "dtype: uint8" << std::endl;
            return 8;
        } else if (this->dtype == "<u2") {
            std::cout << "dtype: uint16" << std::endl;
            return 16;
        } else {
            throw std::runtime_error("Unsupported dtype");
        }
    }

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
                this->data_depth = pre_run();
                this->read_thread = std::thread(&NUMPY::buffer_reading, this);
                break;
            }
            case 1:
            {
                break;
            }
        }
        this->proc_thread = std::thread(&NUMPY::schedule_buffer, this);
        this->starttime  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    };

//-------------------------------------------------------------------------------------------------
// Variables
//-------------------------------------------------------------------------------------------------

    std::vector<int> shape;
    size_t data_offset;
    uint64_t framesize;

    NUMPY(
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
#endif // NUMPY_H