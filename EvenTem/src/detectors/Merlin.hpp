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

#ifndef MERLIN_H
#define MERLIN_H

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

namespace MERLIN_512{
    const int N_CAM = 512;
    const int BUFFER_SIZE = 128;
    const int HEAD_SIZE = 768;
    const int N_BUFFER = 32;
    using PIXEL = uint8_t;
}; 
namespace MERLIN_256{
    const int N_CAM = 256;
    const int BUFFER_SIZE = 128;
    const int HEAD_SIZE = 384;
    const int N_BUFFER = 32;
    using PIXEL = uint8_t;
};


template <int n_cam,int buffer_size, int HEAD_SIZE, int n_buffer, typename pixel>
class MERLIN : public FRAMEBASED<n_cam,buffer_size,HEAD_SIZE,n_buffer,pixel>
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
                read_frame(this->frame_buffer[_buffer_id][_frame_id], !this->first_frame);

                this->first_frame = false;
                ++this->n_frame_filled;

                #ifdef GPRI_OPTION_ENABLED
                if ((this->n_frame_filled%buffer_size == 0) && (this->n_buffer_filled < (n_buffer + this->n_buffer_processed))) 
                {
                    if (this->GPRI_enabled) this->Binned_Tensor_buffer[_buffer_id] = ((torch::from_blob(this->frame_buffer[_buffer_id], {buffer_size,n_cam, n_cam}, torch::TensorOptions().dtype(torch::kUInt8)).to(this->device)).view({buffer_size,n_cam/this->GPRI_detector_bin, this->GPRI_detector_bin, n_cam/this->GPRI_detector_bin, this->GPRI_detector_bin}).sum({2, 4}, /*keepdim=*/false)).to(this->Tensortype);
                }
                #endif
                #ifdef FRAMEBASED_TORCH_ENABLED
                if ((this->n_frame_filled%buffer_size == 0) && (this->n_buffer_filled < (n_buffer + this->n_buffer_processed))) 
                {
                    if (this->frame_torch_enabled) this->Tensor_buffer[_buffer_id] = torch::from_blob(this->frame_buffer[_buffer_id], {buffer_size,n_cam, n_cam}, torch::TensorOptions().dtype(torch::kUInt8)).to(this->device);
                }
                #endif
                if ((this->n_frame_filled%buffer_size == 0) && (this->n_buffer_filled < (n_buffer + this->n_buffer_processed))) ++this->n_buffer_filled;

            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                this->read_wait++;
            }

        }
        this->file.close_file();
    };


    inline void read_head_file()
    {
        this->file.read_data(&head_buffer[0], head_buffer.size());
    };
    
    inline void read_head_socket()
    {
        if (this->socket.read_data(&tcp_buffer[0], tcp_buffer.size()) == -1) perror("Merlin::read_head_socket(): Error reading TCP header from Socket!");
        if (this->socket.read_data(&head_buffer[0], head_buffer.size()) == -1) perror("Merlin::read_head_socket(): Error reading Frame header from Socket!");
    };

    int read_aquisition_header()
    {
        if (this->socket.read_data(&tcp_buffer[0], tcp_buffer.size()) == -1)
        {
            std::cout << "Error on recv() reading TCP-Header" << std::endl;
            return -1;
        }
        int l = decode_tcp_head();
        acq_header.resize(l);
        char *buffer = reinterpret_cast<char *>(&acq_header[0]);
        this->socket.read_data(buffer, l);
        char *p = strtok(buffer, " ");
        while (p)
        {
            acq += std::string(p) + "\n";
            std::cout << p << std::endl; // printing each token
            p = strtok(NULL, "\n");
        }
        this->socket.connection_information = acq;
        return 0;
    }


    int decode_tcp_head()
    {
        rcv.assign(tcp_buffer.cbegin(), tcp_buffer.cend());
        std::string n = rcv.substr(4, 10);

        try
        {
            int ds = stoi(n);
            return ds;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return 0;
        }
    }

    
    inline void read_data_file(char *buffer, int data_size)
    {
        this->file.read_data(buffer, data_size);
    };
    

    inline void read_data_socket(char *buffer, int data_size)
    {
        if (this->socket.read_data(buffer, data_size) == -1)
        {
            perror("Merlin::read_data_socket(): Error reading frame data from Socket!");
        }
    };

    template <typename T>
    inline void convert_binary_to_chars(std::array<T,n_cam*n_cam> &data)
    {
        size_t T_size = static_cast<size_t>(sizeof(T) * 8);
        size_t i_dat = static_cast<size_t>(data.size() / T_size);
        for (size_t i = i_dat - 1; i > 0; i--)
        {
            size_t idx = i * T_size;
            for (size_t j = 0; j < T_size; j++)
            {
                data[idx + j] = static_cast<T>((data[i] >> j) & 1);
            }
        }
    }

    inline bool read_head(bool decode)
    {
        switch (this->mode)
        {
            case 0:
            {
                try
                {
                    read_head_file();
                    if (decode)
                    {
                        decode_head();
                    }
                    return true;
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << '\n';
                    return false;
                }
                break;
            }
            case 1:
            {
                try
                {
                    read_head_socket();
                    if (decode)
                    {
                        decode_head();
                    }
                    return true;
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << '\n';
                    return false;
                }
                break;
            }

        }
    }

    inline void decode_head()
    {
        rcv.assign(head_buffer.cbegin(), head_buffer.cend());
        size_t i = 0;
        head.fill("");
        std::stringstream ss(rcv);
        while (ss.good() && i <= 6)
        {
            std::getline(ss, head[i], ',');
            i++;
        }
        if (i >= 6)
        {
            try
            {
                ds_merlin = stoi(head[4]) * stoi(head[5]);
                this->dtype = head[6];
                std::cout << "dtype: " << this->dtype << std::endl;
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << '\n';
            }
        }
        else
        {
            perror("Frame Header cannot be decoded!");
        }
    }


    inline void read_frame(std::array<pixel,n_cam*n_cam> &data, bool dump_head)
    {
        if (dump_head)
        {
            read_head(false);
        }
        int data_size = static_cast<int>(ds_merlin * sizeof(pixel));
        char *buffer = reinterpret_cast<char *>(&data[0]);
        if (b_binary)
        {
            data_size /= 8;
        }
        switch (this->mode)
        {
            case 0:
            {
                read_data_file(buffer, data_size);
                break;
            }
            case 1:
            {
                read_data_socket(buffer, data_size);
                break;
            }
        }

        if (b_binary)
        {
            convert_binary_to_chars(data);
        }
    };

    inline void init_uv()
    {
        this->u.resize(n_cam);
        this->v.resize(n_cam);

        for (int i = 0; i < n_cam; i++)
        {
            this->v[i] = i;
        }
        if (b_raw)
        // raw format: pixel sequence flipped every 64 bit
        {
            size_t num_per_flip = 64;
            switch (this->data_depth)
            {
            case 1:
                num_per_flip = 64;
                break;
            case 6:
                num_per_flip = 8;
                break;
            case 12:
                num_per_flip = 4;
                break;
            }
            size_t cnt = 0;
            for (size_t i = 0; i < (this->nx / num_per_flip); i++)
            {
                for (size_t j = (num_per_flip * (i + 1)); j > num_per_flip * i; j--)
                {
                    this->u[cnt] = (j - 1);
                    cnt++;
                }
            }
        }
        else
        {
            for (int i = 0; i < n_cam; i++)
            {
                this->u[i] = i;
            }
        }
    };

    int pre_run()
    {
        if (this->mode == 1)
        {
            if (read_aquisition_header() == -1)
            {
                perror("MerlinInterface::pre_run() could not obtain aquisition_header");
                return -1;
            }
        }

        if (read_head(true))
        {
            b_raw = false;
            b_binary = false;

            if (this->dtype == "U08")
            {
                return 8;
            }
            else if (this->dtype == "U16")
            {
                return 16;
            }
            else if (this->dtype == "R64")
            {
                b_raw = true;
                init_uv();
                switch (this->data_depth)
                {
                case 1:
                    b_binary = true;
                    return 8;
                case 6:
                    return 16;
                case 12:
                    return 16;
                default:
                    perror("Merlin::pre_run() supplied data_depth value not valid (must be 1, 6 or 12).");
                    return -1;
                }
            }
            else
            {
                perror("Merlin::pre_run() could not find valid data type description in frame header");
                return -1;
            }
        }
        else
        {
            perror("Merlin::pre_run() could not read frame header");
            return -1;
        }
        return -1;
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
                this->read_thread = std::thread(&MERLIN::buffer_reading, this);
                break;
            }
            case 1:
            {
                //socket connection handled through seperate funtions which have to be called in python binding
                this->data_depth = pre_run();
                this->read_thread = std::thread(&MERLIN::buffer_reading, this);
                break;
            }
        }
        this->proc_thread = std::thread(&MERLIN::schedule_buffer, this);
        this->starttime  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    };


//-------------------------------------------------------------------------------------------------
// Variables
//-------------------------------------------------------------------------------------------------

    std::array<char, HEAD_SIZE> head_buffer;
    std::array<std::string, 8> head;
    std::array<char, 15> tcp_buffer;
    std::string rcv;

    std::string acq;
    std::vector<char> acq_header;

    // Data Properties
    int ds_merlin;
    bool b_raw;
    bool b_binary;

    bool swap_endian;

    MERLIN(
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
#endif // MERLIN_H