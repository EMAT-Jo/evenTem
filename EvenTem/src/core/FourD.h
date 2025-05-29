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

#ifndef FOURD_H
#define FOURD_H

#ifdef _WIN32
#include <io.h>
#pragma warning(disable : 4005 4333 34)
#else
#include <unistd.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>

#include <stdio.h>
#include <cfloat>
#include <vector>
#include <string>
#include <mutex>
#include <future>
#include <thread>
#include <chrono>
#include <algorithm>
#include <ctime>
#include <numeric>
#include <type_traits>
#include <cstdint> 

#include "LiveProcessor.h"
#include "BoundedThreadPool.hpp"
#include "SocketConnector.h"
#include "ProgressMonitor.h"
#include "Cheetah.hpp"
#include "Timepix.hpp"
#include "Advapix.hpp"
#include "Simulated.hpp"
#include "Merlin.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include "H5Cpp.h"

template<int BitDepth>
class FourD : public LiveProcessor
{
private: 
   
public: 

    int bitdepth;
    int deflate_factor;

    std::vector<uint64_t> Dose_image;

    typename std::conditional<BitDepth == 8, std::vector<uint8_t>,
                typename std::conditional<BitDepth == 16, std::vector<uint16_t>,
                    typename std::conditional<BitDepth == 32, std::vector<uint32_t>,
                        void // Handle unsupported bitdepths
                    >::type
                >::type
            >::type chunk_data[2]; 

    size_t det_bin = 1;
    size_t scan_bin = 1;
    size_t diff_pattern_size;
    size_t diff_pattern_length;
    size_t chunksize = 16;

    bool b_save_4D = false;
    bool b_first_run = true;
    bool b_accumulate = false;

    H5::DataSet dataset;
    H5::DataSet dataset4D;

    std::mutex mtx[2]; 

    void run();
    void reset();

    void write_and_clean(int);
    
    void allocate_chunk();

    void init_4D_file_8();
    void init_4D_file_16();
    void init_4D_file_32();
    void init_4D_file();


    void save_dose_image();
    void write_chunk(const std::vector<uint8_t>& chunk , hsize_t index);
    void write_chunk(const std::vector<uint16_t>& chunk , hsize_t index);
    void write_chunk(const std::vector<uint32_t>& chunk , hsize_t index);


    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *p_prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    // H5
    const char* H5FILE_NAME;
    H5::H5File h5file;


    // Constructor
    FourD(std::string f,int repetitions, int bitdepth, int deflate_factor) : LiveProcessor(repetitions),
    bitdepth(bitdepth),
    deflate_factor(deflate_factor)
    {
        f = f + ".hdf5";
        H5FILE_NAME = f.c_str();
        h5file = H5::H5File(H5FILE_NAME, H5F_ACC_TRUNC);

        if (!((bitdepth ==8) || (bitdepth == 16) || (bitdepth == 32)))
        {
            throw std::invalid_argument("bitdepth must be 8, 16 or 32");
        }
        if (deflate_factor < 1 || deflate_factor > 9)
        {
            throw std::invalid_argument("deflate factor must be between 1 and 9");
        }
    };

    // Destructor
    ~FourD()
    { 
    };
};
#endif // !FOURD_H
