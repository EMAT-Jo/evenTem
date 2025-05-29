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

#ifndef ROI_H
#define ROI_H

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

#include "LiveProcessor.h"
#include "BoundedThreadPool.hpp"
#include "SocketConnector.h"
#include "ProgressMonitor.h"
#include "Cheetah.hpp"
#include "Timepix.hpp"
#include "Advapix.hpp"
#include "Simulated.hpp"
#include "Merlin.hpp"
#include "Numpy.hpp"
#include "FrameBased.hpp"
#include "Cheetah_pixeltrig.hpp"
#include "Roi4D.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;



class Roi : public LiveProcessor
{
private: 
   
public: 

    bool b_ROI_4D;
    int det_bin = 1;
    int bitdepth = 8;
    
    std::vector<uint64_t> Roi_diffraction_pattern;
    std::vector<uint64_t> Roi_scan_image;
    std::shared_ptr<Roi4D<uint8_t>> Roi_4D_8;
    std::shared_ptr<Roi4D<uint16_t>> Roi_4D_16;
    std::shared_ptr<Roi4D<uint32_t>> Roi_4D_32;


    std::vector<std::vector<uint64_t>> Roi_diffraction_pattern_stack; 
    std::vector<std::vector<uint64_t>> Roi_scan_image_stack;
    int lower_left[2] = {0, 0};
    int upper_right[2] = {1, 1};
    int L_0;
    int L_1;
    int finish_line;

    bool use_mask = false;
    std::vector<std::vector<int>> roi_mask;
    void set_roi_mask(std::vector<py::array_t<int>> arrays);
    void set_bitdepth(int bitdepth);


    void run();
    void reset();
    void set_roi(int,int,int,int);
    std::array<int,4> get_roi();

    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *p_prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    template <typename T>
    py::array_t<T> create_py_array(const std::vector<std::vector<std::vector<std::vector<T>>>>& data) {
        std::vector<size_t> shape = {data.size(), data[0].size(), data[0][0].size(), data[0][0][0].size()};
        std::vector<T> flat_data;
        for (const auto& dim1 : data) {
            for (const auto& dim2 : dim1) {
                for (const auto& dim3 : dim2) {
                    flat_data.insert(flat_data.end(), dim3.begin(), dim3.end());
                }
            }
        }
        return py::array_t<T>(shape, flat_data.data());
    }

    py::array get_4D() {
        if (bitdepth == 8) {
            return create_py_array<uint8_t>(Roi_4D_8->data);
        } else if (bitdepth == 16) {
            return create_py_array<uint16_t>(Roi_4D_16->data);
        } else if (bitdepth == 32) {
            return create_py_array<uint32_t>(Roi_4D_32->data);
        } else {
            throw std::runtime_error("Unsupported bitdepth");
        }
    }

    // Constructor
    Roi(int repetitions, bool ROI_4D) : LiveProcessor(repetitions),
    b_ROI_4D(ROI_4D)
    {
        Roi_4D_8 = std::make_shared<Roi4D<uint8_t>>();
        Roi_4D_16 = std::make_shared<Roi4D<uint16_t>>();
        Roi_4D_32 = std::make_shared<Roi4D<uint32_t>>();
    };

    // Destructor
    ~Roi(){};
};
#endif // !Roi_H