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

#ifndef VSTEM_H
#define VSTEM_H

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
#include "Cheetah_pixeltrig.hpp"
#include "Timepix.hpp"
#include "Advapix.hpp"
#include "Simulated.hpp"
#include "Merlin.hpp"
#include "Numpy.hpp"
#include "HDF5_DS.hpp"
#include "AtomicWrapper.hpp"
#include "AnnularDetector.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

class vSTEM : public LiveProcessor
{
private: 

public: 
    AnnularDetector detector = AnnularDetector(20, 40);

    std::vector<float> inner_radia = {20};
    std::vector<float> outer_radia = {40};
    std::vector<std::array<float, 2>> offsets = {{0, 0}};

    std::vector<size_t> vSTEM_image;
    std::vector<atomwrapper<int>> atomic_vSTEM_image;
    std::vector<std::vector<size_t>> vSTEM_stack;

    bool allow_torch = false;
    bool allow_cuda = false;

    bool auto_offset = true;

    void run();
    void reset();
    void from_atomic();

    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *p_prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    void set_offsets(std::vector<std::array<float, 2>>);
    std::vector<int> compute_detector();
    std::vector<int> get_detector();
    void set_detector_mask(py::array_t<int> mask);
    std::vector<int> detector_mask;
    bool use_mask = false;

    // Constructor
    vSTEM(int repetitions) : LiveProcessor(repetitions)
    {
    };

    // Destructor
    ~vSTEM(){};
};
#endif // !vSTEM_H