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

#ifndef VAR_H
#define VAR_H

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

#include <pybind11/pybind11.h>
namespace py = pybind11;



class Var : public LiveProcessor
{
private: 
   
    std::mutex stem_mutex;

public: 
    
    std::array<float, 2> offset = {0, 0};
    bool auto_offset = true;

    std::vector<float> Var_image;
    std::vector<size_t> Var_data[2];
    std::atomic<bool> rescale_stem;

    void run();
    void reset();
    void set_offset(std::array<float, 2>);

    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *p_prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    // Constructor
    Var(int repetitions) : LiveProcessor(repetitions)
    {
    };

    // Destructor
    ~Var(){};
};
#endif // !Var_H