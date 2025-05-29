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

#ifndef ELECTRON_H
#define ELECTRON_H

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

#include "dtype_Electron.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

class Electron : public LiveProcessor
{
private: 
   
public: 
    std::ofstream file_electron;

    bool decluster = true;
    uint64_t dtime = 100;
    uint16_t dspace = 6;
    int cluster_range = 256;
    int n_threads = 1;
    std::vector<int> clustersize_histogram = std::vector<int>(50,0);

    int x_crop = 0;
    int y_crop = 0;
    int scan_bin = 1;
    int detector_bin = 1;

    void run();
    void reset();
    void openDatFile();
    void closeDatFile();
    void close();

    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *p_prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    // Constructor
    Electron(int repetitions) : LiveProcessor(repetitions)
    {
    };

    // Destructor
    ~Electron()
    {
    };
};
#endif // !ELECTRON_H