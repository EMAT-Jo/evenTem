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

#ifndef LIVEPROCESSOR_H
#define LIVEPROCESSOR_H

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

#include "BoundedThreadPool.hpp"
#include "SocketConnector.h"
#include "ProgressMonitor.h"
#include "Cheetah.hpp"
#include "Timepix.hpp"
#include "Advapix.hpp"
#include "Simulated.hpp"


namespace CAMERA
{ 
    enum cameras
    {
        ADVAPIX,
        ADVARAW,
        CHEETAH,
        CHEETAH_PIXELTRIG,
        SIMULATED,
        MERLIN,
        NUMPY,
        HDF5,
        DUMMY
    };
}


class LiveProcessor
{protected:
    virtual void run() = 0;
    

public:
    bool event_mode = true;

    bool b_continuous=false;
    bool b_cumulative=true;

    // general
    SocketConnector socket;
    std::string file_path;
    std::string pattern_file;

    ProgressMonitor *p_prog_mon;
    void set_file(std::string filename);
    void set_pattern_file(std::string filename);    
    void set_socket(std::string,int,std::string);
    void set_dwell_time(int);
    void accept_socket();
    void close_socket();

    double progress_percent;
    bool progress_verbose = true;

    // settings
    int mode; // 0:file, 1:tcp
    CAMERA::cameras camera;

    void process_data();
    virtual void line_processor(size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool) = 0;

    // Scan Variables
    int nx; 
    int ny;
    int nxy;
    int n_cam;
    int dt;
    int rep;
    int fr_total;
    int *processor_line = new int;
    int *preprocessor_line  = new int;
    int id_image;

    // Variables for progress and performance
    int n_threads;
    int n_threads_max;
    int queue_size;
    float fr_freq;        // Frequncy per frame
    float fr_count;       // Count all Frames processed in an image
    float fr_count_total; // Count all Frames in a scanning session

    float processing_rate = 0;

    bool rc_quit = false;
    // int max_stall_count = 2147483647; 
    int max_stall_count = 1000000000; 

    // benchmarking
    std::chrono::high_resolution_clock::time_point startime;
    std::chrono::high_resolution_clock::time_point current_time;
    std::chrono::duration<double> elapsed_seconds;
    std::vector<int> reached_pp_id;
    std::vector<float> elapsed_seconds_vec;
        

    // Constructor
    LiveProcessor(int repetitions): 
        socket(),
        file_path(""),
        camera(CAMERA::DUMMY),
        p_prog_mon(nullptr),
        mode(0),
        nx(1024), ny(1024), nxy(0), n_cam(512), dt(0),
        rep(repetitions), fr_total(0),
        n_threads(1), queue_size(64),
        fr_freq(0.0), fr_count(0.0), fr_count_total(0.0)
    {
    };

    // Destructor
    ~LiveProcessor(){};
};
#endif // !LIVEPROCESSOR_H

