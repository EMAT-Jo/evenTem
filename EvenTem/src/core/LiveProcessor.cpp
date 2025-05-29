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

#include "LiveProcessor.h"

void LiveProcessor::process_data()
{
    // Start Thread Pool
    BoundedThreadPool pool;
    if (n_threads > 1)
        pool.init(n_threads, queue_size);

    size_t img_num = 0;
    size_t first_frame = img_num * nxy;
    size_t end_frame = (img_num + 1) * nxy;
    size_t fr_total_u = (size_t)fr_total;

    ProgressMonitor prog_mon(fr_total);
    p_prog_mon = &prog_mon;
    p_prog_mon->verbose = progress_verbose;

    int last_processor_line = 0;
    int stall_count = 0;

    while (*processor_line != -1)
    {
        line_processor(
            img_num, 
            first_frame, 
            end_frame, 
            p_prog_mon, 
            fr_total_u, 
            &pool
        );

        //check for stalling
        // if (last_processor_line == *processor_line)
        // {
        //     stall_count++;

        //     if (stall_count > max_stall_count)
        //     {
        //         std::cerr << " \n The processor seems to have stalled, this usually means that the data does not reach the expected number of probe positions. Check the number of probe positions and dwell time. \n" << std::endl;
        //         rc_quit = true;
        //     }
        // }
        // else
        // {
        //     stall_count = 0;
        //     last_processor_line = *processor_line;
        // }
    }
    p_prog_mon = nullptr;
}


void LiveProcessor::accept_socket()
{
    socket.accept_socket();
}


void LiveProcessor::set_socket(std::string ip, int port,std::string camera_name)
{
    socket.ip = ip;
    socket.port = port;
    mode = 1;
    if (camera_name == "ADVAPIX") 
    {
        camera = CAMERA::ADVAPIX;
        return;
    }

    else if (camera_name == "CHEETAH") 
    {
        socket.socket_type = Socket_type::SERVER;
        camera = CAMERA::CHEETAH;
    }
    else if (camera_name == "CHEETAH_PIXELTRIG")
    {
        socket.socket_type = Socket_type::SERVER;
        camera = CAMERA::CHEETAH_PIXELTRIG;
    }
    else if (camera_name == "MERLIN") 
    {
        socket.socket_type = Socket_type::CLIENT;
        camera = CAMERA::MERLIN;
        event_mode = false;
    }
    else std::cerr << "Unknown camera type: " << camera_name << std::endl;
    socket.connect_socket();
}


void LiveProcessor::close_socket()
{
    socket.close_socket();
}


void LiveProcessor::set_dwell_time(int dwell_time)
{
    dt = dwell_time;
}


void LiveProcessor::set_file(std::string filename)
{

    if (std::filesystem::path(filename).extension() == ".t3p") 
    {
        camera = CAMERA::ADVAPIX;
        n_cam = 256;
    }
    else if (std::filesystem::path(filename).extension() == ".t3r")
    {
        camera = CAMERA::ADVARAW;
        n_cam = 256;
    }
    else if (std::filesystem::path(filename).extension() == ".tpx3")
    {
        camera = CAMERA::CHEETAH;
        n_cam = 512;
    }
    else if (std::filesystem::path(filename).extension() == ".electron") camera = CAMERA::SIMULATED;
    else if (std::filesystem::path(filename).extension() == ".mib") camera = CAMERA::MERLIN;
    else if (std::filesystem::path(filename).extension() == ".npy") camera = CAMERA::NUMPY;
    else if (std::filesystem::path(filename).extension() == ".hdf5") camera = CAMERA::HDF5;
    else camera = CAMERA::DUMMY;

    if ((std::filesystem::path(filename).extension() == ".mib") || (std::filesystem::path(filename).extension() == ".npy") || (std::filesystem::path(filename).extension() == ".hdf5")) event_mode = false;
    else event_mode = true;

    file_path = filename;
    mode = 0;
}

void LiveProcessor::set_pattern_file(std::string  filename)
{
    pattern_file = filename;
    camera = CAMERA::CHEETAH_PIXELTRIG;
    std::cout << "Using pattern file for scan order"<< std::endl;
}
