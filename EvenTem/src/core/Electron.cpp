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

#include "Electron.h"

#ifdef __GNUC__
#define PACK(__Declaration__) __Declaration__ __attribute__((__packed__))
#endif
#ifdef _MSC_VER
#define PACK(__Declaration__) __pragma(pack(push, 1)) __Declaration__ __pragma(pack(pop))
#endif


void Electron::openDatFile()
{
    size_t lastindex =  file_path.find_last_of("."); 
    std::string filename_base = file_path.substr(0, lastindex);   
    file_electron.open(filename_base+".electron", std::ios::out | std::ios::binary);
    if (!(file_electron.is_open()))
        {
            throw std::runtime_error("Error opening dat file!");
        }
}

void Electron::closeDatFile()
{
    if (file_electron.is_open())
    {
        file_electron.close();
    }
}
    


void Electron::run(){
    py::gil_scoped_release release;

    openDatFile();
    reset();
    // Run camera dependent pipeline
    switch (camera)
    {
        case CAMERA::ADVAPIX:
        {   
            using namespace ADVAPIX_ADDITIONAL;
            ADVAPIX<EVENT, BUFFER_SIZE, N_BUFFER> cam(
                nx, 
                ny,  
                dt,
                &b_cumulative,
                rep,
                processor_line,
                preprocessor_line,
                mode,
                file_path,
                socket
            );
            cam.enable_electron(file_electron,decluster,dtime,dspace,cluster_range, x_crop, y_crop, scan_bin, detector_bin, n_threads, &clustersize_histogram);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
        case CAMERA::CHEETAH:
        {

            using namespace CHEETAH_ADDITIONAL;
            CHEETAH<EVENT, BUFFER_SIZE, N_BUFFER> cam(
                nx, 
                ny, 
                dt,
                &b_cumulative,
                rep,
                processor_line,
                preprocessor_line,
                mode,
                file_path,
                socket
            );
            cam.enable_electron(file_electron,decluster,dtime,dspace,cluster_range, x_crop, y_crop, scan_bin, detector_bin, n_threads, &clustersize_histogram);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
    }
    rc_quit = true;
    close();
}

void Electron::reset()
{
    rc_quit = false;
    fr_freq = 0;

    // Initializations
    nxy = nx * ny;
    id_image = 0;
    fr_total = nxy * rep;
    fr_count = 0;

    // Data Processing Progress
    *processor_line = 0;
    *preprocessor_line = 0;

}

void Electron::line_processor(
    size_t &img_num,
    size_t &first_frame,
    size_t &end_frame,
    ProgressMonitor *prog_mon,
    size_t &fr_total_u,
    BoundedThreadPool *pool
)
{
    int idxx = 0;
    // process newly finished lines, if there are any
    if ((int)(prog_mon->fr_count / nx) < *preprocessor_line)
    {
        *processor_line = (int)(prog_mon->fr_count) / nx;
        if (*processor_line%ny==0)
            id_image = *processor_line / ny % 2;
        idxx = (int)(prog_mon->fr_count) % nxy;
        *prog_mon += nx;

        int update_line = idxx / nx;
        if ((prog_mon->report_set) && (update_line)>0)
        {
            fr_freq = prog_mon->fr_freq;
            prog_mon->reset_flags();
        }
    }

    // end of image handler
    if (prog_mon->fr_count >= end_frame)
    {
        if (b_continuous) {
            prog_mon->fr_total += nxy;
            fr_total_u += nxy;
        }
        if (prog_mon->fr_count != fr_total_u)
        {
            img_num++;
            first_frame = img_num * nxy;
            end_frame = (img_num + 1) * nxy;
        }
    }

    // end of recon handler
    if (((prog_mon->fr_count >= fr_total_u) && (!b_continuous)) || rc_quit)
    {
        pool->wait_for_completion();
        p_prog_mon = nullptr;
        b_cumulative = false;
        b_continuous = false;
        *processor_line = -1;
    }
}

void Electron::close()
{
    dtype_Electron electron;
    electron.kx = 0;
    electron.ky = 0;
    electron.rx = 0;
    electron.ry = 0;
    electron.id_image = rep+1;
    file_electron.write((const char *)&electron, sizeof(electron));

    closeDatFile();
}
