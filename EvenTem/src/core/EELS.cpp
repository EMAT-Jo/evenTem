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

 
#include "EELS.h"

void EELS::run(){
    py::gil_scoped_release release;
    reset();
    // Run camera dependent pipeline
    switch (camera)
     {
        case CAMERA::MERLIN:
        {
            if (n_cam == 512){ 
                using namespace MERLIN_512;
                MERLIN<N_CAM,BUFFER_SIZE,HEAD_SIZE,N_BUFFER,PIXEL> cam(
                    nx, 
                    ny, 
                    &b_cumulative,
                    rep,
                    processor_line,
                    preprocessor_line,
                    mode,
                    file_path,
                    socket
                );
                cam.enable_EELS(&EELS_data_stack, &EELS_image_stack);
                cam.run();
                process_data();
                cam.terminate();
                break;
            }
            else if (n_cam == 256){
                using namespace MERLIN_256;
                MERLIN<N_CAM,BUFFER_SIZE,HEAD_SIZE,N_BUFFER,PIXEL> cam(
                    nx, 
                    ny, 
                    &b_cumulative,
                    rep,
                    processor_line,
                    preprocessor_line,
                    mode,
                    file_path,
                    socket
                );
                cam.enable_EELS(&EELS_data_stack, &EELS_image_stack);
                cam.run();
                process_data();
                cam.terminate();
                break;
            }
            else std::runtime_error("Invalid detector size for Merlin");
            break;
        }
    }
    rc_quit = true;
}

void EELS::reset()
{
    rc_quit = false;
    fr_freq = 0;

    // Initializations
    nxy = nx * ny;
    id_image = 0;
    fr_total = nxy * rep;
    fr_count = 0;

    EELS_data.assign(nxy, std::vector<uint64_t>(n_cam, 0));
    EELS_data_stack.assign(rep+1, std::vector<std::vector<uint64_t>>(nxy, std::vector<uint64_t>(n_cam, 0)));

    EELS_image.assign(nxy, 0);
    EELS_image_stack.assign(rep+1, std::vector<uint64_t>(nxy, 0));

    // Data Processing Progress
    *processor_line = 0;
    *preprocessor_line = 0;
}


void EELS::line_processor(
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
        if (*processor_line%ny==0) id_image = *processor_line / ny;
        idxx = (int)(prog_mon->fr_count) % nxy;
        *prog_mon += nx;
    }

    // end of line handler
    int update_line = idxx / nx;
    if ((prog_mon->report_set) && (update_line)>0)
    {
        fr_freq = prog_mon->fr_freq;
        prog_mon->reset_flags();
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
