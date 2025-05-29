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

#include "Var.h"

void Var::run(){
     py::gil_scoped_release release;
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
            cam.enable_var(&Var_data, offset);
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
            cam.enable_var(&Var_data, offset);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
        case CAMERA::SIMULATED:
        {
            using namespace SIMULATED_ADDITIONAL;
            SIMULATED<EVENT, BUFFER_SIZE, N_BUFFER> cam(
                nx, 
                ny, 
                n_cam,
                &b_cumulative,
                rep,
                processor_line,
                preprocessor_line,
                mode,
                file_path,
                socket
            );
            cam.enable_var(&Var_data, offset);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
     }

    rc_quit = true;
 }

void Var::set_offset(std::array<float, 2> _offset)
{
    auto_offset = false;
    offset = _offset;
}



void Var::reset()
{
    rc_quit = false;
    fr_freq = 0;

    std::fill(Var_image.begin(), Var_image.end(), 0);

    // Initializations
    nxy = nx * ny;
    id_image = 0;
    fr_total = nxy * rep;
    fr_count = 0;

    if (auto_offset) offset = {(float)n_cam/2, (float)n_cam/2};

    // Allocate memory for image arrays
    Var_image.assign(nxy, 0);
    for (int i=0; i<2; i++)
    {
        Var_data[i].assign(nxy, 0);
    }

    // Data Processing Progress
    *processor_line = 0;
    *preprocessor_line = 0;

}


void Var::line_processor(
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

        for (size_t i = 0; i < (size_t)nx; i++)
        {
            int idxx_p_i = idxx + i;
            if ((idxx_p_i >= 0) | (nx > 1))
            {
            Var_image[idxx_p_i] = (float)Var_data[id_image][idxx_p_i];
            }
        // end of line handler
        int update_line = idxx / nx; //- kernel.kernel_size*2;
        if ((prog_mon->report_set) && (update_line)>0)
        {
            fr_freq = prog_mon->fr_freq;
            //rescales_recomputes();
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
}

