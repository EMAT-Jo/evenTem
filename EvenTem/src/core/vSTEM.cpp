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

#include "vSTEM.h"

void vSTEM::run(){
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
            if (use_mask) cam.enable_mask_vSTEM(&detector_mask,&vSTEM_stack);
            else if (detector.n_detectors > 1)
            {
                cam.enable_multi_vSTEM(&detector.radia_sqr,&offsets,&vSTEM_stack);
            }
            else  cam.enable_vSTEM(&detector.radia_sqr[0],&offsets[0],&vSTEM_stack);
            // else  cam.enable_atomic_vSTEM(&detector.radia_sqr[0],&offsets[0],&atomic_vSTEM_image);
            cam.run();
            startime = std::chrono::high_resolution_clock::now();
            process_data();
            cam.terminate();
            processing_rate = cam.get_processing_rate();
            // from_atomic();
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

            if (detector.n_detectors > 1)
            {
                cam.enable_multi_vSTEM(&detector.radia_sqr,&offsets,&vSTEM_stack);
            }
            else  cam.enable_vSTEM(&detector.radia_sqr[0],&offsets[0],&vSTEM_stack);
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
            if (detector.n_detectors > 1)
            {
                cam.enable_multi_vSTEM(&detector.radia_sqr,&offsets,&vSTEM_stack);
            }
            else  cam.enable_vSTEM(&detector.radia_sqr[0],&offsets[0],&vSTEM_stack);
            // else  cam.enable_atomic_vSTEM(&detector.radia_sqr[0],&offsets[0],&atomic_vSTEM_image);
            cam.run();
            process_data();
            cam.terminate();
            // from_atomic();
            break;
        }
        case CAMERA::CHEETAH_PIXELTRIG:
        {
            using namespace CHEETAH_ADDITIONAL;
            CHEETAH_pixeltrig<EVENT, BUFFER_SIZE, N_BUFFER> cam(
                nx, 
                ny, 
                &b_cumulative,
                rep,
                processor_line,
                preprocessor_line,
                mode,
                file_path,
                socket,
                pattern_file
            );
            if (detector.n_detectors > 1)
            {
                cam.enable_multi_vSTEM(&detector.radia_sqr,&offsets,&vSTEM_stack);
            }
            else  cam.enable_vSTEM(&detector.radia_sqr[0],&offsets[0],&vSTEM_stack);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
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
                compute_detector();
                cam.enable_vSTEM(&detector.detector_image,&vSTEM_stack, allow_torch);
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
                compute_detector();
                cam.enable_vSTEM(&detector.detector_image,&vSTEM_stack, allow_torch);
                cam.run();
                process_data();
                cam.terminate();
                break;
            }
            else std::runtime_error("Invalid detector size for Merlin");
            break;
        }
        case CAMERA::HDF5:
        {
            using namespace HDF5_ADDITIONAL;
            HDF5<N_CAM,BUFFER_SIZE,HEAD_SIZE,N_BUFFER,PIXEL> cam(
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
            compute_detector();
            cam.enable_vSTEM(&detector.detector_image,&vSTEM_stack, allow_torch);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
        case CAMERA::NUMPY:
        {
            using namespace FRAME_64_ADDITIONAL;
            NUMPY<N_CAM,BUFFER_SIZE,HEAD_SIZE,N_BUFFER,PIXEL> cam(
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
            compute_detector();
            cam.enable_vSTEM(&detector.detector_image,&vSTEM_stack, allow_torch);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
    }
    rc_quit = true;
}

void vSTEM::reset()
{
    rc_quit = false;
    fr_freq = 0;

    std::fill(vSTEM_image.begin(), vSTEM_image.end(), 0);

    // Initializations
    nxy = nx * ny;
    id_image = 0;
    fr_total = nxy * rep;
    fr_count = 0;

    if (auto_offset) offsets = {{(float)n_cam/2, (float)n_cam/2}};

    detector.set_radia(inner_radia, outer_radia);

    vSTEM_image.assign(nxy, 0);
    vSTEM_stack.assign(rep+1, std::vector<size_t>(nxy, 0));
    for (int i = 0; i < nxy; i++)
    {
        atomic_vSTEM_image.push_back(std::atomic<int>(0));
    }

    // Data Processing Progress
    *processor_line = 0;
    *preprocessor_line = 0;
}


void vSTEM::line_processor(
    size_t &img_num,
    size_t &first_frame,
    size_t &end_frame,
    ProgressMonitor *prog_mon,
    size_t &fr_total_u,
    BoundedThreadPool *pool
)
{
    int pp_id = 0;
    // process newly finished lines, if there are any
    if ((int)(prog_mon->fr_count / nx) < *preprocessor_line) 
    {
        *processor_line = (int)(prog_mon->fr_count) / nx;
        if (*processor_line%ny==0) id_image = *processor_line / ny;
        pp_id = (int)(prog_mon->fr_count) % nxy;

        for (size_t i = 0; i < (size_t)nx; i++)
        {
            int idxx_p_i = pp_id + i;
            if (b_cumulative) vSTEM_image[idxx_p_i] += vSTEM_stack[id_image][idxx_p_i];
            if (b_continuous) vSTEM_image[idxx_p_i] = vSTEM_stack[id_image][idxx_p_i];

        }

        *prog_mon += nx;
    }

    // end of line handler
    int update_line = pp_id / nx;
    if ((prog_mon->report_set) && (update_line)>0)
    {
        fr_freq = prog_mon->fr_freq;
        prog_mon->reset_flags();

        current_time = std::chrono::high_resolution_clock::now();
        elapsed_seconds = current_time - startime;
        this->reached_pp_id.push_back(pp_id);
        this->elapsed_seconds_vec.push_back(elapsed_seconds.count());
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

std::vector<int> vSTEM::compute_detector(){
    if (auto_offset) offsets = {{(float)n_cam/2, (float)n_cam/2}};
    detector.set_radia(inner_radia, outer_radia);
    detector.compute_detector(n_cam, n_cam, offsets);
    return detector.detector_image;
}

std::vector<int> vSTEM::get_detector(){
    if (!use_mask) return compute_detector();
    else return detector_mask;
}

void vSTEM::set_detector_mask(py::array_t<int> array){
    py::buffer_info buf_info = array.request();
    if (buf_info.ndim != 1) throw std::runtime_error("Input should be a 1-D array");
    int* data_ptr = static_cast<int*>(buf_info.ptr);
    size_t size = buf_info.size;
    detector_mask = std::vector<int>(data_ptr, data_ptr + size);
    use_mask = true;
};

void vSTEM::set_offsets(std::vector<std::array<float, 2>> _offsets)
{
    auto_offset = false;
    offsets = _offsets;
}

void vSTEM::from_atomic(){
    for (int i = 0; i < nxy; i++)
    {
        vSTEM_image[i] = atomic_vSTEM_image[i].load();
    }
}