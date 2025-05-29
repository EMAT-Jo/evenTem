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

#include "Roi.h"

void Roi::run(){
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
            if (b_ROI_4D)
            {
                if (bitdepth == 8) cam.enable_roi_4D(Roi_4D_8,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 16) cam.enable_roi_4D(Roi_4D_16,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 32) cam.enable_roi_4D(Roi_4D_32,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
            }
            else if (use_mask)
            {
                cam.enable_roi_mask(&roi_mask,&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern);
            }
            else
            {
                cam.enable_roi(&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right);
            }
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
            if (b_ROI_4D)
            {
                if (bitdepth == 8) cam.enable_roi_4D(Roi_4D_8,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 16) cam.enable_roi_4D(Roi_4D_16,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 32) cam.enable_roi_4D(Roi_4D_32,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
            }
            else if (use_mask)
            {
                cam.enable_roi_mask(&roi_mask,&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern);
            }
            else
            {
                cam.enable_roi(&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right);
            }
            cam.run();
            process_data();
            cam.terminate();
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
            if (b_ROI_4D)
            {
                if (bitdepth == 8) cam.enable_roi_4D(Roi_4D_8,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 16) cam.enable_roi_4D(Roi_4D_16,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 32) cam.enable_roi_4D(Roi_4D_32,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
            }
            else if (use_mask)
            {
                cam.enable_roi_mask(&roi_mask,&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern);
            }
            else
            {
                cam.enable_roi(&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right);
            }
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
            if (b_ROI_4D)
            {
                if (bitdepth == 8) cam.enable_roi_4D(Roi_4D_8,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 16) cam.enable_roi_4D(Roi_4D_16,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                else if (bitdepth == 32) cam.enable_roi_4D(Roi_4D_32,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
            }
            else if (use_mask)
            {
                cam.enable_roi_mask(&roi_mask,&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern);
            }
            else
            {
                cam.enable_roi(&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right);
            }
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
                if (b_ROI_4D)
                {
                    // cam.enable_roi_4D(&Roi_4D,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                }
                else
                {
                    cam.enable_roi(&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right);
                }
                cam.run();
                process_data();
                cam.terminate();
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
                if (b_ROI_4D)
                {
                    // cam.enable_roi_4D(&Roi_4D,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
                }
                else
                {
                    cam.enable_roi(&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right);
                }
                cam.run();
                process_data();
                cam.terminate();
            }
            else std::runtime_error("Invalid detector size for Merlin");
            break;
        }
        // case CAMERA::NUMPY:
        // {
        //     if (n_cam == 512) using namespace FRAME_512_ADDITIONAL;
        //     else if (n_cam == 256) using namespace FRAME_256_ADDITIONAL;
        //     else if (n_cam == 128) using namespace FRAME_128_ADDITIONAL;
        //     else if (n_cam == 64) using namespace FRAME_64_ADDITIONAL;
        //     else {
        //         using namespace FRAME_64_ADDITIONAL;
        //         std::runtime_error("Invalid frame size, only 64, 128, 256, 512 are supported");
        //         }
        //     NUMPY<N_CAM,BUFFER_SIZE,HEAD_SIZE,N_BUFFER> cam(
        //         nx, 
        //         ny, 
        //         &b_cumulative,
        //         rep,
        //         processor_line,
        //         preprocessor_line,
        //         mode,
        //         file_path,
        //         socket
        //     );
        //     if (b_ROI_4D)
        //     {
        //         cam.enable_roi_4D(&Roi_4D,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right,det_bin);
        //     }
        //     else
        //     {
        //         cam.enable_roi(&Roi_scan_image_stack,&Roi_diffraction_pattern_stack,&Roi_scan_image,&Roi_diffraction_pattern, lower_left, upper_right);
        //     }
        //     cam.run();
        //     process_data();
        //     cam.terminate();
        //     break;
        // }
    }

    rc_quit = true;
 }


void Roi::set_roi(int x , int y, int width, int height)
{
    lower_left[0] = x;
    lower_left[1] = ny - (y+height); // y is from top to bottom, but the image is from bottom to top, so we need to flip
    upper_right[0] = x + width;
    upper_right[1] = ny-y;

    if (x < 0 || y < 0)
    {
        throw std::invalid_argument("ROI x & y parameters must be positive");
    }
    if (width <= 0 && height <= 0)
    {
        throw std::invalid_argument("ROI width or height parameters must be > 0");
    }
    if (x + width > nx || y + height > ny)
    {
        throw std::invalid_argument("ROI must be within the image dimensions");
    }

    L_0 = upper_right[0] - lower_left[0];
    L_1 = upper_right[1] - lower_left[1];

    finish_line = (y+height+1) + ny*(rep-1);
}

std::array<int,4> Roi::get_roi()
{
    std::array<int,4> roi = {lower_left[0], lower_left[1], upper_right[0], upper_right[1]};
    return roi;
}

void Roi::set_roi_mask(std::vector<py::array_t<int>> arrays){
    for (auto& array : arrays) {
        py::buffer_info buf_info = array.request();
        if (buf_info.ndim != 1) throw std::runtime_error("Input should be a 1-D array");
        int* data_ptr = static_cast<int*>(buf_info.ptr);
        size_t size = buf_info.size;
        roi_mask.push_back(std::vector<int>(data_ptr, data_ptr + size));
    }
    roi_mask.push_back(std::vector<int>(nx*ny,0));
    use_mask = true;
    L_1 = ny;
    L_0 = nx;
    finish_line = ny*rep;
};

void Roi::set_bitdepth(int _bitdepth)
{
    if (_bitdepth == 8) bitdepth = 8;
    else if (_bitdepth == 16) bitdepth = 16;
    else if (_bitdepth == 32) bitdepth = 32;
    else
    {
        throw std::invalid_argument("Invalid bitdepth, only 8 , 16 or 32 are supported");
    }
}

void Roi::reset()
{
    rc_quit = false;
    fr_freq = 0;

    // Initializations
    nxy = nx * ny;
    id_image = 0;
    fr_total = nxy * rep;
    fr_count = 0;

    Roi_diffraction_pattern.assign(n_cam*n_cam, 0);
    Roi_diffraction_pattern_stack.assign(rep+1, std::vector<uint64_t>(n_cam*n_cam, 0));
    Roi_scan_image.assign(L_1*L_0,0);
    Roi_scan_image_stack.assign(rep+1,std::vector<uint64_t>(L_1*L_0,0));

    if (b_ROI_4D)
    {
        if (bitdepth == 8 ) Roi_4D_8->init(L_1,L_0,n_cam,det_bin);
        else if (bitdepth == 16)  Roi_4D_16->init(L_1,L_0,n_cam,det_bin);
        else if (bitdepth == 32)  Roi_4D_32->init(L_1,L_0,n_cam,det_bin);
    }

    // Data Processing Progress
    *processor_line = 0;
    *preprocessor_line = 0;

}


void Roi::line_processor(
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


        // end of line handler
        int update_line = idxx / nx; //- kernel.kernel_size*2;
        if ((prog_mon->report_set) && (update_line)>0)
        {
            fr_freq = prog_mon->fr_freq;
            prog_mon->reset_flags();
        }

        if (*processor_line > finish_line){
            rc_quit = true;
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

