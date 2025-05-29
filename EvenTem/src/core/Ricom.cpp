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

#include "Ricom.h"

// Compute the kernel
void Ricom_kernel::compute_kernel()
{
    float rot_rad = M_PI * rotation / 180;
    float cos_rot = cos(rot_rad);
    float sin_rot = sin(rot_rad);

    k_width_sym = kernel_size * 2 + 1;
    k_area = k_width_sym * k_width_sym;
    kernel_x.assign(k_area, 0);
    kernel_y.assign(k_area, 0);

    for (int iy = 0; iy < k_width_sym; iy++)
    {
        int iy_e = (iy + 1) * k_width_sym - 1;
        for (int ix = 0; ix < k_width_sym; ix++)
        {
            int ix_s = ix - kernel_size;
            int iy_s = iy - kernel_size;
            float d = ix_s * ix_s + iy_s * iy_s;
            int ix_e = k_area - iy_e + ix - 1;

            if (d > 0)
            {
                float ix_sd = (ix_s / d);
                float iy_sd = (iy_s / d);
                kernel_x[ix_e] = cos_rot * ix_sd - sin_rot * iy_sd;
                kernel_y[ix_e] = sin_rot * ix_sd + cos_rot * iy_sd;
            }
            else
            {
                kernel_y[ix_e] = 0;
                kernel_x[ix_e] = 0;
            }
        }
    }

    // Add filter
    if (b_filter)
    {
        compute_filter();
        include_filter();
    }
}

// Compute the filter
void Ricom_kernel::compute_filter()
{
    // kernel_filter.assign(k_area, 0);
    // float lb = pow(kernel_filter_frequency[0], 2);
    // float ub = pow(kernel_filter_frequency[1], 2);

    // for (int iy = 0; iy < k_width_sym; iy++)
    // {
    //     for (int ix = 0; ix < k_width_sym; ix++)
    //     {
    //         float dist = pow(ix - kernel_size, 2) + pow(iy - kernel_size, 2);
    //         int ic = iy * k_width_sym + ix;
    //         if (dist <= ub && dist > lb)
    //         {
    //             kernel_filter[ic] = 1;
    //         }
    //     }
    // }
}

// Applies the filter to the kernel
void Ricom_kernel::include_filter()
{
    
    // FFT2D fft2d(k_width_sym, k_width_sym);
    // std::vector<std::complex<float>> x2c = FFT2D::r2c(kernel_x);
    // std::vector<std::complex<float>> y2c = FFT2D::r2c(kernel_y);
    // fft2d.fft(x2c, x2c);
    // fft2d.fft(y2c, y2c);
    // for (int id = 0; id < k_area; id++)
    // {
    //     if (kernel_filter[id] == 0.0f)
    //     {
    //         x2c[id] = {0, 0};
    //         y2c[id] = {0, 0};
    //     }
    // }
    // fft2d.ifft(x2c, x2c);
    // fft2d.ifft(y2c, y2c);
    // for (int id = 0; id < k_area; id++)
    // {
    //     kernel_x[id] = x2c[id].real();
    //     kernel_y[id] = y2c[id].real();
    // }
    
}

void Ricom_kernel::approximate_frequencies(size_t nx_im)
{
    f_approx.resize(nx_im);
    float f_max = 0;
    float k = kernel_size * 2;
    for (size_t i = 0; i < nx_im; i++)
    {
        float x = 2 * i * M_PI;
        f_approx[i] = (nx_im / x) * (1 - cos(k / 2 * (x / nx_im)));
        if (f_approx[i] > f_max)
        {
            f_max = f_approx[i];
        }
    }
    std::for_each(f_approx.begin(), f_approx.end(), [f_max](float &x)
                  { x /= f_max; });
}


// Entrance function for Ricom_reconstructinon
void Ricom::run()
{
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
            cam.enable_Ricom(&dose_data,&sumx_data,&sumy_data);
            cam.run();
            startime = std::chrono::high_resolution_clock::now();
            process_data();
            cam.terminate();
            processing_rate = cam.get_processing_rate();
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
            cam.enable_Ricom(&dose_data,&sumx_data,&sumy_data);
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
            cam.enable_Ricom(&dose_data,&sumx_data,&sumy_data);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
        case CAMERA::MERLIN:
        {
            if (n_cam == 512) {
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
                cam.enable_Ricom(&comx_image,&comy_image);

                // for (vSTEM *child : vSTEM_children){
                //     setup_child(child);
                //     cam.enable_vSTEM(&(child->get_detector()),&child->vSTEM_data,false,false);
                // }

                cam.run();
                std::thread t_self = std::thread(&Ricom::process_data, this);

                // for (size_t i = 0; i < vSTEM_children.size(); ++i) {threads_for_children[i] = std::thread(&vSTEM::process_data,vSTEM_children[i]);}

                if (t_self.joinable()) t_self.join();
                // for (auto &t : threads_for_children) {if (t && t->joinable()) {t->join();};}

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
                cam.enable_Ricom(&comx_image,&comy_image);

                // for (vSTEM *child : vSTEM_children){
                //     setup_child(child);
                //     cam.enable_vSTEM(&(child->get_detector()),&child->vSTEM_data,false,false);
                // }

                cam.run();
                std::thread t_self = std::thread(&Ricom::process_data, this);

                // for (size_t i = 0; i < vSTEM_children.size(); ++i) {threads_for_children[i] = std::thread(&vSTEM::process_data,vSTEM_children[i]);}

                if (t_self.joinable()) t_self.join();
                // for (auto &t : threads_for_children) {if (t && t->joinable()) {t->join();};}

                cam.terminate();   
                }
            else std::runtime_error("Invalid detector size for Merlin");
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
            cam.enable_Ricom(&comx_image,&comy_image);
            cam.run();
            process_data();
            cam.terminate();
            break;
        }
    }
    rc_quit = true;
}

// Rescales the images according to updated min and max values
// and recomputes the Kernel if settings changed
// inline void Ricom::rescales_recomputes()
// {
//     if (b_recompute_kernel)
//     {
//         kernel.compute_kernel();
//         b_recompute_kernel = false;
//     }
// }

// Compute electric field magnitude
void Ricom::compute_electric_field(std::array<float, 2> &com_xy, size_t id)
{
    float e_mag = std::hypot(com_xy[0] - offset[0], com_xy[1] - offset[1]);
    float e_ang = atan2(com_xy[0] - offset[0], com_xy[1] - offset[1]);
    if (e_mag > e_mag_max)
    {
        e_mag_max = e_mag;
        rescale_e_mag = true;
    }
    if (e_mag < e_mag_min)
    {
        e_mag_min = e_mag;
        rescale_e_mag = true;
    }
    e_field_data[id] = std::polar(e_mag, e_ang);
}

void Ricom::line_processor(
    size_t &img_num,
    size_t &first_frame,
    size_t &end_frame,
    ProgressMonitor *prog_mon,
    size_t &fr_total_u,
    BoundedThreadPool *pool
)
{
    int pp_id = 0;
    if ((int)(prog_mon->fr_count / nx) < *preprocessor_line)
    {
        *processor_line = (int)(prog_mon->fr_count) / nx;
        if (*processor_line%ny==0)
            id_image = *processor_line / ny % 2;
        pp_id = (int)(prog_mon->fr_count) % nxy;
        *prog_mon += nx;

        std::array<float, 2> com_xy = {0.0, 0.0};
        std::array<float, 2> com_xy_sum = {0.0, 0.0};
        for (size_t i = 0; i < (size_t)nx; i++)
        {
            if (event_mode){
                int idxx_p_i = pp_id + i;
                if ((idxx_p_i >= 0) | (nx > 1))
                {
                    if (dose_data[id_image][idxx_p_i] == 0)
                    {
                        comx_image[idxx_p_i] = offset[0];
                        comy_image[idxx_p_i] = offset[1];
                    }
                    else
                    {
                        comx_image[idxx_p_i] = sumx_data[id_image][idxx_p_i] / dose_data[id_image][idxx_p_i];
                        comy_image[idxx_p_i] = sumy_data[id_image][idxx_p_i] / dose_data[id_image][idxx_p_i];
                    }
                    com_xy_sum[0] += comx_image[idxx_p_i];
                    com_xy_sum[1] += comy_image[idxx_p_i];
                }
            }
            else{
                int idxx_p_i = pp_id + i;
                if ((idxx_p_i >= 0) | (nx > 1))
                {
                    com_xy_sum[0] += comx_image[idxx_p_i];
                    com_xy_sum[1] += comy_image[idxx_p_i];
                }
            }
        }

        if (n_threads > 1)
        {
                pool->push_task([=]{ icom_group_classical(pp_id, *processor_line / ny); });
        }
        else
        {
                icom_group_classical(pp_id,*processor_line / ny);
        }

        // end of line handler
        int update_line = pp_id / nx - kernel.kernel_size*2;
        if ((prog_mon->report_set) && (update_line)>0)
        {
            fr_freq = prog_mon->fr_freq;
            // rescales_recomputes(); 
            for (int i = 0; i < 2; i++)
            {
                com_public[i] = com_xy_sum[i] / nx;
                com_xy_sum[i] = 0;
            }
            prog_mon->reset_flags();

            current_time = std::chrono::high_resolution_clock::now();
            elapsed_seconds = current_time - startime;
            this->reached_pp_id.push_back(pp_id);
            this->elapsed_seconds_vec.push_back(elapsed_seconds.count());
        }
    }

    // end of image handler
    if (prog_mon->fr_count >= end_frame)
    {
        if (b_continuous) {
            prog_mon->fr_total += nxy;
            fr_total_u += nxy;
            ricom_image.assign(nxy, 0);
        }
        if (prog_mon->fr_count != fr_total_u)
        {
            img_num++;
            first_frame = img_num * nxy;
            end_frame = (img_num + 1) * nxy;
        }
        if (update_offset)
        {
            offset[0] = com_public[0];
            offset[1] = com_public[1];
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



void Ricom::icom_group_classical(int pp_id, int id_image)
{
    int idk;
    int idc;
    int id;
    int idr_x;
    int idr_delay;
    int k_bias = kernel.k_width_sym * (kernel.kernel_size + 1);

    if (((pp_id / nx - 2*kernel.kernel_size) >= 0))
    {
        for (int iy = -kernel.kernel_size; iy <= kernel.kernel_size; iy++)
        {
            for (int i_line = 0, idr_delay = pp_id - kernel.kernel_size * nx;
                i_line < nx;
                i_line++, idr_delay++)
            {
                idc = idr_delay + iy * nx;
                idr_x = idr_delay % nx;
                idk = (kernel.kernel_size + iy) * kernel.k_width_sym;
                for (int ix = -kernel.kernel_size; ix <= kernel.kernel_size; ix++)
                {
                    if (((idr_x + ix) >= 0) & ((idr_x + ix) < nx))
                    {
                        ricom_image[idr_delay] += ((comx_image[idc + ix] - offset[0]) * -kernel.kernel_x[idk] + (comy_image[idc + ix] - offset[1]) * -kernel.kernel_y[idk]);
                        ricom_image_stack[id_image][idr_delay] += ((comx_image[idc + ix] - offset[0]) * -kernel.kernel_x[idk] + (comy_image[idc + ix] - offset[1]) * -kernel.kernel_y[idk]);
                    }
                    ++idk;
                }
            }
        }
        fr_count = pp_id;
    }

}

void Ricom::reset()
{
    rc_quit = false;
    fr_freq = 0;

    std::fill(ricom_image.begin(), ricom_image.end(), 0);
    std::fill(comx_image.begin(), comx_image.end(), 0);
    std::fill(comy_image.begin(), comy_image.end(), 0);

    // Initializations
    nxy = nx * ny;
    id_image = 0;
    fr_total = nxy * rep;
    fr_count = 0;

    kernel.compute_kernel();

    if (auto_offset) offset = {(float)n_cam/2,(float)n_cam/2};

    // Allocate memory for image arrays
    ricom_image.assign(nxy, 0);
    ricom_image_stack.assign(rep+1, std::vector<float>(nxy, 0));
    comx_image.assign(nxy, 0);
    comy_image.assign(nxy, 0);
    for (int i=0; i<2; i++)
    {
        dose_data[i].assign(nxy, 0);
        sumx_data[i].assign(nxy, 0);
        sumy_data[i].assign(nxy, 0);
    }

    // Data Processing Progress
    *processor_line = 0;
    *preprocessor_line = 0;

}

void Ricom::set_kernel(int kernel_size, int rotation)
{
    kernel.kernel_size = kernel_size;
    kernel.rotation = rotation;
}

void Ricom::set_kernel_filtered(int kernel_size, int rotation, float low_pass, float high_pass)
{
    kernel.kernel_size = kernel_size;
    kernel.rotation = rotation;
    kernel.b_filter = true;
    kernel.kernel_filter_frequency[0] = low_pass;
    kernel.kernel_filter_frequency[1] = high_pass;
}

void Ricom::set_offset(std::array<float, 2> _offset)
{
    auto_offset = false;
    offset = _offset;
}

std::vector<float> Ricom::get_kernel()
{
    kernel.compute_kernel();
    return kernel.kernel_x;
}

// void Ricom::set_masked_com(int radius, std::array<float, 2> offset)
// {
//     masked_com = true;
//     detector.set_radia(0, radius);
//     detector.compute_detector(n_cam, n_cam, offsets);
//     com_mask = detector.detector_image;
// }

void Ricom::add_child(vSTEM* child)
{
    vSTEM_children.push_back(child);
    threads_for_children.emplace_back(std::nullopt);
}

void Ricom::setup_child(vSTEM *child)
{
    child->camera = camera;
    child->nx = nx;
    child->ny = ny;
    child->b_cumulative = b_cumulative;
    child->reset();
    child->processor_line = processor_line;
    child->preprocessor_line = preprocessor_line;
    child->progress_verbose = false;
}
