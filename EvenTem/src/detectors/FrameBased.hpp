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

#ifndef FRAMEBASED_H
#define FRAMEBASED_H

#ifdef __GNUC__
#define PACK(__Declaration__) __Declaration__ __attribute__((__packed__))
#endif

#ifdef _MSC_VER
#define PACK(__Declaration__) __pragma(pack(push, 1)) __Declaration__ __pragma(pack(pop))
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <atomic>
#include <vector>
#include <array>
#include <thread>
#include <chrono>
#include <functional>
#include <algorithm>
#include <sstream>

#include "SocketConnector.h"
#include "FileConnector.h"
#include "BoundedThreadPool.hpp"

#if defined(GPRI_OPTION_ENABLED) || defined(FRAMEBASED_TORCH_ENABLED)
    #include <torch/torch.h>
    #include <cuda_runtime_api.h>
#endif


namespace FRAME_64_ADDITIONAL{
    const int N_CAM = 64;
    const int BUFFER_SIZE = 128;
    const int HEAD_SIZE = 768;
    const int N_BUFFER = 32;
    using PIXEL = uint8_t;
}; 

namespace FRAME_128_ADDITIONAL{
    const int N_CAM = 128;
    const int BUFFER_SIZE = 128;
    const int HEAD_SIZE = 768;
    const int N_BUFFER = 32;
    using PIXEL = uint8_t;
}; 

namespace FRAME_256_ADDITIONAL{
    const int N_CAM = 256;
    const int BUFFER_SIZE = 128;
    const int HEAD_SIZE = 768;
    const int N_BUFFER = 32;
    using PIXEL = uint8_t;
}; 

namespace FRAME_512_ADDITIONAL{
    const int N_CAM = 512;
    const int BUFFER_SIZE = 128;
    const int HEAD_SIZE = 768;
    const int N_BUFFER = 32;
    using PIXEL = uint8_t;
};


template <int n_cam,int buffer_size, int HEAD_SIZE, int n_buffer, typename pixel>
class FRAMEBASED
{
protected:
    // ----------------------------------------------------------------------------------------------- 
    // process methods
    // -----------------------------------------------------------------------------------------------
    #ifdef FRAMEBASED_TORCH_ENABLED
    inline void vstem_torch()
    {
        (*p_stem_data)[0][probe_position] += (torch::mul(Tensor_buffer[buffer_id].index({frame_id}),DetectorTensor)).sum().item().to<uint64_t>();
    };
    #endif
    inline void vstem()
    {
        for (int active_pixel : detector_list) 
        {
            (*p_stem_data)[0][probe_position] += (size_t)frame_buffer[buffer_id][frame_id][active_pixel];
        };
    };

    void pacbed()
    {
        for (int k = 0; k < n_cam*n_cam; k++)
        {
            (*p_pacbed_data)[k] += (size_t)frame_buffer[buffer_id][frame_id][k];
        };
    };

    void com()
    {
        sum_x.assign(n_cam, 0);
        sum_y.assign(n_cam, 0);
        dose = 0;
        COM[0] = 0;
        COM[1] = 0;

        for (int idy = 0; idy < n_cam; idy++)
        {
            sum_x_temp = 0;
            for (int idx = 0; idx < n_cam; idx++)
            {
                current_px = (size_t)frame_buffer[buffer_id][frame_id][idy * n_cam + idx];
                sum_x_temp += current_px;
                sum_y[idx] += current_px;
            }
            sum_x[idy] = sum_x_temp;
            dose += sum_x_temp;
        }
    

        if (dose > 0)
        {
            for (int i = 0; i < n_cam; i++)
            {
                COM[0] += sum_x[i] * v[i];
            }
            for (int i = 0; i < n_cam; i++)
            {
                COM[1] += sum_y[i] * u[i];
            }
            (*p_comx_image)[probe_position] = COM[0] / dose;
            (*p_comy_image)[probe_position] = COM[1] / dose;
        }
    }

    #ifdef GPRI_OPTION_ENABLED
    void GPRI()
    {

    };

    void SparseGPRI()
    {
        for (int k = 0; k < n_cam*n_cam; k++){
            if (frame_buffer[buffer_id][frame_id][k] != 0){
                int _kx = k/n_cam;
                int _ky = k%n_cam;
                
                (*p_k_indices_vec)[probe_position+id_image*GPRI_nxy_scan_bin]->push_back((_kx/GPRI_detector_bin)*GPRI_cam_bin+(_ky/GPRI_detector_bin));
                (*p_N_electrons_map_scangrid)[id_image][probe_position] += 1;
            }
        }
    };
    #endif

    void roi()
    {
        int x = probe_position%nx;
        int y = nx - floor(probe_position/nx);
        if (x >= lower_left[0] && x < upper_right[0] && y > lower_left[1] && y <= upper_right[1])
        {
            for (int k = 0; k < n_cam*n_cam; k++){
                (*p_roi_diffraction_pattern_stack)[id_image][k] += frame_buffer[buffer_id][frame_id][k];
                (*p_roi_scan_image_stack)[id_image][(L_1 - (y-lower_left[1])) * L_0 + (x-lower_left[0])] += frame_buffer[buffer_id][frame_id][k];
                (*p_roi_diffraction_pattern)[k] += frame_buffer[buffer_id][frame_id][k];
                (*p_roi_scan_image)[(L_1 - (y-lower_left[1])) * L_0 + (x-lower_left[0])] += frame_buffer[buffer_id][frame_id][k];
            }
        }
        
    };
    void roi_4D()
    {
    };

    void EELS()
    {
        for (int k = 0; k < n_cam*n_cam; k++)
        {
            (*p_EELS_data_stack)[id_image][probe_position][k%n_cam] += (size_t)frame_buffer[buffer_id][frame_id][k];
            (*p_EELS_image_stack)[id_image][probe_position] += (size_t)frame_buffer[buffer_id][frame_id][k];
        }
    };

    void compress()
    {
        if (probe_position < nxy)
        {
            int _id_chunk = (probe_position/(chunksize*nx))%2;
            std::lock_guard<std::mutex> lock(mtx[_id_chunk]);

            #ifdef FRAMEBASED_TORCH_ENABLED
                temp_Tensor = torch::from_blob(frame_buffer[buffer_id][frame_id].data(), {n_cam, n_cam}, torch::TensorOptions().dtype(torch::kUInt8)).view({n_cam_bin,compress_det_bin,n_cam_bin,compress_det_bin}).sum({1, 3}, /*keepdim=*/false).to(torch::kUInt8);
                (*p_counts_data)[probe_position] += temp_Tensor.sum().item().to<uint64_t>();
                int s = probe_position%(chunksize*nx)*diff_pattern_size;
                auto tensor_data_ptr = temp_Tensor.data_ptr<uint8_t>();
                std::copy(tensor_data_ptr, tensor_data_ptr + diff_pattern_size, (*p_compress_chunk_data_8)[_id_chunk].begin() + s);
            #else
            for (int k = 0; k < n_cam*n_cam; k++)
            {
                (*p_counts_data)[probe_position] += frame_buffer[buffer_id][frame_id][k];
                (*p_compress_chunk_data_8)[_id_chunk][probe_position%(chunksize*nx)*diff_pattern_size + k/(n_cam*compress_det_bin)*n_cam_bin + (k%n_cam)/compress_det_bin] += frame_buffer[buffer_id][frame_id][k];
            }
            #endif
        }
    };


    FileConnector file;
    std::thread read_thread;
    std::thread proc_thread;
    int buffer_id;
    int n_frame_filled=0;
    int n_frame_processed=0;
    int n_buffer_filled=0;
    int n_buffer_processed=0;
    uint8_t id_image = 0 ;
    uint64_t current_line = 0;
    uint64_t probe_position = 0;
    uint64_t probe_position_total = 0;
    uint64_t n_events_processed = 0;
    bool first_frame = true;
    int frame_id;

    inline void schedule_buffer()
    {
        while ((*this->p_processor_line)!=-1)
        {
            if (this->n_buffer_processed < this->n_buffer_filled)
            {
                this->buffer_id = this->n_buffer_processed % n_buffer;
                for (int frm = 0; frm < buffer_size; frm++)
                {
                    this->frame_id = this->n_frame_processed % buffer_size;
                    if ((this->n_frame_processed%nx) != (nx-1)){this->process[0]();}
                    // if ((this->n_frame_processed%nx) != (nx-1)){
                    //     for (int i = 0; i < n_proc; i++) {this->process[i]();}
                    // }

                    ++this->n_frame_processed;
                    this->probe_position = this->n_frame_processed;
                    this->current_line = this->n_frame_processed / ny ;
                    *this->p_preprocessor_line = (int)this->current_line-1;
                    if ((int)this->current_line == ny){*this->p_preprocessor_line = ny;}

                }

                ++this->n_buffer_processed;
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                this->process_wait++;
            }
        }
    };

    
    inline void init_uv()
    {
        u.resize(n_cam);
        v.resize(n_cam);

        for (int i = 0; i < n_cam; i++)
        {
            v[i] = i;
            u[i] = i;
        }
    };

    void reset()
    {
        n_frame_filled = 0;
        n_frame_processed = 0;
        current_line = 0;
    };


public:

//-------------------------------------------------------------------------------------------------
    void enable_vSTEM(std::vector<int> *_p_detector_image,std::vector<std::vector<size_t>> *_p_stem_data, bool _allow_torch = false){
        this->make_detector_list(_p_detector_image);
        p_stem_data = _p_stem_data;
        
        if (!_allow_torch) process.push_back(std::bind(&FRAMEBASED::vstem, this));
        

        #ifdef FRAMEBASED_TORCH_ENABLED
        if (_allow_torch) {
            process.push_back(std::bind(&FRAMEBASED::vstem_torch, this));

            std::cout << "Performing vSTEM using Torch" << std::endl;
            frame_torch_enabled = true;
            // device = at::kCUDA;
            device = at::kCPU;
            DetectorTensor = torch::from_blob((*_p_detector_image).data(), {n_cam, n_cam}, torch::TensorOptions().dtype(torch::kInt32)).to(device).to(torch::kUInt8).requires_grad_(false);

            for (size_t i = 0; i < n_buffer; ++i) {
                Tensor_buffer[i] = torch::zeros({buffer_size, n_cam, n_cam}, torch::kUInt8).to(device);
                Tensor_buffer[i].requires_grad_(false);
            } 
        }
        #endif

        ++n_proc;
    };

    void make_detector_list(std::vector<int> *p_detector_image)
    {
        for (int i = 0; i < n_cam*n_cam; i++)
        {
            if ((*p_detector_image)[i] == 1) detector_list.push_back(i);
        }
    };

    void enable_Pacbed(std::vector<size_t> *_p_pacbed_data){
        p_pacbed_data = _p_pacbed_data;
        process.push_back(std::bind(&FRAMEBASED::pacbed, this));
        ++n_proc;
    };

    void enable_Ricom(std::vector<float> *_p_comx_image,std::vector<float> *_p_comy_image){
        p_comx_image = _p_comx_image;
        p_comy_image = _p_comy_image;
        init_uv();
        process.push_back(std::bind(&FRAMEBASED::com, this));
        ++n_proc;
    };
    
    void enable_roi(std::vector<std::vector<uint64_t>> *_p_roi_scan_image_stack,
                    std::vector<std::vector<uint64_t>> *_p_roi_diffraction_pattern_stack,
                    std::vector<uint64_t> *_p_roi_scan_image,
                    std::vector<uint64_t> *_p_roi_diffraction_pattern,
                    int _lower_left[2] , int _upper_right[2]){

        p_roi_scan_image_stack = _p_roi_scan_image_stack;
        p_roi_diffraction_pattern_stack = _p_roi_diffraction_pattern_stack;
        p_roi_scan_image = _p_roi_scan_image;
        p_roi_diffraction_pattern = _p_roi_diffraction_pattern;
        lower_left[0] = _lower_left[0];
        lower_left[1] = _lower_left[1];
        upper_right[0] = _upper_right[0];
        upper_right[1] = _upper_right[1];
        L_0 = upper_right[0]-lower_left[0];
        L_1 = upper_right[1]-lower_left[1];
        process.push_back(std::bind(&FRAMEBASED::roi, this));
        ++n_proc;
    }
    void enable_roi_4D(std::vector<std::vector<std::vector<std::vector<uint8_t>>>> *_p_roi_4D,std::vector<uint64_t> *_p_roi_scan_image,std::vector<uint64_t> *_p_roi_diffraction_pattern, int _lower_left[2] , int _upper_right[2],int _det_bin){
        p_roi_4D = _p_roi_4D;
        det_bin = _det_bin;
        p_roi_scan_image = _p_roi_scan_image;
        p_roi_diffraction_pattern = _p_roi_diffraction_pattern;
        lower_left[0] = _lower_left[0];
        lower_left[1] = _lower_left[1];
        upper_right[0] = _upper_right[0];
        upper_right[1] = _upper_right[1];
        L_0 = upper_right[0]-lower_left[0];
        L_1 = upper_right[1]-lower_left[1];
        process.push_back(std::bind(&FRAMEBASED::roi_4D, this));
        ++n_proc;
    }

    #ifdef GPRI_OPTION_ENABLED
    void enable_GPRI(std::vector<torch::Tensor> *_p_result_stack,torch::Tensor *_p_G_library,std::vector<std::vector<int32_t>> *_p_scan_index,
                    std::vector<int32_t> *_p_center,std::vector<int32_t> *_p_center_scan,int32_t _interval_R_ratio,int32_t _N_pxl_radius,int detector_bin, 
                    bool _normalize,torch::Device _dev,std::vector<std::vector<uint64_t>> *_p_N_electrons_map_scangrid){

        std::cout << "Performing GPRI in dense mode" << std::endl;

        GPRI_enabled = true;
        p_Result_stack = _p_result_stack;
        p_G_library = _p_G_library;
        p_scan_index = _p_scan_index;
        p_center= _p_center;
        p_center_scan = _p_center_scan;
        interval_R_ratio = _interval_R_ratio;
        N_pxl_radius = _N_pxl_radius;
        GPRI_detector_bin = detector_bin;
        normalize = _normalize;
        device = _dev;
        p_N_electrons_map_scangrid = _p_N_electrons_map_scangrid;
        process.push_back(std::bind(&FRAMEBASED::GPRI, this));
        ++n_proc;

        for (size_t i = 0; i < n_buffer; ++i) {
        Binned_Tensor_buffer[i] = torch::zeros({buffer_size, n_cam/GPRI_detector_bin, n_cam/GPRI_detector_bin}, Tensortype).to(device);
        Binned_Tensor_buffer[i].requires_grad_(false);
        } 

    }

    void enable_SparseGPRI(std::vector<std::shared_ptr<std::vector<int>>> *_p_k_indices_vec,int detector_bin, int scan_bin,
        std::vector<std::vector<uint64_t>> *_p_N_electrons_map_scangrid)
    {
        std::cout << "Performing GPRI in sparse mode" << std::endl;
        sparse_GPRI_enabled = true;
        p_k_indices_vec = _p_k_indices_vec;
        p_N_electrons_map_scangrid = _p_N_electrons_map_scangrid;
        GPRI_detector_bin = detector_bin;
        GPRI_scan_bin = scan_bin;
        GPRI_nxy_scan_bin = nxy/(scan_bin*scan_bin);
        GPRI_cam_bin = n_cam/detector_bin;
        process.push_back(std::bind(&FRAMEBASED::SparseGPRI, this));
        ++n_proc;
    }
    #endif

    void enable_EELS(std::vector<std::vector<std::vector<uint64_t>>> *_p_EELS_data_stack, std::vector<std::vector<uint64_t>> *_p_EELS_image_stack){
        p_EELS_data_stack = _p_EELS_data_stack;
        p_EELS_image_stack = _p_EELS_image_stack;
        process.push_back(std::bind(&FRAMEBASED::EELS, this));
        ++n_proc;
    };

    void enable_compress(std::vector<uint64_t> *_p_counts_data,std::vector<uint8_t> (*_p_compress_chunk_data_8)[2], int _chunksize, int det_bin, std::mutex* _mtx)
    {
        p_compress_chunk_data_8 = _p_compress_chunk_data_8;
        p_counts_data = _p_counts_data;
        process.push_back(std::bind(&FRAMEBASED::compress, this));
        ++n_proc;

        compress_det_bin = det_bin;
        chunksize = _chunksize;
        diff_pattern_size = n_cam/det_bin*n_cam/det_bin;
        diff_pattern_length = n_cam/det_bin;
        n_cam_bin = n_cam/det_bin;
        this->mtx = _mtx;

        #ifdef FRAMEBASED_TORCH_ENABLED
            temp_Tensor = torch::zeros({n_cam_bin, n_cam_bin}, torch::kUInt8).to(torch::kCPU).requires_grad_(false);
        #endif
    }

    void enable_compress(std::vector<uint64_t> *_p_counts_data,std::vector<uint16_t> (*_p_compress_chunk_data_16)[2], int _chunksize, int det_bin, std::mutex* _mtx)
    {
    }

    void enable_compress(std::vector<uint64_t> *_p_counts_data,std::vector<uint32_t> (*_p_compress_chunk_data_32)[2], int _chunksize, int det_bin, std::mutex* _mtx)
    {
    }

//-------------------------------------------------------------------------------------------------

    void terminate()
    {
        while ((!read_thread.joinable()) || (!proc_thread.joinable())) { std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
        read_thread.join();
        proc_thread.join();
        this->endtime = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        uint64_t rate = nxy / ((endtime - starttime) / 1e9);
        uint64_t frametime = ((endtime - starttime) / 1e3) / nxy;
        std::cout << "Processed at " << rate << " fps (" << frametime << " microsec per frame)" << std::endl; 
        std::cout << "read waits: " << read_wait << " process waits: " << process_wait << std::endl;
    };

//-------------------------------------------------------------------------------------------------
// Variables
//-------------------------------------------------------------------------------------------------
    int nx;
    int ny;
    int nxy;
    bool *b_cumulative;
    int repetitions;
    int data_depth;

    bool frame_torch_enabled = false;

    //------------------------
    // VSTEM
    std::vector<std::vector<size_t>> *p_stem_data;
    std::vector<int> detector_list;

    #ifdef FRAMEBASED_TORCH_ENABLED
    torch::Tensor DetectorTensor;
    #endif

    //PACBED
    std::vector<size_t> *p_pacbed_data;

    //RICOM
    std::vector<int> u;
    std::vector<int> v;
    std::vector<float> *p_comx_image;
    std::vector<float> *p_comy_image;
    std::vector<size_t> sum_x;
    std::vector<size_t> sum_y;
    float dose;
    float COM[2];
    size_t sum_x_temp = 0;
    size_t current_px;

    // ROI
    std::vector<std::vector<std::vector<std::vector<uint8_t>>>> *p_roi_4D;
    std::vector<uint64_t> *p_roi_scan_image;
    std::vector<uint64_t> *p_roi_diffraction_pattern;
    std::vector<std::vector<uint64_t>> *p_roi_scan_image_stack;
    std::vector<std::vector<uint64_t>> *p_roi_diffraction_pattern_stack;
    int lower_left[2];
    int upper_right[2];
    int L_0;
    int L_1;
    int det_bin;

    // EELS
    std::vector<std::vector<std::vector<uint64_t>>> *p_EELS_data_stack;
    std::vector<std::vector<uint64_t>> *p_EELS_image_stack;

    // GPRI
    //dense
    bool GPRI_enabled = false;
    #ifdef GPRI_OPTION_ENABLED
    torch::Tensor *p_G_library;
    torch::Tensor *p_Result;
    std::vector<torch::Tensor> *p_Result_stack;
    std::vector<std::vector<int32_t>> *p_scan_index;
    std::vector<int32_t> *p_center;
    std::vector<int32_t> *p_center_scan;
    int32_t interval_R_ratio;
    int32_t N_pxl_radius;
    int32_t position_x;
    int32_t position_y;
    int32_t lim_inf_x;
    int32_t lim_inf_y;
    int32_t lim_sup_x;
    int32_t lim_sup_y;
    bool normalize;
    static const torch::Dtype Tensortype = torch::kFloat32;

    // sparse
    bool sparse_GPRI_enabled = false;
    int GPRI_detector_bin;
    int GPRI_scan_bin; 
    uint16_t GPRI_cam_bin;
    uint64_t GPRI_nxy_scan_bin;
    std::vector<std::shared_ptr<std::vector<int>>> *p_k_indices_vec;
    std::vector<std::vector<uint64_t>> *p_N_electrons_map_scangrid;
    #endif

    #if defined(FRAMEBASED_TORCH_ENABLED) || defined(GPRI_OPTION_ENABLED)
    torch::NoGradGuard no_grad;
    torch::Device device = at::kCPU;
    std::array<torch::Tensor, n_buffer> Tensor_buffer;
    std::array<torch::Tensor, n_buffer> Binned_Tensor_buffer;
    #endif 

    // compress
    std::vector<uint8_t> (*p_compress_chunk_data_8)[2];
    std::vector<uint64_t> *p_counts_data;
    std::mutex* mtx = nullptr; //2 mutexes for 2 chunks
    int compress_det_bin;
    int id_chunk;
    int diff_pattern_size;
    int diff_pattern_length;
    int chunksize;
    int n_cam_bin;
    #ifdef FRAMEBASED_TORCH_ENABLED
    torch::Tensor temp_Tensor;
    #endif

    //-----------------------------------------------------------------------------------------------------------------------

    int *p_processor_line;
    int *p_preprocessor_line;
    int mode;
    std::string file_path;
    SocketConnector socket;

    std::vector<std::function<void()>> process; 
    int n_proc = 0;
    std::vector<std::vector<size_t> (*)[2]> p_images; 
    int n_images = 0;

    std::string dtype;

    // Data Properties
    int ds_merlin;
    bool b_raw;
    bool b_binary;

    typedef std::array<pixel,n_cam*n_cam> frame;

    std::array<frame*, n_buffer> frame_buffer;

    uint64_t starttime;
    uint64_t endtime;

    int read_wait = 0;
    int process_wait = 0;

    FRAMEBASED(
        int &nx,
        int &ny,
        bool *b_cumulative,
        int repetitions,
        int *p_processor_line,
        int *p_preprocessor_line,
        int &mode,
        std::string &file_path,  
        SocketConnector socket
    ) : 
        nx(nx),
        ny(ny),
        b_cumulative(b_cumulative),
        repetitions(repetitions),
        p_processor_line(p_processor_line),
        p_preprocessor_line(p_preprocessor_line),
        mode(mode),
        file_path(file_path), 
        socket(socket)
    {
        nxy = nx*ny;
        n_proc = 0;
        n_images = 0;
        for (size_t i = 0; i < n_buffer; ++i) {frame_buffer[i] = new frame[buffer_size];}
    }

};
#endif //FRAMEBASED_H