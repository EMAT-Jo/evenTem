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

#ifndef TIMEPIX_H
#define TIMEPIX_H

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
#include <cstdlib>

#include "SocketConnector.h"
#include "FileConnector.h"
#include "Declusterer.hpp"
#include "dtype_Electron.hpp"
#include "Logger.hpp"
#include "Roi4D.hpp"
#include "AtomicWrapper.hpp"

template <typename event, int buffer_size, int n_buffer>
class TIMEPIX
{
protected:

    #ifdef GPRI_OPTION_ENABLED
        enum class FunctionType 
        {
            vstem,
            multi_vstem,
            mask_vstem,
            com,
            GPRI,
            count_chunked_8,
            count_chunked_16,
            count_chunked_32,
            pacbed,
            var,
            roi,
            roi_mask,
            roi_4D,
            roi_ToT,
            write_electron,
            write_declusterer_buffer,
            information,
            atomic_vstem
        };
    #else
        enum class FunctionType 
        {
            vstem,
            multi_vstem,
            mask_vstem,
            com,
            count_chunked_8,
            count_chunked_16,
            count_chunked_32,
            pacbed,
            var,
            roi,
            roi_mask,
            roi_4D,
            roi_ToT,
            write_electron,
            write_declusterer_buffer,
            information,
            atomic_vstem
        };
    #endif


    // ----------------------------------------------------------------------------------------------- 
    // process methods
    // -----------------------------------------------------------------------------------------------

    inline void vstem(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        int _d2 = (_kx - x_offset)*(_kx - x_offset) + (_ky - y_offset)*(_ky - y_offset);
        if (_d2 > in_radius_sqr && _d2 <= out_radius_sqr)
        {   
            (*p_stem_data)[_id_image][_probe_position]++;
        }
    };

    inline void atomic_vstem(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        // (*p_atomic_stem_data)[_probe_position]++;
        // atomic_counter++;
    };

    inline void multi_vstem(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        for (int i = 0; i < n_detectors; i++)
        {
            int _d2 = (_kx - offsets[i][0])*(_kx - offsets[i][0]) + (_ky - offsets[i][1])*(_ky - offsets[i][1]);
            if (_d2 >= radia_sqr[i][0] && _d2 <= radia_sqr[i][1])
            {
                (*p_stem_data)[_id_image][_probe_position]++;
            }
        }
    };

    inline void mask_vstem(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
            (*p_stem_data)[_id_image][_probe_position] += detector_mask[_kx*n_cam+_ky];
    };

    inline void com(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        (*p_dose_data)[_id_image%2][_probe_position]++;
        (*p_sumy_data)[_id_image%2][_probe_position] += _ky;
        (*p_sumx_data)[_id_image%2][_probe_position] += _kx;
    };

    inline void com_masked(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        (*p_dose_data)[_id_image%2][_probe_position] += com_mask[_kx*n_cam+_ky];
        (*p_sumy_data)[_id_image%2][_probe_position] += _ky*com_mask[_kx*n_cam+_ky];
        (*p_sumx_data)[_id_image%2][_probe_position] += _kx*com_mask[_kx*n_cam+_ky];
    };

    #ifdef GPRI_OPTION_ENABLED
    inline void GPRI(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        // uint64_t _x_pp = (_probe_position%nx)/GPRI_scan_bin;
        // uint64_t _y_pp = (_probe_position/nx)/GPRI_scan_bin;
        // _probe_position = (_y_pp*(nx/GPRI_scan_bin) + _x_pp);

        (*p_k_indices_vec)[_probe_position+_id_image*GPRI_nxy_scan_bin]->push_back((_kx/GPRI_detector_bin)*GPRI_cam_bin+(_ky/GPRI_detector_bin));
        (*p_N_electrons_map_scangrid)[_id_image][_probe_position] += 1;
    };
    #endif

    inline void count_chunked_8(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        uint64_t _x_pp = _probe_position%nx;
        uint64_t _y_pp = _probe_position/nx;

        uint64_t _bin_probe_position = (_y_pp/fourD_scan_bin)*nx_scan_bin + _x_pp/fourD_scan_bin;

        (*p_counts_data)[_bin_probe_position]++;
        int _id_chunk = (_bin_probe_position/(chunksize_scan_bin*nx_scan_bin))%2;
        std::lock_guard<std::mutex> lock(mtx[_id_chunk]);
        (*p_fourDchunk_data_8)[_id_chunk][((_bin_probe_position)%(chunksize_scan_bin*nx_scan_bin))*diff_pattern_size + (_kx/fourD_det_bin*n_cam/fourD_det_bin+_ky/fourD_det_bin)]++;
    };

    inline void count_chunked_16(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        uint64_t _x_pp = _probe_position%nx;
        uint64_t _y_pp = _probe_position/nx;

        uint64_t _bin_probe_position = (_y_pp/fourD_scan_bin)*nx_scan_bin + _x_pp/fourD_scan_bin;

        (*p_counts_data)[_bin_probe_position]++;
        int _id_chunk = (_bin_probe_position/(chunksize_scan_bin*nx_scan_bin))%2;
        std::lock_guard<std::mutex> lock(mtx[_id_chunk]);
        (*p_fourDchunk_data_16)[_id_chunk][((_bin_probe_position)%(chunksize_scan_bin*nx_scan_bin))*diff_pattern_size + (_kx/fourD_det_bin*n_cam/fourD_det_bin+_ky/fourD_det_bin)]++;
    };

    inline void count_chunked_32(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        uint64_t _x_pp = _probe_position%nx;
        uint64_t _y_pp = _probe_position/nx;

        uint64_t _bin_probe_position = (_y_pp/fourD_scan_bin)*nx_scan_bin + _x_pp/fourD_scan_bin;

        (*p_counts_data)[_bin_probe_position]++;
        int _id_chunk = (_bin_probe_position/(chunksize_scan_bin*nx_scan_bin))%2;
        std::lock_guard<std::mutex> lock(mtx[_id_chunk]);
        (*p_fourDchunk_data_32)[_id_chunk][((_bin_probe_position)%(chunksize_scan_bin*nx_scan_bin))*diff_pattern_size + (_kx/fourD_det_bin*n_cam/fourD_det_bin+_ky/fourD_det_bin)]++;
    };

    inline void pacbed(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        (*p_pacbed_data)[_kx*n_cam+_ky]++;
    };

    inline void var(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        (*p_var_data)[_id_image][_probe_position] += (_kx-offset[0])*(_kx-offset[0])+(_ky-offset[0])*(_ky-offset[0]);
    };

    inline void roi(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        int _x = _probe_position%nx;
        int _y = nx - floor(_probe_position/nx);

        if (_x >= lower_left[0] && _x < upper_right[0] && _y > lower_left[1] && _y <= upper_right[1])
        {
            (*p_roi_diffraction_pattern_stack)[_id_image][_kx*n_cam+_ky]++;
            (*p_roi_scan_image_stack)[_id_image][(L_1 - (_y-lower_left[1])) * L_0 + (_x-lower_left[0])]++;
            (*p_roi_diffraction_pattern)[_kx*n_cam+_ky]++;
            (*p_roi_scan_image)[(L_1 - (_y-lower_left[1])) * L_0 + (_x-lower_left[0])]++;
        } 
    };

    inline void roi_ToT(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        int _x = _probe_position%nx;
        int _y = nx - floor(_probe_position/nx);

        if (_x >= lower_left[0] && _x < upper_right[0] && _y > lower_left[1] && _y <= upper_right[1])
        {
            (*p_roi_diffraction_pattern_stack)[_id_image][_kx*n_cam+_ky] += this->tot;
            (*p_roi_scan_image_stack)[_id_image][(L_1 - (_y-lower_left[1])) * L_0 + (_x-lower_left[0])]++;
            (*p_roi_diffraction_pattern)[_kx*n_cam+_ky] += this->tot;
            (*p_roi_scan_image)[(L_1 - (_y-lower_left[1])) * L_0 + (_x-lower_left[0])]++;
        } 
    };

    inline void roi_mask(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        if (mask_roi[_id_image][_probe_position] == 1)
        {
            (*p_roi_diffraction_pattern_stack)[_id_image][_kx*n_cam+_ky]++;
            (*p_roi_scan_image_stack)[_id_image][_probe_position]++;
            (*p_roi_diffraction_pattern)[_kx*n_cam+_ky] += 1;
            (*p_roi_scan_image)[_probe_position]++;
        } 
    };
  
    inline void roi_4D(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        int _x = _probe_position%nx;
        int _y = nx - floor(_probe_position/nx);

        if (_x >= lower_left[0] && _x < upper_right[0] && _y > lower_left[1] && _y <= upper_right[1] )
        {
            (*p_roi_diffraction_pattern)[_kx*n_cam+_ky]++;
            (*p_roi_scan_image)[(L_1 - (_y-lower_left[1])) * L_0 + (_x-lower_left[0])]++;
            p_roi_4D->increment(L_1 - (_y-lower_left[1]), _x-lower_left[0], _kx/det_bin, _ky/det_bin);
        } 
    };

    inline void write_electron(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        // FIX ELECTRON NOT LOCAL VAR NOW
        electron.kx = _kx/det_bin_electron;
        electron.ky = _ky/det_bin_electron;
        electron.rx = (_probe_position%nx)/scan_bin_electron;
        electron.ry = (_probe_position/nx)/scan_bin_electron;
        electron.id_image = _id_image;

        if (electron.rx < (x_crop/scan_bin_electron) && electron.ry < (y_crop/scan_bin_electron))
        {
            (*p_file).write((const char *)&electron, sizeof(electron));
        }
    };

    inline void write_declusterer_buffer(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image, uint64_t _toa, uint16_t _tot)
    {
        uint16_t _rx = _probe_position%nx;
        uint16_t _ry = _probe_position/nx;
        declusterer.buffer[declusterer.buffer_id_filling]->push_back({_kx, _ky, _rx, _ry, _id_image,_toa,_tot});
    };


    inline void information(uint64_t _probe_position, uint16_t _kx, uint16_t _ky, uint16_t _id_image)
    {
        (*p_information_image)[_probe_position] += -log2((*p_probability_distribution)[_kx*n_cam+_ky]);
        (*p_count_image)[_probe_position]++;
    };

protected:
    
    FileConnector file;
    std::thread read_thread;
    std::thread proc_thread;
    int n_buf = n_buffer;
    int n_buffer_filled=0;
    int n_buffer_processed=0;
    uint16_t id_image = 0 ;
    uint64_t current_line = 0;
    uint64_t probe_position = 0;
    uint64_t probe_position_total = 0;
    uint16_t kx;
    uint16_t ky;
    uint64_t n_events_processed = 0;

    inline void read_file()
    {
        int buffer_id;
        while ((!this->repetitions_reached) && (*p_processor_line!=-1))
        {
            if (n_buffer_filled < (n_buffer + n_buffer_processed)) 
            {
                buffer_id = n_buffer_filled % n_buffer; 
                file.read_data((char *)&(buffer[buffer_id]), sizeof(buffer[buffer_id]));
                ++n_buffer_filled;
            }
            else
            {
                read_wait++;
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
        file.close_file();
    };

    inline void read_socket()
    {
        int buffer_id;
        while (*p_processor_line != -1)
        {
            if (n_buffer_filled < (n_buffer + n_buffer_processed))
            {
                buffer_id = n_buffer_filled % n_buffer;
                socket.read_data((char *)&(buffer[buffer_id]), sizeof(buffer[buffer_id]));
                ++n_buffer_filled;
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
       
    };

    void flush_image(int id_img)
    {
        // //flushes intermidiate images like com data
        // for (int i_image = 0; i_image < n_images; i_image++)
        // {
        //     std::fill(
        //         (*(p_images[i_image]))[id_img].begin(),
        //         (*(p_images[i_image]))[id_img].end(), 0);
        // }
    };

    void reset()
    {
        n_buffer_filled = 0;
        n_buffer_processed = 0;
        current_line = 0;
    };

public:
//-------------------------------------------------------------------------------------------------
    void enable_multi_vSTEM(std::vector<std::array<float, 2>> *_p_radia_sqr,std::vector<std::array<float, 2>> *_p_offsets,std::vector<std::vector<size_t>> *_p_stem_data)
    {
        p_stem_data = _p_stem_data;
        radia_sqr = (*_p_radia_sqr);
        offsets = *_p_offsets;
        n_detectors = radia_sqr.size();
        process.push_back(std::bind(&TIMEPIX::multi_vstem, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::multi_vstem;
        ++n_proc;
    }

    void enable_vSTEM(std::array<float, 2> *_p_radius_sqr,std::array<float, 2> *_p_offset,std::vector<std::vector<size_t>> *_p_stem_data)
    {
        p_stem_data = _p_stem_data;
        in_radius_sqr = (int)(*_p_radius_sqr)[0];
        out_radius_sqr = (int)(*_p_radius_sqr)[1];
        x_offset = (int)(*_p_offset)[0];
        y_offset = (int)(*_p_offset)[1];
        process.push_back(std::bind(&TIMEPIX::vstem, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::vstem;
        ++n_proc;
    }

    void enable_atomic_vSTEM(std::array<float, 2> *_p_radius_sqr,std::array<float, 2> *_p_offset,std::vector<atomwrapper<int>> *_p_atomic_stem_data)
    {
        p_atomic_stem_data = _p_atomic_stem_data;
        in_radius_sqr = (int)(*_p_radius_sqr)[0];
        out_radius_sqr = (int)(*_p_radius_sqr)[1];
        x_offset = (int)(*_p_offset)[0];
        y_offset = (int)(*_p_offset)[1];
        functionType = FunctionType::atomic_vstem;
        process.push_back(std::bind(&TIMEPIX::atomic_vstem, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        ++n_proc;
    }

    void enable_mask_vSTEM(std::vector<int> *_p_detector_mask,std::vector<std::vector<size_t>> *_p_stem_data)
    {
        p_stem_data = _p_stem_data;
        detector_mask = *_p_detector_mask;
        process.push_back(std::bind(&TIMEPIX::mask_vstem, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::mask_vstem;
        ++n_proc;
    }

    void enable_Ricom(std::vector<size_t> (*_p_dose_data)[2],std::vector<size_t> (*_p_sumx_data)[2],std::vector<size_t> (*_p_sumy_data)[2])
    {
        p_dose_data = _p_dose_data;
        p_sumx_data = _p_sumx_data;
        p_sumy_data = _p_sumy_data;

        process.push_back(std::bind(&TIMEPIX::com, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::com;
        p_images.push_back(p_sumx_data);
        p_images.push_back(p_sumy_data);
        p_images.push_back(p_dose_data);
        ++n_proc; 
        n_images += 3;
    }

    void enable_Ricom_masked(std::vector<int> *_p_com_mask,std::vector<size_t> (*_p_dose_data)[2],std::vector<size_t> (*_p_sumx_data)[2],std::vector<size_t> (*_p_sumy_data)[2])
    {
        p_dose_data = _p_dose_data;
        p_sumx_data = _p_sumx_data;
        p_sumy_data = _p_sumy_data;

        com_mask = *_p_com_mask;
        process.push_back(std::bind(&TIMEPIX::com_masked, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::com;
        p_images.push_back(p_sumx_data);
        p_images.push_back(p_sumy_data);
        p_images.push_back(p_dose_data);
        ++n_proc; 
        n_images += 3;
    }


    void enable_FourD(std::vector<uint64_t> *_p_counts_data, std::vector<uint8_t> (*_p_fourD_data)[2], size_t det_bin, size_t scan_bin, size_t _chunksize, std::mutex* _mtx)
    {
        p_counts_data = _p_counts_data;
        p_fourDchunk_data_8 = _p_fourD_data;
        process.push_back(std::bind(&TIMEPIX::count_chunked_8, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::count_chunked_8;
        ++n_proc;

        fourD_det_bin = det_bin;
        fourD_scan_bin = scan_bin;
        chunksize_scan_bin = _chunksize/scan_bin;
        n_cam_det_bin = n_cam/det_bin;
        nx_scan_bin = nx/scan_bin;
        diff_pattern_size = n_cam/fourD_det_bin*n_cam/fourD_det_bin;
        this->mtx = _mtx;
    }
    void enable_FourD(std::vector<uint64_t> *_p_counts_data, std::vector<uint16_t> (*_p_fourD_data)[2], size_t det_bin, size_t scan_bin, size_t _chunksize, std::mutex* _mtx)
    {
        p_counts_data = _p_counts_data;
        p_fourDchunk_data_16 = _p_fourD_data;
        process.push_back(std::bind(&TIMEPIX::count_chunked_16, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::count_chunked_16;
        ++n_proc;

        fourD_det_bin = det_bin;
        fourD_scan_bin = scan_bin;
        chunksize_scan_bin = _chunksize/scan_bin;
        n_cam_det_bin = n_cam/det_bin;
        nx_scan_bin = nx/scan_bin;
        diff_pattern_size = n_cam/fourD_det_bin*n_cam/fourD_det_bin;
        this->mtx = _mtx;
    }
    void enable_FourD(std::vector<uint64_t> *_p_counts_data, std::vector<uint32_t> (*_p_fourD_data)[2], size_t det_bin, size_t scan_bin, size_t _chunksize, std::mutex* _mtx)
    {
        p_counts_data = _p_counts_data;
        p_fourDchunk_data_32 = _p_fourD_data;
        process.push_back(std::bind(&TIMEPIX::count_chunked_32, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::count_chunked_32;
        ++n_proc;

        fourD_det_bin = det_bin;
        fourD_scan_bin = scan_bin;
        chunksize_scan_bin = _chunksize/scan_bin;
        n_cam_det_bin = n_cam/det_bin;
        nx_scan_bin = nx/scan_bin;
        diff_pattern_size = n_cam/fourD_det_bin*n_cam/fourD_det_bin;
        this->mtx = _mtx;
    }

    void enable_Pacbed(std::vector<size_t> *_p_pacbed_data)
    {
        p_pacbed_data = _p_pacbed_data;
        process.push_back(std::bind(&TIMEPIX::pacbed, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::pacbed;
        ++n_proc;
    }

    void enable_var(std::vector<size_t> (*_p_var_data)[2],std::array<float, 2> _offset)
    {
        p_var_data = _p_var_data;
        offset = _offset;
        process.push_back(std::bind(&TIMEPIX::var, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        p_images.push_back(p_var_data);
        functionType = FunctionType::var;
        ++n_proc;
        ++n_images;
    }

    void enable_roi(std::vector<std::vector<uint64_t>> *_p_roi_scan_image_stack,std::vector<std::vector<uint64_t>> *_p_roi_diffraction_pattern_stack,
    std::vector<uint64_t> *_p_roi_scan_image,std::vector<uint64_t> *_p_roi_diffraction_pattern,int _lower_left[2] , int _upper_right[2])
    {
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
        process.push_back(std::bind(&TIMEPIX::roi, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::roi;
        ++n_proc;
        // this->b_tot = true;
    }

    void enable_roi_mask(std::vector<std::vector<int>> *_p_roi_mask,std::vector<std::vector<uint64_t>> *_p_roi_scan_image_stack,std::vector<std::vector<uint64_t>> *_p_roi_diffraction_pattern_stack,
    std::vector<uint64_t> *_p_roi_scan_image,std::vector<uint64_t> *_p_roi_diffraction_pattern)
    {
        p_roi_scan_image_stack = _p_roi_scan_image_stack;
        p_roi_diffraction_pattern_stack = _p_roi_diffraction_pattern_stack;
        p_roi_scan_image = _p_roi_scan_image;
        p_roi_diffraction_pattern = _p_roi_diffraction_pattern;
        mask_roi = *_p_roi_mask;

        process.push_back(std::bind(&TIMEPIX::roi_mask, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::roi_mask;
        ++n_proc;

    }

    template <typename T>
    void enable_roi_4D(std::shared_ptr<Roi4D<T>> _p_roi_4D,std::vector<uint64_t> *_p_roi_scan_image,
    std::vector<uint64_t> *_p_roi_diffraction_pattern, int _lower_left[2] , int _upper_right[2],int _det_bin)
    {
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
        process.push_back(std::bind(&TIMEPIX::roi_4D, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::roi_4D;
        ++n_proc;
    }

    void enable_electron(std::ofstream& _p_file, bool _decluster, uint64_t _dtime, uint16_t _dspace, int _cluster_range, int _x_crop, int _y_crop,
    int _scan_bin_electron, int _det_bin_electron, int _n_threads, std::vector<int> *_p_clustersize_histogram)
    {
        decluster = _decluster;
        if (decluster) declusterer.init(_dtime, _dspace, _cluster_range,_x_crop,_y_crop,_scan_bin_electron,_det_bin_electron,_p_file,_n_threads,_p_clustersize_histogram);
        p_file = &_p_file;
        x_crop = _x_crop;
        y_crop = _y_crop;
        scan_bin_electron = _scan_bin_electron;
        det_bin_electron = _det_bin_electron;
        if (decluster) functionType = FunctionType::write_declusterer_buffer;
        else functionType = FunctionType::write_electron;
        if (decluster) decluster_thread = std::thread(&Declusterer::run, &declusterer);
        ++n_proc;        
    } 

    #ifdef GPRI_OPTION_ENABLED
    void enable_GPRI(std::vector<std::shared_ptr<std::vector<int>>> *_p_k_indices_vec,int detector_bin, int scan_bin,
    std::vector<std::vector<uint64_t>> *_p_N_electrons_map_scangrid)
    {
        GPRI_enabled = true;
        p_k_indices_vec = _p_k_indices_vec;
        GPRI_detector_bin = detector_bin;
        GPRI_scan_bin = scan_bin;
        GPRI_nxy_scan_bin = nxy/(scan_bin*scan_bin);
        GPRI_cam_bin = n_cam/detector_bin;
        p_N_electrons_map_scangrid = _p_N_electrons_map_scangrid;
    
        process.push_back(std::bind(&TIMEPIX::GPRI, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::GPRI;
        ++n_proc;
    }
    #endif

    void enable_information(std::vector<float> *_p_information_image, std::vector<float> *_p_probability_distribution, std::vector<float> *_p_count_image)
    {
        p_information_image = _p_information_image;
        p_probability_distribution = _p_probability_distribution;
        p_count_image = _p_count_image;
        process.push_back(std::bind(&TIMEPIX::information, this,std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
        functionType = FunctionType::information;
        ++n_proc;
    }

    //-------------------------------------------------------------------------------------------------

    void terminate()
    {
        while ((!read_thread.joinable()) || (!proc_thread.joinable())) { std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
        read_thread.join();
        proc_thread.join();
        if (decluster)
        {
            declusterer.still_reading = false;
            while (declusterer.still_processing || declusterer.still_writing) {std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
            declusterer.terminate();
            decluster_thread.join();
        }
        this->endtime = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        processing_rate = n_events_processed / ((endtime - starttime) / 1e9);
        std::cout << n_events_processed << " events processed at " << processing_rate/1000000 << " M events/s " << std::endl; 
        if (decluster) std::cout << declusterer.n_electrons_kept << " electrons kept" << std::endl;
        std::cout << "reading waited " << read_wait << " times, processing waited " << process_wait << " times"<< std::endl;
        std::cout << "atomic counter: " << atomic_counter << std::endl;
        // std::cout << "avg BM time: " << BM_duration_sum/BM_count << " us for " << BM_count << " calls" << std::endl;
    };

    float get_processing_rate(){
        return processing_rate;
    }

    //-------------------------------------------------------------------------------------------------
    // Variables
    //-------------------------------------------------------------------------------------------------
    int nx;
    int ny;
    int nxy;
    int n_cam;
    bool *b_cumulative;
    // bool b_continuous;
    int repetitions;
    bool b_tot = false;
    bool repetitions_reached = false;

    uint64_t toa;
    uint16_t tot;

    float processing_rate;

    //------------------------
    // VSTEM
    std::vector<atomwrapper<int>> *p_atomic_stem_data;
    std::vector<std::vector<size_t>> *p_stem_data;
    int d2;
    int in_radius_sqr;
    int out_radius_sqr;
    int x_offset;
    int y_offset;

    // masked VSTEM
    std::vector<int> detector_mask;

    // multi VSTEM
    int n_detectors;
    std::vector<std::array<float, 2>> radia_sqr;
    std::vector<std::array<float, 2>> offsets;

    // ricom 
    std::vector<size_t> (*p_dose_data)[2];
    std::vector<size_t> (*p_sumx_data)[2];
    std::vector<size_t> (*p_sumy_data)[2];
    std::vector<uint64_t> *p_counts_data;
    std::vector<int> com_mask;

    // GPRI
    #ifdef GPRI_OPTION_ENABLED
    bool GPRI_enabled = false;
    int GPRI_detector_bin;
    int GPRI_scan_bin; 
    uint16_t GPRI_cam_bin;
    uint64_t GPRI_nxy_scan_bin;
    std::vector<std::shared_ptr<std::vector<int>>> *p_k_indices_vec;
    std::vector<std::vector<uint64_t>> *p_N_electrons_map_scangrid;
    #endif


    // FourD
    std::vector<uint8_t> (*p_fourDchunk_data_8)[2];
    std::vector<uint16_t> (*p_fourDchunk_data_16)[2];
    std::vector<uint32_t> (*p_fourDchunk_data_32)[2];
    std::mutex* mtx = nullptr; //2 mutexes for 2 chunks
    size_t id_chunk;
    size_t fourD_det_bin;
    size_t fourD_scan_bin;
    size_t n_cam_det_bin;
    size_t diff_pattern_size;
    int x_pp;
    int y_pp;
    int bin_probe_position;
    int chunksize_scan_bin;
    int nx_scan_bin;

    // PACBED
    std::vector<size_t> *p_pacbed_data;

    // Variance
    std::array<float, 2> offset;
    std::vector<size_t> (*p_var_data)[2];

    // ROI
    std::shared_ptr<Roi4DBase> p_roi_4D;
    std::vector<uint64_t> *p_roi_scan_image;
    std::vector<uint64_t> *p_roi_diffraction_pattern;
    std::vector<std::vector<uint64_t>> *p_roi_scan_image_stack;
    std::vector<std::vector<uint64_t>> *p_roi_diffraction_pattern_stack;
    int lower_left[2];
    int upper_right[2];
    int L_0;
    int L_1;
    int det_bin; 
    std::vector<std::vector<int>> mask_roi;

    //Information
    std::vector<float> *p_information_image;
    std::vector<float> *p_probability_distribution;
    std::vector<float> *p_count_image;

    //Electron
    std::ofstream* p_file;
    dtype_Electron electron;
    int x_crop;
    int y_crop;
    int scan_bin_electron;
    int det_bin_electron;

    // declustering
    bool decluster = false;
    Declusterer declusterer;
    std::thread decluster_thread;

    // atomic counter
    std::atomic<uint64_t> atomic_counter = 0;
    //-----------------------------------------------------------------------------------------------------------------------

    int *p_processor_line;
    int *p_preprocessor_line;
    int mode;
    std::string file_path;
    SocketConnector socket;

    //std::array<std::array<event, buffer_size>, n_buffer> buffer; // allocated on stack -> limited by 2gig stack frame 
    std::array<event, buffer_size> *buffer = new std::array<event, buffer_size>[n_buffer]; // dynamic allocation on heap

    std::vector<std::function<void(uint64_t, uint16_t, uint16_t, uint16_t)>> process; 
    TIMEPIX<event, buffer_size, n_buffer>::FunctionType functionType;

    int n_proc = 0;
    std::vector<std::vector<size_t> (*)[2]> p_images; 
    int n_images = 0;

    uint64_t starttime;
    uint64_t endtime;

    //bemchmarking
    std::chrono::high_resolution_clock::time_point BM_start;
    std::chrono::high_resolution_clock::time_point BM_stop;
    float BM_duration_sum = 0;
    int BM_count = 0;

    float BM_duration_sum_2 = 0;
    int BM_count_2 = 0;

    int read_wait = 0;
    int process_wait = 0;

    TIMEPIX(
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
    }
};
#endif // TIMEPIX_H