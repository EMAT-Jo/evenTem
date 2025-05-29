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

#ifndef RICOM_H
#define RICOM_H

#ifdef _WIN32
#include <io.h>
#pragma warning(disable : 4005 4333 34)
#else
#include <unistd.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>

#include <stdio.h>
#include <complex>
#include <cfloat>
#include <vector>
#include <string> 
#include <mutex>
#include <future>
#include <thread>
// #include <fftw3.h>
#include <chrono>
#include <algorithm>
#include <ctime>

#include "LiveProcessor.h"
#include "BoundedThreadPool.hpp"
// #include "fft2d.hpp"
#include "SocketConnector.h"
#include "ProgressMonitor.h"
#include "Cheetah.hpp"
#include "Timepix.hpp"
#include "Advapix.hpp"
#include "Merlin.hpp"
#include "Numpy.hpp"
#include "AnnularDetector.hpp"

#include "vSTEM.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace chc = std::chrono;


class Ricom_kernel
{
public:
    // Properties
    uint16_t nx_cam;
    uint16_t ny_cam;
    int kernel_size = 5;
    bool b_filter;
    std::array<int, 2> kernel_filter_frequency;
    int k_width_sym;
    int k_area;
    float rotation;
    std::vector<float> kernel_x;
    std::vector<float> kernel_y;
    std::vector<float> kernel_filter;
    std::vector<float> f_approx;
    // Methods
    void compute_kernel();
    void compute_filter();
    void include_filter();
    std::vector<int> fftshift_map(int x, int y); 
    // Constructor
    Ricom_kernel() : kernel_size(5),
                     b_filter(false),
                     kernel_filter_frequency{1, 4},
                     k_width_sym(0),
                     k_area(0),
                     rotation(0.0),
                     kernel_x(),
                     kernel_y(),
                     kernel_filter(),
                     f_approx()
    {
        compute_kernel();
    };
    void approximate_frequencies(size_t n_im);
    // Destructor
    ~Ricom_kernel(){};
};

class Ricom : public LiveProcessor
{
private:
    // Electric field magnitude
    float e_mag_max = -FLT_MAX;
    float e_mag_min = FLT_MAX;

    // inline void rescales_recomputes();

    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    void icom_group_classical(int idxx, int id_image);

    inline void compute_electric_field(std::array<float, 2> &p_com_xy, size_t id);

public:

    bool update_offset = true;

    bool b_e_mag = false;

    // bool b_recompute_detector;
    bool b_recompute_kernel = false;
    bool masked_com;

    Ricom_kernel kernel;
    std::array<float, 2> offset = {0,0};
    std::array<float, 2> com_public = {0, 0};
    std::vector<float> comx_image;
    std::vector<float> comy_image;
    std::vector<float> ricom_image;
    // std::vector<std::atomic<float>> ricom_image;
    std::vector<std::vector<float>> ricom_image_stack;
    
    std::vector<size_t> dose_data[2];
    std::vector<size_t> sumx_data[2];
    std::vector<size_t> sumy_data[2];
    std::vector<std::complex<float>> e_field_data;

    
    std::atomic<bool> rescale_ricom = false;
    std::atomic<bool> rescale_e_mag = false;

    bool auto_offset = true;

    void run();
    void reset();

    //pybind
    void set_kernel(int kernel_size, int rotation);
    void set_kernel_filtered(int kernel_size, int rotation,float low_pass, float high_pass);
    void set_offset(std::array<float, 2>);
    std::vector<float> get_kernel();

    // void set_masked_com(int radius, std::array<float, 2> offset);
    // std::vector<int> com_mask;
    // AnnularDetector detector = AnnularDetector(0, 100);
    // std::vector<std::array<float, 2>> offsets = {{256, 256}};


    // Children
    std::vector<vSTEM*> vSTEM_children;
    std::vector<std::optional<std::thread>> threads_for_children;
    void setup_child(vSTEM* child);
    void add_child(vSTEM* child);


    // Constructor
    Ricom(int repetitions) : LiveProcessor(repetitions)
    {
        n_threads_max = std::thread::hardware_concurrency();
    }

    //Destructor
    ~Ricom(){};
    };

#endif // __RICOM_H__