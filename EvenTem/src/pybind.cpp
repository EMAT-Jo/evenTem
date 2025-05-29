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

#include <iostream>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include "Ricom.h"
#include "LiveProcessor.h"
#include "vSTEM.h"
#include "Pacbed.h"
#include "Var.h"
#include "Roi.h"
#include "Electron.h"
#include "EELS.h"
#include "FourD.h"

#ifdef GPRI_OPTION_ENABLED
        #include "GPRI.h"
#endif

namespace py = pybind11;

#ifdef EVENTEM
#define MODULE_NAME eventem
#endif
#ifdef EVENTEMTORCH
#define MODULE_NAME eventemTorch
#endif


PYBIND11_MODULE(MODULE_NAME, m) {

        py::class_<LiveProcessor>(m, "LiveProcessor")
        .def_readwrite("nx", &LiveProcessor::nx)
        .def_readwrite("ny", &LiveProcessor::ny)
        .def_readwrite("dt", &LiveProcessor::dt)
        .def_readwrite("detector_size", &LiveProcessor::n_cam)
        .def("set_socket", &LiveProcessor::set_socket)
        .def("accept_socket", &LiveProcessor::accept_socket)
        .def("close_socket", &LiveProcessor::close_socket)
        .def_readwrite("n_threads", &LiveProcessor::n_threads)
        .def_readwrite("file_path", &LiveProcessor::file_path)
        .def_readwrite("repetitions", &LiveProcessor::rep)
        .def("set_dwell_time", &LiveProcessor::set_dwell_time)
        .def_readonly("progress", &LiveProcessor::progress_percent)
        .def_readwrite("b_cumulative", &LiveProcessor::b_cumulative)   
        .def_readwrite("b_continuous", &LiveProcessor::b_continuous)
        .def_readwrite("rc_quit", &LiveProcessor::rc_quit)
        .def_readonly("elapsed_seconds", &LiveProcessor::elapsed_seconds_vec)
        .def_readonly("reached_pp_id", &LiveProcessor::reached_pp_id)
        .def("set_pattern_file", &LiveProcessor::set_pattern_file)
        .def_readonly("processing_rate", &LiveProcessor::processing_rate)
        .def("set_file", &LiveProcessor::set_file);


        py::class_<Ricom,LiveProcessor>(m, "Ricom")
        .def(py::init<int>(),py::arg("repetitions"))
        .def("add_child", py::overload_cast<vSTEM*>(&Ricom::add_child))
        .def("set_kernel", &Ricom::set_kernel)
        .def("set_kernel_filtered", &Ricom::set_kernel_filtered)
        .def("run", &Ricom::run)
        .def_readonly("offset", &Ricom::offset)
        .def("set_offset",&Ricom::set_offset)
        .def("get_kernel", &Ricom::get_kernel)
        .def_readonly("comx_image", &Ricom::comx_image)
        .def_readonly("comy_image", &Ricom::comy_image)
        .def_readonly("ricom_stack", &Ricom::ricom_image_stack)
        .def_readonly("ricom_image", &Ricom::ricom_image);


        py::class_<vSTEM,LiveProcessor>(m, "vSTEM")
        .def(py::init<int>(),py::arg("repetitions"))
        .def_readwrite("inner_radia", &vSTEM::inner_radia)
        .def_readwrite("outer_radia", &vSTEM::outer_radia)
        .def_readonly("offsets", &vSTEM::offsets)
        .def("set_offsets", &vSTEM::set_offsets)
        .def("get_detector", &vSTEM::get_detector,py::return_value_policy::copy)
        .def("set_detector_mask", &vSTEM::set_detector_mask)
        .def("run", &vSTEM::run)
        .def_readwrite("allow_torch", &vSTEM::allow_torch)
        .def_readwrite("allow_cuda", &vSTEM::allow_cuda)        
        .def_readonly("vSTEM_stack", &vSTEM::vSTEM_stack)
        .def_readonly("vSTEM_image", &vSTEM::vSTEM_image);

                
        py::class_<FourD<8>, LiveProcessor>(m, "FourD8")
        .def(py::init<const std::string&, int, int,int>(), py::arg("output_filename"), py::arg("repetitions"), py::arg("bitdepth"),py::arg("compression_factor"))
        .def("run", &FourD<8>::run)
        .def("allocate_chunk", &FourD<8>::allocate_chunk)
        .def("save_dose_image", &FourD<8>::save_dose_image)
        .def("init_4D_file", &FourD<8>::init_4D_file)
        .def_readwrite("det_bin", &FourD<8>::det_bin)
        .def_readwrite("scan_bin", &FourD<8>::scan_bin)
        .def_readwrite("chunksize", &FourD<8>::chunksize)
        .def_readonly("Dose_image", &FourD<8>::Dose_image);


        py::class_<FourD<16>, LiveProcessor>(m, "FourD16")
        .def(py::init<const std::string&, int, int,int>(), py::arg("output_filename"), py::arg("repetitions"), py::arg("bitdepth"),py::arg("compression_factor"))
        .def("run", &FourD<16>::run)
        .def("allocate_chunk", &FourD<16>::allocate_chunk)
        .def("save_dose_image", &FourD<16>::save_dose_image)
        .def("init_4D_file", &FourD<16>::init_4D_file)
        .def_readwrite("det_bin", &FourD<16>::det_bin)
        .def_readwrite("scan_bin", &FourD<16>::scan_bin)
        .def_readwrite("chunksize", &FourD<16>::chunksize)
        .def_readonly("Dose_image", &FourD<16>::Dose_image);


        py::class_<FourD<32>, LiveProcessor>(m, "FourD32")
        .def(py::init<const std::string&, int, int,int>(), py::arg("output_filename"), py::arg("repetitions"), py::arg("bitdepth"),py::arg("compression_factor"))
        .def("run", &FourD<32>::run)
        .def("allocate_chunk", &FourD<32>::allocate_chunk)
        .def("save_dose_image", &FourD<32>::save_dose_image)
        .def("init_4D_file", &FourD<32>::init_4D_file)
        .def_readwrite("det_bin", &FourD<32>::det_bin)
        .def_readwrite("scan_bin", &FourD<32>::scan_bin)
        .def_readwrite("chunksize", &FourD<32>::chunksize)
        .def_readonly("Dose_image", &FourD<32>::Dose_image);


        #ifdef GPRI_OPTION_ENABLED
                py::class_<GPRI,LiveProcessor>(m, "GPRI")
                .def(py::init<int,std::string&,bool,bool>(),py::arg("repetitions"), py::arg("path_to_library"),py::arg("allow_cuda"),py::arg("allow_mps"))
                .def("add_child", py::overload_cast<vSTEM*>(&GPRI::add_child))
                .def("get_GPRI_result", &GPRI::get_GPRI_result)
                .def("get_GPRI_stack_result", &GPRI::get_GPRI_stack_result)
                .def("set_SparseFrame", &GPRI::set_SparseFrame)
                .def_readwrite("scan_index", &GPRI::scan_index)
                .def_readwrite("interval_R_ratio", &GPRI::interval_R_ratio)
                .def_readwrite("center", &GPRI::center)
                .def_readwrite("center_scan", &GPRI::center_scan) 
                .def_readwrite("N_pxl_radius", &GPRI::N_pxl_radius)
                .def_readwrite("detector_bin", &GPRI::detector_bin)
                .def_readwrite("scan_bin", &GPRI::scan_bin)
                .def_readonly("result_shape", &GPRI::result_shape)
                .def_readonly("N_electrons_map_scangrid", &GPRI::N_electrons_map_scangrid)
                .def_readwrite("normalize", &GPRI::normalize)   
                .def_readonly("N_image_modes", &GPRI::image_modes)
                .def_readonly("kernelsizes", &GPRI::kernelsizes)
                .def("run", &GPRI::run);
        #endif
        
        py::class_<Pacbed,LiveProcessor>(m, "Pacbed")
        .def(py::init<int>(),py::arg("repetitions"))
        .def("run", &Pacbed::run)
        .def_readonly("Pacbed_image", &Pacbed::Pacbed_image);

        py::class_<Roi,LiveProcessor>(m, "Roi")
        .def(py::init<int,bool>(), py::arg("repetitions"), py::arg("extract_4D"))
        .def("run", &Roi::run)
        .def_readonly("Roi_scan_image", &Roi::Roi_scan_image)
        .def_readonly("Roi_diffraction_pattern", &Roi::Roi_diffraction_pattern)
        .def_readonly("Roi_scan_image_stack", &Roi::Roi_scan_image_stack)
        .def_readonly("Roi_diffraction_pattern_stack", &Roi::Roi_diffraction_pattern_stack)
        .def("get_4D", &Roi::get_4D)
        .def_readwrite("det_bin", &Roi::det_bin)
        .def("get_roi", &Roi::get_roi)
        .def("set_roi_mask", &Roi::set_roi_mask)
        .def("set_bitdepth", &Roi::set_bitdepth)
        .def("set_roi", &Roi::set_roi,py::arg("x"),py::arg("y"),py::arg("width"),py::arg("height"));

        py::class_<Var,LiveProcessor>(m, "Var")
        .def(py::init<int>(),py::arg("repetitions"))
        .def("run", &Var::run)
        .def("set_offset", &Var::set_offset)
        .def_readonly("offset", &Var::offset)
        .def_readonly("Var_image", &Var::Var_image);

        py::class_<Electron,LiveProcessor>(m, "Electron")
        .def(py::init<int>(),py::arg("repetitions"))
        .def("run", &Electron::run)
        .def_readwrite("decluster", &Electron::decluster)
        .def_readwrite("dtime", &Electron::dtime)
        .def_readwrite("dspace", &Electron::dspace)
        .def_readwrite("cluster_range", &Electron::cluster_range)
        .def_readwrite("n_threads", &Electron::n_threads)
        .def_readwrite("x_crop", &Electron::x_crop)
        .def_readwrite("y_crop", &Electron::y_crop)
        .def_readwrite("scan_bin", &Electron::scan_bin)
        .def_readwrite("detector_bin", &Electron::detector_bin)
        .def_readwrite("clustersize_histogram", &Electron::clustersize_histogram);

        py::class_<EELS,LiveProcessor>(m, "EELS")
        .def(py::init<int>(),py::arg("repetitions"))
        .def("run", &EELS::run)
        .def_readonly("EELS_data", &EELS::EELS_data);

}
