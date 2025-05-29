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

#include "FourD.h"

template<int BitDepth>
void FourD<BitDepth>::run(){
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
           
            cam.enable_FourD(&Dose_image, &chunk_data, det_bin,scan_bin, chunksize,mtx);
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

               
            cam.enable_FourD(&Dose_image, &chunk_data, det_bin,scan_bin, chunksize,mtx); 
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
            cam.enable_FourD(&Dose_image, &chunk_data, det_bin,scan_bin, chunksize,mtx); 
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
                cam.enable_compress(&Dose_image,&chunk_data, chunksize,det_bin,mtx);
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
                cam.enable_compress(&Dose_image,&chunk_data,chunksize,det_bin,mtx);
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
    h5file.close();
}

template<int BitDepth>
void FourD<BitDepth>::allocate_chunk()
{
    switch (camera)
     {
        case CAMERA::ADVAPIX:
        {   
            diff_pattern_size = 256/det_bin*256/det_bin;
            diff_pattern_length = 256/det_bin;
            break;
        }
        case CAMERA::CHEETAH:
        {
            diff_pattern_size = 512/det_bin*512/det_bin;
            diff_pattern_length = 512/det_bin;
            break;
        }
        case CAMERA::MERLIN:
        {
            diff_pattern_size = n_cam/det_bin*n_cam/det_bin;
            diff_pattern_length = n_cam/det_bin;
            break;
        }
    }
    std:: cout << chunksize/scan_bin << " " << nx/scan_bin << " " << diff_pattern_size << std::endl;
    double bits = chunksize/scan_bin*nx/scan_bin*diff_pattern_size*bitdepth;
    std::cout << "Chunk size: " << bits/8000000. << "MB" << std::endl;

    chunk_data[0].assign(chunksize/scan_bin*nx/scan_bin*diff_pattern_size,0);
    chunk_data[1].assign(chunksize/scan_bin*nx/scan_bin*diff_pattern_size,0);

}

template<int BitDepth>
void FourD<BitDepth>::reset()
{
    rc_quit = false;
    fr_freq = 0;

    // Initializations
    nxy = nx * ny;
    id_image = 0;
    fr_total = nxy * rep;
    fr_count = 0;

    // Allocate memory for image arrays
    if (b_first_run) {
        Dose_image.assign(nxy/(scan_bin*scan_bin), 0);
        b_first_run = false;
    }

    // Data Processing Progress
    *processor_line = 0;
    *preprocessor_line = 0;

}


template<int BitDepth>
void FourD<BitDepth>::line_processor(
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

        if (b_save_4D && *processor_line%chunksize == 0 && *processor_line > 0){
            write_and_clean(*processor_line);
        }

        progress_percent = prog_mon->progress_percent;


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

        if (b_save_4D){
            write_and_clean(nx);
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

template<int BitDepth>
void FourD<BitDepth>::write_and_clean(int line)
{
    int chunk_id = (line-chunksize)/(chunksize);
    int mod_chunk_id = chunk_id%2;
    std::lock_guard<std::mutex> lock(mtx[mod_chunk_id]);
    write_chunk(chunk_data[mod_chunk_id], chunk_id);
    chunk_data[mod_chunk_id].assign(chunksize/scan_bin*nx/scan_bin*diff_pattern_size,0);
       
}

template<int BitDepth>
void FourD<BitDepth>::init_4D_file_32() {

    const char* __DATASET_NAME = "shape";
    hsize_t __dims[1] = {4};
    H5::DataSpace __dataspace(1, __dims);
    H5::DataSet __dataset = h5file.createDataSet(__DATASET_NAME,H5::PredType::NATIVE_UINT16, __dataspace);
    hsize_t shape[4] = {ny/scan_bin, nx/scan_bin, diff_pattern_length, diff_pattern_length};
    __dataset.write(shape, H5::PredType::NATIVE_UINT16);

    b_save_4D = true;

    const char* DATASET_NAME_4 = "4D";
    
    hsize_t dims_4[4] = {static_cast<hsize_t>(nx/scan_bin),static_cast<hsize_t>(ny/scan_bin),diff_pattern_length,diff_pattern_length};
    H5::DataSpace dataspace_4(4, dims_4);

    H5::DSetCreatPropList prop_4;
    hsize_t chunk_dims_4[4] = {chunksize/scan_bin,nx/scan_bin,diff_pattern_length, diff_pattern_length};
    prop_4.setChunk(4, chunk_dims_4);
    prop_4.setDeflate(deflate_factor); // 0-9, 9 is maximum compression (slower), 0 is no compression
    dataset4D = h5file.createDataSet(DATASET_NAME_4,H5::PredType::NATIVE_UINT32, dataspace_4,prop_4);
    }


template<int BitDepth>
void FourD<BitDepth>::init_4D_file_16() {

    const char* __DATASET_NAME = "shape";
    hsize_t __dims[1] = {4};
    H5::DataSpace __dataspace(1, __dims);
    H5::DataSet __dataset = h5file.createDataSet(__DATASET_NAME,H5::PredType::NATIVE_UINT16, __dataspace);
    hsize_t shape[4] = {ny/scan_bin, nx/scan_bin, diff_pattern_length, diff_pattern_length};
    __dataset.write(shape, H5::PredType::NATIVE_UINT16);

    b_save_4D = true;

    const char* DATASET_NAME_4 = "4D";
    
    hsize_t dims_4[4] = {static_cast<hsize_t>(nx/scan_bin),static_cast<hsize_t>(ny/scan_bin),diff_pattern_length,diff_pattern_length};
    H5::DataSpace dataspace_4(4, dims_4);

    H5::DSetCreatPropList prop_4;
    hsize_t chunk_dims_4[4] = {chunksize/scan_bin,nx/scan_bin,diff_pattern_length, diff_pattern_length};
    prop_4.setChunk(4, chunk_dims_4);
    prop_4.setDeflate(deflate_factor); // 0-9, 9 is maximum compression (slower), 0 is no compression
    dataset4D = h5file.createDataSet(DATASET_NAME_4,H5::PredType::NATIVE_UINT16, dataspace_4,prop_4);
    }
 
template<int BitDepth>
void FourD<BitDepth>::init_4D_file_8() {

    try{
        const char* __DATASET_NAME = "shape";
        hsize_t __dims[1] = {4};
        H5::DataSpace __dataspace(1, __dims);
        H5::DataSet __dataset = h5file.createDataSet(__DATASET_NAME,H5::PredType::NATIVE_UINT16, __dataspace);
        hsize_t shape[4] = {ny/scan_bin, nx/scan_bin, diff_pattern_length, diff_pattern_length};
        __dataset.write(shape, H5::PredType::NATIVE_UINT16);

        b_save_4D = true;

        const char* DATASET_NAME_4 = "4D";

        hsize_t dims_4[4] = {static_cast<hsize_t>(nx/scan_bin),static_cast<hsize_t>(ny/scan_bin),diff_pattern_length,diff_pattern_length};
        H5::DataSpace dataspace_4(4, dims_4);

        H5::DSetCreatPropList prop_4;
        hsize_t chunk_dims_4[4] = {chunksize/scan_bin,nx/scan_bin,diff_pattern_length, diff_pattern_length};
        // hsize_t chunk_dims_4[4] = {chunksize/scan_bin,chunksize/scan_bin,diff_pattern_length, diff_pattern_length};
        prop_4.setChunk(4, chunk_dims_4);
        prop_4.setDeflate(deflate_factor); // 0-9, 9 is maximum compression (slower), 0 is no compression
        dataset4D = h5file.createDataSet(DATASET_NAME_4,H5::PredType::NATIVE_UINT8, dataspace_4,prop_4);
    }
    catch (const H5::FileIException& e) {
        std::cerr << "File error: " << e.getDetailMsg() << std::endl;
        throw;} 
    catch (const H5::DataSetIException& e) {
        std::cerr << "Dataset error: " << e.getDetailMsg() << std::endl;
        throw;} 
    catch (const H5::DataSpaceIException& e) {
        std::cerr << "Dataspace error: " << e.getDetailMsg() << std::endl;
        throw;}
    catch (const H5::PropListIException& e) {
        std::cerr << "Property list error: " << e.getDetailMsg() << std::endl;
        throw;}
    catch (const std::exception& e) {
        std::cerr << "Standard exception: " << e.what() << std::endl;
        throw;} 
    catch (...) {
    std::cerr << "Unknown error occurred during HDF5 file initialization." << std::endl;
    throw;}
}

template<int BitDepth>
void FourD<BitDepth>::init_4D_file()
{
    switch (bitdepth)
    {
    case 8:
        init_4D_file_8();
        break;
    case 16:
        init_4D_file_16();
        break;
    case 32:
        init_4D_file_32();
        break;
    }
}

template<int BitDepth>
void FourD<BitDepth>::write_chunk(const std::vector<uint32_t>& chunk, hsize_t index) {
    hsize_t offset4D[4];
    offset4D[1] = 0;
    offset4D[0] = (index*chunksize)/scan_bin;
    offset4D[2] = 0;
    offset4D[3] = 0;

    if (offset4D[1] < ny/scan_bin - chunksize/scan_bin)
    {
        hsize_t chunk_dims4D[4] = {chunksize/scan_bin,nx/scan_bin,diff_pattern_length, diff_pattern_length};
        H5::DataSpace filespace4D = dataset4D.getSpace();
        filespace4D.selectHyperslab(H5S_SELECT_SET, chunk_dims4D, offset4D);
        H5::DataSpace memspace4D(4, chunk_dims4D);
        dataset4D.write(chunk.data(), H5::PredType::NATIVE_UINT32, memspace4D, filespace4D);
    }
}

template<int BitDepth>
void FourD<BitDepth>::write_chunk(const std::vector<uint16_t>& chunk, hsize_t index) {
    hsize_t offset4D[4];
    offset4D[1] = 0;
    offset4D[0] = (index*chunksize)/scan_bin;
    offset4D[2] = 0;
    offset4D[3] = 0;

    if (offset4D[1] < ny/scan_bin - chunksize/scan_bin)
    {
        hsize_t chunk_dims4D[4] = {chunksize/scan_bin,nx/scan_bin,diff_pattern_length, diff_pattern_length};
        H5::DataSpace filespace4D = dataset4D.getSpace();
        filespace4D.selectHyperslab(H5S_SELECT_SET, chunk_dims4D, offset4D);
        H5::DataSpace memspace4D(4, chunk_dims4D);
        dataset4D.write(chunk.data(), H5::PredType::NATIVE_UINT16, memspace4D, filespace4D);
    }
}

template<int BitDepth>
void FourD<BitDepth>::write_chunk(const std::vector<uint8_t>& chunk, hsize_t index) {
    hsize_t offset4D[4];
    offset4D[1] = 0;
    offset4D[0] = (index*chunksize)/scan_bin;
    offset4D[2] = 0;
    offset4D[3] = 0;

    if (offset4D[1] < ny/scan_bin - chunksize/scan_bin)
    {
        hsize_t chunk_dims4D[4] = {chunksize/scan_bin,nx/scan_bin,diff_pattern_length, diff_pattern_length};
        // hsize_t chunk_dims4D[4] = {chunksize/scan_bin,chunksize/scan_bin,diff_pattern_length, diff_pattern_length};
        H5::DataSpace filespace4D = dataset4D.getSpace();
        filespace4D.selectHyperslab(H5S_SELECT_SET, chunk_dims4D, offset4D);
        H5::DataSpace memspace4D(4, chunk_dims4D);
        dataset4D.write(chunk.data(), H5::PredType::NATIVE_UINT8, memspace4D, filespace4D);
    }
}

template<int BitDepth>
void FourD<BitDepth>::save_dose_image(){

    const char*  DATASET_NAME = "dose_image";
    
    hsize_t dims[2] = {nx/scan_bin , ny/scan_bin};
    H5::DataSpace dataspace(2, dims);
    H5::DataSet dataset_ = h5file.createDataSet(DATASET_NAME, H5::PredType::NATIVE_UINT64, dataspace);
    dataset_.write(Dose_image.data(), H5::PredType::NATIVE_UINT64);
}


template class FourD<8>;
template class FourD<16>;
template class FourD<32>;