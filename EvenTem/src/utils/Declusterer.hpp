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

#ifndef DECLUSTERER_HPP
#define DECLUSTERER_HPP

#include <cmath>
#include <atomic>
#include <vector>
#include <array>
#include <thread>
#include <chrono>
#include <functional>
#include <algorithm>
#include <cstdlib>

#include "dtype_Electron.hpp"
#include "BoundedThreadPool.hpp"

#include "Logger.hpp"


#pragma pack(push, 1)
struct cluster_event
{
    uint16_t kx;
    uint16_t ky;
    uint16_t rx;
    uint16_t ry;
    uint16_t id_image;
    uint64_t toa;
    uint16_t tot;
};
#pragma pack(pop)

class Declusterer
{
private:

    static const int n_buffer = 128;

    uint64_t dtime;
    uint16_t dspace; 
    int cluster_range;

    int x_crop;
    int y_crop;
    int scan_bin;
    int det_bin;

    dtype_Electron electron;
    std::ofstream* p_file;

    std::thread decluster_thread;
    std::thread write_thread;

    BoundedThreadPool *pool = new BoundedThreadPool;

    void decluster(int _buffer_id)
    {
        int upper_lim;
        int clustersize = 1;
        int _bs = buffer[_buffer_id]->size();
        std::vector<bool> used(_bs, false);
        std::shared_ptr<std::vector<int>> lcl_keep = std::make_shared<std::vector<int>>();

        for (int i = 0; i < _bs; ++i) 
        {
            if (!used[i]){
                upper_lim = std::min(i + cluster_range, _bs);

                for (int j = i+1; j < upper_lim; ++j) 
                {
                    if (!used[j]){
                        bool x_condition = ((*buffer[_buffer_id])[i].kx > (*buffer[_buffer_id])[j].kx) ? 
                                        ((*buffer[_buffer_id])[i].kx - (*buffer[_buffer_id])[j].kx <= dspace) : 
                                        ((*buffer[_buffer_id])[j].kx - (*buffer[_buffer_id])[i].kx <= dspace);

                        bool y_condition = ((*buffer[_buffer_id])[i].ky > (*buffer[_buffer_id])[j].ky) ? 
                                        ((*buffer[_buffer_id])[i].ky - (*buffer[_buffer_id])[j].ky <= dspace) : 
                                        ((*buffer[_buffer_id])[j].ky - (*buffer[_buffer_id])[i].ky <= dspace);

                        bool t_condition = ((*buffer[_buffer_id])[i].toa > (*buffer[_buffer_id])[j].toa) ? 
                                        ((*buffer[_buffer_id])[i].toa - (*buffer[_buffer_id])[j].toa <= dtime) : 
                                        ((*buffer[_buffer_id])[j].toa - (*buffer[_buffer_id])[i].toa <= dtime);
                                        
                        if (x_condition && y_condition && t_condition)
                        {
                            used[j] = true;
                            clustersize++;
                        }
                    }
                }
                if (clustersize < max_clustersize) (*p_clustersize_histogram)[clustersize]++;
                clustersize = 1;
                lcl_keep->push_back(i);
            }
            used[i] = true;
        }

        n_electrons_kept += lcl_keep->size();
        ++this->n_buffer_declustered;
        keep[_buffer_id] = lcl_keep;
    };

    
    void write_to_file(int _buffer_id)
    {
        for (int i: *keep[_buffer_id]){
            electron.kx = (*buffer[_buffer_id])[i].kx/det_bin;
            electron.ky = (*buffer[_buffer_id])[i].ky/det_bin;
            electron.rx = (*buffer[_buffer_id])[i].rx/scan_bin;
            electron.ry = (*buffer[_buffer_id])[i].ry/scan_bin;
            electron.id_image = (*buffer[_buffer_id])[i].id_image;

            if (electron.rx < (x_crop/scan_bin) && electron.ry < (y_crop/scan_bin))
            {
                (*p_file).write((const char *)&electron, sizeof(electron));
            }
        }
    };

    void schedule_declustering()
    {
        while ((still_reading) || (n_buffer_declustered < n_buffer_filled))
        {
            if (n_buffer_declustered < n_buffer_filled)
            {
                // pool->push_task([=]{decluster(n_buffer_declustered % n_buffer);}); 
                #ifdef LOG 
                    Logger::getInstance().log("Declustering buffer " + std::to_string(n_buffer_declustered) + " of " + std::to_string(n_buffer_filled));
                #endif
                decluster(n_buffer_declustered % n_buffer);   
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::microseconds(10));
            }
        }
        still_processing = false; 
    };

    void schedule_writing()
    {
        int _buffer_id;
        while (still_processing || n_buffer_written < n_buffer_declustered)
        {
            if (n_buffer_written < n_buffer_declustered)
            {
                _buffer_id = n_buffer_written % n_buffer;
                write_to_file(_buffer_id);
                buffer[_buffer_id]->clear();
                #ifdef LOG 
                    Logger::getInstance().log("clearing buffer " + std::to_string(n_buffer_written) + " of " + std::to_string(n_buffer_filled));
                #endif
                ++n_buffer_written;
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::microseconds(10));
            }
        }
        still_writing = false;
    };


public:
    std::vector<cluster_event> *buffer[n_buffer];
    std::shared_ptr<std::vector<int>> keep[n_buffer];
    
    bool still_reading = true;
    bool still_processing = true;
    bool still_writing = true;
    int n_buffer_filled = 0;
    int buffer_id_filling = 0;
    std::atomic<int> n_buffer_declustered = 0;
    int n_buffer_written = 0;
    std::atomic<int> n_electrons_kept = 0;
    std::vector<int> *p_clustersize_histogram;
    int max_clustersize = 0;

    void init(uint64_t dtime, uint16_t dspace, int cluster_range,int x_crop,int y_crop,int scan_bin, int det_bin, std::ofstream& _p_file, int n_threads, std::vector<int> *_p_clustersize_histogram)
    {
        this->dtime = dtime;
        this->dspace = dspace;
        this->cluster_range = cluster_range;
        this->p_file = &_p_file;
        this->scan_bin = scan_bin;
        this->det_bin = det_bin;
        this->x_crop = x_crop;
        this->y_crop = y_crop;
        this->p_clustersize_histogram = _p_clustersize_histogram;
        max_clustersize = p_clustersize_histogram->size();

        std::cout << "Declustering param: dtime = " << dtime << ", dspace = " << dspace << ", cluster_range = " << cluster_range << std::endl;

        // pool->init(n_threads,n_threads);
    };

    void set_buffer_read()
    {
        ++n_buffer_filled;
        buffer_id_filling = n_buffer_filled % n_buffer;

        while (n_buffer_filled - n_buffer_written >= n_buffer)
        {
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
    };

    void run()
    {
        decluster_thread = std::thread(&Declusterer::schedule_declustering, this);
        write_thread = std::thread(&Declusterer::schedule_writing, this);
    };

    void terminate()
    {
        while (!write_thread.joinable() || !decluster_thread.joinable()) {std::this_thread::sleep_for(std::chrono::milliseconds(1));}
        write_thread.join();
        decluster_thread.join();
    };

    Declusterer(){
        for (int i = 0; i < n_buffer; ++i) {
            buffer[i] = new std::vector<cluster_event>();
            keep[i] = std::make_shared<std::vector<int>>();
        }
    }
};

#endif // DECLUSTERER_HPP


 // void decluster(int _buffer_id)
    // {
    //     int lower_lim;
    //     int upper_lim;
    //     int len_kx = buffer[_buffer_id]->size();
    //     std::vector<bool> visited;
    //     std::vector<std::vector<int>> adjacency_list = std::vector<std::vector<int>>();
    //     adjacency_list.assign(len_kx, std::vector<int>());
    //     visited.assign(len_kx, false);

    //     // Construct adjacency list
    //     for (int i = 0; i < len_kx; ++i) 
    //     {
    //         lower_lim = std::max(i - cluster_range, 0);
    //         upper_lim = std::min(i + cluster_range, len_kx);

    //         for (int j = lower_lim; j < upper_lim; ++j) 
    //         {
    //             if (std::abs((*buffer[_buffer_id])[i].kx - (*buffer[_buffer_id])[j].kx) <= dspace && std::abs((*buffer[_buffer_id])[i].ky - (*buffer[_buffer_id])[j].ky) <= dspace && std::llabs(static_cast<long long>((*buffer[_buffer_id])[i].toa - (*buffer[_buffer_id])[j].toa)) <= dtime)
    //             {
    //                 adjacency_list[i].push_back(j);
    //                 adjacency_list[j].push_back(i);
    //             }
    //         }
    //     }

    //     for (int i = 0; i < len_kx; ++i)
    //     {
    //         if (!visited[i])
    //         {
    //             visited[i] = true;
    //             for (int neighbor : adjacency_list[i]) visited[neighbor] = true;
    //             keep[_buffer_id]->push_back(i);
    //         }
    //     }

    //     n_electrons_kept += keep[_buffer_id]->size();
    //     ++this->n_buffer_declustered;
    //     std::cout << "Buffer " << _buffer_id << " declustered" << std::endl;


    //     // for (int i = 0; i < len_kx; ++i)
    //     // {
    //     //     if (!visited[i])
    //     //     {
    //     //         visited[i] = true;
    //     //         TOAs.push_back((*buffer[buffer_id])[i].toa);
    //     //         TOTs.push_back((*buffer[buffer_id])[i].tot);
    //     //         INDEXs.push_back(i);
    //     //         for (int neighbor : adjacency_list[i])
    //     //         { 
    //     //             visited[neighbor] = true;
    //     //             TOAs.push_back((*buffer[buffer_id])[neighbor].toa);
    //     //             TOTs.push_back((*buffer[buffer_id])[neighbor].tot);
    //     //             INDEXs.push_back(neighbor);
    //     //         }
    
    //             // // centerElementIter16 = std::max_element(TOTs.begin(), TOTs.end());
    //             // // centerIndex = std::distance(TOTs.begin(), centerElementIter16);
    //             // // centerElementIter64 = std::min_element(TOAs.begin(), TOAs.end());
    //             // // centerIndex = std::distance(TOAs.begin(), centerElementIter64);
    //             // // keep.push_back(INDEXs[centerIndex]);
    //             // // keep.push_back(i);
    //             // keep.push_back(INDEXs[0]);

    //             // TOAs.clear();
    //             // TOTs.clear();
    //             // INDEXs.clear();
    //     //     }
    //     // }
    // };

// void decluster_alt()
// {
//     len_kx = kx_vec.size();
//     adjacency_list.assign(len_kx, std::set<int>());
//     visited.assign(len_kx, false);
    
//     // Construct adjacency list
//     for (int i = 0; i < len_kx; ++i) 
//     {
//         lower_lim = std::max(i - cluster_range, 0);
//         upper_lim = std::min(i + cluster_range, len_kx);

//         for (int j = lower_lim; j < upper_lim; ++j) 
//         {
//             if (std::abs(kx_vec[i] - kx_vec[j]) <= dspace && std::abs(ky_vec[i] - ky_vec[j]) <= dspace && std::llabs(static_cast<long long>(toa_vec[i] - toa_vec[j])) <= dtime) 
//             {
//                 adjacency_list[i].insert(j);
//                 adjacency_list[j].insert(i);
//             }
//         }
//     }
    

//     for (int i = 0; i < len_kx; ++i)
//     {
//         if (!visited[i])
//         {
//             visited[i] = true;
//             KXs.push_back(kx_vec[i]);
//             KYs.push_back(ky_vec[i]);
//             RXs.push_back(rx_vec[i]);
//             RYs.push_back(ry_vec[i]);
//             TOAs.push_back(toa_vec[i]);
//             TOTs.push_back(tot_vec[i]);
//             for (int neighbor : adjacency_list[i])
//             { 
//                 visited[neighbor] = true;
//                 KXs.push_back(kx_vec[neighbor]);
//                 KYs.push_back(ky_vec[neighbor]);
//                 RXs.push_back(rx_vec[neighbor]);
//                 RYs.push_back(ry_vec[neighbor]);
//                 TOAs.push_back(toa_vec[neighbor]);
//                 TOTs.push_back(tot_vec[neighbor]);
//             }

//         // centerElementIter16 = std::max_element(TOTs.begin(), TOTs.end());
//         // centerIndex = std::distance(TOTs.begin(), centerElementIter16);
//         centerIndex = 0;
//         kx_keep.push_back(KXs[centerIndex]);
//         ky_keep.push_back(KYs[centerIndex]);
//         rx_keep.push_back(RXs[centerIndex]);
//         ry_keep.push_back(RYs[centerIndex]);
        

//         TOAs.clear();
//         TOTs.clear();
//         KXs.clear();
//         KYs.clear();
//         RXs.clear();
//         RYs.clear();
//         }
//     }
// };