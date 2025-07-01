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

#ifndef EELS_H
#define EELS_H

#include "LiveProcessor.h"


class EELS : public LiveProcessor
{
private: 

public: 

    std::vector<std::vector<uint64_t>> EELS_data;
    std::vector<uint64_t> EELS_image;
    std::vector<std::vector<std::vector<uint64_t>>> EELS_data_stack; 
    std::vector<std::vector<uint64_t>> EELS_image_stack;

    void run();
    void reset();

    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *p_prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    // Constructor
    EELS(int repetitions) : LiveProcessor(repetitions)
    {
    };

    // Destructor
    ~EELS(){};
};
#endif // !EELS_H