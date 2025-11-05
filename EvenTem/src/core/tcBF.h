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

#ifndef TCBF_H
#define TCBF_H

#include "LiveProcessor.h"

class tcBF : public LiveProcessor
{
private: 

public: 
    std::vector<size_t> BF_image;
    std::vector<std::vector<size_t>> tcBF_stack;

    bool allow_torch = false;
    bool allow_cuda = false;

    void run();
    void reset();

    int N_angles;

    void line_processor(
        size_t &img_num,
        size_t &first_frame,
        size_t &end_frame,
        ProgressMonitor *p_prog_mon,
        size_t &fr_total_u,
        BoundedThreadPool *pool
    );

    void set_detector_mask(py::array_t<int> mask);
    std::vector<int> detector_mask;

    // Constructor
    tcBF(int repetitions) : LiveProcessor(repetitions)
    {
    };

    // Destructor
    ~tcBF(){};
};
#endif // !tcBF_H