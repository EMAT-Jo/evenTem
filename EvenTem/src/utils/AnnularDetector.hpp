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


#ifndef AnnularDetector_hpp
#define AnnularDetector_hpp

#include <vector>

class AnnularDetector
{
public:
    // Properties
    std::vector<std::array<float, 2>> radia = {{0, 0}};
    std::vector<std::array<float, 2>> radia_sqr = {{0, 0}};
    std::vector<int> detector_image;
    int n_detectors = 1;

    // Methods
    void compute_detector(int nx_cam, int ny_cam, std::vector<std::array<float, 2>> &offsets)
    {
        float d2;
        detector_image.assign(nx_cam * ny_cam, 0);
        for (int iy = 0; iy < ny_cam; iy++)
        {
            for (int ix = 0; ix < nx_cam; ix++)
            {
                for (int i = 0; i < radia.size(); i++)
                {
                    d2 = pow((float)ix - offsets[i][0], 2) + pow((float)iy - offsets[i][1], 2);
                    if (d2 >= radia_sqr[i][0] && d2 <= radia_sqr[i][1])
                    {
                        detector_image[iy * nx_cam + ix] = 1;
                    }
                }
            }
        }
    }

    void set_radia(std::vector<float> r_inner, std::vector<float> r_outer)
    {
        n_detectors= r_inner.size();
        radia.resize(n_detectors);
        radia_sqr.resize( n_detectors);
        for (int i = 0; i <  n_detectors; i++)
        {
            radia[i][0] = r_inner[i];
            radia[i][1] = r_outer[i];
            radia_sqr[i][0] = pow(r_inner[i], 2);
            radia_sqr[i][1] = pow(r_outer[i], 2);
        }
    }

    // Constructor
    AnnularDetector(float r_inner, float r_outer){
        radia[0][0] = r_inner;
        radia[0][1] = r_outer;
        radia_sqr[0][0] = pow(r_inner, 2);
        radia_sqr[0][1] = pow(r_outer, 2);
    };
    // Destructor
    ~AnnularDetector(){};
};

#endif // AnnularDetector_hpp