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

#ifndef ROI4D_HPP
#define ROI4D_HPP

#include <vector>

class Roi4DBase {
public:
    virtual ~Roi4DBase() = default;
    virtual void init(int L_1, int L_0, int n_cam, int det_bin) = 0;
    virtual void* getData() = 0; 
    virtual void increment(int,int,int,int) = 0; 

};

template <typename T>
class Roi4D : public Roi4DBase {
public:
    std::vector<std::vector<std::vector<std::vector<T>>>> data;
    void init(int L_1, int L_0, int n_cam, int det_bin) {
        data.assign(L_1, std::vector<std::vector<std::vector<T>>>(
            L_0, std::vector<std::vector<T>>(n_cam / det_bin, std::vector<T>(n_cam / det_bin, 0))));
    }
    void* getData() override {
        return &data;
    }
    void increment(int x, int y, int kx, int ky) override {
        data[x][y][kx][ky]++;
    }
};

#endif // ROI4D_HPP