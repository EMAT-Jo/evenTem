#!/bin/bash
cd BuildTools
rm -rf UNIXbuildTorch
mkdir UNIXbuildTorch
cd UNIXbuildTorch
export CUDACXX=/usr/local/cuda-11.8/bin/nvcc
export __LIBTORCH_VERSION="2.4.0"
export __CUDA_VERSION="cu118" 
# export CUDACXX=/usr/local/cuda-12.4/bin/nvcc
# export __LIBTORCH_VERSION="2.5.1"
# export __CUDA_VERSION="cu124" 
cmake -DGPRI_OPTION=ON -D_TORCH=OFF -DPIXET=OFF -DLOG=OFF ..
make

