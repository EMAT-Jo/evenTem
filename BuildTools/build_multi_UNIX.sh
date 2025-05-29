#!/bin/bash
cd BuildTools
mkdir build_all_py
cd build_all_py 
mkdir build_39
cd build_39
source /home/arno/anaconda3/bin/activate py39
cmake -DGPRI_OPTION=OFF -D_TORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
make -s
cd ..
mkdir build_310
cd build_310
source /home/arno/anaconda3/bin/activate py310
cmake -DGPRI_OPTION=OFF -D_TORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
make -s
cd ..
mkdir build_311
cd build_311
source /home/arno/anaconda3/bin/activate py311
cmake -DGPRI_OPTION=OFF -D_TORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
make -s
cd ..
mkdir build_312
cd build_312
source /home/arno/anaconda3/bin/activate py312
cmake -DGPRI_OPTION=OFF -D_TORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
make -s



