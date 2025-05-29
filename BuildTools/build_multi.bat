cd BuildTools
mkdir build_all_py
cd build_all_py 
mkdir build_39
cd build_39 
call conda activate py39
call cmake -DGPRI_OPTION=OFF -DTORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
call msbuild eventem.vcxproj /p:Configuration=Release /v:minimal
cd .. 
mkdir build_310
cd build_310
call conda activate py310
call cmake -DGPRI_OPTION=OFF -DTORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
call msbuild eventem.vcxproj /p:Configuration=Release /v:minimal
cd ..
mkdir build_311
cd build_311
call conda activate py311
call cmake -DGPRI_OPTION=OFF -DTORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
call msbuild eventem.vcxproj /p:Configuration=Release /v:minimal
cd ..
mkdir build_312
cd build_312
call conda activate py312
call cmake -DGPRI_OPTION=OFF -DTORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ../..
call msbuild eventem.vcxproj /p:Configuration=Release /v:minimal

