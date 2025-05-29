cd BuildTools
rmdir /S /Q GNUbuild
mkdir GNUbuild
cd GNUbuild
cmake -G "MinGW Makefiles" -DCMAKE_TOOLCHAIN_FILE=../toolchain.cmake -DGPRI_OPTION=ON -DTORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ..
mingw32-make 
