cd BuildTools
rmdir /S /Q build_MSVC
mkdir build_MSVC
cd build_MSVC
cmake -DGPRI_OPTION=OFF -DTORCH=OFF -DPIXET=OFF -DLOG=OFF -DCMAKE_CUDA_FLAGS="--allow-unsupported-compiler" ..
msbuild eventem.vcxproj /p:Configuration=Release
copy Release\*.pyd "..\include" 
copy Release\*.pyd "..\..\EvenTem\binaries"