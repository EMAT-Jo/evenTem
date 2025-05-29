# /* Copyright (C) 2025 Thomas Friedrich, Chu-Ping Yu, Arno Annys
#  * University of Antwerp - All Rights Reserved. 
#  * You may use, distribute and modify
#  * this code under the terms of the GPL3 license.
#  * You should have received a copy of the GPL3 license with
#  * this file. If not, please visit: 
#  * https://www.gnu.org/licenses/gpl-3.0.en.html
#  * 
#  * Authors: 
#  *   Thomas Friedrich <>
#  *   Chu-Ping Yu <>
#  *   Arno Annys <arno.annys@uantwerpen.be>
#  */

cmake_minimum_required(VERSION 3.16.0)
if(USE_CUDA)
    project(eventem LANGUAGES C CXX CUDA)
else()
    project(eventem LANGUAGES C CXX)
endif()
set(CMAKE_CXX_STANDARD 17)

option(GPRI_OPTION "option for enabling or disabling GPRI" OFF)
option(TORCH "option for enabling or disabling Torch use for frame processing (besides GPRI)" OFF)
option(PIXET "option for enabling or disabling Pixet" OFF)
option(LOG "option for enabling or disabling debug logging" OFF)
set(CMAKE_BUILD_TYPE Release)

set(SOURCES
    ../EvenTem/src/utils/Declusterer.hpp
    ../EvenTem/src/utils/SocketConnector.cpp
    ../EvenTem/src/utils/SocketConnector.h
    ../EvenTem/src/utils/FileConnector.cpp
    ../EvenTem/src/utils/FileConnector.h
    ../EvenTem/src/utils/ProgressMonitor.cpp
    ../EvenTem/src/utils/ProgressMonitor.h
    ../EvenTem/src/utils/Logger.hpp
    ../EvenTem/src/utils/Roi4D.hpp
    ../EvenTem/src/utils/AtomicWrapper.hpp
    ../EvenTem/src/utils/AnnularDetector.hpp
    ../EvenTem/src/core/LiveProcessor.cpp
    ../EvenTem/src/core/LiveProcessor.h
    ../EvenTem/src/core/Ricom.cpp 
    ../EvenTem/src/core/Ricom.h
    ../EvenTem/src/core/vSTEM.h
    ../EvenTem/src/core/vSTEM.cpp
    ../EvenTem/src/core/Pacbed.h
    ../EvenTem/src/core/Pacbed.cpp
    ../EvenTem/src/core/Var.h
    ../EvenTem/src/core/Var.cpp
    ../EvenTem/src/core/Roi.h
    ../EvenTem/src/core/Roi.cpp
    ../EvenTem/src/core/Electron.h 
    ../EvenTem/src/core/Electron.cpp
    ../EvenTem/src/core/EELS.h
    ../EvenTem/src/core/EELS.cpp
    ../EvenTem/src/core/FourD.h
     ../EvenTem/src/core/FourD.cpp
    ../EvenTem/src/detectors/Cheetah.hpp
    ../EvenTem/src/detectors/Cheetah_pixeltrig.hpp
    ../EvenTem/src/detectors/Advapix.hpp
    ../EvenTem/src/detectors/Timepix.hpp
    ../EvenTem/src/detectors/Simulated.hpp
    ../EvenTem/src/detectors/Merlin.hpp
    ../EvenTem/src/detectors/HDF5_DS.hpp
    ../EvenTem/src/pybind.cpp
)

#----------------------------------------------------------------------------------------------------------------
if (WIN32)
    find_package(pybind11 REQUIRED)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

    # set(FFTW3_DIR ./FFTW)
    # link_directories(${FFTW3_DIR})

    # add_definitions(-DH5_OPTION_ENABLED)
    set(HDF5_ROOT "C:/Program Files/HDF_Group/HDF5/1.14.3")
    find_package(HDF5 1.14.3 REQUIRED COMPONENTS CXX)

    set(BUILD_DIR ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
    message(STATUS "Build directory: ${BUILD_DIR}")
    add_definitions(-DNOMINMAX)
    add_definitions(-DEVENTEM)

    if (LOG)
        add_definitions(-DDBG_LOG)
        message(STATUS "Logging enabled")
    endif()

    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        MESSAGE(STATUS "Using ${CMAKE_CXX_COMPILER_ID} Compiler!")
        set(CMAKE_CXX_FLAGS "-Wall -Wextra -pthread -std=c++17")
        set(CMAKE_CXX_FLAGS_DEBUG "-g")
        set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive")
    endif()

    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
	MESSAGE(STATUS "Using ${CMAKE_CXX_COMPILER_ID} Compiler!")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /nologo")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GA")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O3")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MT")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qpar")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /favor:INTEL64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /LTCG:PGInstrument")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /LTCG:PGUpdate")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /LTCG:PGOptimize")
    endif()

    if (PIXET)
        list(APPEND SOURCES ../EvenTem/src/detectors/Advapix.cpp)
        message(STATUS "Pixet enabled")
        add_definitions(-DPIXET_ENABLED)
    endif()

    if (GPRI_OPTION)
        add_definitions(-DGPRI_OPTION_ENABLED)
        message(STATUS "GPRI support enabled")
    endif()  

    if (TORCH)
        add_definitions(-DFRAMEBASED_TORCH_ENABLED)
        message(STATUS "Torch support enabled")
    endif()

    if (GPRI_OPTION OR TORCH)
        set(__LIBTORCH_VERSION "2.6.0") 
        set(__CUDA_VERSION "cu126") 

        message(STATUS "CUDA support enabled")
        set(CMAKE_PREFIX_PATH "c:/libtorch${__LIBTORCH_VERSION}+${__CUDA_VERSION}/libtorch") 
        set(Torch_DIR "c:/libtorch${__LIBTORCH_VERSION}+${__CUDA_VERSION}/libtorch/share/cmake/Torch")  

        find_package(CUDA REQUIRED)
        find_package(Torch REQUIRED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}") 
    endif()

    if(GPRI_OPTION)
        enable_language(CUDA)
        list(INSERT SOURCES 0 ../EvenTem/src/core/GPRI.h ../EvenTem/src/core/GPRI.cpp ../EvenTem/src/core/GPRI.cu)
    endif()

    pybind11_add_module(eventem ${SOURCES}) 

    if(GPRI_OPTION OR TORCH)
        target_include_directories(eventem PRIVATE ${TORCH_INCLUDE_DIRS} ${CUDA_INCLUDE_DIRS})
        target_link_libraries(eventem PRIVATE ${TORCH_LIBRARIES} ${CUDA_LIBRARIES})
    endif()

    if (PIXET)
        target_link_directories(eventem PRIVATE include)
        target_link_libraries(eventem PRIVATE pxcore.lib)
    endif()

    # target_include_directories(eventem PRIVATE include ../EvenTem/src/core ../EvenTem/src/utils ../EvenTem/src/detectors)
    target_include_directories(eventem PRIVATE include ../EvenTem/src/core ../EvenTem/src/utils ../EvenTem/src/detectors ${Python3_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS}) #${FFTW3_DIR}
    # target_link_libraries(eventem PRIVATE pybind11::module ${CMAKE_DL_LIBS} ws2_32)
    target_link_libraries(eventem PRIVATE ${Python3_LIBRARIES} ${CMAKE_DL_LIBS} ${HDF5_LIBRARIES} ws2_32) #libfftw3f-3

endif (WIN32)
#----------------------------------------------------------------------------------------------------------------
if (APPLE)
    find_package(PkgConfig)
    find_package(pybind11 REQUIRED)
    find_package(HDF5 REQUIRED COMPONENTS CXX)
    
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -pthread")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

    if (GPRI_OPTION)
        pybind11_add_module(eventemTorch ${SOURCES})
        target_link_libraries(eventemTorch PRIVATE pybind11::module ${CMAKE_DL_LIBS})
        target_include_directories(eventemTorch  PRIVATE include ../EvenTem/src/core ../EvenTem/src/utils ../EvenTem/src/detectors)
    else()
        pybind11_add_module(eventem ${SOURCES})
        target_link_libraries(eventem PRIVATE pybind11::module ${CMAKE_DL_LIBS})
        target_include_directories(eventem  PRIVATE include ../EvenTem/src/core ../EvenTem/src/utils ../EvenTem/src/detectors)
    endif()

    if (GPRI_OPTION)
        add_definitions(-DEVENTEMTORCH)
        add_definitions(-DGPRI_OPTION_ENABLED)
        message(STATUS "GPRI support enabled")
        set(CMAKE_PREFIX_PATH "/Users/arnoannys/Programming/libtorch")
        set(Torch_DIR "/Users/arnoannys/Programming/libtorch/share/cmake/Torch")
        find_package(Torch REQUIRED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

        list(APPEND SOURCES ../EvenTem/src/core/GPRI.h ../EvenTem/src/core/GPRI.cpp)
        target_link_libraries(eventemTorch PRIVATE ${TORCH_LIBRARIES})
        target_include_directories(eventemTorch  PRIVATE ${TORCH_INCLUDE_DIRS})
    else()
        add_definitions(-DEVENTEM)
    endif()

    target_link_libraries(eventem PRIVATE ${HDF5_LIBRARIES} pybind11::module ${CMAKE_DL_LIBS})
    target_include_directories(eventem  PRIVATE include ../EvenTem/src/core ../EvenTem/src/utils ../EvenTem/src/detectors ${HDF5_INCLUDE_DIRS})
    

endif (APPLE) 
#----------------------------------------------------------------------------------------------------------------
if (UNIX AND NOT APPLE)
    set(pybind11_ROOT "/usr")
    find_package(PkgConfig)
    find_package(pybind11 REQUIRED)
    add_definitions(-DNOMINMAX)

    set(CMAKE_CUDA_COMPILER "/usr/local/cuda-11.8/bin/nvcc")

    message(STATUS "HDF5 support enabled")
    set(HDF5_ROOT "/usr")
    find_package(HDF5 REQUIRED COMPONENTS CXX)

    if (TORCH)
        add_definitions(-DTORCH_ENABLED)
        message(STATUS "Torch support enabled")
    endif()

    set(CMAKE_CXX_FLAGS "-w -Wall -Wextra -pthread")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

    if (GPRI_OPTION)
        add_definitions(-DTORCH_ENABLED)
        add_definitions(-DGPRI_OPTION_ENABLED)
        message(STATUS "GPRI support enabled")

        if (USE_CUDA)
            add_definitions(-DUSE_CUDA)
            set(CMAKE_CUDA_HOST_COMPILER "/usr/bin/g++")
            message(STATUS "CUDA support enabled")
            set(CMAKE_PREFIX_PATH "/home/arno/libtorch$ENV{__LIBTORCH_VERSION}+$ENV{__CUDA_VERSION}/libtorch") 
            set(Torch_DIR "/home/arno/libtorch$ENV{__LIBTORCH_VERSION}+$ENV{__CUDA_VERSION}/libtorch/share/cmake/Torch") 
        else()
            message(STATUS "CUDA support disabled")
            set(CMAKE_PREFIX_PATH "/home/arno/libtorch$ENV{__LIBTORCH_VERSION}+cpu/libtorch") 
            set(Torch_DIR "/home/arno/libtorch$ENV{__LIBTORCH_VERSION}+cpu/libtorch/share/cmake/Torch") 
        endif()
        
        find_package(Torch REQUIRED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")  
    endif()  


    if(GPRI_OPTION)
        add_definitions(-DEVENTEMTORCH)

        if (USE_CUDA)
            enable_language(CUDA)

            list(INSERT SOURCES 0 ../EvenTem/src/core/GPRI.h ../EvenTem/src/core/GPRI.cpp ../EvenTem/src/core/GPRI.cu)
        endif()

    else()
        add_definitions(-DEVENTEM)
    endif()


    if (GPRI_OPTION AND USE_CUDA)
        pybind11_add_module(eventemTorch ${SOURCES})
        target_include_directories(eventemTorch PRIVATE include ../EvenTem/src/core ../EvenTem/src/utils ../EvenTem/src/detectors)
        target_link_libraries(eventemTorch PRIVATE pybind11::module)
    else()
        pybind11_add_module(eventem ${SOURCES})
        target_include_directories(eventem PRIVATE include ../EvenTem/src/core ../EvenTem/src/utils ../EvenTem/src/detectors)
        target_link_libraries(eventem PRIVATE pybind11::module)
    endif()
        
    target_include_directories(eventem PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(eventem PRIVATE ${HDF5_LIBRARIES})
    
    if(GPRI_OPTION)
        target_include_directories(eventemTorch PRIVATE ${TORCH_INCLUDE_DIRS})
        target_link_libraries(eventemTorch PRIVATE ${TORCH_LIBRARIES})

        message(STATUS "Torch libraries: ${TORCH_LIBRARIES}")
        message(STATUS "Torch include directories: ${TORCH_INCLUDE_DIRS}")
    endif()

endif (UNIX AND NOT APPLE)
