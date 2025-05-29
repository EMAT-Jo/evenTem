#!/bin/bash
cd BuildTools
rm -rf UNIXbuild
mkdir UNIXbuild
cd UNIXbuild
cmake -DGPRI_OPTION=OFF -D_TORCH=OFF -DPIXET=OFF -DLOG=OFF ..
make -s

