#!/bin/bash

mkdir build
cd build
cmake ..
#cmake -D CMAKE_BUILD_TYPE:STRING=DEBUG ..
make
