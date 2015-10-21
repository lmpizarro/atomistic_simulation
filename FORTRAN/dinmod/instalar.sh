#!/bin/bash

mkdir build
cd build
#cmake -Doptimize=ON ..
cmake -Ddebug=OFF -Dopenmp=ON  -Dverbose=OFF -Dprofile=OFF -Doptimize=ON ..
make 

# Para pasarle opciones de compilación al cmake:
# cmake -D[opc]=ON ..
# donde [opc] pueden ser:
# debug
# profile
# optimize
# openmp
# verbose

# Por default, las opciones del falg son -cpp -std=f2003 -fbacktrace
# Todo esto se puede modificar desde el archivo CMakelists.txt
