#!/bin/bash

mkdir build
cd build
#cmake ..
cmake -Dverbose=OFF  ..
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
