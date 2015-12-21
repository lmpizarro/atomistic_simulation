#!/bin/bash

# descomentar y setear adecuadamente, ver abajo SETEO DE COMPILADOR
# export FC=/opt/intel/composer_xe_2013.5.192/bin/intel64/ifort

mkdir build
cd build
#cmake -Doptimize=ON ..
cmake -Ddebug=ON -Dopenmp=ON -Dverbose=OFF -Dprofile=OFF -Doptimize=OFF ..
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

######  SETEO DE COMPILADOR #######
# Para usar otro compilador distinto de gfortran se
# debe inicializar la variable de entorno FC con el 
# path al ejecutable del compilador que se quiere usar
# por ejemplo en sheldon se usa 
# export FC=/opt/intel/composer_xe_2013.5.192/bin/intel64/ifort
