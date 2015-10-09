# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

cmake_minimum_required(VERSION 2.8)

# Nombre del proyecto
project(ising)
# Lenguaje de programacion
enable_language (Fortran)
# Nombre del paquete ¿para qué?
set(PACKAGE ising)
# EL lugar donde se guardaran los ejecutables (dentro de build)
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )

#============================================================================
# FFLAGS depend on the compiler
#============================================================================
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
#  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  message("Fortran compiler: " ${CMAKE_Fortran_COMPILER})
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#============================================================================

#============================================================================
# CCPIO ARCHIVOS SECUNDARIOS EL LUGAR DONDE VOY A EJECUTAR EL PROGRAMA
#============================================================================

file (COPY "${PROJECT_SOURCE_DIR}/src/parametros.dat" DESTINATION ${EXECUTABLE_OUTPUT_PATH}) 
file (COPY "${PROJECT_SOURCE_DIR}/scripts/corridas_paralelo.py" DESTINATION ${EXECUTABLE_OUTPUT_PATH})              
file (COPY "${PROJECT_SOURCE_DIR}/scripts/corridas.py" DESTINATION ${EXECUTABLE_OUTPUT_PATH})        
file (COPY "${PROJECT_SOURCE_DIR}/scripts/parametros.py" DESTINATION ${EXECUTABLE_OUTPUT_PATH})        


#============================================================================
# LISTA DE LOS DIRECTORIOS A INCLUIR
#============================================================================
add_subdirectory(src)
add_subdirectory(libs)
