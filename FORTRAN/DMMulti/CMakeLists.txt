#Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

cmake_minimum_required(VERSION 2.8)

# Nombre del proyecto
project(DMMulti)
# Lenguaje de programacion
enable_language (Fortran)
# EL lugar donde se guardaran los ejecutables (dentro de build)
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
# Estas dos líneas permiten agregar el archivo de preprocesamiento
configure_file ("${PROJECT_SOURCE_DIR}/src/control.h" "${PROJECT_BINARY_DIR}/libs/control.h")
include_directories("${PROJECT_BINARY_DIR}/libs")

#===============================================================================
# OPCIONES A TRAVÉS DE CMAKE ( cmake -Dprofile=ON , por ejemplo)
#===============================================================================

option(openmp   "Enable shared-memory parallelism with OpenMP"   OFF)
option(profile  "Compile with profiling flags"                   OFF)
option(debug    "Compile with debug flags"                       OFF)
option(optimize "Turn on all compiler optimization flags"        OFF)
option(verbose  "Create verbose Makefiles"                       OFF)

if (verbose)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

#===============================================================================
# DETECTA LA VARIBLA FC Y HABILITA EL USO DE MPI 
#===============================================================================

set(MPI_ENABLED FALSE)
if($ENV{FC} MATCHES "mpi[^/]*$")
  message("-- Detected MPI wrapper: $ENV{FC}")
  add_definitions(-DMPI)
  set(MPI_ENABLED TRUE)
endif()
  
#============================================================================
# A su vez, se agregan los flags despendiendo de las opciones especificadas
# para la compilación ( cmake -D[opt] = ON)

message("-- Compilador utilizado: ${CMAKE_Fortran_COMPILER}")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # GNU Fortran compiler options
  set(f90flags  "-cpp -std=f2003 -fbacktrace")
  #  set(ldflags   "-lfftw3") # Se agregó directamente en src/CMakeLists.txt
  if(debug)
    set(f90flags  "-g -Wall -pedantic -fbounds-check -ffpe-trap=invalid,overflow,underflow ${f90flags}")
    set(ldflags   "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-fopenmp ${f90flags}")
    set(ldflags   "-fopenmp ${ldflags}")
  endif()
  
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Intel Fortran compiler options
  set(f90flags  "-fpp -std08 -assume byterecl -traceback")
  if(debug)
    set(f90flags "-g -warn -ftrapuv -fp-stack-check -check all -fpe0 ${f90flags}")
    set(ldflags  "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-openmp ${f90flags}")
    set(ldflags   "-openmp ${ldflags}")
  endif()
endif()

#============================================================================
# CCPIO ARCHIVOS SECUNDARIOS EL LUGAR DONDE VOY A EJECUTAR EL PROGRAMA
#============================================================================
# Cosas que se podrían hacer con INSTALL
file (COPY "${PROJECT_SOURCE_DIR}/src/input.par" DESTINATION ${EXECUTABLE_OUTPUT_PATH}) 
# Se guarda en la variable scripts_py todos los archivos python de la carpeta scripts
file (GLOB scripts_py "${PROJECT_SOURCE_DIR}/scripts/*.py")
# Se copian los scripts de python a la carpeta del ejecutable
file (COPY ${scripts_py} DESTINATION ${EXECUTABLE_OUTPUT_PATH})

#============================================================================
# LISTA DE LOS DIRECTORIOS A INCLUIR
#============================================================================
add_subdirectory(src)
add_subdirectory(libs)

