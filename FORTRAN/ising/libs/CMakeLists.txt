# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

set(libs_LIBS_SRCS
        globales.f90
        ziggurat.f90
        estadistica.f90
	io_parametros.f90
	isingmods.f90
	strings.f90
	usozig.f90 
    )


add_library(libs ${libs_LIBS_SRCS})
