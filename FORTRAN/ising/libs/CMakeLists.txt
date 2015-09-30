# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

set(lib_LIB_SRCS
        ziggurat.f90
        globales.f90
        estadistica.f90
	io_parametros.f90
	isingmods.f90
	strings.f90
	usozig.f90 
    )


add_library(lib ${lib_LIB_SRCS})
