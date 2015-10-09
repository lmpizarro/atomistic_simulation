# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

set(libs_LIBS_SRCS
        globales.f90
        ziggurat.f90
        types.f90
        constants.f90
        dinmods.f90
	io_parametros.f90
    )


add_library(libs ${libs_LIBS_SRCS})
