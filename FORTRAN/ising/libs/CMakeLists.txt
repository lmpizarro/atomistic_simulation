# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

set(lib_LIB_SRCS
        mc.f90
        parameters_io.f90
        ziggurat.f90
        utils.f90
        globals.f90
        constants.f90
        types.f90
    )

add_library(lib ${lib_LIB_SRCS})
