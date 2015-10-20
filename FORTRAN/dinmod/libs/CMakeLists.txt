# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

file (GLOB libs_LIBS_SRCS "*.f90")

add_library(libs ${libs_LIBS_SRCS})

set_target_properties(libs PROPERTIES
	COMPILE_FLAGS "${f90flags}"
	LINK_FLAGS "${ldflags}")

