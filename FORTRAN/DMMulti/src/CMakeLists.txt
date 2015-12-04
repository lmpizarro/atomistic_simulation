# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

include_directories(${CMAKE_BINARY_DIR}/libs)

add_executable(dmmulti main.f90)
# Línea modificada para poder linkear la biblioteca fftw3
target_link_libraries(dmmulti libs -lfftw3)
#target_link_libraries(dmmulti libs)
set_target_properties(dmmulti PROPERTIES
	COMPILE_FLAGS "${f90flags}"
	LINK_FLAGS "${ldflags}")
