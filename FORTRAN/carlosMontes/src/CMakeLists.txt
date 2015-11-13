# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

include_directories(${CMAKE_BINARY_DIR}/libs)

add_executable(mc_pote mc_pote.f90)
target_link_libraries(mc_pote libs)
set_target_properties(mc_pote PROPERTIES
	COMPILE_FLAGS "${f90flags}"
	LINK_FLAGS "${ldflags}")

