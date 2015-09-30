# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

include_directories(${CMAKE_BINARY_DIR}/lib)

add_executable(ising ising.f90)

target_link_libraries(ising lib)
