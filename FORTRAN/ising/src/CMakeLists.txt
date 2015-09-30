# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

include_directories(${CMAKE_BINARY_DIR}/libs)

add_executable(ising ising.f90)
add_executable(crea_histograma crea_histograma.f90)

target_link_libraries(ising libs)
target_link_libraries(crea_histograma libs)
