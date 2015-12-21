# Top-level CMakeLists.txt for example
###
### Process this file with cmake to produce Makefile
###

include_directories(${CMAKE_BINARY_DIR}/libs)

add_executable(dmmulti main.f90)
add_executable(modo_vib modo_vib.f90)
add_executable(autocorr autocorrelaciones.f90)
# Línea modificada para poder linkear la biblioteca fftw3
target_link_libraries(dmmulti libs -L/home/pizarro/local/lib -lfftw3)
target_link_libraries(modo_vib libs -L/home/pizarro/local/lib -lfftw3)
target_link_libraries(autocorr libs -L/home/pizarro/local/lib -lfftw3)
# usar estas lineas si fftw3 no se encuentra en su lugar  
#target_link_libraries(dmmulti libs -L/home/pizarro/local/lib -lfftw3)
#target_link_libraries(modo_vib libs -L/home/pizarro/local/lib -lfftw3)
#target_link_libraries(autocorr libs -L/home/pizarro/local/lib -lfftw3)
#target_link_libraries(dmmulti libs)
set_target_properties(dmmulti PROPERTIES
	COMPILE_FLAGS "${f90flags}"
	LINK_FLAGS "${ldflags}")
set_target_properties(modo_vib PROPERTIES
	COMPILE_FLAGS "${f90flags}"
	LINK_FLAGS "${ldflags}")
set_target_properties(autocorr PROPERTIES
	COMPILE_FLAGS "${f90flags}"
	LINK_FLAGS "${ldflags}")

