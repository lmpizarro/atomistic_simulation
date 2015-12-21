program freq_vib
  use modos_vib
  use types 
  use FFTW3
  use io_parametros

  !  power spectral density function fft
  ! http://www.cygres.com/OcnPageE/Glosry/SpecE.html
  ! http://www.mathworks.com/help/signal/ug/psd-estimate-using-fft.html
  ! http://www.ece.umd.edu/~tretter/commlab/c6713slides/ch4.pdf pag. 29
  implicit none

  logical ::  leido
  character(32) :: nombre
  real(dp), allocatable ::  v(:,:)

  nombre = "./veloc_fac_1.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_fac_1.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)


  nombre = "./veloc_fac_2.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_fac_2.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)

  nombre = "./veloc_fac_3.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_fac_3.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)

  nombre = "./veloc_ver_1.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_ver_1.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)

  deallocate(v)

endprogram freq_vib
