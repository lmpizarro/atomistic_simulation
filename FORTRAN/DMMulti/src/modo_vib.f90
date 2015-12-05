program freq_vib
  use modos_vib
  use types 
  use FFTW3
  use io_parametros

  implicit none

  logical ::  leido
  character(32) :: nombre
  real(dp), allocatable ::  v(:,:)
  real(dp), allocatable ::  corr(:,:)

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

  allocate(corr(1:3, 1:2*size(v(1,:)) + 1))

  call calcula_autocorr_vel_3D (v, corr)


  call escribe_en_columnas (corr, "autocorr.dat", 512.0_dp)

  deallocate(v)

endprogram freq_vib
