program autocorr
  use modos_vib
  use types 
  use FFTW3
  use io_parametros

  implicit none

  logical ::  leido
  character(32) :: nombre
  real(dp), allocatable ::  v(:,:)
  real(dp), allocatable ::  corr(:,:)
  integer :: ns

  nombre = "./veloc_fac_1.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_fac_1.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)

  ns = size(v(1,:))

  allocate(corr(1:3, 1:ns))
  corr = 0
  call calcula_autocorr_vel_3D (v, corr)
  corr(1,:) = corr(1,:) + corr(2,:) + corr(3,:)
  call escribe_en_columnas (corr(1,1:ns), "autocorr_fac_1.dat", 512.0_dp)


  nombre = "./veloc_fac_2.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_fac_2.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)

  corr = 0
  call calcula_autocorr_vel_3D (v, corr)
  corr(1,:) = corr(1,:) + corr(2,:) + corr(3,:)
  call escribe_en_columnas (corr(1,1:ns), "autocorr_fac_2.dat", 512.0_dp)

  nombre = "./veloc_fac_3.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_fac_3.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)

  corr = 0
  call calcula_autocorr_vel_3D (v, corr)
  corr(1,:) = corr(1,:) + corr(2,:) + corr(3,:)
  call escribe_en_columnas (corr(1,1:ns), "autocorr_fac_3.dat", 512.0_dp)


  nombre = "./veloc_ver_1.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)
  nombre = "./modo_ver_1.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)

  corr = 0
  call calcula_autocorr_vel_3D (v, corr)
  corr(1,:) = corr(1,:) + corr(2,:) + corr(3,:)

  call escribe_en_columnas (corr(1,1:ns), "autocorr_ver_1.dat", 512.0_dp)

  deallocate(v)
  deallocate(corr)



end program autocorr 
