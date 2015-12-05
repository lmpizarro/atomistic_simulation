program freq_vib

  use modos_vib
  use types 

  implicit none

  logical ::  leido
  character(32) :: nombre
  real(dp), allocatable ::  v(:,:)
  integer :: i

  nombre = "./veloc_ver_1.dat"
  nombre = trim(nombre)
  call lee_velocidades(nombre, v, leido)


  nombre = "./modo_ver_1.dat"
  nombre = trim(nombre)
  call modos_posicion(v, nombre)


  deallocate(v)

endprogram freq_vib
