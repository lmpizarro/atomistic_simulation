module globales
  
  implicit none

  integer, parameter :: dp=kind(0.d0)
  ! gT: Temperatura de la corrida
  ! gDt: Paso de tiempo de la corrida
  ! gL: longitud de un lado del cubo
  real(dp) :: gT, gDt, gL
  ! gNpart: cantidad de partículas del sistema
  ! gNtime: cantidad de pasos de  tiempo
  integer :: gNpart, gNtime
  !real(dp),  dimension(:), allocatable:: 

end module globales
