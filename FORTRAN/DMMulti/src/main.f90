!
!
!
!
!
program main
  use types
  use globales

  use inic_fin

  use mediciones

  implicit none

  call inicializacion()

  call integracion()

  print *, "llama fin"
  call finalizacion()
endprogram main
