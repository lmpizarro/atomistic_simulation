program main

  !use types
  !use globales
  use inic_fin,        only: inicializacion, finalizacion 
  !use mediciones
  use integra,         only: integracion

  implicit none

  call inicializacion()

  call integracion()

  print *, "llama fin"
  call finalizacion()

endprogram main
