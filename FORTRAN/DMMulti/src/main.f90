program main

  use inic_fin,        only: inicializacion, finalizacion 
  use integra,         only: integracion

  implicit none

  call inicializacion()

  call integracion()

  call finalizacion()

endprogram main
