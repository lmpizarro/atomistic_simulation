program main

  use inic_fin,        only: inicializacion, finalizacion 
  use integra,         only: integracion
  use types,           only: dp
  use utils,           only: wtime

  implicit none

  real(dp)   :: wt ! Variable auxiliar para calcular tiempos de ejecución

  wt = wtime()

  call inicializacion()

  call integracion()

  call finalizacion()

  write(*,*) '********************************'
  write(*,*) 'Tiempo transcurrido: ', wtime()-wt
  write(*,*) '********************************'

endprogram main
