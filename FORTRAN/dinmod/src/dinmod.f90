program main_dimod 

  use globales,         only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM, gNmed
  use dinmods,          only: inicializacion, finalizacion, integracion 
  use utils,            only: wtime
  use types,            only: dp

  implicit none

  real(dp)   :: wt ! Variable auxiliar para calcular tiempos de ejecución

  ! Inicializa parámetros 
  call inicializacion()
    
  wt = wtime()

  ! Integración de la dinámica
  call integracion()

  print *, 'Tiempo transcurrido: ', wtime()-wt

  ! Finalización del programa
  call finalizacion()

end program main_dimod
