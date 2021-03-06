program main_dimod 
 
    use omp_lib

  use inic_fin,         only: inicializacion, finalizacion 
  use utils,            only: wtime
  use types,            only: dp
  use integra,          only: integracion

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
