program main_dimod 

#include "control.h"

  use inic_fin,         only: inicializacion, finalizacion 
  use utils,            only: wtime
  use types,            only: dp
  use integra,          only: integracion

  implicit none

  real(dp)   :: wt ! Variable auxiliar para calcular tiempos de ejecución

print *, 'Estoy en src'
#if THERM == 1
print *, 'Se definió THERM == 1'
#endif
#if THERM == 0
print *, 'Se definió THERM == 0'
#endif

  ! Inicializa parámetros 
  call inicializacion()
    
  wt = wtime()

  ! Integración de la dinámica
  call integracion()

  print *, 'Tiempo transcurrido: ', wtime()-wt

  ! Finalización del programa
  call finalizacion()

end program main_dimod
