program main_ising 
    
  use isingmods,           only: read_parameters, inicializacion, calcula_EM,   &
                                 metropolis, finalizacion, hace_estadistica 

  implicit none

  ! Lee los datos necesario
  call read_parameters()

  ! Inicializo la matriz de spines con la condicion de contorno
  call inicializacion()
 
  ! Calcula la energía y magnetización
  call calcula_EM()

  ! Ejecuta el algoritmo de metrópolis
  call metropolis()
  
  ! Finaliza programa
  call finalizacion()

end program main_ising
