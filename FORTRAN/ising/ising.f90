program main_ising 
    
  use globales,            only: dp
  use isingmods,           only: read_parameters, inicializacion, calcula_EM,   &
                                 metropolis, finalizacion, hace_estadistica 
  use estadistica,         only: running_mean_var, vector_mean_var

  implicit none

  real(dp) :: x,y,w,z
  
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

  ! Calcula varlor medio y varianza de los resultados
  call hace_estadistica()

end program main_ising
