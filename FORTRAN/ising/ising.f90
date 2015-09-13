program main_ising 
    
  use globales,            only: dp
  use isingmods,           only: read_parameters, inicializacion, calcula_EM,   &
                                 metropolis, finalizacion 
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

 ! call running_mean_var(x,y,'energia.dat')
 !print *, x,sqrt(y), w, sqrt(z)
  call vector_mean_var(x,y,'energia.dat')
  call vector_mean_var(w,z,'magneti.dat')
  print *, x,sqrt(y), w, sqrt(z)

end program main_ising
