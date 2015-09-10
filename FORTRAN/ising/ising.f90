program main_ising 
    
  use isingmods 
  use globals

  implicit none
    
  ! Lee los datos necesario
  call read_parameters()

  ! Inicializo la matriz de spines con la condicion de contorno
  call inicializacion()
 
  ! Calcula la energía y magnetización
  call calcula_EM(Eng,Mag,RED)
  
  print *, 'Valores iniciales'
  print *, 'Energía = ',  Eng, 'Magnetización = ',  Mag 

  ! Ejecuta el algoritmo de metrópolis
  call metropolis()
  
  ! Finaliza programa
  call finalizacion()

end program main_ising
