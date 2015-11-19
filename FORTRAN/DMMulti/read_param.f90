program read_param
  implicit none

  integer :: Nespecies, i
  real :: v1, v2, v3

  REAL, ALLOCATABLE :: lj_param(:,:), percent(:)
  real :: Lado_caja, densidad, Temperatura
  real :: v_temp

  OPEN(UNIT=10,FILE='input.par',STATUS='UNKNOWN')

  ! lee la cantidad de especies
  read(10,*) Nespecies

  ALLOCATE (lj_param(1:Nespecies,1:3))
  allocate(percent(1:Nespecies))

  ! lee los parámetros de lj de cada especie
  do i=1,Nespecies
    read(10,*) lj_param(i,1), lj_param(i,2), lj_param(i,3)
  enddo
 
  ! lee el tamaño de la caja y la densidad y temperatura
  read(10,*) Lado_caja, densidad, Temperatura

  ! lee los porcentajes de Nespecies - 1
  do i=1,Nespecies - 1
    read(10,*) percent(i)
  enddo
  close(10)

  ! determina el porcentaje de la especie restante
  v_temp = 0.0
  do i=1, Nespecies - 1
    v_temp = v_temp + percent(i) 
  enddo
  percent(Nespecies) = 1 - v_temp




  print *, "Nespecies", Nespecies
  print *, "epsilon           sigma           masa"
  do i=1, Nespecies
     print *, lj_param(i,1), lj_param(i,2), lj_param(i,3)
  enddo
  print *, "Lado_caja: ", Lado_caja, "densidad", densidad, "temperatura: ",&
     Temperatura

  do i=1, Nespecies
     print *, "Especie: ", i, "porcentaje: " , percent(i) 
  enddo

endprogram read_param
