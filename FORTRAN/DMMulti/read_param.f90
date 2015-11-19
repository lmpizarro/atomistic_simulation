!
!
!
!
!
program read_param
  use types
  use globales

  implicit none

  integer :: i

  real(dp) :: v_temp

  OPEN(UNIT=10,FILE='input.par',STATUS='UNKNOWN')

  ! lee la cantidad de especies
  read(10,*) gNespecies

  ! Alloca memoria de acuerdo a la cantidad de especies
  ALLOCATE (gLj_param(1:gNespecies,1:3))
  allocate(gPercent(1:gNespecies))
  allocate(gNp(1:gNespecies))

  ! lee los parámetros de lj de cada especie
  do i=1,gNespecies
    read(10,*) glj_param(i,1), glj_param(i,2), glj_param(i,3)
  enddo
 
  ! lee el tamaño de la caja, la densidad, temperatura, 
  ! paso de tiempo, cantidad de pasos de tiempo
  read(10,*) gLado_caja, gRho, gTemperatura, gDt, gNtime, gNmed

  ! lee los porcentajes de Nespecies - 1
  do i=1,gNespecies - 1
    read(10,*) gpercent(i)
  enddo
  close(10)

  ! determina el porcentaje de la especie restante
  v_temp = 0.0
  do i=1, gNespecies - 1
    v_temp = v_temp + gpercent(i) 
  enddo
  gpercent(gNespecies) = 1 - v_temp

  gNPart = int (gRho * gLado_caja ** 3)

  ! calcula la cantidad de partículas de cada especie
  do i=1, gNespecies
    gNp(i) = gNPart * gpercent(i) 
  enddo
  ! Ajuste de la cantidad total de partículas
  ! La conversión a enteros pierde partículas
  gNPart = sum(gNp)

  call print_gvars()


endprogram read_param
