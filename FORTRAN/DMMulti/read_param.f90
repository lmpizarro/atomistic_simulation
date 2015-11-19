!
!
!
!
!
program read_param
  implicit none

  integer :: N_total_part, Nespecies, i
  real :: v1, v2, v3

  ! guarda los parámetros del modelo lj de
  ! cada especie
  REAL, ALLOCATABLE :: lj_param(:,:) 
  ! guarda el porcentaje de cada especie
  REAL, ALLOCATABLE :: percent(:) 
  ! guarda la cantidad de  partículas de cada especie
  Integer, ALLOCATABLE :: Np(:)
  real :: Lado_caja, densidad, Temperatura
  real :: v_temp

  OPEN(UNIT=10,FILE='input.par',STATUS='UNKNOWN')

  ! lee la cantidad de especies
  read(10,*) Nespecies

  ! Alloca memoria de acuerdo a la cantidad de especies
  ALLOCATE (lj_param(1:Nespecies,1:3))
  allocate(percent(1:Nespecies))
  allocate(Np(1:Nespecies))

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

  N_total_part = int (densidad * Lado_caja ** 3)

  ! calcula la cantidad de partículas de cada especie
  do i=1, Nespecies
    Np(i) = N_total_part * percent(i) 
  enddo
  ! Ajuste de la cantidad total de partículas
  ! La conversión a enteros pierde partículas
  N_total_part = sum(Np)

  call print_info()

contains

  subroutine print_info()
    print *, "Nespecies", Nespecies
    print *, "epsilon           sigma           masa"
    do i=1, Nespecies
       print *, lj_param(i,1), lj_param(i,2), lj_param(i,3)
    enddo
    print *, "Lado_caja: ", Lado_caja, "densidad", densidad, "temperatura: ",&
        Temperatura

    print *, "Cantidad Total de Partículas: ", N_total_part
    do i=1, Nespecies
       print *, "Especie: ", i, "porcentaje: " , percent(i), "cant de particulas: ", Np(i) 
    enddo
  endsubroutine
endprogram read_param
