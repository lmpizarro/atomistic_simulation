! 
!
!
!
!
module read_param
  use types
  use globales

  implicit none

contains
  !
  ! Subroutina que lee los parámetros del problema
  !
  subroutine leer_parametros()
    integer :: i

    OPEN(UNIT=10,FILE='input.par',STATUS='UNKNOWN')

    ! lee la cantidad de especies
    read(10,*) gNespecies

    ! Alloca memoria de acuerdo a la cantidad de especies
    ALLOCATE (gLj_param(1:gNespecies,1:3))

    ! lee los parámetros de lj de cada especie
    do i=1,gNespecies
      read(10,*) glj_param(i,1), glj_param(i,2), glj_param(i,3)
    enddo
 
    ! lee el tamaño de la caja, la densidad, temperatura, 
    ! paso de tiempo, cantidad de pasos de tiempo
    read(10,*) gLado_caja, gRho, gTemperatura, gDt, gNtime, gNmed


    allocate(gPercent(1:gNespecies))

    ! lee los porcentajes de Nespecies - 1
    do i=1,gNespecies - 1
      read(10,*) gpercent(i)
    enddo
    close(10)
          
  endsubroutine leer_parametros

endmodule read_param
