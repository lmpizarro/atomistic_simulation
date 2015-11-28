module read_param

  use types
  use globales

  implicit none

contains
  !===============================================================================
  ! LEE LOS DATOS DEL PROBLEMA EN EL ARCHIVO DE ENTRADA 
  !===============================================================================
  ! Se lee el archivo 'input.par' que continene todos los datos del problema
  subroutine leer_parametros()

    integer :: i

    open(unit=10,file='input.par',status='UNKNOWN')

    ! lee la cantidad de especies
    read(10,*) gNespecies, gLiqSol

    ! Alloca memoria de acuerdo a la cantidad de especies
    allocate(gLj_param(1:gNespecies,1:3))

    ! lee los parámetros de lj de cada especie
    do i=1,gNespecies
      read(10,*) glj_param(i,1), glj_param(i,2), glj_param(i,3)
    enddo

    if (gLiqSol .eq. 0) then
      ! lee el tamaño de la caja, la densidad, temperatura, 
      ! paso de tiempo, cantidad de pasos de tiempo
      read(10,*) gLado_caja, gRho, gTemperatura, gDt, gNtime, gNmed, gGamma

    else if ( gLiqSol .eq. 1) then
      read(10,*) gLado_caja, gPeriodos, gCubicStructure, gTemperatura, gDt, gNtime, gNmed, gGamma
    else
      print *, "no implementado"
      stop 212121212
    endif

    allocate(gPercent(1:gNespecies))

    ! lee los porcentajes de Nespecies - 1
    do i=1,gNespecies - 1
      read(10,*) gpercent(i)
    enddo

    close(10)

  endsubroutine leer_parametros

endmodule read_param
