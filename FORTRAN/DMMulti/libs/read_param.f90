!===============================================================================
! MODULO PARA LECTURA DE PARAMETROS DEL PROBLEMA 
!===============================================================================
!
module read_param

  use types,           only: dp
  use globales

  implicit none

contains
  !===============================================================================
  ! LEE LOS DATOS DEL PROBLEMA EN EL ARCHIVO DE ENTRADA 
  !===============================================================================
  ! Se lee el archivo 'input.par' que continene todos los datos del problema
  subroutine leer_parametros()

    integer :: i

    open(unit=10,file='input.par',status='unknown')

    ! lee la cantidad de especies
    read(10,*) gNespecies, gLiqSol

    ! Alloca memoria de acuerdo a la cantidad de especies
    ! TODO tal vez esto vaya en la inicialización
    allocate(gLj_param(1:gNespecies,1:3))
    allocate(gMasa(1:gNespecies))
    ! lee los parámetros de lj de cada especie
    ! epsilon, sigma, masa
    do i=1,gNespecies
      read(10,*) glj_param(i,1), glj_param(i,2), glj_param(i,3)
      gMasa(i) = glj_param(i,3)
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

    !TODO No entiendo para qué esto
    ! ¿Por qué no se define directamente con el N de cada una?
    allocate(gPercent(1:gNespecies))

    ! lee los porcentajes de Nespecies - 1
    do i=1,gNespecies - 1
      read(10,*) gPercent(i)
    enddo

    close(10)

  endsubroutine leer_parametros

endmodule read_param
