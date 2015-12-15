!===============================================================================
! MODULO PARA LECTURA DE PARAMETROS DEL PROBLEMA 
!===============================================================================
!
module read_param

#include "control.h"

  use types,           only: dp
  use globales

  implicit none

contains
  !===============================================================================
  ! LEE LOS DATOS DEL PROBLEMA EN EL ARCHIVO DE ENTRADA 
  !===============================================================================
  ! Se lee el archivo 'input.par' que continene todos los datos del problema
  subroutine leer_parametros()
    logical :: es
    integer :: i

    inquire(file='./input.par',exist=es)

    if(es .eqv. .false.) then
       stop 1
    endif
    open(unit=10,file='./input.par',status='old')

    ! lee la cantidad de especies
    read(10,*) gNespecies, gLiqSol

    write(*,'(a)') ''
    write(*,'(a)')      '************ PARAMETROS LEIDOS **************'
    write(*,'(a,I8)') '************ Número de especies    = ' , gNespecies

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
      read(10,*) gLado_caja, gRho
      write(*,'(a)')      '************ Sistema  líquido ***************'
      write(*,'(a,F8.3)') '************ Lado del cubo         = ' , gLado_caja
      write(*,'(a,F8.3)') '************ Densidad partículas   = ' , gRho

    else if ( gLiqSol .eq. 1) then
      read(10,*) gLado_caja, gPeriodos, gCubicStructure
      write(*,'(a)')      '************ Sistema  sólido ****************'
      if (gCubicStructure .eq. 2) then
        write(*,'(a)')    '************ Estructura cristalina =      FCC'
      end if
      write(*,'(a,F8.3)') '************ Lado del cubo         = ' , gLado_caja
      write(*,'(a,I8)')   '************ Periodos cel. unidad  = ' , gPeriodos
      
    else
      print *, "no implementado"
      stop 3
    end if

    read(10,*) gTemperatura, gDt, gNtime, gNmed, gGamma

    write(*,'(a,F8.3)') '************ Temperatura           = ' , gTemperatura
    write(*,'(a,E8.2)') '************ Paso temporal dt      = ' , gDt
    write(*,'(a,I8)')   '************ Número de pasos dt    = ' , gNtime
    write(*,'(a,I8)')   '************ Pasos mediciones      = ' , gNmed
    write(*,'(a)')      '*********************************************'

    !TODO No entiendo para qué esto
    ! ¿Por qué no se define directamente con el N de cada una?
    allocate(gPercent(1:gNespecies))

    ! lee los porcentajes de Nespecies - 1
    do i=1,gNespecies - 1
      read(10,*) gPercent(i)
    enddo

    close(10)

#if THERM == 0
    write(*,'(a)') ''
    write(*,'(a)')      '******** SIMULACION A E CONSTANTE ***********'
    write(*,'(a)')      '*********************************************'
#elif THERM == 1
    write(*,'(a)') ''
    write(*,'(a)')      '******** SIMULACION A T CONSTANTE ***********'
    write(*,'(a)')      '********* TERMOSTATO DE LANGEVIN ************'
    write(*,'(a,F8.3)') '************ Gamma                 = ' , gGamma
    write(*,'(a)')      '*********************************************'
#endif

  endsubroutine leer_parametros

endmodule read_param
