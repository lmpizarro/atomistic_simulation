module integra

#include "control.h"

  use types,          only: dp
  use globales,       only: gt, gL,gDt, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil, gM, & 
                            gNmed, gPot, gKin, gVir 
  use mediciones,     only: calcula_kin, calcula_pres, calcula_fuerza, calcula_temp 
  use io_parametros,  only: escribe_en_columnas 
implicit none

private

public      :: integracion, cpc_vec

contains

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO - VELOCITY VERLET
  !===============================================================================
  ! Integra las ecuaciones dinámicas con el algoritmo de Velocity-Verlet

  subroutine integracion()

    real(dp)                                :: Pres    ! Presión instantánea
    real(dp)                                :: Temp    ! Presión instantánea
    real(dp), dimension(:,:), allocatable   :: Eng_t   ! Energía en función del tiempo
    real(dp), dimension(:), allocatable     :: Pres_t  ! Presión en función del tiempo
    real(dp), dimension(:), allocatable     :: Temp_t  ! Temperatura en función del tiempo
#ifdef CONTROL_TEMP
    real(dp), dimension(:,:), allocatable   :: Vel_t   ! Temperatura en función del tiempo
#endif
    integer    :: i, j
    integer    :: Kmed                                 ! Cantidad de puntos medidos

    ! Se define la cantidad de puntos que se van a medir
    Kmed = int(gNtime/abs(gNmed)) + 1              ! Se agrega +1 para poner el inicial
    allocate( Eng_t(1:3,1:Kmed), Pres_t(1:Kmed), Temp_t(1:Kmed) )
#ifdef CONTROL_TEMP
    allocate( Vel_t(1:3,1:Kmed) )
#endif

    ! El primer punto son los valores iniciales
    call calcula_pres(Pres)
    Pres_t(1) = Pres
    Eng_t(:,1) = (/ gPot, gKin, gPot + gKin/) 
    ! Imprime en pantalla info
    write(*,*) '********************************************'
    print *, '* Comienza integracion temporal (Vel-Verlet)'
    print *, '* Energias por partícula al comienzo de la integración'
    call print_info(Pres)

! -------------------------------------------------------------------------------------------
! COMIENZA EL LOOP PRINCIPAL DE INTEGRACION
! ------------------------------------------------------------------------------------------
    j = 2                                              ! Contador para el loop de mediciones
                                                       ! (j=1 lo usé para el valor inicial)
    do i = 1, gNtime 
      gR = gR + gDt*gV + 0.5_dp * gF * gDt**2 / gM     ! gR(t+dt)
      ! Aplica condiciones peródicas de contorno
      call cpc_vec()
      gV =          gV + 0.5_dp * gF * gDt / gM        ! gV(t+0.5dt) 
      call calcula_fuerza()                            ! Calcula fuerzas y potencial
      gV =          gV + 0.5_dp * gF * gDt / gM        ! gV(t+dt)
      ! Se realizan las mediciones
      if (mod(i,gNmed) == 0) then
        ! Energia cinetica
        call calcula_kin()
        ! Temperatura
        call calcula_temp(Temp)
        ! Presión
        call calcula_pres(Pres)
        ! Escribe energía total
        Eng_t(:,j) = (/gPot, gKin, gPot + gKin/) 
        Pres_t(j)  = Pres
        Temp_t(j)  = Temp
#ifdef CONTROL_TEMP
        ! Guardo las velocidades de una partícula arbitraria
        Vel_t(:,j) = gV(:,15)
#endif
        ! Escribe posiciones de las partículas
        !call escribe_trayectoria(gR,i)
        j = j + 1                                      ! Actualiza contador mediciones
      end if
    end do
! -------------------------------------------------------------------------------------------
! FIN DEL LOOP PRINCIPAL DE INTEGRACION
! -------------------------------------------------------------------------------------------

    ! Escritura de las magnitudes en función del tiempo
    if (gNmed > 0) then
      ! Guarda las energías en un archivo [Potencial, Cinética, Total]
      call escribe_en_columnas(Eng_t,'energias.dat',gNmed*gDt)

      ! Guarda la presion en un archivo
      call escribe_en_columnas(Pres_t,'presion.dat',gNmed*gDt)
      
      ! Guarda la temperatura en un archivo
      call escribe_en_columnas(Temp_t,'temperatura.dat',gNmed*gDt)

#ifdef CONTROL_TEMP  
      ! Guarda la velocidad en un archivo [v_x  v_y  v_z]
      call escribe_en_columnas(Vel_t,'velocidades_control_T.dat',gNmed*gDt)
      ! Libera memoria 
      deallocate( Vel_t )
#endif
    end if
    
    ! Se imprime en pantalla los resultados finales
    print *, '* Energias por partícula al final de la integración'
    call print_info(Pres)
    print *, '* Fin de la integracion temporal'
    ! Se libera memoria
    deallocate( Eng_t, Pres_t, Temp_t )

contains

    subroutine print_info(presion)
    ! Subrutina para imprimir en pantalla resultados de interes
 
      real(dp), intent(in)        :: presion     ! Presión instantánea

      print *, 'Pot=' , gPot/gNpart, 'Kin=', gKin/gNpart, 'Tot=', (gPot+gKin)/gNpart
      print *, 'Presion= ' , presion
     
    end subroutine print_info

  end subroutine integracion

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================

  subroutine cpc_vec()

    gR = gR - gL*floor(gR/gL)

  end subroutine

end module integra
