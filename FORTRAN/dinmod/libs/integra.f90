module integra

#include "control.h"

  use types,          only: dp
  use globales,       only: gt, gL,gDt, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil, gM, & 
                            gNmed, gPot, gKin, gVir 
  use mediciones,     only: calcula_kin, calcula_pres, calcula_fuerza, calcula_temp 
  use io_parametros,  only: escribe_en_columnas 
  use procesamiento,  only: hace_estadistica 

! Si graba la posición de las partículas
#ifdef GRABA_TRAYECTORIA
  use io_parametros,  only: escribe_trayectoria
#endif

! Si calcula la correlación entre pares
#ifdef CORR_PAR
  use globales,       only: gNhist, gDbin, gRho, gCorr_par, gNgr 
  use constants,      only: PI
#endif

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
    integer                                 :: i, j
    integer                                 :: Kmed    ! Cantidad de puntos medidos
#ifdef CORR_PAR
    real(dp), dimension(1:2,1:gNhist)       :: cor_par ! Función g(r)
#endif

    ! Se define la cantidad de puntos que se van a medir
    Kmed = int(gNtime/abs(gNmed)) + 1              ! Se agrega +1 para poner el inicial
    allocate( Eng_t(1:3,1:Kmed), Pres_t(1:Kmed), Temp_t(1:Kmed) )
#ifdef CONTROL_TEMP
    allocate( Vel_t(1:3,1:Kmed) )
#endif
#ifdef GRABA_TRAYECTORIA
    call escribe_trayectoria(gR,.TRUE.)
#endif

    ! El primer punto son los valores iniciales
    Eng_t(:,1) = (/ gPot, gKin, gPot + gKin/) 
    call calcula_pres(Pres)
    Pres_t(1) = Pres
    call calcula_temp(Temp)
    Temp_t(1) = Temp
#ifdef CONTROL_TEMP
    ! Guardo la velocidade inicial de una partícula arbitraria
    Vel_t(:,1) = gV(:,15)
#endif   
    ! Imprime en pantalla info
    write(*,*) '********************************************'
    write(*,*) '* Comienza integracion temporal (Vel-Verlet)'
    write(*,*) 
    write(*,*) '* Energias por partícula al comienzo de la integración'
    call print_info(Pres,Temp)
    write(*,*) '********************************************'

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
        ! Guarda magnitudes
        Eng_t(:,j) = (/gPot, gKin, gPot + gKin/) 
        Pres_t(j)  = Pres
        Temp_t(j)  = Temp
#ifdef CONTROL_TEMP
        ! Guardo las velocidades de una partícula arbitraria
        Vel_t(:,j) = gV(:,15)
#endif
#ifdef GRABA_TRAYECTORIA
        call escribe_trayectoria(gR,.FALSE.)
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
      ! Guarda las energías por partícula  en un archivo [Potencial, Cinética, Total]
      call escribe_en_columnas(Eng_t/gNpart,'energias.dat',gNmed*gDt)

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
#ifdef CORR_PAR
      ! Calcula y guarda la función de correlación de pares g(r)
      call calcula_gr(cor_par)
      call escribe_en_columnas(cor_par,'gr.dat',0.0_dp)
#endif
    end if

    ! Se imprime en pantalla los resultados finales
    write(*,*) '* Energias por partícula al final de la integración'
    call print_info(Pres,Temp)

    ! Se calculan valores medios de presión y temperatura y se escriben a archivo
    call hace_estadistica(Pres_t, Temp_t, Eng_t/gNpart)

    write(*,*) '********************************************'
    write(*,*) '* Fin de la integracion temporal'
    write(*,*) '********************************************'
    ! Se libera memoria
    deallocate( Eng_t, Pres_t, Temp_t )

contains

#ifdef CORR_PAR
    subroutine calcula_gr(grnor)
    ! Subrutina para normalizar y calcular la función g(r)     
    ! ver pg. 86 del Frenkel

      real(dp), dimension(:,:),intent(out) :: grnor  ! Función g(r)
      real(dp)                             :: dvol   ! Diferencie de volumen entre bines
      real(dp)                             :: nid    ! Parte de gas ideal en dvol
      integer                              :: i

      do i = 1, gNhist
        grnor(1,i) = gDbin * (i + 0.5)                ! Posición centrada de cada bin
        dvol  = ( (i+1)**3 - i**3 ) * gDbin**3        ! Diferencia de vol entre los bin (i+1) e i
        nid   = (4.0_dp/3.0_dp) * PI * dvol * gRho    ! Parte del gas ideal en dvol 

        grnor(2,i) = real(gCorr_par(i),kind=dp) / (gNgr*gNpart*nid)      ! Función g(r) normalizada
      end do

    end subroutine calcula_gr
#endif

    subroutine print_info(presion,temperatura)
    ! Subrutina para imprimir en pantalla resultados de interes
 
      real(dp), intent(in)        :: presion     ! Presión instantánea
      real(dp), intent(in)        :: temperatura ! Temperatura instantánea

      write(*,100) gPot/gNpart, gKin/gNpart, (gPot+gKin)/gNpart
      100 format(1X,'Potencial = ', E14.7, 5X, 'Cinética = ', E14.7, 5X, 'Total = ', E14.7) 
      write(*,*)
      write(*,200)  presion, temperatura
      200 format(1X,'Presion = ' ,E14.7,5X,'Temperatura = ',E14.7)

    end subroutine print_info

  end subroutine integracion

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================

  subroutine cpc_vec()

    gR = gR - gL*floor(gR/gL)

  end subroutine

end module integra
