module integra

#include "control.h"

  use types,             only: dp
  use globales
  use mediciones,        only: calcula_pres, calcula_temp, calcula_kin, calcula_fuerza
  use utils,             only: hace_estadistica
#ifdef MODOS_VIB
  use mediciones,        only: acumula_velocidades_equivalentes,  acumula_velocidades_posicion,&
#endif

  use io_parametros,     only: escribe_trayectoria, escribe_en_columnas
#ifdef CORR_PAR
  use mediciones,        only: normaliza_gr
#endif

  implicit none

  private

  public       :: integracion, integracion_min, cpc_vec

contains 

  !===============================================================================
  ! integración de las ecuaciones de movimiento - velocity verlet
  !===============================================================================
  ! integra las ecuaciones dinámicas con el algoritmo de velocity-verlet
  subroutine integracion()

    real(dp)                              :: pres    ! presión instantánea
    real(dp)                              :: temp    ! presión instantánea
    integer                               :: Kmed    ! Cantidad de puntos medidos
    real(dp), dimension(:,:), allocatable :: Eng_t   ! energía en función del tiempo
    real(dp), dimension(:), allocatable   :: Pres_t  ! presión en función del tiempo
    real(dp), dimension(:), allocatable   :: Temp_t  ! temperatura en función del tiempo
#ifdef CORR_PAR
    real(dp), dimension(5,gNhist)         :: cor_par ! Función de correlación de pares    
#endif

    character(50), parameter :: nombre='trayectoria.vtf' ! Nombre archivo para guardar trayectorias
    integer :: i, j, k

    ! Se define la cantidad de puntos que se van a medir
    Kmed = int(gNtime/abs(gNmed)) + 1              ! Se agrega +1 para poner el inicial
    gKmed = Kmed
    allocate( Eng_t(1:3,1:Kmed), Pres_t(1:Kmed), Temp_t(1:Kmed))

#ifdef MODOS_VIB
    ! para acumular velcidades en función del tiempo para el calculo 
    ! de correlaciones
    allocate (gCorrVfac_1(1:3, 1: gNtime))
    allocate (gCorrVfac_2(1:3, 1: gNtime))
    allocate (gCorrVver_1(1:3, 1: gNtime))
    allocate (gCorrVfac_3(1:3, 1: gNtime))
    ! inicializacion para modos vibracion pos indep
    gNmodosVibra = 1 
    ! inicializacion para modos vibracion pos equiv
    gNCorrVfac_1 = 1
    gNCorrVfac_2 = 1
    gNCorrVfac_3 = 1
    gNCorrVver_1 = 1 
#endif

#ifdef GRABA_TRAYECTORIA
    call escribe_trayectoria(gR,nombre,.TRUE.)
#endif

    ! El primer punto son los valores iniciales
    Eng_t(:,1) = (/ gPot, gKin, gPot + gKin/) 
    call calcula_pres(pres)
    Pres_t(1) = pres
    call calcula_temp(temp)
    Temp_t(1) = temp

    ! -------------------------------------------------------------------------------------------
    ! COMIENZA EL LOOP PRINCIPAL DE INTEGRACION
    ! ------------------------------------------------------------------------------------------
    write(*,*) '********************************************'
    write(*,*) '* Comienza la integracion temporal'
    write(*,*) '********************************************'

    j = 2                                              ! Contador para el loop de mediciones
    do i = 1, gNtime 
      do k = 1, gNpart
        gR(:,k) = gR(:,k) + gDt * gV(:,k) + 0.5_dp *  gF(:,k) * gDt ** 2 / gMasa(gIndice_elemento(k))   ! gR(t+dt)
      end do 
     ! Aplica condiciones periódicas de contorno
      call cpc_vec()
      do k = 1, gNpart
        gV(:,k) = gV(:,k) + 0.5_dp * gF(:,k) * gDt / gMasa(gIndice_elemento(k))     ! gV(t+0.5dt) 
      end do
      call calcula_fuerza()                            ! Calcula fuerzas y potencial
      do k = 1, gNpart
        gV(:,k) = gV(:,k) + 0.5_dp * gF(:,k) * gDt / gMasa(gIndice_elemento(k))     ! gV(t+dt)
      end do
      !
      ! para calcular modos de vibracion sobre posiciones
      ! independientes
      !
#ifdef MODOS_VIB
      call acumula_velocidades_posicion_simple (gV, 32)
#endif 


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
        Pres_t(j)  = pres
        Temp_t(j)  = temp

#ifdef GRABA_TRAYECTORIA
    call escribe_trayectoria(gR,nombre,.FALSE.)
#endif
        j = j + 1                                      ! Actualiza contador mediciones
      end if
    end do
  ! ----------------------------------------------------------------------------
  ! FIN DEL LOOP PRINCIPAL DE INTEGRACION
  ! ----------------------------------------------------------------------------
  write(*,*) '********************************************'
  write(*,*) '* Fin de la integracion temporal'
  write(*,*) '********************************************'

#ifdef MODOS_VIB
    call escribe_en_columnas(gCorrVfac_1,'veloc_fac_1.dat',gNmed*gDt)
    call escribe_en_columnas(gCorrVfac_2,'veloc_fac_2.dat',gNmed*gDt)
    call escribe_en_columnas(gCorrVfac_3,'veloc_fac_3.dat',gNmed*gDt)
    call escribe_en_columnas(gCorrVver_1,'veloc_ver_1.dat',gNmed*gDt)
#endif

  ! Escritura de las magnitudes en funci..n del tiempo
  if (gNmed > 0) then
    ! Guarda las energías por partícula  en un archivo [Potencial, Cinética, Total]
    call escribe_en_columnas(Eng_t/gNpart,'energias.dat',gNmed*gDt)

    ! Guarda la presion en un archivo
    call escribe_en_columnas(Pres_t,'presion.dat',gNmed*gDt)
    
    ! Guarda la temperatura en un archivo
    call escribe_en_columnas(Temp_t,'temperatura.dat',gNmed*gDt)
  end if

#ifdef CORR_PAR
    ! Calcula y guarda la función de correlación de pares g(r)
    call normaliza_gr(cor_par)
    ! Como se escribe también r, el dr no es necesario. Se pasa 0.0 para evitar errores
    call escribe_en_columnas(cor_par,'gr.dat',0.0_dp)
#endif
 
  ! Se imprime en pantalla los resultados finales
  write(*,*) '* Energias por partícula al final de la integración'
  call print_info(Pres,Temp)
  
  ! Se calculan valores medios de presión y temperatura y se escriben a archivo
  call hace_estadistica(Pres_t, Temp_t, Eng_t/gNpart)
  
  ! Se libera memoria
  deallocate( Eng_t, Pres_t, Temp_t)

  contains

    subroutine print_info(presion,temperatura)
    ! Subrutina para imprimir en pantalla resultados de interes
 
      real(dp), intent(in)        :: presion     ! Presi..n instant..nea
      real(dp), intent(in)        :: temperatura ! Temperatura instant..nea

      write(*,100) gPot/gNpart, gKin/gNpart, (gPot+gKin)/gNpart
      100 format(1X,'Potencial = ', E14.7, 5X, 'Cinética = ', E14.7, 5X, 'Total = ', E14.7) 
      write(*,*)
      write(*,200)  presion, temperatura
      200 format(1X,'Presión = ' ,E14.7,5X,'Temperatura = ',E14.7)
      write(*,*) '********************************************'

    end subroutine print_info

  end subroutine integracion

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO - MINIMIZACIÓN ENERGÍA
  !===============================================================================
  ! Subrutina de integración de las ecuaciones de movimiento para minimizar energía
  subroutine integracion_min()

    real(dp), dimension(gNtime+1) :: Eng_t   ! Energía en función del tiempo
    integer :: i, j
    ! Nombre del archivo para guardar la trayectoria
    character(50), parameter   :: nombre='trayectoria_pot_minimiza.vtf'
     
    ! Escribe encabezado y primer punto de la trayectoria
    call escribe_trayectoria(gR,nombre,.TRUE.) 
    ! El primer punto es la energía inicial
    Eng_t(1) = gPot

    do i = 1, gNtime 
      do j = 1, gNpart
        gR(:,j) = gR(:,j) + 0.5_dp * gF(:,j) * (gDt**2) / gMasa(gIndice_elemento(j))
      end do
      ! Aplica condiciones periódicas de contorno
      call cpc_vec()   

      ! Esta subrutine abre y cierra un archivo. Se puede optimizar haciéndolo acá.
      call escribe_trayectoria(gR,nombre,.FALSE.)    
      ! Calcula fuerza y energía
      call calcula_fuerza()
      ! Escribe energía potencial en vector
      Eng_t(i+1) = gPot
    end do

    call escribe_en_columnas(Eng_t/gNpart,'energia_pot_minimizada.dat',gDt)

  end subroutine integracion_min

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================
  subroutine cpc_vec()

    real(dp), dimension(3,gNpart) :: tmp     ! Variable temporal
    
    !gR = gR - gLado_caja*floor(gR/gLado_caja)
    
    ! Lo escribo de esta forma porque de lo contrario da error de compilación
    ! en el cluster (dice ser un bug de gfortran)
    tmp = gLado_caja*floor(gR/gLado_caja)

    gR = gR - tmp

  end subroutine cpc_vec

end module integra
