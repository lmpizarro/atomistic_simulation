module inic_fin 

#include "control.h"

  use types,      only: dp
  use globales,   only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil, gM, & 
                        gNmed, gRc2, gPot_cut, gRho, gVol, gPot, gKin, gVir 
  use mediciones, only: calcula_kin, calcula_pres, calcula_fuerza, calcula_temp 
  use utils,      only: write_array3D_lin  ! No se usa. Se pone para que cmake procese dependencias
  use ziggurat
  use usozig
  use io_parametros,  only: escribe_trayectoria, escribe_estados, lee_estados, &
                            read_parameters, escribe_en_columnas
  use integra,        only: cpc_vec

! Si se calcula la función g(r)
#ifdef CORR_PAR
  use globales,      only: gCorr_par, gNhist, gNgr, gDbin 
#endif

! Si se utiliza openmp
#ifdef _OPENMP
  use omp_lib
  use utils,      only: init_openmp
#endif

! Si se utiliza mpi
#ifdef MPI
  use mpi
#endif

  implicit none

  private

  public :: inicializacion,  finalizacion
  
contains

  !===============================================================================
  ! INICIALIZA SISTEMA de CALCULO
  !===============================================================================

  subroutine inicializacion()

    real(dp)  :: Pres   ! Presión instantánea
    real(dp)  :: Temp   ! Temperatura instantánea
    logical   :: leido  ! Indica si se leyo bien el archivo con la configuracion inicial

    ! Inicializa parametros para el calculo en paralelo con OpenMP
    ! Por ahora nada en particular.
#ifdef _OPENMP
    call init_openmp()
#endif

    ! Lee los datos necesario del archivo parametros.dat
    ! Lee todas las variables globales
    call read_parameters()
    
    allocate(gR(3,gNpart))
    allocate(gV(3,gNpart))
    allocate(gF(3,gNpart))

    ! Calcula volumen y densidad
    gVol = gL**3
    gRho = gNpart / gVol

    write(*,'(a)') 
    write(*,'(a)')      '*********** PARAMETROS DERIVADOS ************'
    write(*,'(a,F8.3)') '************ Volumen               = ' , gVol 
    write(*,'(a,F8.4)') '************ Densidad              = ' , gRho
    write(*,'(a)')      '*********************************************'

    ! Se fija si se grabarán los datos
    if (gNmed == 0) then
      ! Pensando en la función módulo, con esto se evita que mod(j,gMed)
      ! se haga 1 en algún momento
      gNmed = gNtime + 1  
    end if

    ! Inicial generador de número aleatorios
    call inic_zig()
    ! Define el radio de corte y el potencial desplazado
    call corta_desplaza_pote()
#ifdef CORR_PAR
    call inic_gr()
#endif   
    ! Trata de leer el archivo con la configuracion inicial
    call lee_estados(gR,gV,leido)
    ! Si no se leyo, define posicion y velocidades 

   if (leido .eqv. .FALSE. ) then
      write(*,*) '* Se generan posiciones de r y v aleatorias'
      ! Define posiciones iniciales
      call inicia_posicion_rn()
      ! Busca el mínimo de energía
      write(*,*) '* Se busca un mínimo de energía potencial'
      call integracion_min()
      ! DKinefine velocidades iniciales
      call vel_inic()
    end if
    ! Calcula energía cinética inicial 
    call calcula_kin()
    ! Calcula la fuerza inicial
    call calcula_fuerza()
    ! Calcula la presión inicial
    call calcula_pres(Pres)
    ! Calcula la temperatura inicial
    call calcula_temp(Temp)

#ifdef CORR_PAR
    call inic_gr()
#endif

    write(*,*) '* Valores iniciales por partícula'
    write(*,100) gPot/gNpart, gKin/gNpart, (gPot+gKin)/gNpart
    100 format(1X,'Potencial = ', E14.7, 5X, 'Cinética = ', E14.7, 5X, 'Total = ', E14.7) 
    write(*,*)
    write(*,200)  Pres, Temp 
    200 format(1X,'Presion = ' ,E14.7,5X,'Temperatura = ',E14.7)

  end subroutine inicializacion 

  !===============================================================================
  ! RADIO DE CORTE 
  !===============================================================================
  ! Define el radio de corte y calcula el desplazamiento del potencial de L-J 

  subroutine corta_desplaza_pote()   

    ! Radio de corte fijo
    gRc2 = (4.0_dp)**2                
    ! Potencial de L-J evaluado en el radio de corte
    gPot_cut = 4.0_dp * ( 1.0_dp / (gRc2**6) - 1.0_dp / (gRc2**3) ) 
    ! Imprime en pantalla
    write(*,'(a)')      ''
    write(*,'(a)')      '*************** RADIO DE CORTE **************'
    write(*,'(a,F9.4)') '************ Rc    = ' , sqrt(gRc2)
    write(*,'(a,F9.5)') '************ V(Rc) = ' , gPot_cut
    write(*,'(a)')      '*********************************************'
    write(*,*)

  end subroutine corta_desplaza_pote

#ifdef CORR_PAR
  !===============================================================================
  ! INICIALIZA PARAMETROS PARA LA g(r) 
  !===============================================================================
  ! Inicialización de variables para galcular la g(r)
  
  subroutine inic_gr()
   
    ! Este if es porque la rutina de inicialización llama a la de fuerza para minimizar energía
    ! Se están haciendo los cálculos de forma repetida, pero es mejor que modificar y agregar
    ! parámetros de entrada al loop de fuerza. Esto sucede sólo si se inicializan las partículas
    ! de forma aleatoria, de lo contrario no se llama a la rutina de minimización de energía. 
    if ( .not. allocated(gCorr_par) ) then
      gNhist = 400                   ! Número de bines
      allocate(gCorr_par(1:gNhist))

      gNgr      = 0                  ! Contador para saber cuántas veces se acumuló la g(r)
      gDbin     = gL / (2 * gNhist)  ! Ancho del bin
      gCorr_par = 0                  ! Inicializo la g(r) sin normalizar
    
      ! Imprime en pantalla
      write(*,'(a)')      ''
      write(*,'(a)')      '************* CALCULO DE LA g(r) ************'
      write(*,'(a,I0)') '************ Nhist    = ' , gNhist
      write(*,'(a)')      '*********************************************'
      write(*,*)
    else
      ! Vuelvo a inicializar. Los datos guardados correspondían a la minimización de energía y no
      ! tienen ningún significado físico.
      gNgr      = 0 
      gCorr_par = 0
    end if

  end subroutine inic_gr
#endif
 
  !===============================================================================
  ! VELOCIDADES INICIALES 
  !===============================================================================
  ! Subrutina para inicializar las velocidades del problema 

  subroutine vel_inic()   

    integer    :: i, j

    ! Asigna a cada componente de la velocidad una distribución gaussiana N(0,1)
    do i = 1, gNpart
      do j = 1, 3
        gV(j,i) = rnor()
      end do
    end do
    ! Define la desviación estandar sqrt(kT/m) en unidades adimensionales
    gV = sqrt( gT ) * gV

  end subroutine vel_inic

  !===============================================================================
  ! INICIALIZA Posicion en Red periódica cúbica simple
  !===============================================================================

  subroutine inicia_posicion_cs()

    integer  :: i, j, k, l, nl
    real(dp) :: rx, ry, rz

    nl = gNpart ** (1.0_dp/3.0_dp) 

    rx = 0.0_dp
    ry = 0.0_dp
    rz = 0.0_dp
    j = 1
    i = 1
    k = 1
    do l = 1, gNpart
      gR(1, l) = (rx + uni() - 0.5_dp) * gL / nl 
      gR(2, l) = (ry + uni() - 0.5_dp) * gL / nl 
      gR(3, l) = (rz + uni() - 0.5_dp) * gL / nl 

      call cpc_vec()
           
      j = j + 1
      k = k + 1
      rx = rx + 1.0_dp
      if ( mod(j, nl) .eq. 0) then
        rx = 0.0_dp
        ry = ry  + 1.0_dp
      end if        
      if ( mod(k, nl ** 2) .eq. 0) then
        ry = 0.0_dp
        rx = 0.0_dp
        rz = rz + 1.0_dp
       end if        
     enddo     
  
  end subroutine inicia_posicion_cs

  !===============================================================================
  ! INICIALIZA Posicion aleatoria dentro de la caga 
  !===============================================================================

  subroutine inicia_posicion_rn()
 
    integer :: i, j

    do i = 1, gNpart
      do j= 1, 3
        gR(j, i) = uni() * gL 
      end do
   end do     

  end subroutine inicia_posicion_rn

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO - MINIMIZACIÓN ENERGÍA
  !===============================================================================
  ! Subrutina de integración de las ecuaciones de movimiento para minimizar energía
  ! Es el Problema 3 de la Guia_2a

  subroutine integracion_min()

    real(dp), dimension(gNtime+1)   :: Eng_t   ! Energía en función del tiempo
    integer    :: i

    ! Escribe encabezado y primer punto de la trayectoria
    ! call escribe_trayectoria(gR,.FALSE.) 
    ! El primer punto es la energía inicial
    Eng_t(1) = gPot

    do i = 1, gNtime 
      gR = gR + 0.5_dp * gF * gDt**2 / gM
      ! Aplica condiciones peródicas de contorno
      call cpc_vec()    
      ! Esta subrutine abre y cierra un archivo. Se puede optimizar haciéndolo acá.
      ! call escribe_trayectoria(gR,.FALSE.)    
      ! Calcula fuerza y energía
      call calcula_fuerza()
      ! Escribe energía potencial en vector
      Eng_t(i+1) = gPot
    end do

    ! Guarda la energía potencial en un archivo
    call escribe_en_columnas(Eng_t/gNpart,'energia_pot_min.dat',gDt)
  
    write(*,*) '* Energía potencial minimizada: ' , gPot

  end subroutine integracion_min

  !===============================================================================
  ! FINALIZA PARAMETROS
  !===============================================================================

  subroutine finalizacion()

    ! Finaliza el generador de número aleatorios
    call fin_zig()

    ! Escribe la última configuración de posición y velocidad en un archivo
    call escribe_estados(gR,gV)

    ! Libera memoria
    deallocate(gR,gV,gF)

#ifdef CORR_PAR
    deallocate(gCorr_par)
#endif

  endsubroutine finalizacion


end module inic_fin 
