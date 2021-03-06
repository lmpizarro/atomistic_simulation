!===============================================================================
! DEFINICION DE VARIABLES GLOBALES Y SUBRUTINAS ASOCIADAS 
!===============================================================================
!
module globales

#include "control.h"

  use types,       only: dp

  implicit none

  !=============================================================================
  ! PARAMETROS FIJOS DEL SISTEMA 
  !=============================================================================
  ! Se definen a través del archivo de entrada 'input.par'. Luego no cambian 

  real(dp) :: gTemperatura   ! Temperatura del sistema
  real(dp) :: gDt            ! Paso de tiempo de la corrida
  real(dp) :: gLado_caja     ! Longitud de un lado del cubo
  integer  :: gNpart         ! Cantidad de partículas del sistema
  integer  :: gNtime         ! Cantidad de pasos de tiempo
  integer  :: gNmed          ! Cantidad de pasos entre mediciones
  real(dp) :: gGamma         ! Parámetro para el termostato de Langevin  
  
  ! --- Parámetros asociados la cantidad de especies
  real(dp), allocatable :: gLj_param(:,:)  ! Parámetros del modelo lj de cada especie
  real(dp), allocatable :: gPercent(:)     ! Porcentaje de cada especie
  integer, allocatable  :: gNp(:)          ! Cantidad de partículas de cada especie
  integer               :: gNespecies      ! Cantidad de especies en el sistema
  integer               :: gLiqSol         ! =1 si es un sólido. =0 otra cosa 

  ! Define la cantidad de veces que se repite el cubo basico
  ! aplica a z, y, x
  integer :: gPeriodos
  ! Define el tipo de estrucura cubica 
  ! 0: simple, 
  ! 1: centrada en el cuerpo, 
  ! 2: centrada en las caras
  integer :: gCubicStructure


  ! --- Parametros asociados al radio de corte
  real(dp), allocatable :: gRc2(:,:)  ! Cuadrado del radio de corte
  real(dp), allocatable :: gPot_Cut(:,:)  ! Potencial L-J en el radio de corte
  ! --- Parámetros derivados
  real(dp) :: gVol       ! Volumen del cubo
  real(dp) :: gRho       ! Densidad numero de particulas N/V
 
  !=============================================================================
  ! VARIABLES DEL SISTEMA 
  !=============================================================================

  real(dp),  dimension(:,:), allocatable  :: gR    ! Posicion de las partículas
  real(dp),  dimension(:,:), allocatable  :: gV    ! Velocidades de las partículas
  real(dp),  dimension(:,:), allocatable  :: gF    ! Fuerza sobre las partículas
  real(dp)                                :: gPot  ! Energía potencial del sistema
  real(dp)                                :: gKin  ! Energia cinetica del sistema
  real(dp)                                :: gVir  ! Cálculo termino del virial
  integer, allocatable :: gIndice_elemento(:)      ! Indice de partícula

  !============================================================================
  ! Para la medicion de correlaciones de velocidad
  !============================================================================
  ! Se hace sobre posiciones equivalentes por simetría por especie
  !
  real(dp),  dimension(:,:), allocatable  :: gCorrVfac_1 
  real(dp),  dimension(:,:), allocatable  :: gCorrVfac_2 
  real(dp),  dimension(:,:), allocatable  :: gCorrVfac_3 
  real(dp),  dimension(:,:), allocatable  :: gCorrVver_1 
    
  integer :: gNmodosVibra ! cantidad de velocidades acumuladas para modos de
                          ! vibracion 
  integer :: gNCorrVfac_1 ! cantidad de particula en las caras 1
  integer :: gNCorrVfac_2 ! idem cara 2 
  integer :: gNCorrVfac_3 ! idem cara 3
  integer :: gNCorrVver_1 ! cantidad  de particula en los vertices 

  integer ::   gKmed   ! Se define la cantidad de puntos que se van a medir
  ! Magnitudes utilizadas para describir muchas especies
  real(dp),  dimension(:,:), allocatable  :: gCombSigma    ! Matriz con sigma
  real(dp),  dimension(:,:), allocatable  :: gCombEpsilon  ! Matriz con epsilon
  real(dp),  dimension(:), allocatable    :: gMasa         ! Vector de masas

#ifdef CORR_PAR
  ! g(r) sin normalizar (por eso definida como entera)  
  ! TODO  Poner con el prepro, hace falta ver la subrutina distribucion_radial
  integer, dimension(:,:), allocatable   :: gCorr_par    
  integer           :: gNhist   ! Cantidad de bines de la g(r)
  integer           :: gNgr     ! Contador para saber cuántas veces se acumuló la g(r)
  real(dp)          :: gDbin    ! Ancho del bin de la g(r)    
#endif

contains
  
  !=============================================================================
  ! TODO AGREGAR DESCRIPCION DE LO QUE HACE 
  ! ¿por qué se pone acá esto?
  !=============================================================================

  subroutine print_gvars()

    integer :: i      
    
    write( *,600) "Nespecies", gNespecies

    write( *, 500)  "epsilon           sigma           masa"
    do i=1, gNespecies
       write( *, 800) glj_param(i,1), glj_param(i,2), glj_param(i,3)
    enddo

    write (*,400) "Lado_caja: ", gLado_caja, " densidad", gRho, &
            " Temperatura: ", gTemperatura

    write (*,300) "paso de tiempo: ", gDt, " cantidad de pasos: ",&
            gNtime, " paso de medida ", gNmed, "gamma langevin: ", gGamma

    print *, "N Partículas Total: ", gNPart
    print *, "Especie: porcentaje: N particulas: " 

    do  i = 1, gNespecies
      write(*,700)  i, gpercent(i), gNp(i) 
    enddo 


    300 format (a, F10.3, a, I10, a, I10, a, F10.3)
    400 format (a, F10.3, a, F10.3, a, F10.3)
    500 format (a)
    600 format (a, I3)
    700 format (I3, F15.3, I10)
    800 format (F12.3, F12.3, F12.3)

  endsubroutine print_gvars

  !=============================================================================
  ! INICIALIZA VARIABLES GLOBALES SI SE TRABAJA CON POSICIONES ALEATORIAS 
  !=============================================================================
  
  subroutine inicializar_globales_random()
    
    integer  :: i
    real(dp) :: v_temp

    ! determina el volumen de la caja
    gVol = gLado_caja ** 3

    ! determina el porcentaje de la especie restante
    v_temp = 0.0
    do i=1, gNespecies - 1
      v_temp = v_temp + gpercent(i) 
    enddo

    gpercent(gNespecies) = 1 - v_temp

    gNPart = int ( gLado_caja ** 3 * gRho)

    allocate(gNp(1:gNespecies))
    ! calcula la cantidad de partículas de cada especie
    do i=1, gNespecies
      gNp(i) = gNPart * gpercent(i) 
    enddo

    ! Ajuste de la cantidad total de partículas
    ! La conversión a enteros pierde partículas
    gNPart = sum(gNp)
    allocate(gIndice_elemento(1:gNPart))

    allocate(gR(1:3, 1:gNpart))
    allocate(gF(1:3, 1:gNpart))
    allocate(gV(1:3, 1:gNpart))

    call print_gvars()

  endsubroutine inicializar_globales_random

  !=============================================================================
  ! INICIALIZA VARIABLES GLOBALES SI SE TRABAJA CON UNA RED FCC 
  !=============================================================================

  subroutine inicializar_globales_fcc()

    integer  :: i
    real(dp) :: v_temp

    ! define la cantidad de veces que se repite el cubo basico
    ! aplica a z, y, x
    !!  gPeriodos
    ! define el tipo de estrucura cubica 
    ! 0: simple, 
    ! 1: centrada en el cuerpo, 
    ! 2: centrada en las caras
    !! gCubicStructure
   
    ! para calcular la cantidad de particulas
    ! se considera un cubo de lado iguales
    ! 4 es la cantidad de atomos por primitiva fcc
    gNpart = 4 * gPeriodos ** 3

    allocate(gIndice_elemento(1:gNPart))
    allocate(gR(1:3, 1:gNpart))
    allocate(gF(1:3, 1:gNpart))
    allocate(gV(1:3, 1:gNpart))

    ! calcula la densidad
    gRho = gNpart / (gLado_caja ** 3)

    ! determina el volumen de la caja
    gVol = gLado_caja ** 3

    ! determina el porcentaje de la especie restante
    v_temp = 0.0
    do i=1, gNespecies - 1
      v_temp = v_temp + gpercent(i) 
    enddo

    gpercent(gNespecies) = 1 - v_temp

    allocate(gNp(1:gNespecies))
    ! calcula la cantidad de partículas de cada especie
    do i=1, gNespecies
      gNp(i) = gNPart * gpercent(i) 
    enddo

    !print *, "fcc ", gNpart, gLado_caja, gRho, gVol

  endsubroutine inicializar_globales_fcc

  !=============================================================================
  ! LIBERA MEMORIA DE LAS VARIABLES GLOBALES 
  !=============================================================================

  subroutine finalizar_globales()
   
    deallocate(gR,gV,gF, gNp, gLj_param, gPercent)
    deallocate(gCombSigma, gCombEpsilon)
    deallocate(gRc2, gPot_Cut)

  endsubroutine finalizar_globales

endmodule globales
