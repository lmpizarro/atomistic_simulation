! definición de variables globales que usa
! el programa
!
!
!
module globales
  use types, only: dp

  implicit none

  ! gT: Temperatura de la corrida
  ! gDt: Paso de tiempo de la corrida
  ! gL: longitud de un lado del cubo
  !real(dp) :: gT, gDt, gL
  real(dp) :: gLado_caja, gTemperatura, gDt
  ! gNpart: cantidad de partículas del sistema
  ! gNtime: cantidad de pasos de  tiempo
  ! gNmed : cantidad de pasos entre mediciones
  integer :: gNpart, gNtime, gNmed

  ! guarda los parámetros del modelo lj de
  ! cada especie
  ! epsilon_2 sigma_2 masa_2
  REAL(dp), ALLOCATABLE :: gLj_param(:,:) 
  ! guarda el porcentaje de cada especie
  REAL(dp), ALLOCATABLE :: gPercent(:) 
  ! guarda la cantidad de  partículas de cada especie
  Integer, ALLOCATABLE :: gNp(:)
  ! guarda un indice de particula
  Integer, ALLOCATABLE :: gIndice_elemento(:)

  integer :: gNespecies, gLiqSol

  ! define la cantidad de veces que se repite el cubo basico
  ! aplica a z, y, x
  integer :: gPeriodos
  ! define el tipo de estrucura cubica 
  ! 0: simple, 
  ! 1: centrada en el cuerpo, 
  ! 2: centrada en las caras
  integer :: gCubicStructure

  ! gR:  posicion de la particula
  ! gF:  fuerza entre particulas particula
  ! gV:  velocidad de la particula
  real(dp),  dimension(:,:), allocatable  :: gR, gF, gV

  !
  !
  real(dp),  dimension(:,:), allocatable  :: gCombSigma
  real(dp),  dimension(:,:), allocatable  :: gCombEpsilon

  real(dp), allocatable :: gRc2(:,:)  ! Cuadrado del radio de corte
  real(dp), allocatable :: gPot_Cut(:,:)  ! Potencial L-J en el radio de corte
  real(dp) :: gVol       ! Volumen del cubo
  real(dp) :: gRho       ! Densidad numero de particulas N/V
 
  real(dp)        :: gPot       ! Energía potencial del sistema
  real(dp)        :: gKin       ! Energia cinetica del sistema
  real(dp)        :: gVir       ! Cálculo del virial para la presión 

  real(dp)        :: gGamma     ! Parámetro para el termostato de Langevin  
contains

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


  subroutine inicializar_globales_random()
    integer :: i,j, inic, fin
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
    ALLOCATE(gIndice_elemento(1:gNPart))

    !! esto vale para sistema binario
    !! TODO: revisar para sistema multielemento
    gIndice_elemento = 2
    
    inic = 1
    do i=1, gNespecies
      fin = gNp(i) + inic - 1
      do j=inic, fin 
        gIndice_elemento(j) = i
      enddo
      inic = fin + 1
    enddo  

    allocate(gR(1:3, 1:gNpart))
    allocate(gF(1:3, 1:gNpart))
    allocate(gV(1:3, 1:gNpart))

    call print_gvars()

  endsubroutine inicializar_globales_random

  subroutine inicializar_globales_fcc()
    integer :: i
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

    print *, "fcc ", gNpart, gLado_caja, gRho, gVol

  endsubroutine inicializar_globales_fcc

  subroutine finalizar_globales()
    deallocate(gR,gV,gF, gNp, gLj_param, gPercent)
    deallocate(gCombSigma, gCombEpsilon)
    deallocate(gRc2, gPot_Cut)
  endsubroutine finalizar_globales

endmodule globales
