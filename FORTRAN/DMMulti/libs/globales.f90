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
  real(dp) :: gLado_caja, gRho, gTemperatura, gDt
  ! gNpart: cantidad de partículas del sistema
  ! gNtime: cantidad de pasos de  tiempo
  ! gNmed : cantidad de pasos entre mediciones
  integer :: gNpart, gNtime, gNmed

  ! guarda los parámetros del modelo lj de
  ! cada especie
  REAL(dp), ALLOCATABLE :: gLj_param(:,:) 
  ! guarda el porcentaje de cada especie
  REAL(dp), ALLOCATABLE :: gPercent(:) 
  ! guarda la cantidad de  partículas de cada especie
  Integer, ALLOCATABLE :: gNp(:)

  integer :: gNespecies

  ! gR:  posicion de la particula
  ! gF:  fuerza entre particulas particula
  ! gV:  velocidad de la particula
  real(dp),  dimension(:,:), allocatable  :: gR, gF, gV

  !
  !
  real(dp),  dimension(:,:), allocatable  :: gCombSigma
  real(dp),  dimension(:,:), allocatable  :: gCombEpsilon

  real(dp), allocatable :: gRc2(:) 
  real(dp), allocatable :: gPot_Cut(:) 
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
            gNtime, " paso de medida", gNmed



    print *, "N Partículas Total: ", gNPart
    print *, "Especie: porcentaje: N particulas: " 

    do  i = 1, gNespecies
      write(*,700)  i, gpercent(i), gNp(i) 
    enddo  

    300 format (a, F10.3, a, I10, a, I10)
    400 format (a, F10.3, a, F10.3, a, F10.3)
    500 format (a)
    600 format (a, I3)
    700 format (I3, F15.3, I10)
    800 format (F12.3, F12.3, F12.3)
  endsubroutine print_gvars


  subroutine inicializar_globales()
    integer :: i
    real(dp) :: v_temp
    ! determina el porcentaje de la especie restante
    v_temp = 0.0
    do i=1, gNespecies - 1
      v_temp = v_temp + gpercent(i) 
    enddo
    gpercent(gNespecies) = 1 - v_temp

    gNPart = int (gRho * gLado_caja ** 3)

    allocate(gNp(1:gNespecies))
    ! calcula la cantidad de partículas de cada especie
    do i=1, gNespecies
      gNp(i) = gNPart * gpercent(i) 
    enddo
    ! Ajuste de la cantidad total de partículas
    ! La conversión a enteros pierde partículas
    gNPart = sum(gNp)


    allocate(gR(1:3, 1:gNpart))
    allocate(gF(1:3, 1:gNpart))
    allocate(gV(1:3, 1:gNpart))

    call print_gvars()

  endsubroutine inicializar_globales

  subroutine finalizar_globales()
    deallocate(gR,gV,gF, gNp, gLj_param, gPercent)
    deallocate(gCombSigma, gCombEpsilon)
    deallocate(gRc2, gPot_Cut)
  endsubroutine finalizar_globales

endmodule globales
