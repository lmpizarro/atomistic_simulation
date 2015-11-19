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


contains

  subroutine print_gvars()
    integer :: i      
    print *, "Nespecies", gNespecies
    print *, "epsilon           sigma           masa"
    do i=1, gNespecies
       print *, glj_param(i,1), glj_param(i,2), glj_param(i,3)
    enddo
    print *, "Lado_caja: ", gLado_caja, "densidad", gRho, "temperatura: ",&
        gTemperatura
    print *, "paso de tiempo: ", gDt, "cantidad de pasos: ", gNtime, "paso de&
    medida", gNmed

    print *, "Cantidad Total de Partículas: ", gNPart
    do i=1, gNespecies
       print *, "Especie: ", i, "porcentaje: " , gpercent(i), "cant de particulas: ", gNp(i) 
    enddo
  endsubroutine

endmodule globales
