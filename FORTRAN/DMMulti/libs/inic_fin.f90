module inic_fin 
  use types,      only: dp
  use globales
  use ziggurat
  use usozig
  use read_param
  use combinacion
  use mediciones
  use integra

implicit none
contains

  !===============================================================================
  ! INICIALIZA SISTEMA de CALCULO
  !===============================================================================
  subroutine inicializacion()
    integer :: j 
    real(dp)  :: insPres   ! Presión instantánea
    real(dp)  :: insTemp   ! Temperatura instantánea

    call inic_zig()
    call leer_parametros()

    call inicializar_globales()
    call comb_sigma()
    call comb_epsilon()

    ! como se calcularon los sigmas
    ! calculamos los rc y potencial en rc
    call corta_desplaza_pote()

    call inicia_posicion_rn()

    print *, "pos inic"
    do j=1,gNpart 
      print *, gR(1,j) , gR(2,j), gR(3,j)
    enddo  

    call calcula_fuerza()

    print *, "energia potencial inicial: ", gPot

    call integracion_min()

    print *, "pos final after min"
    do j=1,gNpart 
      print *, gR(1,j) , gR(2,j), gR(3,j)
    enddo  

    print *, "energia potencial after min: ", gPot

    call vel_inic()

    call calcula_kin ()

    write (*,100) "gKin inicial: ", gKin

    call calcula_fuerza()

    ! Calcula la presión inicial
    call calcula_pres(insPres)
    ! Calcula la temperatura inicial
    call calcula_temp(insTemp)
    write (*,100) "insTemp inicial: ", insTemp
    write (*,100) "insPres inicial: ", insPres

    100 format (a, F20.3)
  endsubroutine inicializacion

  !===============================================================================
  ! FINALIZA SISTEMA de CALCULO
  !===============================================================================
  subroutine finalizacion()
    call finalizar_globales()
    call fin_zig()
  endsubroutine finalizacion


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
    gV = sqrt( gTemperatura ) * gV

  end subroutine vel_inic


  !===============================================================================
  ! INICIALIZA Posicion aleatoria dentro de la caga 
  !===============================================================================
  subroutine inicia_posicion_rn()
 
    integer :: i, j

    do i = 1, gNpart
      do j= 1, 3
        gR(j, i) = uni() * gLado_caja
      end do
   end do     

  end subroutine inicia_posicion_rn

end module inic_fin 
