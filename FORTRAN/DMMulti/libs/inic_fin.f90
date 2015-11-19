module inic_fin 
  use types,      only: dp
  use globales
  use ziggurat
  use usozig
  use read_param
  use combinacion
  use mediciones


contains

  !===============================================================================
  ! INICIALIZA SISTEMA de CALCULO
  !===============================================================================
  subroutine inicializacion()
    call leer_parametros()
    call inicializar_globales()
    call comb_sigma()
    call comb_epsilon()

    call inic_zig()

    call inicia_posicion_rn()
    call vel_inic()

    call calcula_kin ()


    write (*,100) "gKin: ", gKin

    100 format (a, F15.3)
  endsubroutine inicializacion

  !===============================================================================
  ! FINALIZA SISTEMA de CALCULO
  !===============================================================================
  subroutine finalizacion()
    call finalizar_globales()
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
