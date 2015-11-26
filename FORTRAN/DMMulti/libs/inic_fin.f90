module inic_fin 
#include "control.h"
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


    call comb_sigma()
    call comb_epsilon()

    ! como se calcularon los sigmas
    ! calculamos los rc y potencial en rc
    call corta_desplaza_pote()

    if (gLiqSol .eq. 0) then 
       call inicializar_globales_random()
       call inicia_posicion_rn()
    else if ( gLiqSol .eq. 1) then
      print *, "llama a posiciones en  cubica"
      print *, "parametros cubica periodos: ", gPeriodos, "tipo: ", gCubicStructure
      if (gCubicStructure .eq. 0) then
        print *, "llama a cubica simple"
        print *, "no implementado"
        stop 212121212
      else if (gCubicStructure .eq. 1) then 
        print *, "llama a cubica centrada en el cuerpo"
        print *, "no implementado"
        stop 212121212
      else if (gCubicStructure .eq. 2) then 
        print *, "llama a cubica centrada en las caras"
        call inicializar_globales_fcc()
        call inicia_posicion_fcc(gPeriodos)
      endif
    else
      print *, "no implementado"
      stop 212121212
    endif        

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

  !===============================================================================
  ! INICIALIZA Posicion fcc dentro de la caga 
  !===============================================================================
  ! param_red: parametro de red, la caja es cubica(todos los lados iguales) 
  ! n_pred: indica la cantidad de veces que se repite el cubo primitivo
  ! tanto en z como en x y la variable y
  !
  subroutine inicia_posicion_fcc(n_pred)
 
    integer :: i, j, k, n_pred, i_gr
    print *, "inicia posicion red fcc"

    if (gPercent(1).ne. 0.5 ) then
      stop 1212121212
    endif

    print *, "cincuenta"

    ! inicia para un binario
    gIndice_elemento = 2
    i_gr = 1
    do i = 1, n_pred 
      do j= 1,  n_pred
        do k= 1,  n_pred
           gR(1, i_gr) = i - 1 
           gR(2, i_gr) = j - 1
           gR(3, i_gr) = k - 1
           i_gr = i_gr + 1
           gR(1, i_gr) = i - 1 
           gR(2, i_gr) = j - 0.5
           gR(3, i_gr) = k - 0.5
           gIndice_elemento(i_gr) = 1
           i_gr = i_gr + 1
           gR(1, i_gr) = i - 0.5 
           gR(2, i_gr) = j - 1 
           gR(3, i_gr) = k - 0.5
           gIndice_elemento(i_gr) = 1
           i_gr = i_gr + 1
           gR(1, i_gr) = i - 0.5 
           gR(2, i_gr) = j - 0.5 
           gR(3, i_gr) = k - 1 
           i_gr = i_gr + 1
           ! dividir por 4
           ! colocar dentro de la caja
           !  gLado_caja * gPeriodos / (gPeriodos + 1)
        enddo
      enddo
   enddo

   gR = gLado_caja * gR / gPeriodos

#if DEBUG == 1
   do i=1, gNpart
      print *, "punto ", gR(1, i), gR(2, i), gR(3, i), gIndice_elemento(i)
   enddo
#endif

  end subroutine inicia_posicion_fcc


end module inic_fin 
