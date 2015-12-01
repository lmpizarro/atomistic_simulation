!===============================================================================
! INICIALIZACION Y FINALIZACION DEL PROGRAMA 
!===============================================================================
!
module inic_fin 

#include "control.h"

  use types,            only: dp
  use globales
  use ziggurat
  use usozig,           only: inic_zig, fin_zig
  use read_param,       only: leer_parametros
  use combinacion,      only: comb_sigma, comb_epsilon, corta_desplaza_pote
  use mediciones,       only: calcula_fuerza, calcula_kin, calcula_pres, calcula_temp
  use integra

  implicit none
  
  private

  public  :: inicializacion, finalizacion

contains

  !===============================================================================
  ! INICIALIZA SISTEMA de CALCULO
  !===============================================================================

  subroutine inicializacion()

    integer   :: j 
    real(dp)  :: Pres   ! Presión instantánea
    real(dp)  :: Temp   ! Temperatura instantánea

    call inic_zig()
    call leer_parametros()

    ! Combina los sigma de cada especie ii
    ! para obtener el de interacción entre especies ij
    call comb_sigma()
    ! Combina los epsilon de cada especie ii
    ! para obtener el de interacción entre especies ij
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

   ! print *, "pos inic"
   ! do j=1,gNpart 
   !   print *, gR(1,j) , gR(2,j), gR(3,j)
   ! enddo  

    call calcula_fuerza()
    
    print *, "energia potencial inicial: ", gPot
    write(*,*) '***************************************'

   !call integracion_min()

   ! print *, "pos final after min"
   ! do j=1,gNpart 
   !   print *, gR(1,j) , gR(2,j), gR(3,j)
   ! enddo  

  !print *, "energia potencial after min: ", gPot
  !
  ! Inicializa las velocidades de todas las partículas
    call vel_inic()

    call calcula_kin()

 !   write (*,100) "gKin inicial: ", gKin

    call calcula_fuerza()

    ! Calcula la presión inicial
    call calcula_pres(Pres)
    ! Calcula la temperatura inicial
    call calcula_temp(Temp)
!    write (*,100) "insTemp inicial: ", Temp
!    write (*,100) "insPres inicial: ", Pres

!    100 format (a, F20.3)

    write(*,*) '* Valores iniciales por partícula'
    write(*,100) gPot/gNpart, gKin/gNpart, (gPot+gKin)/gNpart
    100 format(1X,'Potencial = ', E14.7, 5X, 'Cinética = ', E14.7, 5X, 'Total = ', E14.7) 
    write(*,*)
    write(*,200)  Pres, Temp 
    200 format(1X,'Presion = ' ,E14.7,5X,'Temperatura = ',E14.7)

   
    write(*,*) '* Finaliza subrutina de inicializacións'
    write(*,*) '***************************************'

  end subroutine inicializacion

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
  ! TODO: Está mal, hay que cambiarlo
  subroutine vel_inic()   

    integer    :: i, j

    ! Asigna a cada componente de la velocidad una distribución gaussiana N(0,1)
    do i = 1, gNpart
      do j = 1, 3
        gV(j,i) = 0.0_dp !rnor()
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
  ! INICIALIZA POSICIONES EN UNA RED FCC  
  !===============================================================================
  ! La caja es cubica (todos los lados son iguales) 
  ! Se inicializa con una mezcla fija de especies 1 y 2

  subroutine inicia_posicion_fcc(n_pred)
  
    ! Esto es en verdad una variable gloabal, se puede omitir el parámetro. 
    integer, intent(in)    :: n_pred    ! Periodicidad de la red, # veces que se
                                        ! repite la red en las direcciones x, y, z
    ! De define las coordenadas de los cuatro puntos de la celda unidad para
    ! la red FCC. Se trabaja con una celda unidad de tamaño unidad.
    real(dp), dimension(3,4), parameter :: ucell =  &  ! Celda unidad para la FCC
            reshape ( (/ 0.0_dp , 0.0_dp , 0.0_dp , &  ! Posición 1 de la ucell
                         0.5_dp , 0.5_dp , 0.0_dp , &  ! Posición 2 de la ucell
                         0.0_dp , 0.5_dp , 0.5_dp , &  ! Posición 3 de la ucell
                         0.5_dp , 0.0_dp , 0.5_dp   &  ! Posición 4 de la ucell
                       /) , (/3,4/) )
    integer :: m       ! Auxiliar para mapear todos los puntos de la red en gR
    integer :: i_gr    ! Indice que recorre las cuatro posiciones de la celda unidad
    integer :: i, j, k ! Indices para trasladar la celda unidad en x, y, z
        
    ! --- Para debug
    print *, "inicia posicion red fcc"
    if (gPercent(1).ne. 0.5 ) then
      stop 1212121212
    endif
    ! --- Fin debug

    ! Se inicializa con todos los elementos de especie 1
    gIndice_elemento = 1
    ! Comienza la construcción del a red FCC
    m = 0
    do i = 0, n_pred - 1       ! Traslada en la dirección x
      do j = 0,  n_pred - 1    ! Traslada en la dirección y
        do k = 0,  n_pred - 1  ! Traslada en la dirección z
          do i_gr = 1, 4       ! Recorre los cuatro componentes de la celda unidad
            gR(1, i_gr+m) = ucell(1,i_gr) + real(i,dp) 
            gR(2, i_gr+m) = ucell(2,i_gr) + real(j,dp)
            gR(3, i_gr+m) = ucell(3,i_gr) + real(k,dp)
            ! En las posiciones 2 y 3 de la ucell hay partículas de especie 2
            if ( (i_gr == 2) .or. (i_gr == 3) ) then
              gIndice_elemento(i_gr+m) = 2
            end if  
          end do
          m = m + 4   ! Se actualiza para trabajar con las 4 próximas posiciones
        enddo
      enddo
   enddo

   ! Se pasa a las unidades del problema.
   ! 1) Se divide por la cantidad de repeticiones de la celda unidad para
   !    obtener una caja de lado unidad
   ! 2) Se multiplica por la longitud deseada de la caja
   !
   gR = gLado_caja * gR / n_pred 

#if DEBUG == 1
   do i=1, gNpart
      print *, "punto ", gR(1, i), gR(2, i), gR(3, i), gIndice_elemento(i)
   enddo
#endif

  end subroutine inicia_posicion_fcc

  subroutine inicia_posicion_fcc_random(n_pred)
 
    integer, intent(in)    :: n_pred     ! Periodicidad de la red, # veces que se
                                        ! repite la red en las direcciones x, y e z
    integer :: i, j, k, i_gr

    real (dp) :: tmp

    print *, "inicia posicion red fcc"

    if (gPercent(1).ne. 0.5 ) then
      stop 1212121212
    endif

    tmp = uni()
    ! inicia para un binario el arreglo de indice
    ! de elementos a 2
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

   gR = gLado_caja * gR / n_pred

   ! para un sistema binario cambio la naturaleza de la 
   ! partícula para la mitad de los elementos en forma aleatoria
   do i=1, gNpart / 2
      gIndice_elemento(int(gNpart * uni())) = 1
   enddo
   ! En una futura vesión debería ser una rutina independiente

#if DEBUG == 1
   do i=1, gNpart
      print *, "punto ", gR(1, i), gR(2, i), gR(3, i), gIndice_elemento(i)
   enddo
#endif

  end subroutine inicia_posicion_fcc_random

end module inic_fin 
