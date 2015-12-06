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
  use io_parametros,    only: lee_estados, escribe_estados
  use utils,            only: write_array3D_lin ! No se usa. Se pone para que cmake procese dependencias

! Si se utiliza openmp
#ifdef _OPENMP
  use omp_lib
  use utils,      only: init_openmp
#endif

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
    logical   :: leido  ! Flag para saber si se leyó el archivo de estados

  ! Inicializa parametros para el calculo en paralelo con OpenMP
  ! Por ahora nada en particular.
#ifdef _OPENMP
  call init_openmp()
#endif

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
#ifdef CORR_PAR
    ! Inicializa las variables asociadas para el cálculo de la g(r)
    call inic_gr()
#endif  
 
     if (gLiqSol .eq. 0) then 
        call inicializar_globales_random()
      else if ( gLiqSol .eq. 1) then
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
        endif
      else
        print *, "no implementado"
        stop 212121212
      endif  
    ! Trata de leer el archivo con la configuración inicial
    call lee_estados(gR,gV,gIndice_elemento,leido)

    !
    ! Si no se leyó el archivo de estados, lo define
    !
    if (leido .eqv. .FALSE. ) then
    !  write(*,*) '* Se generan posiciones y velocidades'
    ! Bloque para inicializar algunas variables globales
    ! y posiciones dependiendo de los datos de entrada
    !
    if (gLiqSol .eq. 0) then 
        ! call inicializar_globales_random()
        call inicia_posicion_rn()
      else if ( gLiqSol .eq. 1) then
        print *, "llama a posiciones en  cubica"
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
        !  call inicializar_globales_fcc()
          call inicia_posicion_fcc(gPeriodos)
        endif
      else
        print *, "no implementado"
        stop 212121212
      endif        
      ! Inicializa las velocidades 
      call vel_inic()
    end if

    !
    !--- En este punto, todas las variables y parámetros quedan definidos 
    !

    call calcula_fuerza()
    
    print *, "energia potencial inicial: ", gPot
    write(*,*) '********************************************'

   !call integracion_min()
    ! Calcula energía cinética
    call calcula_kin()
    ! Calcula fuerza
    call calcula_fuerza()
    ! Calcula la presión inicial
    call calcula_pres(Pres)
    ! Calcula la temperatura inicial
    call calcula_temp(Temp)

#ifdef CORR_PAR
    ! Se vuelve a llamar la inicialización para resetear el contador. Como la g(r)
    ! está implementada en el loop de fuerzas, cada vez que esa subrutina es llamada
    ! actualiza los valores aunque no sean de interés. Ver en la subrutina inic_gr()
    ! que hace cosas distintas dependiendo si gCorr_par() está alocada o no.
    call inic_gr()
#endif
    ! Escribe información en pantalla
    write(*,*) '* Valores iniciales por partícula'
    write(*,100) gPot/gNpart, gKin/gNpart, (gPot+gKin)/gNpart
    100 format(1X,'Potencial = ', E14.7, 5X, 'Cinética = ', E14.7, 5X, 'Total = ', E14.7) 
    write(*,*)
    write(*,200)  Pres, Temp 
    200 format(1X,'Presion = ' ,E14.7,5X,'Temperatura = ',E14.7)

   
    write(*,*) '* Finaliza subrutina de inicialización'
    write(*,*) '********************************************'

  end subroutine inicializacion

  !===============================================================================
  ! FINALIZA SISTEMA de CALCULO
  !===============================================================================

  subroutine finalizacion()

    ! Escribe la última configuración de x, v e indice de elemento en el
    ! archivo 'estado.dat'
    call escribe_estados(gR,gV,gIndice_elemento)
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
      ! Se tiene en cuenta la masa, pues sigma=sqrt(kt/m) 
      gV(:,i) = gV(:,i)*sqrt(1/gMasa(gIndice_elemento(i)))
    end do
    ! Define la desviación estandar sqrt(kT/m) en unidades adimensionales
    gV = sqrt( gTemperatura ) * gV

  end subroutine vel_inic

  !===============================================================================
  ! INICIALIZA Posicion aleatoria dentro de la caga 
  !===============================================================================

  subroutine inicia_posicion_rn()
 
    integer :: i, j, inic, fin

    do i = 1, gNpart
      do j= 1, 3
        gR(j, i) = uni() * gLado_caja
      end do
   end do     

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
           i_gr = i_gr + 1
           gR(1, i_gr) = i - 0.5 
           gR(2, i_gr) = j - 1 
           gR(3, i_gr) = k - 0.5
           i_gr = i_gr + 1
           gR(1, i_gr) = i - 0.5 
           gR(2, i_gr) = j - 0.5 
           gR(3, i_gr) = k - 1 
           i_gr = i_gr + 1
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
      allocate(gCorr_par(1:3,1:gNhist))

      gNgr      = 0                  ! Contador para saber cuántas veces se acumuló la g(r)
      gDbin     = gLado_caja / (2 * gNhist)  ! Ancho del bin
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


end module inic_fin 
