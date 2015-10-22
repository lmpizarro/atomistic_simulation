module dinmods 

  use types,      only: dp
  use globales,   only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil, gM, & 
                        gNmed, gRc2, gPot_cut, gRho, gVol, gPot, gKin, gVir 
  use utils,      only: write_array3D_lin, init_openmp
  use mediciones, only: calcula_kin, calcula_pres, calcula_fuerza 
  use constants,  only: PI
  use ziggurat
  use usozig
  use io_parametros,  only: escribe_trayectoria, escribe_estados, lee_estados, &
                            read_parameters

! Si se utiliza openmp
#ifdef _OPENMP
  use omp_lib
#endif

! Si se utiliza mpi
#ifdef MPI
  use mpi
#endif

  implicit none

  private

  public :: inicializacion, inicia_posicion_cs, finalizacion, cpc, &
            integracion_min, integracion, inicia_posicion_rn
  
  !===============================================================================
  ! VARIABLE PROPIAS DEL MÓDULO
  !===============================================================================

  real(dp)        :: Pres      ! Presión instantánea

contains

  !===============================================================================
  ! INICIALIZA SISTEMA de CALCULO
  !===============================================================================

  subroutine inicializacion()

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

    write(*,'(a)') ''
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

    print *, '* Valores iniciales por partícula'
    print *, 'Pot=', gPot/gNpart, 'Kin=' , gKin/gNpart, 'Tot=', (gPot+gKin)/gNPart
    print *, 'Presion= ' , Pres
    print *, '* Se termina la inicialización de parámetros'

  end subroutine inicializacion 

  !===============================================================================
  ! RADIO DE CORTE 
  !===============================================================================
  ! Define el radio de corte y calcula el desplazamiento del potencial de L-J 

  subroutine corta_desplaza_pote()   

    ! Radio de corte fijo
    gRc2 = (2.5_dp)**2                
    ! Potencial de L-J evaluado en el radio de corte
    gPot_cut = 4.0_dp*gEpsil*( 1.0_dp/(gRc2**6) - 1.0_dp/(gRc2**3) ) 
    ! Imprime en pantalla
    write(*,'(a)')      ''
    write(*,'(a)')      '*************** RADIO DE CORTE **************'
    write(*,'(a,F9.4)') '************ Rc    = ' , sqrt(gRc2)
    write(*,'(a,F9.5)') '************ V(Rc) = ' , gPot_cut
    write(*,'(a)')      '*********************************************'

  end subroutine corta_desplaza_pote
  
  !===============================================================================
  ! VELOCIDADES INICIALES 
  !===============================================================================
  ! Subrutina para inicializar las velocidades del problema 

  subroutine vel_inic()   

    integer    :: i, j

!    ! Muestreo una dirección aleatoria uniforme 
!    ! Este loop genera una velocidad de módulo unidad               
!    do i = 1, gNpart
!      phi     = 2.0_dp * PI * uni()
!      gV(3,i) = 1.0_dp - 2.0_dp * uni()
!      gV(1,i) = sqrt( 1.0_dp - gV(3,i)**2 ) * cos(phi)
!      gV(2,i) = gV(1,i) * tan(phi)
!    end do
!    ! Le especifico el módulo
!    ! Se debe poner una Gaussiana
!    gV = 5.0_dp*gV 

    ! Asigna a cada componente de la velocidad una distribución gaussiana N(0,1)
    do i = 1, gNpart
      do j = 1, 3
        gV(j,i) = rnor()
      end do
    end do
    ! Define la desviación estandar sqrt(kT/m) en unidades adimensionales
    gV = sqrt( gT * gEpsil / gM) * gV

  end subroutine vel_inic

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================

  subroutine cpc(l)
  
    integer, intent(in) :: l

    if (gR(1,l) .lt. 0) then
      gR(1,l) = gR(1,l) + gL
    endif        

    if (gR(2,l) .lt. 0) then
     gR(2,l) = gR(2,l) + gL
    endif        

    if (gR(3,l) .lt. 0) then
     gR(3,l) = gR(3,l) + gL
    endif

    if (gR(1,l) .gt. gL) then
     gR(1,l) = gR(1,l) - gL
    endif        

    if (gR(2,l) .gt. gL) then
     gR(2,l) = gR(2,l) - gL
    endif        

    if (gR(3,l) .gt. gL) then
     gR(3,l) = gR(3,l) - gL
    endif

  endsubroutine cpc

  ! Verson vectorizada

  subroutine cpc_vec()

    gR = gR - gL*floor(gR/gL)

  end subroutine

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

      call cpc(l)
           
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

    ! El primer punto es la energía inicial
    Eng_t(1) = gPot

    !call write_array3D_lin(gR)

    do i = 1, gNtime 
      gR = gR + 0.5_dp * gF * gDt**2 / gM
      ! Aplica condiciones peródicas de contorno
      call cpc_vec()    
      ! Esta subrutine abre y cierra un archivo. Se puede optimizar haciéndolo acá.
      call escribe_trayectoria(gR,i)    
      ! Calcula fuerza y energía
      call calcula_fuerza()
      ! Escribe energía potencial en vector
      Eng_t(i+1) = gPot
    end do

    ! Guarda la energía potencial en un archivo
    open(unit=10,file='./energia_pot_min.dat',status='unknown')
    !write(10,'(F10.4)') gDt
    write(10,'(E16.9)') Eng_t
    close(10)
  
    !call write_array3D_lin(gR)
    write(*,*) '* Energía minimizada: ' , gPot

  end subroutine integracion_min

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO - VELOCITY VERLET
  !===============================================================================
  ! Integra las ecuaciones dinámicas con el algoritmo de Velocity-Verlet

  subroutine integracion()
!TODO Ver cómo escbirir energías y presiones. Ahora sólo anda si Nmed=1
    real(dp), dimension(gNtime+1)   :: Eng_t   ! Energía en función del tiempo
    real(dp), dimension(gNtime+1)   :: Pres_t  ! Presión en función del tiempo
    integer    :: i

    Pres_t(1) = Pres
    ! El primer punto es la energía inicial
    Eng_t(1) = gPot + gKin
    write(*,*) '********************************************'
    print *, '* Comienza integracion temporal (Vel-Verlet)'
    print *, '* Energias por partícula al comienzo de la integración'
    print *, 'Pot=' , gPot/gNpart, 'Kin=', gKin/gNpart, 'Tot=', (gPot+gKin)/gNpart

    do i = 1, gNtime 
      gR = gR + gDt*gV + 0.5_dp * gF * gDt**2 / gM     ! gR(t+dt)
      ! Aplica condiciones peródicas de contorno
      call cpc_vec()
      gV =          gV + 0.5_dp * gF * gDt / gM        ! gV(t+0.5dt) 
      call calcula_fuerza()                            ! Calcula fuerzas y potencial
      gV =          gV + 0.5_dp * gF * gDt / gM        ! gV(t+dt)
      ! Se realizan las mediciones
      if (mod(i,gNmed) == 0) then
        ! Energia cinetica
        call calcula_kin()
        ! Presión
        call calcula_pres(Pres)
        ! Escribe energía total
        Eng_t(i+1) = gPot + gKin
        Pres_t(i+1) = Pres
        ! Escribe posiciones de las partículas
        !call escribe_trayectoria(gR,i)
      end if
    end do

    ! Guarda la energía potencial en un archivo
    open(unit=10,file='./energia_tot.dat',status='unknown')
    !write(10,'(F10.4)') gDt
    write(10,'(E16.9)') Eng_t/gNpart
    close(10)

    ! Guarda la presion en un archivo
    open(unit=20,file='./presion.dat',status='unknown')
    !write(10,'(F10.4)') gDt
    write(20,'(E16.9)') Pres_t 
    close(20)
    
    ! En caso de que no se hagan mediciones, se calculan los valores finales
    call calcula_kin()
    call calcula_pres(Pres)
    ! Se imprime en pantalla los resultados finales
    print *, '* Energias por partícula al final de la integración'
    print *, 'Pot=' , gPot/gNpart, 'Kin=', gKin/gNpart, 'Tot=', (gPot+gKin)/gNpart
    print *, 'Presion= ' , Pres
    print *, '* Fin de la integracion temporal'

  end subroutine integracion

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

  endsubroutine finalizacion


end module dinmods 
