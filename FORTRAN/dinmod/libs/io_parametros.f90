module io_parametros

#include "control.h"

  use types,     only: dp
  use globales,  only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM, gNmed, &
                       gGamma

  implicit none

  ! ------------------------------------------------------------------------------
  ! Se crea un wrapper para poder escribir vectores en columnas independientemente
  ! de su rango. Debo duplicar subrutinas.

  interface escribe_en_columnas
    module procedure escribe_en_columnas_1D
    module procedure escribe_en_columnas_2D
  end interface escribe_en_columnas
  ! -----------------------------------------------------------------------------

  private

  public  :: read_parameters, escribe_trayectoria, escribe_estados, lee_estados, &
             escribe_en_columnas

contains

!===============================================================================
! LEE LOS DATOS DE ENTRADA DEl PROBLEMA 
!===============================================================================

  subroutine read_parameters()
    logical :: es
 
    inquire(file='./parametros.dat',exist=es)
    if(es) then
      open(unit=10,file='./parametros.dat',status='old')
      read(10,*) gT, gNpart, gL, gDt, gNtime, gSigma, gEpsil, gM, gNmed
#if THERM == 1
      read(10,*) gGamma
#endif
      close(10)
    else
      print *, "Usando parametros por default"
      gT = 293.0_dp
      gNpart = 1000
      gL = 10.0_dp
      gDt = 0.1_dp
      gNtime = 100
      gSigma = 1.0_dp
      gEpsil = 1.0_dp
      gM = 1.0_dp
    end if
   
    write(*,'(a)') ''
    write(*,'(a)')      '************ PARAMETROS LEIDOS **************'
    write(*,'(a,F8.3)') '************ Temperatura           = ' , gT 
    write(*,'(a,I8)')   '************ Número de partículas  = ' , gNpart 
    write(*,'(a,F8.3)') '************ Masa de la partícula  = ' , gM
    write(*,'(a,F8.3)') '************ Lado del cubo         = ' , gL 
    write(*,'(a,E8.2)') '************ Paso temporal dt      = ' , gDt 
    write(*,'(a,I8)')   '************ Número de pasos dt    = ' , gNtime 
    write(*,'(a,I8)')   '************ Pasos mediciones      = ' , gNmed 
    write(*,'(a)')      '************ Potencial de L-J' 
    write(*,'(a,F8.4)') '************ Epsilon               = ' , gEpsil
    write(*,'(a,F8.4)') '************ Sigma                 = ' , gSigma
    write(*,'(a)')      '*********************************************'

#if THERM == 0
    write(*,'(a)') ''
    write(*,'(a)')      '******** SIMULACION A E CONSTANTE ***********'
    write(*,'(a)')      '*********************************************'
#elif THERM == 1
    write(*,'(a)') ''
    write(*,'(a)')      '******** SIMULACION A T CONSTANTE ***********'
    write(*,'(a)')      '********* TERMOSTATO DE LANGEVIN ************'
    write(*,'(a,F8.3)') '************ Gamma                 = ' , gGamma
    write(*,'(a)')      '*********************************************'
#endif
#ifdef CONTROL_TEMP
    write(*,'(a)') ''
    write(*,'(a)')      '****** COMPILADO CON VERIFICACIÓN DE T ******'
    write(*,'(a)')      '**** Se guardan valores de [vx vy vz]  ******'
    write(*,'(a)')      '**** para una partícula arbitraria     ******'
    write(*,'(a)')      '*********************************************'
#endif
  end subroutine read_parameters

!===============================================================================
! ESBRIBE POSICIONES PARA VISUALIZAR TRAYECTORIAS 
!===============================================================================

  subroutine escribe_trayectoria(r,primera)
    ! Se puede optimizar abriendo y cerrando el archivo en el momento de usarla

    real(dp), dimension(3,gNpart), intent(in)    :: r    ! Posición
    logical, intent(in)                          :: primera
    integer                                      :: i,j

    if ( primera .eqv. .TRUE. ) then
      open(unit=20, file='./trayectoria.vtf',status='unknown')
    else
      open(unit=20, file ='./trayectoria.vtf',status='unknown',position='append')
    end if

    ! Escribe el encabezado del archivo
    if (primera .eqv. .TRUE.) then
      write(20,*) '### Trayectorias ###'
      write(20,*) 'atom default   radius 0.1 name P'
      write(20,'(A,I0)') ' atom 0:',gNpart-1 
      write(20,*) 'timstep'
      write(20,*) ''
    end if
    ! Escribe las coordenadas de cada partícula a un dado tiempo
    do j = 1, gNpart
      write(20,200) (r(i,j), i=1,3)
      200 format (3(F10.5))
    end do
    write(20,*) 'timestep'    
    write(20,*) ''

    close(20)

  end subroutine 

!===============================================================================
! ESBRIBE UN VECTOR 1D A UN ARCHIVO EN COLUMNAS
!===============================================================================
! Escribe al vector con el siguiente formato:
! #puntos  dt
! x(1)
! x(2)
! ...

  subroutine escribe_en_columnas_1D(x,nombre,dt)

    real(dp), dimension(:), intent(in)      :: x       ! Vector que se escribirá
    character(*), intent(in)                :: nombre  ! Nombre del archivo
    real(dp), intent(in)                    :: dt      ! dt entre dos puntos

    
    open(unit=20, file=nombre, status='unknown')

    write(20, *) size(x,1), dt          ! Primer linea con # puntos y dt
    write(20,200) x                     ! Vector 
    200 format (1(E16.9))

    close(20)

  end subroutine escribe_en_columnas_1D 

!===============================================================================
! ESBRIBE UN VECTOR 2D A UN ARCHIVO EN COLUMNAS
!===============================================================================
! Escribe al vector con el siguiente formato:
! #puntos  dt
! x(1,1) x(2,1) x(3,1)
! x(1,2) x(2,2) x(3,2)
! ...  

  subroutine escribe_en_columnas_2D(x,nombre,dt)

    real(dp), dimension(:,:), intent(in)    :: x       ! Vector que se escribirá
    character(*), intent(in)                :: nombre  ! Nombre del archivo
    real(dp), intent(in)                    :: dt      ! dt entre dos puntos
    integer                                 :: npuntos ! # puntos del vector
    integer                                 :: ncompo  ! # compunentes del vector
    integer                                 :: i,j
    character(50)                           :: fmt1
    
    npuntos = size(x,2)
    ncompo  = size(x,1)

    open(unit=20, file=nombre, status='unknown')

    if (ncompo==2) then
      fmt1 = '(2(E16.9,2X))'
    else if (ncompo==3) then
      fmt1 = '(3(E16.9,2X))'
    end if

    write(20, *) npuntos, dt          ! Primer linea con # puntos y dt
    do j = 1, npuntos
      write(20,fmt1)  ( x(i,j) , i = 1, ncompo )   
    end do
    !200 format (3(E16.9,2X))

    close(20)

  end subroutine escribe_en_columnas_2D 

!===============================================================================
! ESBRIBE ESTADOS DE POSICION Y VELOCIDAD
!===============================================================================
! Graba los vectores posición y velocidad en un archivo
! Se pretende utilizar esta información para inicializar una nueva corrida
! Se escribe el archivo 'estado.dat' y el contenido tiene el formato
! N    (Número total de partículas)
! x(1)  y(1)  z(1)  v_x(1)  v_y(1)  v_z(1)
!                ......
! x(N)  y(N)  z(N)  v_x(N)  v_y(N)  v_z(N)

  subroutine escribe_estados(r,v)

    real(dp), dimension(3,gNpart)    :: r, v
    integer                          :: i, j

    open(unit=30, file='./estados.dat', status='replace')
    ! Escribe la cantidad de partículas en la primer línea
    write(30,'(I8)')  gNpart
    ! Escribe posiciones y velocidades de las partículas
    do i = 1, gNpart
     write(30,100) (r(j,i),j=1,3) , (v(j,i),j=1,3)
     100 format (6(E24.17,3X))
    end do 

    close(30)

  end subroutine escribe_estados 

!===============================================================================
! LEE ESTADOS DE POSICION Y VELOCIDAD
!===============================================================================
! Lee los vectores posición y velocidad 
! Primero lee la primera lína para ver si coinciden el número de partículas del
! archivo con las del problema que se está corriendo 

  subroutine lee_estados(r,v,exito)

    real(dp), dimension(3,gNpart),intent(out)   :: r, v  ! Vectores pos. y vel.
    logical, intent(out)                        :: exito ! Controla si se leyo bien
    integer          :: N     ! Número de partículas leidas    
    integer          :: i, j
    logical          :: es

    exito = .FALSE.

    inquire(file='./estados.dat', exist=es)

    if (es) then
      open(unit=20,file='./estados.dat',status='old')
      ! Lee la cantidad de partículas existentes en el archivo
      read(20,'(I8)') N
      if ( N==gNpart ) then
        ! Lee el resto del archivo
        do i = 1, gNpart
            read(20,100) (r(j,i),j=1,3) , (v(j,i),j=1,3)
            100 format (6(E24.117,3X))
        end do 
        write(*,*) '* Archivo de configuracion inicial leido correctamente'
        exito = .TRUE.
        close(20)
      else
        ! No coinciden las cantidades de partículas
        write(*,*) '* La cantidad de datos leidos no coincide con las partículas actuales'
      end if
    else
      ! No existe el archivo
      write(*,*) '* No existe archivo de estados iniciales'
    end if

  end subroutine lee_estados 

end module io_parametros
