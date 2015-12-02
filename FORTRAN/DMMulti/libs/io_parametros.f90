module io_parametros

#include "control.h"

  use types,     only: dp
  use globales

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

  public  :: escribe_trayectoria, escribe_estados, lee_estados, escribe_en_columnas

contains

  subroutine escribe_multiple_atom_vtf(r,nombre,primera)
    ! Se puede optimizar abriendo y cerrando el archivo en el momento de usarla

    real(dp), dimension(3,gNpart), intent(in)    :: r    ! Posición
    character(*), intent(in)                     :: nombre  ! Nombre del archivo
    logical, intent(in)                          :: primera
    integer                                      :: i,j

    if ( primera .eqv. .TRUE. ) then
      open(unit=20, file=nombre, status='unknown')
    else
      open(unit=20, file=nombre, status='unknown',position='append')
    end if

    ! Escribe el encabezado del archivo
    if (primera .eqv. .TRUE.) then
      write(20,*) '### Trayectorias ###'
      ! Atomo tipo 1 el default
      write(20,*) 'atom default   radius 1.0 name P'
      !
      ! aca hay que declarar al o los otros atomos de la estructura
      !
      
      ! buscar en el arreglo de atomos todos los atomos
      ! que no son de tipo 1
      write(20,"(a)", advance='no') 'atom '
      do i=1, gNpart
        if (gIndice_elemento(i) .eq. 2) then
           write (20, "(i4, a)", advance='no') i - 1, ", "
        endif 
      enddo
      write(20,"(a)", advance='no') 'radius 1.1 name S'

      write(20,*) 'timestep'
      write(20,'(A,1X,3(F4.1,1X))') 'pbc', gLado_caja, gLado_caja, gLado_caja
      write(20,*) ''
    end if

  endsubroutine escribe_multiple_atom_vtf

!===============================================================================
! ESBRIBE POSICIONES PARA VISUALIZAR TRAYECTORIAS 
!===============================================================================

  subroutine escribe_trayectoria(r,nombre,primera)
    ! Se puede optimizar abriendo y cerrando el archivo en el momento de usarla

    real(dp), dimension(3,gNpart), intent(in)    :: r    ! Posición
    character(*), intent(in)                     :: nombre  ! Nombre del archivo
    logical, intent(in)                          :: primera
    integer                                      :: i,j

    if ( primera .eqv. .TRUE. ) then
      open(unit=20, file=nombre, status='unknown')
    else
      open(unit=20, file=nombre, status='unknown',position='append')
    end if

    ! Escribe el encabezado del archivo
    if (primera .eqv. .TRUE.) then
      write(20,*) '### Trayectorias ###'
      ! Para escribir una sola especie, es más compacto esto:
      !write(20,*) 'atom default   radius 0.1 name P'
      ! write(20,'(A,I0)') ' atom 0:',gNpart-1 
      ! Para escribir más de una especie:
      do i = 1, gNpart
        write(20,100) ' atom ', i-1, ' radius ', gCombSigma(gIndice_elemento(i),gIndice_elemento(i)), ' type ', gIndice_elemento(i) 
      end do
      100 format (A,I0,A,F5.2,A,I0)
      write(20,'(A,1X,3(F4.1,1X))') 'pbc', gLado_caja, gLado_caja, gLado_caja
      write(20,*) ''
    end if
    ! Escribe las coordenadas de cada partícula a un dado tiempo
    write(20,*) 'timestep'    
    do j = 1, gNpart
      write(20,200) (r(i,j), i=1,3)
      200 format (3(F10.5))
    end do
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
! Graba los vectores posición, velocidad y tipo de partícula en un archivo
! Se pretende utilizar esta información para inicializar una nueva corrida
! Se escribe el archivo 'estado.dat' y el contenido tiene el formato
! N    (Número total de partículas)
! x(1)  y(1)  z(1)  v_x(1)  v_y(1)  v_z(1) tipo_particula(i)
!                       .................
! x(N)  y(N)  z(N)  v_x(N)  v_y(N)  v_z(N  tipo_particula(N))

  subroutine escribe_estados(r,v,tipo)

    real(dp), intent(in), dimension(3,gNpart)    :: r, v  ! Posición y velocidad
    integer,  intent(in), dimension(gNpart)      :: tipo  ! Tipo de partícula
    integer   :: i, j

    open(unit=30, file='./estados.dat', status='replace')
    ! Escribe la cantidad de partículas en la primer línea
    write(30,'(I8)')  gNpart
    ! Escribe posiciones y velocidades de las partículas
    do i = 1, gNpart
     write(30,100) (r(j,i),j=1,3) , (v(j,i),j=1,3) , tipo(i)
     100 format (6(E24.17,3X),I3)
    end do 

    close(30)

  end subroutine escribe_estados 

!===============================================================================
! LEE ESTADOS DE POSICION, VELOCIDAD Y TIPO DE PARTICULA
!===============================================================================
! Lee los vectores posición, velocidad y tipo de partícula
! Primero lee la primera lína para ver si coinciden el número de partículas del
! archivo con las del problema que se está corriendo 

  subroutine lee_estados(r,v,tipo,exito)

    real(dp), intent(out), dimension(3,gNpart) :: r, v  ! Vectores pos. y vel.
    logical, intent(out)                       :: exito ! Controla si se leyo bien
    integer, intent(out), dimension(gNpart)    :: tipo  ! Tipo de partícula
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
        ! Quedan definidos los tamaños de los vectores
        do i = 1, gNpart
            read(20,100) (r(j,i),j=1,3) , (v(j,i),j=1,3), tipo(i)
            100 format (6(E24.17,3X),I3)
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
