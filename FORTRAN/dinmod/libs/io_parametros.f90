module io_parametros

  use types,     only: dp
  use globales,  only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM, gNmed

  implicit none

 ! public  :: read_parameters, escribe_trayectoria

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
    write(*,'(a,F8.5)') '************ Paso temporal dt      = ' , gDt 
    write(*,'(a,I8)')   '************ Número de pasos dt    = ' , gNtime 
    write(*,'(a,I8)')   '************ Pasos mediciones      = ' , gNmed 
    write(*,'(a)')      '************ Potencial de L-J' 
    write(*,'(a,F8.4)') '************ Epsilon               = ' , gEpsil
    write(*,'(a,F8.4)') '************ Sigma                 = ' , gSigma
    write(*,'(a)')      '*********************************************'

    
  end subroutine read_parameters

!===============================================================================
! ESBRIBE POSICIONES PARA VISUALIZAR TRAYECTORIAS 
!===============================================================================

  subroutine escribe_trayectoria(r,k)
    ! Se puede optimizar abriendo y cerrando el archivo en el momento de usarla

    real(dp), dimension(3,gNpart), intent(in)    :: r    ! Posición
    integer                                      :: i,j
    integer, intent(in)                          :: k

    if ( k==1 ) then
      open(unit=20, file='./trayectoria.vtf',status='unknown')
    else
      open(unit=20, file ='./trayectoria.vtf',status='unknown',position='append')
    end if

    ! Escribe el encabezado del archivo
    if (k == 1) then
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
            read(20,*) (r(j,i),j=1,3) , (v(j,i),j=1,3)
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
