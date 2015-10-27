module integra

  use types,          only: dp
  use globales,       only: gt, gL,gDt, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil, gM, & 
                            gNmed, gPot, gKin, gVir 
  use mediciones,     only: calcula_kin, calcula_pres, calcula_fuerza, calcula_temp 

implicit none

private

public      :: integracion, cpc_vec

contains

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO - VELOCITY VERLET
  !===============================================================================
  ! Integra las ecuaciones dinámicas con el algoritmo de Velocity-Verlet

  subroutine integracion()

    real(dp)                                :: Pres    ! Presión instantánea
    real(dp)                                :: Temp    ! Presión instantánea
    real(dp), dimension(:,:), allocatable   :: Eng_t   ! Energía en función del tiempo
    real(dp), dimension(:), allocatable     :: Pres_t  ! Presión en función del tiempo
    real(dp), dimension(:), allocatable     :: Temp_t  ! Temperatura en función del tiempo
    integer    :: i, j
    integer    :: Kmed                                 ! Cantidad de puntos medidos

    ! Se define la cantidad de puntos que se van a medir
    Kmed = int(gNtime/abs(gNmed)) + 1              ! Se agrega +1 para poner el inicial
    allocate( Eng_t(1:3,1:Kmed), Pres_t(1:Kmed), Temp_t(1:Kmed) )

    ! El primer punto son los valores iniciales
    call calcula_pres(Pres)
    Pres_t(1) = Pres
    Eng_t(:,1) = (/ gPot, gKin, gPot + gKin/) 
    ! Imprime en pantalla info
    write(*,*) '********************************************'
    print *, '* Comienza integracion temporal (Vel-Verlet)'
    print *, '* Energias por partícula al comienzo de la integración'
    call print_info(Pres)

    j = 1                                              ! Contador para el loop de mediciones
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
        ! Temperatura
        call calcula_temp(Temp)
        ! Presión
        call calcula_pres(Pres)
        ! Escribe energía total
        Eng_t(:,j) = (/gPot, gKin, gPot + gKin/) 
        Pres_t(j)  = Pres
        Temp_t(j)  = Temp
        ! Escribe posiciones de las partículas
        !call escribe_trayectoria(gR,i)
        j = j + 1                                      ! Actualiza contador mediciones
      end if
    end do

    if (gNmed > 0) then
      ! Guarda la energía potencial en un archivo
      open(unit=10,file='./energia_tot.dat',status='unknown')
      !write(10,'(F10.4)') gDt
      do j= 1,Kmed
        write(10,'(4(E16.9,2X))') (/ (j-1)*gDt, ( Eng_t(i,j)/gNpart, i=1,3) /)
      end do
      close(10)

      ! Guarda la presion en un archivo
      open(unit=20,file='./presion.dat',status='unknown')
      !write(10,'(F10.4)') gDt
      write(20,'(E16.9)') Pres_t 
      close(20)
      
      ! Guarda la temperatura en un archivo
      open(unit=30,file='./temperatura.dat',status='unknown')
      !write(10,'(F10.4)') gDt
      write(30,'(E16.9)') Temp_t 
      close(30)

    end if
    
    ! Se imprime en pantalla los resultados finales
    print *, '* Energias por partícula al final de la integración'
    call print_info(Pres)
    print *, '* Fin de la integracion temporal'
    ! Se libera memoria
    deallocate( Eng_t, Pres_t, Temp_t )

contains

    subroutine print_info(presion)
    ! Subrutina para imprimir en pantalla resultados de interes
 
      real(dp), intent(in)        :: presion     ! Presión instantánea

      print *, 'Pot=' , gPot/gNpart, 'Kin=', gKin/gNpart, 'Tot=', (gPot+gKin)/gNpart
      print *, 'Presion= ' , presion
     
    end subroutine print_info

  end subroutine integracion

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================

  subroutine cpc_vec()

    gR = gR - gL*floor(gR/gL)

  end subroutine

end module integra
