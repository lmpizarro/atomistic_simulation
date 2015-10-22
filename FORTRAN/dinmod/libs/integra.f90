module integra

  use types,          only: dp
  use globales,       only: gt, gL,gDt, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil, gM, & 
                            gNmed, gPot, gKin, gVir 
  use mediciones,     only: calcula_kin, calcula_pres, calcula_fuerza 

implicit none

private

public      :: integracion, cpc_vec

contains

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO - VELOCITY VERLET
  !===============================================================================
  ! Integra las ecuaciones dinámicas con el algoritmo de Velocity-Verlet

  subroutine integracion()

    real(dp)        :: Pres      ! Presión instantánea

!TODO Ver cómo escbirir energías y presiones. Ahora sólo anda si Nmed=1
    real(dp), dimension(3,gNtime+1)   :: Eng_t   ! Energía en función del tiempo
    real(dp), dimension(gNtime+1)   :: Pres_t  ! Presión en función del tiempo
    integer    :: i, j

    Pres_t(1) = Pres
    ! El primer punto es la energía inicial
    Eng_t(:,1) = (/ gPot, gKin, gPot + gKin/) 
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
        Eng_t(:,i+1) = (/gPot, gKin, gPot + gKin/) 
        Pres_t(i+1) = Pres
        ! Escribe posiciones de las partículas
        !call escribe_trayectoria(gR,i)
      end if
    end do

    ! Guarda la energía potencial en un archivo
    open(unit=10,file='./energia_tot.dat',status='unknown')
    !write(10,'(F10.4)') gDt
    do j= 1,gNtime+1
      write(10,'(4(E16.9,2X))') (/ (j-1)*gDt, ( Eng_t(i,j)/gNpart, i=1,3) /)
    end do
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
  ! Condiciones períodicas de contorno
  !===============================================================================

  subroutine cpc_vec()

    gR = gR - gL*floor(gR/gL)

  end subroutine

end module integra
