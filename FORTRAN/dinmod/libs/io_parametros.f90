module io_parametros

  use types,     only: dp
  use globales,  only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM

  implicit none

  public  :: read_parameters, escribe_trayectoria

contains

!===============================================================================
! LEE LOS DATOS DE ENTRADA DEl PROBLEMA 
!===============================================================================

  subroutine read_parameters()
    logical :: es
 
    inquire(file='./parametros.dat',exist=es)
    if(es) then
        open(unit=10,file='./parametros.dat',status='old')
        read(10,*) gT, gNpart, gL, gDt, gNtime, gSigma, gEpsil, gM
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

    write(*, 700)  gT, gNpart, gL, gDt, gNtime, gSigma, gEpsil, gM
    700 format (F8.3 I7 F8.3 F8.3 I7 F8.3 F8.3 F8.3)
    
  end subroutine read_parameters

!===============================================================================
! ESBRIBE POSICIONES PARA VISUALIZAR TRAYECTORIAS 
!===============================================================================

  subroutine escribe_trayectoria(r)
    ! Se puede optimizar abriendo y cerrando el archivo en el momento de usarla

    real(dp), dimension(3,gNpart), intent(in)    :: r    ! Posición
    integer                                      :: i,j

    open(unit=20, file ='./trayectoria.vtf',status='unknown',position='append')

    do j = 1, gNpart
      write(20,200) (r(i,j), i=1,3)
      200 format (3(F10.5))
    end do
      write(20,*) 'timestep'    
      write(20,*) ''
    close(20)

  end subroutine


end module io_parametros
