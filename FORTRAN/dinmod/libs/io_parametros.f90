module io_parametros
    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil

    implicit none

    public  :: read_parameters
contains

!===============================================================================
! LEE LOS DATOS DE ENTRADA DEl PROBLEMA 
!===============================================================================

  subroutine read_parameters()
        logical :: es
     
        inquire(file='./parametros.dat',exist=es)
        if(es) then
            open(unit=10,file='./parametros.dat',status='old')
            read(10,*) gT, gNpart, gL, gDt, gNtime, gSigma, gEpsil
            close(10)
        else
            print *, "Usando parametros por default"
            gT = 293.0_dp
            gNpart = 1000
            gL = 10.0_dp
            gDt = .1_dp
            gNtime = 100
            gSigma = 1.0_dp
            gEpsil = 1.0_dp
        end if
        write(*, 700)  gT, gNpart, gL, gDt, gNtime, gSigma, gEpsil

        700 format (F8.3 I7 F8.3 F8.3 I7 F8.3 F8.3)
  end subroutine read_parameters


end module io_parametros
