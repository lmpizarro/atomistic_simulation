module io_parametros
    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, sigma, epsil

    implicit none

    public  :: read_parameters
contains

!===============================================================================
! LEE LOS DATOS DE ENTRADA DEl PROBLEMA 
!===============================================================================

  subroutine read_parameters()
        logical :: es
     
        inquire(file='./parameters.dat',exist=es)
        if(es) then
            open(unit=10,file='./parameters.dat',status='old')
            read(10,*) gT, gNpart, gL, gDt, gNtime, sigma, epsil
            close(10)
        else
            print *, "Usando parametros por default"
            gT = 293
            gNpart = 1000
            gL = 10
            gDt = .1
            gNtime = 100
            sigma = 1.0
            epsil = 1.0
        end if
        print *,"* read_parameters ", gT, gNpart, gL, gDt, gNtime, sigma, epsil

  end subroutine read_parameters


end module io_parametros
