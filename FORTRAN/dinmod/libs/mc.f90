module mc
    use ziggurat
    use usozig
    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil,&
        gAbs_ener
    use constants
    use potenciales

    implicit none

    private        
    public :: metropolis


contains

    subroutine metropolis ()
        real(dp) :: rx, ry, rz
        real(dp) :: n_energy, deltaE
        real(dp):: pr, beta
        integer :: iPart, i

        !beta = 1.0 / (gT * K_BOLTZMANN)
        beta = 1.0 / (gT * kb)

        allocate(gAbs_ener( gNtime))

        n_energy = poten_lj() 

        do i=1, gNtime
            iPart = rand_int(gNpart)
            ! guardo el valor original de la energia
            rx = gR(1, iPart) 
            ry = gR(2, iPart)
            rz = gR(3, iPart)

            gR(1, iPart) = gR(1, iPart) + .1*gSigma * uni()
            gR(1, iPart) = gR(2, iPart) + .1*gSigma * uni()
            gR(1, iPart) = gR(3, iPart) + .1*gSigma * uni()
       
            call cpc(iPart)

            deltaE = poten_lj() - n_energy

            if (deltaE .lt. 0) then
                n_energy = n_energy + deltaE    
            else
                pr = exp(-beta * deltaE)
                if (uni() .lt. pr) then
                    n_energy = n_energy + deltaE    
                else
                    gR(1, iPart) = rx
                    gR(2, iPart) = ry
                    gR(3, iPart) = rz
                endif        
            endif        
            gAbs_ener(i) = n_energy
        enddo
        print *, gAbs_ener
    endsubroutine metropolis

endmodule mc
