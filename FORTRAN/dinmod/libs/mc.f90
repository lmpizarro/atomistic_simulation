module mc
    use ziggurat
    use usozig
    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, sigma, epsil
    use dinmods

    implicit none

    private        
    public :: metropolis


contains

    subroutine metropolis ()
        real(dp) :: rx, ry, rz
        real(dp) :: n_energy
        integer :: iPart

        n_energy = potencial() 

        do i=1, gNtime
            iPart = rand_int(gNpart)
            ! guardo el valor original de la energia
            rx = gR(1, iPart) 
            ry = gR(2, iPart)
            rz = gR(3, iPart)

            gR(1, iPart) = gR(1, iPart) + sigma * uni()
            gR(1, iPart) = gR(2, iPart) + sigma * uni()
            gR(1, iPart) = gR(3, iPart) + sigma * uni()

            deltaE = potencial() - n_energy

            if (deltaE .lt. 0) then
                n_energy = n_energy + deltaE    
                print *, "valor aceptado", n_energy
            else
                gR(1, iPart) = rx
                gR(2, iPart) = ry
                gR(3, iPart) = rz
            endif        
        enddo
    endsubroutine metropolis

endmodule mc
