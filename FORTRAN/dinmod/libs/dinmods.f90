module dinmods 
    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, sigma, epsil


    implicit none

    private
    public :: inicializacion, potencial
contains
    !===============================================================================
    ! INICIALIZA SISTEMA de CALCULO
    !===============================================================================
    subroutine inicializacion()
        allocate(gR(3,gNpart))
        allocate(gV(3,gNpart))
        allocate(gF(3,gNpart))

    end subroutine inicializacion 

    ! calculo del potencial de pag 18 del allen-tildesley
    function potencial () result(v)
        integer :: i, j     
        real(dp) :: v, rxi, ryi, rzi
        real(dp) :: rxij, ryij, rzij, rijsq
        real(dp) :: sr2, sr6, sr12

        v = 0

        do i = 1, gNpart - 1

            rxi = gR(1, i)
            ryi = gR(2, i)
            rzi = gR(3, i)

            do j = i + 1, gNpart
                rxij = rxi - gR(1,j)
                ryij = ryi - gR(2,j)
                rzij = rzi - gR(3,j)
    
                rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
                sr2 = sigma / rijsq
                sr6 = sr2 **3
                sr12 = sr6 ** 2
                v = v + sr12 - sr6
            enddo
        enddo
        v = 4.0 * epsil * v
    endfunction potencial

end module dinmods 
