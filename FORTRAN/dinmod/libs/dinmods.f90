module dinmods 
    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, sigma, epsil
    use utils, only: write_array3D_lin
    use ziggurat


    implicit none

    private
    public :: inicializacion, potencial, inicia_posicion_cs

contains
    !===============================================================================
    ! INICIALIZA SISTEMA de CALCULO
    !===============================================================================
    subroutine inicializacion()
        allocate(gR(3,gNpart))
        allocate(gV(3,gNpart))
        allocate(gF(3,gNpart))

    end subroutine inicializacion 

   
    !===============================================================================
    ! Condiciones períodicas de contorno
    !===============================================================================
    subroutine cpc(nl, l)
            integer, intent(in) :: l, nl
            if (gR(1,l) .lt. 0) then
                    gR(1,l) = gR(1,l) + gL
            endif        

            if (gR(2,l) .lt. 0) then
                    gR(2,l) = gR(2,l) + gL
            endif        

            if (gR(3,l) .lt. 0) then
                    gR(3,l) = gR(3,l) + gL
            endif

            if (gR(1,l) .gt. nL) then
                    gR(1,l) = gR(1,l) - nL
            endif        

            if (gR(2,l) .gt. nL) then
                    gR(2,l) = gR(2,l) - nL
            endif        

            if (gR(3,l) .gt. nL) then
                    gR(3,l) = gR(3,l) - nL
            endif

    endsubroutine cpc

    !===============================================================================
    ! INICIALIZA Posicion en Red periódica cúbica simple
    !===============================================================================
    subroutine inicia_posicion_cs()
        integer :: i, j, k, l, nl
        real(dp) :: rx, ry, rz

        nl = gNpart ** (1.0/3.0) 

        print *, "nl", nl

        rx = 0.0
        ry = 0.0
        rz = 0.0
        j = 1
        i = 1
        k = 1
        do l = 1, gNpart
            gR(1, l) = rx + uni() - 0.5 
            gR(2, l) = ry + uni() - 0.5 
            gR(3, l) = rz + uni() - 0.5 

            call cpc(nl, l)
           
            j = j + 1
            k = k + 1
            rx = rx + 1
            if ( mod(j, nl) .eq. 0) then
                rx = 0.0
                ry = ry  + 1
            end if        
            if ( mod(k, nl ** 2) .eq. 0) then
                ry = 0.0
                rx = 0.0
                rz = rz + 1
            end if        
        enddo     
        call write_array3D_lin (gR)
        print *, "potencial", potencial()
    end subroutine inicia_posicion_cs

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
