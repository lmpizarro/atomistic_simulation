module dinmods 

    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, sigma, epsil
    use utils, only: write_array3D_lin
    use ziggurat
    use usozig

    implicit none

    private

    public :: inicializacion, inicia_posicion_cs, finalizacion, cpc, fuerza

contains

    !===============================================================================
    ! INICIALIZA SISTEMA de CALCULO
    !===============================================================================
    subroutine inicializacion()

        allocate(gR(3,gNpart))
        allocate(gV(3,gNpart))
        allocate(gF(3,gNpart))

        call inic_zig()

    end subroutine inicializacion 

   
    !===============================================================================
    ! Condiciones períodicas de contorno
    !===============================================================================

    subroutine cpc( l)
            integer, intent(in) :: l

            if (gR(1,l) .lt. 0) then
                    gR(1,l) = gR(1,l) + gL
            endif        

            if (gR(2,l) .lt. 0) then
                    gR(2,l) = gR(2,l) + gL
            endif        

            if (gR(3,l) .lt. 0) then
                    gR(3,l) = gR(3,l) + gL
            endif

            if (gR(1,l) .gt. gL) then
                    gR(1,l) = gR(1,l) - gL
            endif        

            if (gR(2,l) .gt. gL) then
                    gR(2,l) = gR(2,l) - gL
            endif        

            if (gR(3,l) .gt. gL) then
                    gR(3,l) = gR(3,l) - gL
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
            gR(1, l) = (rx + uni() - 0.5) * gL / nl 
            gR(2, l) = (ry + uni() - 0.5) * gL / nl 
            gR(3, l) = (rz + uni() - 0.5) * gL / nl 

            call cpc( l)
           
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
    end subroutine inicia_posicion_cs

    !===============================================================================
    ! CALCULA LA FUERZA  
    !===============================================================================

    subroutine fuerza

    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: Fij       ! Módulo fuerza entre partículas i y j
    real(dp)                :: Pot       ! Energía potencial
    real(dp)                :: rc2       ! Cuadrado del radio de corte
    real(dp)                :: pot_cut   ! Potencial L-J en el radio de corte
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    integer                 :: i,j       

    rc2 = (2.5_dp)**2                ! Definición provisoria del radio de corte
    pot_cut = 4*epsil*( 1/(rc2**6) - 1/(rc2**3) ) 

      ! Se van a acumular las fuerzas. Se comienza poniendo todas a cero.
      gF  = 0.0_dp
      Pot = 0.0_dp
      ! Paso a trabajar distancias en unidades de sigma
      gR = gR/sigma
      gL = gL/sigma

      do i = 1, gNpart - 1       
        do j = i+1, gNpart
          rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
          ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
          ! Siempre en distancias relativas de sigma
          rij_vec = rij_vec - gL*nint(rij_vec/gL)
          r2ij   = dot_product( rij_vec , rij_vec )    ! Cuadrado de la distancia
          if ( r2ij < rc2 ) then               
            r2in = 1/r2ij                              ! Inversa al cuadrado
            r6in = r2in**3                             ! Inversa a la sexta
            Fij     = r2in * r6in * (r6in - 0.5_dp)    ! Fuerza entre partículas
            gF(:,i) = gF(:,i) + Fij * rij_vec          ! Contribución a la partícula i
            gF(:,j) = gF(:,j) - Fij * rij_vec          ! Contribucion a la partícula j
            Pot     = Pot + r6in * ( r6in - 1.0_dp)    ! Energía potencial
          end if
        end do
      end do

      ! Constantes que faltaban en la energía
      Fij = 48 * epsil * Fij                
      ! Constantes que faltaban en el potencial
      ! Agrego el desplazamiento del potencial considerando la cantidad de
      ! pares con que se obtuvo la energía potencial N(N-1)/2
      Pot =  4 * epsil * Pot - gNpart*(gNpart-1)*pot_cut/2  
 
      ! Se vuelven a pasar a las coordenadas absolutas
      gR = gR*sigma
      gL = gL*sigma

      write(*,'(A,2X,3(E15.5,3X))')  'Sumatoria de fuerzas:' , sum(gF,2)
      print *, 'Potencial: ', Pot


    end subroutine fuerza 

    !===============================================================================
    ! FINALIZA PARAMETROS
    !===============================================================================
    subroutine finalizacion()

      call fin_zig()

    endsubroutine finalizacion


end module dinmods 
