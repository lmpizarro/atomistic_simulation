module dinmods 

  use types,      only: dp
  use globales,   only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil, gM
  use utils,      only: write_array3D_lin
  use ziggurat
  use usozig
  use io_parametros,  only: escribe_trayectoria

  implicit none

  private


  public :: inicializacion, inicia_posicion_cs, finalizacion, cpc, fuerza, &
            integracion_min, inicia_posicion_rn
  
  !===============================================================================
  ! VARIABLE PROPIAS DEL MÓDULO
  !===============================================================================

  real(dp)        :: Pot     ! Energía potencial del sistema
  real(dp)        :: rc2       ! Cuadrado del radio de corte
  real(dp)        :: pot_cut   ! Potencial L-J en el radio de corte
  
contains

  !===============================================================================
  ! INICIALIZA SISTEMA de CALCULO
  !===============================================================================
  subroutine inicializacion()

    allocate(gR(3,gNpart))
    allocate(gV(3,gNpart))
    allocate(gF(3,gNpart))

    ! Inicial generador de número aleatorios
    call inic_zig()
    ! Define el radio de corte y el potencial desplazado
    call corta_desplaza_pote()

  end subroutine inicializacion 

  !===============================================================================
  ! RADIO DE CORTE 
  !===============================================================================
 
  subroutine corta_desplaza_pote()   

    ! Radio de corte fijo
    rc2 = (2.5_dp)**2                
    ! Potencial de L-J evaluado en el radio de corte
    pot_cut = 4.0_dp*gEpsil*( 1.0_dp/(rc2**6) - 1.0_dp/(rc2**3) ) 
    
  end subroutine corta_desplaza_pote

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================

  subroutine cpc(l)
  
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

  ! Verson vectorizada

  subroutine cpc_vec()

    gR = gR - gL*floor(gR/gL)

  end subroutine

  !===============================================================================
  ! INICIALIZA Posicion en Red periódica cúbica simple
  !===============================================================================

  subroutine inicia_posicion_cs()

    integer  :: i, j, k, l, nl
    real(dp) :: rx, ry, rz

    nl = gNpart ** (1.0_dp/3.0_dp) 

    rx = 0.0_dp
    ry = 0.0_dp
    rz = 0.0_dp
    j = 1
    i = 1
    k = 1
    do l = 1, gNpart
      gR(1, l) = (rx + uni() - 0.5_dp) * gL / nl 
      gR(2, l) = (ry + uni() - 0.5_dp) * gL / nl 
      gR(3, l) = (rz + uni() - 0.5_dp) * gL / nl 

      call cpc(l)
           
      j = j + 1
      k = k + 1
      rx = rx + 1.0_dp
      if ( mod(j, nl) .eq. 0) then
        rx = 0.0_dp
        ry = ry  + 1.0_dp
      end if        
      if ( mod(k, nl ** 2) .eq. 0) then
        ry = 0.0_dp
        rx = 0.0_dp
        rz = rz + 1.0_dp
       end if        
     enddo     
  
  end subroutine inicia_posicion_cs

  !===============================================================================
  ! INICIALIZA Posicion aleatoria dentro de la caga 
  !===============================================================================

  subroutine inicia_posicion_rn()
 
    integer :: l

    do l = 1, gNpart
      gR(1, l) = uni() * gL 
      gR(2, l) = uni() * gL 
      gR(3, l) = uni() * gL 

       call cpc(l)
   enddo     
    
  end subroutine inicia_posicion_rn

  !===============================================================================
  ! CALCULA LA FUERZA  
  !===============================================================================

  subroutine fuerza()

    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: Fij       ! Módulo fuerza entre partículas i y j
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    integer                 :: i,j       

    ! Se van a acumular las fuerzas. Se comienza poniendo todas a cero.
    gF  = 0.0_dp
    Pot = 0.0_dp
    ! Paso a trabajar distancias en unidades de sigma
    gR = gR/gSigma
    gL = gL/gSigma

    do i = 1, gNpart - 1       
      do j = i+1, gNpart
        rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
        ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
        ! Siempre en distancias relativas de sigma
        rij_vec = rij_vec - gL*nint(rij_vec/gL)
        r2ij   = dot_product( rij_vec , rij_vec )    ! Cuadrado de la distancia
        if ( r2ij < rc2 ) then               
          r2in = 1.0_dp/r2ij                         ! Inversa al cuadrado
          r6in = r2in**3                             ! Inversa a la sexta
          Fij     = r2in * r6in * (r6in - 0.5_dp)    ! Fuerza entre partículas
          gF(:,i) = gF(:,i) + Fij * rij_vec          ! Contribución a la partícula i
          gF(:,j) = gF(:,j) - Fij * rij_vec          ! Contribucion a la partícula j
          Pot     = Pot + r6in * ( r6in - 1.0_dp)    ! Energía potencial
        end if
      end do
    end do

    ! Constantes que faltaban en la energía
    Fij = 48.0_dp * gEpsil * Fij                
    ! Constantes que faltaban en el potencial
    ! Agrego el desplazamiento del potencial considerando la cantidad de
    ! pares con que se obtuvo la energía potencial N(N-1)/2
    Pot =  4.0_dp * gEpsil * Pot - gNpart*(gNpart-1)*pot_cut/2.0_dp  

    ! Se vuelven a pasar a las coordenadas absolutas
    gR = gR*gSigma
    gL = gL*gSigma

!    write(*,'(A,2X,3(E15.5,3X))')  'Sumatoria de fuerzas:' , sum(gF,2)
!    print *, 'Potencial: ', Pot

  end subroutine fuerza 

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO  
  !===============================================================================

  subroutine integracion_min()
  ! Subrutina de integración de las ecuaciones de movimiento para minimizar energía
  ! Es el Problema 3 de la Guia_2a

    real(dp), dimension(gNtime+1)   :: Eng_t   ! Energía en función del tiempo
    integer    :: i

    ! El primer punto es la energía inicial
    Eng_t(1) = Pot

    call write_array3D_lin(gR)

    call cpc_vec()

    do i = 1, gNtime 
      gR = gR + 0.5_dp * gF * gDT**2 / gM
      call cpc_vec()    
     ! Esta subrutine abre y cierra un archivo. Se puede optimizar haciéndolo acá.
      call escribe_trayectoria(gR)    
      call fuerza()
      Eng_t(i+1) = Pot
    end do

    open(unit=10,file='./energia.dat',status='unknown')
    !write(10,'(F10.4)') gDt
    write(10,'(E15.5)') Eng_t
    close(10)
  
    call write_array3D_lin(gR)
    print *, 'Energía' , Pot

  end subroutine integracion_min

  !===============================================================================
  ! FINALIZA PARAMETROS
  !===============================================================================

  subroutine finalizacion()

    call fin_zig()

    deallocate(gR,gV,gF)

  endsubroutine finalizacion


end module dinmods 
