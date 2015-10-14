module potenciales 
  use types, only: dp
  use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil

  implicit none

  private
  public :: poten_lj, delta_poten_lj, poten_lj_vec, delta_poten_lj_vec

  !===============================================================================
  ! VARIABLE PROPIAS DEL MÓDULO
  !===============================================================================

  real(dp)        :: rc2       ! Cuadrado del radio de corte
  real(dp)        :: pot_cut   ! Potencial L-J en el radio de corte
 
contains

  !===============================================================================
  ! RADIO DE CORTE 
  !===============================================================================
  ! Define el radio de corte y calcula el desplazamiento del potencial de L-J 

  subroutine corta_desplaza_pote()   

    ! Radio de corte fijo
    rc2 = (2.5_dp)**2                
    ! Potencial de L-J evaluado en el radio de corte
    pot_cut = 4.0_dp*gEpsil*( 1.0_dp/(rc2**6) - 1.0_dp/(rc2**3) ) 
    
  end subroutine corta_desplaza_pote
 
  ! calculo del potencial de pag 18 del allen-tildesley
  function poten_lj () result(v)
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
          sr2 = gSigma / rijsq
          sr6 = sr2 ** 3
          sr12 = sr6 ** 2
          v = v + sr12 - sr6
        enddo
      enddo
      v = 4.0 * gEpsil * v
  endfunction poten_lj


  function poten_lj_vec () result(Pot)
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    real(dp)                :: Pot
    integer                 :: i,j       

    call corta_desplaza_pote()

    ! Paso a trabajar distancias en unidades de sigma
    gR = gR/gSigma
    gL = gL/gSigma

    Pot = 0.0_dp

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
          Pot     = Pot + r6in * ( r6in - 1.0_dp)    ! Energía potencial
        end if
      enddo
    enddo

    ! Constantes que faltaban en el potencial
    ! Agrego el desplazamiento del potencial considerando la cantidad de
    ! pares con que se obtuvo la energía potencial N(N-1)/2
    Pot =  4.0_dp * gEpsil * Pot - gNpart*(gNpart-1)*pot_cut/2.0_dp  

    ! Se vuelven a pasar a las coordenadas absolutas
    gR = gR*gSigma
    gL = gL*gSigma

  endfunction poten_lj_vec

  function kernel_pot (j, ri_vec) result (Pot)
    real(dp), intent(in), dimension(3)  :: ri_vec
    integer, intent(in) :: j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    real(dp)                :: Pot
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j

    Pot = 0

    rij_vec = ri_vec - gR(:,j)               ! Distancia vectorial
    ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
    ! Siempre en distancias relativas de sigma
    rij_vec = rij_vec - gL*nint(rij_vec/gL)
    r2ij   = dot_product( rij_vec , rij_vec )    ! Cuadrado de la distancia
    if ( r2ij < rc2 ) then               
       r2in = 1.0_dp/r2ij                         ! Inversa al cuadrado
       r6in = r2in**3                             ! Inversa a la sexta
       Pot = r6in * ( r6in - 1.0_dp)    ! Energía potencial
    end if


  endfunction kernel_pot


  function delta_poten_lj_vec  (i, ri_vec) result(deltaPot)
    ! la particula con la que se hace le calculo
    integer, intent(in) :: i    
    ! la posición original de la partícula 
    real(dp), intent(in), dimension(3)  :: ri_vec   
    integer :: j
    real(dp) :: deltaPot
    real(dp) :: vo, vn

    deltaPot = 0.0_dp

    vo = 0.0_dp
    vn = 0.0_dp
    do j = 1, i - 1
      vo     = vo + kernel_pot(j, ri_vec)   ! Energía potencial
      vn     = vn + kernel_pot(j, gR(:,i))   ! Energía potencial
    enddo

    do j = i + 1, gNpart 
      vo     = vo + kernel_pot(j, ri_vec)   ! Energía potencial
      vn     = vn + kernel_pot(j, gR(:,i))   ! Energía potencial
    enddo


    deltaPot =  4.0_dp * gEpsil * (vn - vo)   

  endfunction delta_poten_lj_vec


  ! i: particula i-esima
  ! rx, ry, rz: vieja posicion de la particula
  function delta_poten_lj (i, rxi, ryi, rzi) result(v)
    real(dp) :: v, vo, vn, rxi, ryi, rzi
    real(dp) :: rxij, ryij, rzij, rijsq
    real(dp) :: sr2, sr6, sr12
    integer :: i, j

    vo = 0.0_dp

    do j = 1, i - 1
      rxij = rxi - gR(1, j)
      ryij = rxi - gR(2, j)
      rzij = rxi - gR(3, j)

      rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
      sr2 = gSigma / rijsq
      sr6 = sr2 ** 3
      sr12 = sr6 ** 2
      vo = vo + sr12 - sr6

    enddo

    do j = i + 1, gNpart 
      rxij = rxi - gR(1, j)
      ryij = ryi -  gR(2, j)
      rzij = rzi -  gR(3, j)

      rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
      sr2 = gSigma / rijsq
      sr6 = sr2 **3
      sr12 = sr6 ** 2
      vo = vo + sr12 - sr6
    enddo

    vn = 0.0_dp
    rxi = gR(1, i)
    ryi = gR(2, i)
    rzi = gR(3, i)

    do j = 1, i - 1
      rxij = rxi - gR(1, j)
      ryij = ryi - gR(2, j)
      rzij = rzi - gR(3, j)

      rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
      sr2 = gSigma / rijsq
      sr6 = sr2 **3
      sr12 = sr6 ** 2
      vn = vn + sr12 - sr6

    enddo

    do j = i + 1, gNpart 
      rxij = rxi - gR(1, j)
      ryij = ryi -  gR(2, j)
      rzij = rzi -  gR(3, j)

      rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
      sr2 = gSigma / rijsq
      sr6 = sr2 **3
      sr12 = sr6 ** 2
      vn = vn + sr12 - sr6
    enddo

    v = vn - vo

  endfunction delta_poten_lj

  function cuadrado () result(v)
    real(dp) :: v
  endfunction cuadrado

  function esfera_rigida () result(v)
    real(dp) :: v
  endfunction esfera_rigida

end module potenciales 
