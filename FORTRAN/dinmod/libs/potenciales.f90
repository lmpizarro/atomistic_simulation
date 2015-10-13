module potenciales 
  use types, only: dp
  use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV, gSigma, gEpsil

  implicit none

  private
  public :: poten_lj, delta_poten_lj

contains

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
          sr6 = sr2 **3
          sr12 = sr6 ** 2
          v = v + sr12 - sr6
        enddo
      enddo
      v = 4.0 * gEpsil * v
  endfunction poten_lj


  ! i: particula i-esima
  ! rx, ry, rz: vieja posicion de la particula
  function delta_poten_lj (i, rxi, ryi, rzi) result(v)
    real(dp) :: v, vo, vn, rxi, ryi, rzi
    real(dp) :: rxij, ryij, rzij, rijsq
    real(dp) :: sr2, sr6, sr12
    integer :: i, j

    vo = 0

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
      rxij = gR(1, j)
      ryij = gR(2, j)
      rzij = gR(3, j)

      rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
      sr2 = gSigma / rijsq
      sr6 = sr2 **3
      sr12 = sr6 ** 2
      vo = vo + sr12 - sr6
    enddo

    vn = 0
    rxi = gR(1, i)
    ryi = gR(2, i)
    rzi = gR(3, i)

    do j = 1, i - 1
      rxij = rxi - gR(1, j)
      ryij = rxi - gR(2, j)
      rzij = rxi - gR(3, j)

      rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
      sr2 = gSigma / rijsq
      sr6 = sr2 **3
      sr12 = sr6 ** 2
      vn = vn + sr12 - sr6

    enddo

    do j = i + 1, gNpart 
      rxij = gR(1, j)
      ryij = gR(2, j)
      rzij = gR(3, j)

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
