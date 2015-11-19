! subroutinas que calculan los parámetros
! de interacción de la mezcla
!
!
!
module combinacion
  use types, only: dp
  use globales

  implicit none

  integer :: i, j

contains

  subroutine comb_epsilon()
    allocate (gCombEpsilon(1:gNespecies,1:gNespecies))

    ! asigna la diagonal
    do i=1, gNespecies
      gCombEpsilon(i,i) = gLj_param(i,1)
    enddo

    ! calcula los elementos fuera de la diagonal
    do i =1, gNespecies
      do j=i+1, gNespecies
         gCombEpsilon(i,j) = sqrt(gLj_param(i,1) *  gLj_param(j,1))
         gCombEpsilon(j,i) = gCombEpsilon(i,j) 
      enddo
    enddo
   
    gCombEpsilon = gCombEpsilon / gLj_param (1,1)

    ! imprime la matriz por pantalla
    print *, "gEpsilon"
    do i = 1, gNespecies
      write (*, 700) (gCombEpsilon (i,j), j = 1,gNespecies)
    enddo

    700 format (100g15.5)
  endsubroutine comb_epsilon

  subroutine comb_sigma()
    allocate (gCombSigma(1:gNespecies,1:gNespecies))

    ! asigna la diagonal
    do i=1, gNespecies
      gCombSigma(i,i) = gLj_param(i,2)
    enddo

    ! calcula los elementos fuera de la diagonal
    do i =1, gNespecies
      do j=i+1, gNespecies
         gCombSigma(i,j) =  (gLj_param(i,2) + gLj_param(j,2)) / 2
         gCombSigma(j,i) = gCombSigma(i,j)
      enddo
    enddo
 
    gCombSigma = gCombSigma / gLj_param(1,2)

    ! imprime la matriz por pantalla
    print *, "gSigma"
    do i = 1, gNespecies
      write (*, 700) (gCombSigma (i,j), j = 1,gNespecies)
    enddo

    ! como se calcularon los sigmas
    ! calculamos los rc y potencial en rc
    call corta_desplaza_pote()

    700 format (100g15.5)
  endsubroutine comb_sigma

  ! calcula el radio de corte y el potencial
  ! en ese lugar
  subroutine corta_desplaza_pote ()
    allocate (gRc2(1:gNespecies))
    allocate (gPot_Cut(1:gNespecies))

    ! calcula para cada una de las especies 
    ! el radio de corte y el potencial de corte
    ! de acuerdo al cálculo en dinmod
    print *, "Rc**2       V(Rc)"
    do i = 1, gNespecies
       gRc2(i) = (4 * gCombSigma(i,i) ) ** 2
       gPot_Cut(i) = 4.0_dp * ( 1.0_dp / (gRc2(i)**6) - 1.0_dp / (gRc2(i)**3) ) 
       write (*, 700) gRc2(i), gPot_Cut(i)
    enddo
    700 format (F10.3, E15.3)
  endsubroutine corta_desplaza_pote

endmodule combinacion
