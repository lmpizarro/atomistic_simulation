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

    ! imprime la matriz por pantalla
    print *, "gSigma"
    do i = 1, gNespecies
      write (*, 700) (gCombSigma (i,j), j = 1,gNespecies)
    enddo

    700 format (100g15.5)


  endsubroutine comb_sigma

endmodule combinacion
