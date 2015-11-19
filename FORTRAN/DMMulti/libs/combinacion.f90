!
!
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
    do i=1, gNespecies
      gCombEpsilon(i,i) = gLj_param(i,1)
    enddo

    do i =1, gNespecies
      do j=i+1, gNespecies
         gCombEpsilon(i,j) = gLj_param(i,1) *  gLj_param(j,1)
         gCombEpsilon(j,i) = gCombEpsilon(i,j) 
      enddo
    enddo
  endsubroutine comb_epsilon

  subroutine comb_sigma()
    allocate (gCombSigma(1:gNespecies,1:gNespecies))
    do i=1, gNespecies
      gCombSigma(i,i) = gLj_param(i,2)
    enddo
    do i =1, gNespecies
      do j=i+1, gNespecies
         gCombSigma(i,j) =  (gLj_param(i,2) + gLj_param(j,2)) / 2
         gCombSigma(j,i) = gCombSigma(i,j)
      enddo
    enddo

  endsubroutine comb_sigma

endmodule combinacion
