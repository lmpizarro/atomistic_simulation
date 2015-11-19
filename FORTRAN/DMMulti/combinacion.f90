!
!
!
!
!
module combinacion
  use types, only: dp
  use globales

  implicit none

contains

  subroutine comb_epsilon()
    ALLOCATE (gCombEpsilon(1:gNespecies,1:gNespecies))
  endsubroutine comb_epsilon

  subroutine comb_sigma()
    ALLOCATE (gCombSigma(1:gNespecies,1:gNespecies))
  endsubroutine comb_sigma

endmodule combinacion
