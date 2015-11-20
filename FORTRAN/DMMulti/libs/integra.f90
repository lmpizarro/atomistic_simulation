module integra
  use types,          only: dp
  use globales

  implicit none
contains 

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================

  subroutine cpc_vec()

    real(dp), dimension(3,gNpart) :: temp     ! Variable temporal
    
    !gR = gR - gL*floor(gR/gL)
    
    ! Lo escribo de esta forma porque de lo contrario da error de compilación
    ! en el cluster (dice ser un bug de gfortran)
    temp = gLado_caja*floor(gR/gLado_caja)
    gR = gR - temp

  end subroutine

end module integra
