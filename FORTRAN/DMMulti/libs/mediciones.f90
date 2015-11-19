module mediciones
  use types,      only: dp
  use globales
  !===============================================================================
  ! CALCULA ENERGIA CINETICA 
  !===============================================================================
  ! Calcula la anergia cinetica total del sistema
contains
  subroutine calcula_kin()
    
    real(dp), allocatable   :: v2(:)    ! Vector con la velocidad cuadratica
    integer :: i

    allocate(v2(1:gNpart))

    do i = 1, gNpart
      v2(i) =  dot_product(gV(:,i),gV(:,i))  
    end do

    ! TODO: calcular de acuerdo a la masa de las partículas
    gKin = 0.5_dp * gM * sum( v2 )

  endsubroutine calcula_kin

end module mediciones
