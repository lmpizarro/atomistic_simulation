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
    integer :: i, inic, fin
    real(dp) :: masa


    allocate(v2(1:gNpart))

    do i = 1, gNpart
      v2(i) =  dot_product(gV(:,i),gV(:,i))  
    end do

    gKin = 0.0_dp
    masa = 0.0_dp 
    inic = 1
    do i=1, gNespecies
      fin = gNp(i) + inic
      masa = gLj_param(i, 3)
      gKin = gKin + masa * sum(v2(inic:fin))
    enddo

    gKin = 0.5_dp * gKin

  endsubroutine calcula_kin

end module mediciones
