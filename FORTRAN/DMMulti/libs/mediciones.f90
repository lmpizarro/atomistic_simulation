module mediciones
  use types,      only: dp
  use globales
  !===============================================================================
  ! CALCULA ENERGIA CINETICA 
  !===============================================================================
  ! Calcula la anergia cinetica total del sistema
  implicit none

  integer :: i, j, ke
  integer :: inic, fin
contains

  subroutine calcula_fuerza()
    integer :: inic_i, fin_i      
    ! calcula las fuerzas de interacción 
    ! de elementos iguales

    ! loop sobre los distintos elementos
    inic = 1
    do ke=1, gNespecies
      fin = gNp(i) + inic
      do i= inic, fin - 1
        do j= inic + 1, fin
        !! calcula la fuerza
        enddo
      enddo
      inic = fin + 1
    enddo


    ! calculo de las fuerzas de interacción de  
    ! elementos distintos
    inic = 1
    do ke=1, gNespecies - 1
      fin = gNp(ke) + inic - 1
      print *, inic, fin
      inic_i = fin + 1 
      do j = ke, gNespecies - 1
        fin_i = inic_i + gNp(j + 1) - 1
        print *, "loop interno:", inic_i, fin_i
        inic_i = fin_i + 1
      enddo
      inic = fin + 1
    enddo
  endsubroutine calcula_fuerza

  subroutine calcula_kin()
    
    real(dp), allocatable   :: v2(:)    ! Vector con la velocidad cuadratica
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
      inic = fin + 1
    enddo

    gKin = 0.5_dp * gKin

  endsubroutine calcula_kin

end module mediciones
