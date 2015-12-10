module matrices

  implicit none

  private
  
  public  :: iepsilon 

  real (kind=8), allocatable :: gEpsilon(:,:) 

contains

  !
  !
  subroutine iepsilon (epsilon, gNelements)
    integer, intent(in) :: gNelements       
    real (kind=8),intent(in), dimension(0:gNelements - 1 ,0: gNelements - 1 ):: epsilon
!f2py intent(in) :: epsilon 
    integer :: i, j

    allocate (gEpsilon(1:gNelements,1:gNelements ))

    do i=0, gNelements - 1
      do j=0, gNelements - 1
         gEpsilon(i+1, j + 1)  = epsilon(i,j)
      enddo  
    enddo  
    call lepsilon()
  endsubroutine iepsilon
  !
  !
  subroutine lepsilon ()
    integer :: i, j, size_

    size_ = size(gEpsilon, 1)
  
    do i=1, size_ 
      do j=1, size_ 
        write(*, 100, advance="no")   gEpsilon(i, j)
      enddo
        write(*, *) ""
    enddo

    100 format (F10.3)
  endsubroutine lepsilon
 
endmodule matrices
