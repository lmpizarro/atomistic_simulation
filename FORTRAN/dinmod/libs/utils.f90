module utils
    use types, only :dp
    use globales, only : gR

    implicit none

    private

    public :: write_pos, write_array3D_lin

    contains

    !
    ! Para grabar las posiciones ?? 
    !
    subroutine write_pos()
    endsubroutine


    !================================================================================
    ! imprime con formato por pantalla un array 3D en lineal 
    !================================================================================
    subroutine write_array3D_lin (b)
        integer :: i, j, k, n, m, l
        real(dp), dimension(:,:), intent(in) :: b

        n = size(b,1)
        m = size(b,2)

        print *, "n m ", n, m


        !do i = 1, n
        !    write(*,700) (b(i,j), j = 1,m)
        !enddo

        do i = 1, m
            write(*,700) b(1,i), b(2,i), b(3,i)
        enddo

        700 format (F7.3 F7.3 F7.3)
    endsubroutine


end module utils
