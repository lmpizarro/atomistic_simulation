module utils
  
  use types, only :dp
  use globales, only : gR

  implicit none

  private

  public :: write_pos, write_array3D_lin, wtime

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

    integer  :: i, j, k, n, m, l
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

  !================================================================================
  ! FUNCION PARA MEDIR EL WALL-TIME
  !================================================================================

  function wtime ( )

  !*****************************************************************************80
  !           https://people.sc.fsu.edu/~jburkardt/f_src/wtime/wtime.f90
  !! WTIME returns a reading of the wall clock time.
  !
  !  Discussion:
  !
  !    To get the elapsed wall clock time, call WTIME before and after a given
  !    operation, and subtract the first reading from the second.
  !
  !    This function is meant to suggest the similar routines:
  !
  !      "omp_get_wtime ( )" in OpenMP,
  !      "MPI_Wtime ( )" in MPI,
  !      and "tic" and "toc" in MATLAB.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    27 April 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, real ( kind = 8 ) WTIME, the wall clock reading, in seconds.
  !
    implicit none

    integer ( kind = 4 ) clock_max
    integer ( kind = 4 ) clock_rate
    integer ( kind = 4 ) clock_reading
    real ( kind = dp ) wtime

    call system_clock ( clock_reading, clock_rate, clock_max )

    wtime = real ( clock_reading, kind = dp ) &
          / real ( clock_rate, kind = dp )

    return
  end


end module utils
