module utils
  
  use types, only :dp
  use globales, only : gR

! Si se utiliza openmp
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  private

  public :: write_array3D_lin, wtime
#ifdef _OPENMP
  init_openmp
#endif

contains

  !================================================================================
  ! imprime con formato por pantalla un array 3D en lineal 
  !================================================================================

  subroutine write_array3D_lin (b)

    integer  :: i, n, m
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

    700 format (F7.3 , F7.3 , F7.3)

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
  end function

  !================================================================================
  ! INICIALIZA PARAMETROS DE OPENMP
  !================================================================================
  ! Subrutina para inicializar parametros en el caso de haber compilado con OpenMP

#ifdef _OPENMP
  subroutine init_openmp()

    integer     :: num_proc         ! Cantidad de procesadores disponibles
    integer     :: num_threads      ! Cantidad de threads disponibles (especificados)

    num_proc    = omp_get_num_procs()
    num_threads = omp_get_max_threads()

    write(*,'(a)') ''
    write(*,'(a)')      '************ PARAMETROS OPENMP **************'
    write(*,'(a,I8)')   '********* Procesadores disponibles = ' , num_proc 
    write(*,'(a,I8)')   '********* Threads disponibles      = ' , num_threads
    write(*,'(a)')      '*********************************************'

  end subroutine init_openmp
#endif

end module utils
