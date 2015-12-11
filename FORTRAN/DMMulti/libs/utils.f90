module utils
  
  use types,        only : dp
  use globales     
  use estadistica,  only: vector_mean_var_dir

! Si se utiliza openmp
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  private

  public :: write_array3D_lin, wtime, hace_estadistica
#ifdef _OPENMP
  public ::  init_openmp
!  public :: write_array3D_lin, wtime, init_openmp
!#else
!  public :: write_array3D_lin, wtime
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

#ifdef _OPENMP
  !================================================================================
  ! INICIALIZA PARAMETROS DE OPENMP
  !================================================================================
  ! Subrutina para inicializar parametros en el caso de haber compilado con OpenMP

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

  !===============================================================================
  ! ESTADISTICA SOBRE VECTORES Y ESCRITURA DE RESULTADOS  
  !===============================================================================
  ! Se calculan valores medios y desviacions estándares
  ! Se escriben los resultados al archivo <val_medios.dat> 

  subroutine hace_estadistica(x,y,z)

    real(dp), dimension(:), intent(in)      :: x, y  ! Datos en 1D para hacer estadística
    real(dp), dimension(:,:), intent(in)    :: z     ! Datos en 2D para hacer estadística

    real(dp)                         :: x_mean, y_mean ! Valor medio de los datos
    real(dp)                         :: x_std, y_std   ! Desviación estandar de los datos
    real(dp), dimension(1:size(z,1)) :: z_mean, z_std  ! para 2D
    integer                          :: i

    ! Calcula valores medios para el primer parámetro (presion) 
    write(*,*) '* Se hace estadistica con los valores de presión'
    call vector_mean_var_dir(x_mean,x_std,x)
    x_std = sqrt(x_std)

    ! Calcula valores medios para el segundo parámetro (temperatura) 
    write(*,*) '* Se hace estadistica con los valores de temperatura'
    call vector_mean_var_dir(y_mean,y_std,y)
    y_std = sqrt(y_std)
    
    ! Parche para calcular estadistica por columnas en una matriz
    ! Hacerlo más prolijo en algún momento
    write(*,*) '* Se hace estadistica con los valores de energías'
    z_mean = sum(z,2) / size(z,2) 

    do i = 1, size(z,1)
      z_std(i) =  sum( (z(i,:)-z_mean(i))**2 ) / (size(z,2) -1) 
    end do
    z_std = sqrt(z_std)
    ! Fin parchado 

    !-------------------------------------------------------------------------------------
    ! ESCRITURA DE DATOS
    !------------------------------------------------------------------------------------

    ! Se guardan los resultados en un archivo
    write(*,*) '* Se guardan valores medios y desv. estandar en <val_medios.dat>'

    open(unit=10, file='val_medios.dat', status='unknown')
    ! Primera parte con datos de la corrida   
    write(10,100)   'Temp', '#_part', 'Lado_cubo', 'Dens', 'dt', '#_pasos',  '#_pasos_med.'
    100 format (2X,7(A,8X))
    write(10,200) gTemperatura, gNpart, gLado_caja, float(gPeriodos), &
                   gDt, float(gNtime), float(gNmed)
    200 format (F6.2,10X,I6,9X,F11.5,5X,F11.5,5X,3(E9.2,5X))
    write(10,*)

    ! Segunda linea con los valores medios y desviaciones estándars
    write(10,300) 'mean(p)', 'std(p)' , 'mean(T)', 'std(T)'
    300 format (4X,4(A,12X))
    write(10,400) x_mean, x_std, y_mean, y_std
    400 format (4(E14.7,4X))
    write(10,*)

    ! Tercera linea con los valores medios y desviaciones estandars
    write(10,500) 'mean(E_pot)', 'std(E_pot)', 'mean(E_kin)', 'std(E_kin)', 'mean(E_tot)', 'std(E_kin)'
    500 format (2X,6(A,8X))
    write(10,600) (  (/z_mean(i), z_std(i)/)  , i=1,3)
    600 format (6(E14.7,4X))

    close(10)

  end subroutine hace_estadistica



end module utils
