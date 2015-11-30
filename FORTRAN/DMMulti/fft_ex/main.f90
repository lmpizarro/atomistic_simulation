! compilar: gfortran -c main.f90
! linkear: gfortran -o main main.o -lfftw3
! 
! genera main ejecutable que corre en linux debian
!
!

module FFTW3
  use, intrinsic :: iso_c_binding
  include '/usr/include/fftw3.f03'
end module

program main
  use FFTW3      
  implicit none

endprogram main


subroutine fft_array_real (arr)
endsubroutine fft_array_real 

subroutine test02 ()
  use FFTW3      
  implicit none

  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ), parameter :: nc = 501
  real ( kind = 8 ) in(n), in2(n)
  real ( kind = 8 ) f, dt
  integer ( kind = 4 ) i
  integer ( kind = 8 ) plan_forward
  complex ( kind = 8 ) out(nc)
  real ( kind = 8 ) out_real(nc), out_img(nc) 


  f = 1.0
  dt = .1
  call  gen_sin (n, f, dt, in)
  f = 0.5
  call  gen_sin (n, f, dt, in2)

  in = in + in2



  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input Data:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    !write ( *, '(i4, 2g14.6)' ) i, in(i)
  end do

  call dfftw_plan_dft_r2c_1d_ ( plan_forward, n, in, out, FFTW_ESTIMATE )



  call dfftw_execute_ ( plan_forward )

  !
  ! para determinar la fase atan2(aimag(z),realpart(z))
  !
  out_real = REALPART (out)
  out_img = AIMAG (out)

  out = out_img ** 2  + out_real ** 2


  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output FFT Coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, nc
    write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)
  end do



endsubroutine test02

subroutine test01 ()
  use FFTW3      
  implicit none

  type(C_PTR) :: plan
  complex(C_DOUBLE_COMPLEX), dimension(1024,1000) :: in, out

  plan = fftw_plan_dft_2d(1000,1024, in,out, FFTW_FORWARD,FFTW_ESTIMATE)
  call fftw_execute_dft(plan, in, out)


endsubroutine test01
!
! gen_sin: genera una onda senoidal
! n: cantidad de elementos
! f: frecuencia
! dt: delta de tiempo 
! r: array
subroutine gen_sin (n, f, dt, r)
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) i

  real ( kind = 8 ) r(n)
  real ( kind = 8) dt, f
  real ( kind = 8) pi

  pi = 4.0 * atan(1.0)
 
  do i =1, n
    r(i) = sin(2*PI*f*(i -1)* dt) 
  enddo

endsubroutine gen_sin 
