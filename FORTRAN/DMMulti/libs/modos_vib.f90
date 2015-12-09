module FFTW3
  use, intrinsic :: iso_c_binding
  include '/usr/include/fftw3.f03'
  ! También podría estar en '/usr/local/include/fftw3.f03
end module

module modos_vib

#include "control.h"
  use FFTW3      

  ! forma de llamar fftw3
  ! ver: http://fftw.org/doc/Calling-FFTW-from-Modern-Fortran.html
  !type(C_PTR) :: plan
  !complex(C_DOUBLE_COMPLEX), dimension(1024,1000) :: in, out
  !plan = fftw_plan_dft_2d(1000,1024, in,out, FFTW_FORWARD,FFTW_ESTIMATE)

  !call fftw_execute_dft(plan, in, out)


  use globales,   only: gKmed, gV, gF, gR, gNpart, gKin, gNespecies,&
                        gCorrVfac_1, gCorrVfac_2, gCorrVfac_3, gCorrVver_1,&
                        gNCorrVfac_1, gNCorrVfac_2, gNCorrVfac_3, gNCorrVver_1,&
                        gNmodosVibra, gRho, gVir, gVol, gGamma, gDt,&
                        gPot, gTemperatura, gLado_caja, &
                        gIndice_elemento, gMasa, gPot_cut, gRc2, gCombSigma,&
                        gCombEpsilon


  use  io_parametros, only: escribe_en_columnas 
  use types
  implicit none

  private

  public :: lee_velocidades, modos_posicion,&
            calcula_modos_vibracion, calcula_modos_vibracion_vel, &
            calcula_autocorr_vel_3D, calcula_autocorr_vel_3D_1D
contains


  subroutine modos_posicion(v, nombre)
    real(dp), allocatable ::  freqs(:,:) 
    real(dp), dimension(:,:), intent(in) :: v
    character(*), intent(in)                     :: nombre  ! Nombre del archivo
    integer :: n

    n=size(v(1,:))

  print *, "calculo de modos de vibracion ", n 
    allocate (freqs(1:3, 1:(n/2) + 1 ))
    call calcula_modos_vibracion_vel(v, freqs)
    call escribe_en_columnas(freqs,nombre,n*1.0_dp)
    ! grabar resultados que estan en Corr_tr
    deallocate (freqs)
  end subroutine modos_posicion


  subroutine lee_velocidades(nombre, v,exito)
    real(dp), intent(out), allocatable ::  v(:,:)  ! Vectores pos. y vel.
    character(*), intent(in)                     :: nombre  ! Nombre del archivo
    logical, intent(out)                       :: exito ! Controla si se leyo bien
    integer          :: N     ! Número de velocidades leidas    
    integer          :: i, j
    logical          :: es

    exito = .FALSE.

    inquire(file=nombre, exist=es)
    if (es) then
      open(unit=20, file=nombre, status='unknown')
      read(20,'(I12)') N
      print *, "Cantidad de velocidades a leer ", N, " del archivo ", nombre
      allocate (v(1:3,1:N))
        do i = 1, N
            read(20,100) (v(j,i),j=1,3)
            100 format (3(E16.9,2X))
        end do 
        write(*,*) '* Archivo de velocidades  leido correctamente'
        exito = .TRUE.
      close(20)
    else
      ! No existe el archivo
      write(*,*) '* No existe archivo de estados iniciales'
    end if

  endsubroutine lee_velocidades

  !
  ! Viene un vector 3D de N, devuelve un vector c de N/2 + 1
  !
  subroutine calcula_modos_vibracion_vel (v, c)
    real(dp), dimension(:,:), intent(in) :: v
    real(dp), dimension(:,:), intent(inout) :: c
    complex(dp), dimension(:), allocatable :: tmp
    integer :: n, nc
    real(kind=8), dimension(:), allocatable :: in

    n = size(c(1,:))
    nc = size(v(1,:))


    allocate( tmp(1:n ))
    allocate(in(1:nc))

    print *,  "size 1 ", size(v(1,:)), size(c(1,:))
    in = v(1,:)
    call calcula_modos_vibracion (in, tmp(1:n))
    c(1,1:n) = REAL(tmp(1:n)) ** 2 + AIMAG(tmp(1:n)) ** 2
    in = v(2,:)
    call calcula_modos_vibracion (in, tmp(1:n))
    c(2,1:n) = REAL(tmp(1:n)) ** 2 + AIMAG(tmp(1:n)) ** 2
    in = v(3,:)
    call calcula_modos_vibracion (in, tmp(1:n))
    c(3,1:n) = REAL(tmp(1:n)) ** 2 + AIMAG(tmp(1:n)) ** 2

    deallocate (tmp)
    deallocate(in)
  endsubroutine calcula_modos_vibracion_vel 

  !
  ! Viene un vector 1 D (in) de N, 
  ! devuelve un vector 1D (c) de N/2 + 1
  !
  subroutine calcula_modos_vibracion (in, out)
    real(dp), dimension(:), intent(in) :: in
    complex(dp), dimension(:), intent(inout) :: out
    integer :: n, nc
    integer ( kind = 8 ) plan_forward

    n = size(in)
    nc = n / 2  + 1 
    write (*,'(a)') " calculo de modos de vibracion 1D " 
    write (*,'(a , I8, a, I8)') "size in ", n, "size out ", nc

    ! out tiene tamaño nc
    call dfftw_plan_dft_r2c_1d_ (plan_forward, n, in, out, FFTW_ESTIMATE)
    call dfftw_execute_ (plan_forward)
 
  endsubroutine calcula_modos_vibracion

  !
  ! calcula la auto_correlacion de un vector 1D
  ! con fftw3
  ! el tamaño de vel es del tipo 2 ** n, n entero
  ! los ultimos n elementos en 0
  ! de acuerdo a Allen
  subroutine calcula_autocorr (v, c)
    real(dp), dimension(:), intent(in) :: v
    real(dp), dimension(:), intent(inout) :: c
    integer :: n, nc, i
    complex (kind=8), allocatable :: out_forwd(:), in_backw(:)
    integer ( kind = 8 ) plan_forward
    integer ( kind = 8 ) plan_backward
    real(kind = 8), allocatable :: in(:)
    real(kind = 8), allocatable :: out_backw(:)

    n = size(v)
    nc = n / 2 + 1

    allocate (out_forwd(1:nc))
    allocate (out_backw(1:n))
    allocate (in_backw(1:nc))
    allocate(in(1:n))

    in(1:n) = v(1:n)

    print *, "v out_backw salida out_forwd ", n, nc, size(c), size(out_forwd)

    ! out tiene tamaño nc
    call dfftw_plan_dft_r2c_1d_ (plan_forward, n, in, out_forwd, FFTW_ESTIMATE)

    call dfftw_execute_ (plan_forward)
    

    ! paso c) página 190 del Allen
    !in_backw = REAL(out_forwd) ** 2 + AIMAG(out_forwd) ** 2
    in_backw = out_forwd * conjg(out_forwd)

    ! paso d) aplicar una fft inversa  a c
    call dfftw_plan_dft_c2r_1d_ ( plan_backward, n, in_backw, out_backw, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_backward )


    c = out_backw(1:n/2) / (DBLE(SIZE(in)))

    DO i=2,nc +1
      c(i-1)=c(i-1)/DBLE(nc-(i-1))
    END DO

    deallocate (out_forwd)
    deallocate (in_backw)

  endsubroutine calcula_autocorr

  ! calcula la correlacion de velocidad de un vector temporal 3D
  subroutine calcula_autocorr_vel_3D (v, corr)
    real(dp), dimension(:,:), intent(in) :: v
    real(dp), allocatable, intent(inout) :: corr(:,:)
    real(dp), allocatable :: v2(:)
    integer :: n

    ! viene el resultado de la operacion
    n =  size(v(1,:))
    allocate (v2(1: 2 * n))
    v2 = 0
    corr = 0

    v2(1:n) = v(1,:)
    call calcula_autocorr (v2, corr(1,:))
    v2 = 0
    v2(1:n) = v(2,:)
    call calcula_autocorr (v2, corr(2,:))
    v2 = 0
    v2(1:n) = v(3,:)
    call calcula_autocorr (v2, corr(3,:))
    
 
    deallocate (v2)
  endsubroutine calcula_autocorr_vel_3D 

  ! calcula la correlacion de velocidad de un vector temporal 3D
  subroutine calcula_autocorr_vel_3D_1D (v, corr)
    real(dp), dimension(:,:), intent(in) :: v
    real(dp), allocatable, intent(inout) :: corr(:)
    real(dp), allocatable :: v2(:)
    integer :: n

    ! viene el resultado de la operacion
    n =  size(v(1,:))
    allocate (v2(1: 2 * n))
    v2 = 0
    corr = 0

    !v2(1:n) = sqrt(v(1,1:n) ** 2 + v(2,1:n) ** 2 + v(3,1:n) ** 2)
    v2(1:n) = v(1,1:n) + v(2,1:n) + v(3,1:n) 
    call calcula_autocorr (v2, corr(:))

    
 
    deallocate (v2)
  endsubroutine calcula_autocorr_vel_3D_1D 




#ifdef MODOS_VIB_EQUIVALENTES   
  subroutine calcula_autocorr_vel_3D_b (v, c)
    ! Se puede pensar de otra manera
    ! acumular en un vector 1D cada una de las 3dimensiones
    ! del array de velocidades
    real(dp), dimension(:,:), intent(in) :: v
    real(dp), dimension(:), intent(inout) :: c
    real(dp), allocatable :: corr(:)
    real(dp), allocatable :: v2(:)

    ! viene el resultado de la operacion
    allocate (corr(2 * size(v)))
    allocate (v2(2 * size(v)))
    v2 = 0
    corr = 0

    v2(1:size(v)) = (v(1,:) + v(2,:) + v(3,:)) / 3
    call calcula_autocorr (v2(:), corr(:))

    deallocate (corr)
    deallocate (v2)

  endsubroutine calcula_autocorr_vel_3D_b 

#endif

endmodule modos_vib
