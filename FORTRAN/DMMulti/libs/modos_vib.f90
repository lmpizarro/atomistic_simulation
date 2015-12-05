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
                        gNmodosVibra, gCorr_par, gRho, gVir, gVol, gGamma, gDt,&
                        gPot, gTemperatura, gLado_caja, gDbin, gNgr, &
                        gIndice_elemento, gMasa, gPot_cut, gRc2, gCombSigma,&
                        gCombEpsilon


  use  io_parametros, only: escribe_en_columnas 
  use types
  implicit none

  private

  public :: modos_equivalentes, lee_velocidades, modos_posicion,&
            calcula_modos_vibracion, calcula_modos_vibracion_vel, &
            calcula_autocorr_vel_3D


contains

  subroutine modos_equivalentes ()
#ifdef MODOS_VIB_EQUIVALENTES   
  !#############################################################################
  !               Calculo de los modos de vibracion de velocidad
  !#############################################################################

  ! calcular autocorrelacion para particula tipo 1 en vertice
  ! gNCorrVver_2 
  print *, "gNcorr*", gNCorrVver_1, gNCorrVfac_3, gNCorrVfac_1, gNCorrVfac_2 
  if (gNCorrVfac_3 .gt. 1) then
    ! para la autocorrelacion hay que allocar 2 * (gNCorrVver_2-1)
    ! para los modos de vibracion hay que allocar gNCorrVver_2 - 1
    allocate (Modos_vibra(1:3, 1:(gNCorrVfac_3)))
    print *, "Calcula auto corre v particula tipo 2 en vertice"
    !call calcula_autocorr_vel_3D_b (gCorrVver_2(:,1:(gNCorrVver_2 - 1)), Corr_tr)
    call calcula_modos_vibracion_vel(gCorrVfac_3(:,1:(gNCorrVfac_3)),&
      Modos_vibra)
    ! grabar resultados que estan en Corr_tr
    print *, "corr_tr", Corr_tr
    call escribe_en_columnas(Modos_vibra,'corr_vver_2.dat',gNmed*gDt)
    deallocate (Modos_vibra)
  endif 

  if (gNCorrVver_1 .gt. 1) then
    ! calcular autocorrelacion para particula tipo 2 en vertice
    ! gNCorrVver_1 
    allocate (Modos_vibra(1:3, 1:(gNCorrVver_1)))
    print *, "Calcula auto corre v particula tipo 1 en vertice"
    !call calcula_autocorr_vel_3D_b (gCorrVver_1(:,1:gNCorrVver_1 - 1), Corr_tr)
    call calcula_modos_vibracion_vel(gCorrVver_1(:,1:(gNCorrVver_1)),&
      Modos_vibra)
    call escribe_en_columnas(Modos_vibra,'corr_vver_1.dat',gNmed*gDt)
    ! grabar resultados que estan en Corr_tr
    deallocate (Modos_vibra)
  endif

  if (gNCorrVfac_1 .gt. 1) then
    ! calcular autocorrelacion para particula tipo 1 en cara
    ! gNCorrVfac_1
    print *, "Calcula auto corre v particula tipo 1 en cara"
    allocate (Modos_vibra(1:3, 1:(gNCorrVfac_1)))
    !call calcula_autocorr_vel_3D_b (gCorrVfac_1(:,1:(gNCorrVfac_1 - 1)), Corr_tr)
    call calcula_modos_vibracion_vel(gCorrVfac_1(:,1:(gNCorrVfac_1)),&
      Modos_vibra)
    call escribe_en_columnas(Modos_vibra,'corr_vfac_1.dat',gNmed*gDt)
    ! grabar resultados que estan en Corr_tr
    deallocate (Modos_vibra)
  endif

  if (gNCorrVfac_2 .gt. 1) then
    ! calcular autocorrelacion para particula tipo 1 en cara
    ! calcular autocorrelacion para particula tipo 2 en cara
    ! gNCorrVfac_2
    print *, "Calcula auto corre v particula tipo 2 en cara"
    allocate (Modos_vibra(1:3, 1:(gNCorrVfac_2)))
    !call calcula_autocorr_vel_3D_b (gCorrVfac_2(:,1:gNCorrVfac_2 - 1), Corr_tr)
    call calcula_modos_vibracion_vel(gCorrVfac_2(:,1:(gNCorrVfac_2)),&
      Modos_vibra)
    call escribe_en_columnas(Modos_vibra,'corr_vfac_2.dat',gNmed*gDt)
    ! grabar resultados que estan en Corr_tr
    deallocate (Modos_vibra)
  endif  
#endif

  end subroutine modos_equivalentes


#ifdef MODOS_VIB
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

#endif 


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
    integer :: n, nc
    complex (dp), allocatable :: out(:), tmp(:)
    integer ( kind = 8 ) plan_forward
    integer ( kind = 8 ) plan_backward
    real(kind = 8), allocatable :: in(:)

    n = size(v)
    nc = n / 2 + 1

    allocate (out(1:nc))
    allocate (tmp(1:nc))
    allocate(in(1:n))

    in(1:n) = v(1:n)

    print *, "calcula_autocorr_fft n nc", n, nc, size(in), size(c), size(out)

    ! out tiene tamaño nc
    call dfftw_plan_dft_r2c_1d_ (plan_forward, n, in, out, FFTW_ESTIMATE)
    call dfftw_execute_ (plan_forward)
    

    ! paso c) página 190 del Allen
    tmp = REAL(out) ** 2 + AIMAG(out) ** 2

    print *, "luego de fft3w"
    ! paso d) aplicar una fft inversa  a c
    call dfftw_plan_dft_c2r_1d_ ( plan_backward, n, tmp, c, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_backward )

    print *, c

    deallocate (out)
    deallocate (tmp)

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
    v2(1:n) = v(2,:)
    call calcula_autocorr (v2, corr(2,:))
    v2(1:n) = v(3,:)
    call calcula_autocorr (v2, corr(3,:))
    
 
    deallocate (v2)
  endsubroutine calcula_autocorr_vel_3D 


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
