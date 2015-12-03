module FFTW3
  use, intrinsic :: iso_c_binding
  include '/usr/include/fftw3.f03'
  ! También podría estar en '/usr/local/include/fftw3.f03
end module

module mediciones

#include "control.h"

  use types,      only: dp
  use globales
  use FFTW3      

  ! forma de llamar fftw3
  ! ver: http://fftw.org/doc/Calling-FFTW-from-Modern-Fortran.html
  !type(C_PTR) :: plan
  !complex(C_DOUBLE_COMPLEX), dimension(1024,1000) :: in, out
  !plan = fftw_plan_dft_2d(1000,1024, in,out, FFTW_FORWARD,FFTW_ESTIMATE)

  !call fftw_execute_dft(plan, in, out)

  use ziggurat
!#include  "mpif.h"
  !use mpi

  implicit none

  private

  public   :: calcula_fuerza, calcula_pres, calcula_kin, calcula_temp
#ifdef LUIS
  public   ::   acumula_velocidades_equivalentes, calcula_corr_vel_3D_b
#endif

contains
 
  !===============================================================================
  ! CALCULA g(r)
  !===============================================================================

  subroutine radial_distribution (switch, nhis)
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    integer :: i, j, switch, ig, i_el, j_el
    integer :: nhis ! numero total de bines

    if (switch.eq.0) then !inicio
      gNgr = 0
      gDbin = gLado_caja / (2*nhis)
      ! alloca una matriz con espacio suficiente para
      ! alojar todas las correlaciones espaciales
      ! en las filas de la matriz estan las relaciones
      ! 11 12 13 14 .. 22 23 24 .. 33 34 .. 44
      allocate (gCorr_par(1:gNespecies * (gNespecies - 1) / 2 + gNespecies , 1:nhis))
      do i=1, nhis
        gCorr_par = 0
      enddo  
    else if (switch.eq.1) then ! muestreo       
      do i = 1, gNpart -1
        do j = 1, gNpart
         
          rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
          rij_vec = rij_vec - gLado_caja*nint(rij_vec/gLado_caja)
          r2ij   = sqrt(dot_product( rij_vec , rij_vec ))          ! Cuadrado de la distancia
          if (r2ij .lt. gLado_caja / 2) then

            ! determina sobre que elementos se toma
            ! la distancia
            i_el = gIndice_elemento(i)
            j_el = gIndice_elemento(j)

            ig = int(r2ij/gDbin)
            if (i_el .ne. 1) then
                i_el = 2 * i_el
            endif        
            gCorr_par(i_el + j_el -1 , ig) = gCorr_par(i_el + j_el - 1, ig) + 2
          endif 
        enddo
      enddo
    else if (switch.eq.2) then ! resultados
     ! normalizacion de la gdr       
     gCorr_par = gCorr_par / maxval(gCorr_par)
    endif        
  endsubroutine radial_distribution

!
! Se puso la subrutina de fuerza igual a la de dinmod
! Se puede cambiar entre una u otra con el prepro (PABLO==0 no tiene en cuenta
! lo que yo hice)
!
#if PABLO == 0
  !===============================================================================
  ! KERNEL FUERZA 
  !===============================================================================

  subroutine kernel_fuerza(i, j, i_inter, j_inter)

    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: Fij       ! Módulo fuerza entre partículas i y j
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6

    real(dp) :: rc2       ! Radio de corte para la interacción ij
    real(dp) :: cut4      ! Cuarta parte del potencial en r_c
    integer, intent(in) :: i, j, i_inter, j_inter

    ! Cuarta parte del potencial evaluado en el radio de corte
    cut4 = gPot_cut(i_inter, j_inter) / 4.0_dp
    ! El radio de corte es siempore 2.5 por el sigma de cada interacción
    ! Lo comparo luego con la distancia dividida por el sigma de cada interacción
    rc2  = gRc2(i_inter, j_inter) 


    rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
    ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
    ! Siempre en distancias relativas de sigma
    rij_vec = rij_vec - gLado_caja*nint(rij_vec/gLado_caja)
    ! Divido por el sigma de cada interacción
    ! Después sólo lo utilizo para calcular fuerza y potencial, siempre
    ! aparece la distancia dividida por el sigma.
    rij_vec = rij_vec / gCombSigma(i_inter,j_inter) 

    r2ij   = dot_product( rij_vec , rij_vec )  ! Cuadrado de la distancia
    if ( r2ij < rc2) then                
      r2in = 1.0_dp / r2ij                            ! Inversa al cuadrado
      r6in = r2in**3                                  ! Inversa a la sexta
      ! Se debe dividir por el sigma que proviene de la derivación del potencial
      Fij     = r2in * r6in * gCombEpsilon(i_inter, j_inter) * &
                (r6in - 0.5_dp) / gCombSigma(i_inter,j_inter)     ! Fuerza entre partículas
      gF(:,i) = gF(:,i) + Fij * rij_vec               ! Contribución a la partícula i
      gF(:,j) = gF(:,j) - Fij * rij_vec               ! Contribucion a la partícula j
      gPot    = gPot +  r6in * ( r6in - 1.0_dp) *  &
                gCombEpsilon(i_inter, j_inter) - cut4 ! Energía potencial
      gVir    = gVir + Fij * r2ij                     ! Término del virial para la presión
                                                      ! pg 48 de Allen W=-1/3 sum(r dv/dr)
    end if  ! Termina if del radio de corte
  endsubroutine kernel_fuerza

  !===============================================================================
  ! CALCULA LA FUERZA 
  !===============================================================================

  subroutine calcula_fuerza()

    integer :: i, j, i_inter, j_inter

    gF    = 0.0_dp
    gPot  = 0.0_dp
    gVir  = 0.0_dp

    do i=1, gNpart - 1
       do j=i + 1, gNpart
          i_inter = gIndice_elemento(i)
          j_inter = gIndice_elemento(j)
          !print *, "---->", i, j, i_inter, j_inter 
          call kernel_fuerza (i, j, i_inter, j_inter)
       end do
    end do

    ! Constantes que faltaban en la energía
    gF = 48.0_dp * gF
    ! Constantes que faltaban en el potencial
    gPot =  4.0_dp * gPot
    ! Se agregan las constantes que faltan para el término del virial
    gVir = 48.0_dp * gVir / 3.0_dp
! Si se utiliza el termostato de Langevin
#if THERM == 1
    call fuerza_langevin()
#endif

  endsubroutine calcula_fuerza

#else /* PABLO distinto de 0 */
  !===============================================================================
  ! CALCULA LA FUERZA (DINMOD)  
  !===============================================================================
  ! Calcula la fuerza ya la energía potencial del potencial de L-J

  subroutine calcula_fuerza()

    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: Fij       ! Módulo fuerza entre partículas i y j
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    integer                 :: i,j
    real(dp)                :: cut4      ! Cuarta parte del potencial en r_c

    ! Por cuestiones de eficiencia. Se evita hacer una multiplicación dentro
    ! del loop anidado
    cut4 = gPot_cut / 4.0_dp

!$omp parallel workshare
    ! Se van a acumular las fuerzas. Se comienza poniendo todas a cero.
    gF    = 0.0_dp
    gPot  = 0.0_dp
    gVir  = 0.0_dp
!$omp end parallel workshare

! Se escribe dos loops distintos dependiendo de si se compila el programa con
! OPENMP o no. La razón es que para paralelizar correctamente el loop de fuerza
! se debe cambiar la forma de recorrer las interacciones i,j en vez de j<i

#ifdef _OPENMP
! Si se compila con OPENMP
! Se usa un loop corriendo sobre todas los ij. Se calcula la fuerza sólo una vez y se
! divide el potencial por 1/2 para cada partícula

!$omp parallel &
!$omp shared (gNpart, gR, gL, gRc2, gF, cut4 ) &
!$omp private (i, j, rij_vec, r2ij, r2in, r6in, Fij)

!$omp do schedule(static,5) reduction( + : gPot, gVir)
! El static es casi irrelevante en este loop, porque no hay un desbalance de carga
! significativo. Serviría si se recorre el loop con i<j

     do i = 1, gNpart
      do j = 1, gNpart
        if (i /= j) then
          rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
          ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
          ! Siempre en distancias relativas de sigma
          rij_vec = rij_vec - gL*nint(rij_vec/gL)
          r2ij   = dot_product( rij_vec , rij_vec )         ! Cuadrado de la distancia
          if ( r2ij < gRc2 ) then
            r2in = 1.0_dp/r2ij                              ! Inversa al cuadrado
            r6in = r2in**3                                  ! Inversa a la sexta
            Fij     = r2in * r6in * (r6in - 0.5_dp)         ! Fuerza entre partículas
            gF(:,i) = gF(:,i) + Fij * rij_vec               ! Contribución a la partícula i
            gPot    = gPot + r6in * ( r6in - 1.0_dp) - cut4 ! Energía potencial
            gVir    = gVir + Fij * r2ij                     ! Término del virial para la presión
                                                            ! pg 48 de Allen W=-1/3 sum(r dv/dr)
          end if  ! Termina if del radio de corte
        end if   ! Termina if de i /= j
      end do
    end do

!$omp end do
!$omp end parallel
    gPot = 0.5_dp * gPot   ! En este loop se cuentan dos veces las interacciones
    gVir = 0.5_dp * gVir    ! En este loop se cuentan dos veces las interacciones

#else /* Si no se compila con OPENMP */
! Si no se compila con OPENMP
! Se usa un loop corriendo sólo sobre los j<i. Se asigna la fuerza a dos partículas
! con signo contrario. Se calcula el potencial por cada interacción (sin repetir)

    do i = 1, gNpart - 1
      do j = i+1, gNpart
        rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
        ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
        ! Siempre en distancias relativas de sigma
        rij_vec = rij_vec - gL*nint(rij_vec/gL)
        r2ij   = dot_product( rij_vec , rij_vec )          ! Cuadrado de la distancia
        if ( r2ij < gRc2 ) then
          r2in = 1.0_dp/r2ij                               ! Inversa al cuadrado
          r6in = r2in**3                                   ! Inversa a la sexta
          Fij     = r2in * r6in * (r6in - 0.5_dp)          ! Fuerza entre partículas
          gF(:,i) = gF(:,i) + Fij * rij_vec                ! Contribución a la partícula i
          gF(:,j) = gF(:,j) - Fij * rij_vec                ! Contribucion a la partícula j
          gPot    = gPot + r6in * ( r6in - 1.0_dp) - cut4  ! Energía potencial
          gVir    = gVir + Fij * r2ij                      ! Término del virial para la presión
                                                           ! pg 48 de Allen W=-1/3 sum(r dv/dr)
        end if  ! Termina if del radio de corte
      end do
    end do

#endif /* Fin _OPENMP */

!$omp parallel workshare
    ! Constantes que faltaban en la energía
    gF = 48.0_dp * gF
    ! Constantes que faltaban en el potencial
    gPot =  4.0_dp * gPot
    ! Se agregan las constantes que faltan para el término del virial
    gVir = 48.0_dp*gVir / 3.0_dp
!$omp end parallel workshare

! Si se utiliza el termostato de Langevin
#if THERM == 1
    call fuerza_langevin()
#endif

  end subroutine calcula_fuerza

#endif /* Fin PABLO */

#if THERM == 1
  !===============================================================================
  ! TERMOSTATO DE LANGEVIN
  !===============================================================================
  ! Agrega la parte estocástica al a fuerza de cada partícula
  ! TODO No está andando bien
  subroutine fuerza_langevin()

    real(dp) :: ma    ! Masa de la partícula
    integer  :: i, j

    do i = 1, gNpart
      ma = gMasa(gIndice_elemento(i))
      do j = 1, 3
        gF(j,i) = gF(j,i) - gGamma * gV(j,i) + &
                  sqrt(2.0_dp*gTemperatura*gGamma*ma/gDt) * rnor()
      end do
    end do

  end subroutine fuerza_langevin
#endif
 
  !===============================================================================
  ! CALCULA  PRESION 
  !===============================================================================
  ! Calcula la presion en base al teorema del virial (Ver 2.4 del Allen)
  ! TODO: Revisar si está bien

  subroutine calcula_pres(presion)

    real(dp), intent(out)    :: presion
    real(dp)                 :: tempe

   ! Calcula la temperatura instantánea
   call calcula_temp(tempe)
   ! Calcula la presión con la temperatura instantánea
   presion = gRho * tempe + gVir / gVol 

  end subroutine calcula_pres

  !===============================================================================
  ! CALCULA ENERGIA CINETICA 
  !===============================================================================
  ! Calcula la anergia cinetica total del sistema

  subroutine calcula_kin()
    
    real(dp), dimension(gNpart)   :: mv2    ! Vector con la masa* velocidad cuadratica
    integer ::  i

    do i = 1, gNpart
      mv2(i) =  gMasa(gIndice_elemento(i))*dot_product(gV(:,i),gV(:,i))
    end do

    gKin = 0.5_dp * sum(mv2) 

  endsubroutine calcula_kin

  !===============================================================================
  ! CALCULA TEMPERATURA INSTANTANEA 
  !===============================================================================
  ! Calcula la anergia cinetica total del sistema

  subroutine calcula_temp(Temp_ins)
    
    real(dp), intent(out)      :: Temp_ins ! Temperatura instantánea

    Temp_ins = 2.0_dp * gKin / (3*gNpart-3)

  end subroutine calcula_temp

#ifdef LUIS
  !############################################################################!
  ! calculo de autocorrelaciones 
  ! valido para sistemas binarios
  !############################################################################!
  ! se llama  cuando se tiene un estado de velocidades a analizar
  ! antes de llamar hay que inicializar gNCorr* y gCorr* a cero
  subroutine acumula_velocidades_equivalentes ()
    integer :: i, j
    ! recorrer el array de velocidades
    ! y asigna al array de velocidades 
    ! equivalentes
    do i=1,gNpart / 4
       j = 4 * (i - 1)
       ! vertice 
       if (gIndice_elemento(j + 1) .eq. 1) then 
         gCorrVver_1(:, gNCorrVver_1) = gV(:, j + 1)
         gNCorrVver_1 = gNCorrVver_1 + 1 
       else 
         gCorrVver_2(:, gNCorrVver_2) = gV(:, j + 1)
         gNCorrVver_2 = gNCorrVver_2 + 1 
       endif  
             
       ! cara 1
       if (gIndice_elemento(j + 2) .eq. 1) then 
         gCorrVfac_1(:, gNCorrVfac_1) = gV(:, j + 1)
         gNCorrVfac_1 = gNCorrVfac_1  + 1 
       else 
         gCorrVfac_2(:, gNCorrVfac_2) = gV(:, j + 1)
         gNCorrVfac_2 = gNCorrVfac_2 + 1 
       endif 

       ! cara 2
       if (gIndice_elemento(j + 3) .eq. 1) then 
         gCorrVfac_1(:, gNCorrVfac_1) = gV(:, j + 1)
         gNCorrVfac_1 = gNCorrVfac_1 + 1 
       else 
         gCorrVfac_2(:, gNCorrVfac_2) =  gV(:, j + 1)
         gNCorrVfac_2 = gNCorrVfac_2 + 1 
       endif 

       ! cara 3
       if (gIndice_elemento(j + 4) .eq. 1) then 
         gCorrVfac_1(:, gNCorrVfac_1) = gV(:, j + 1)
         gNCorrVfac_1 = gNCorrVfac_1 + 1 
       else 
         gCorrVfac_2(:, gNCorrVfac_2) = gV(:, j + 1)
         gNCorrVfac_2 = gNCorrVfac_2 + 1 
       endif 
    enddo
  endsubroutine acumula_velocidades_equivalentes

  !
  ! calcula la auto_correlacion de un vector 1D
  ! con fftw3
  ! el tamaño de vel es del tipo 2 ** n, n entero
  ! los ultimos n elementos en 0
  ! de acuerdo a Allen
  subroutine calcula_corr_vel (in, c)
    real(dp), dimension(:) :: in, c
    integer :: n, nc
    complex (dp), allocatable :: out(:)
    integer ( kind = 8 ) plan_forward
    integer ( kind = 8 ) plan_backward
    real (dp), allocatable :: v2(:)
    real (dp), allocatable :: in2(:)

    n = size(in)
    nc = n / 2 + 1

    allocate (out(1:nc))

    ! out tiene tamaño nc
    call dfftw_plan_dft_r2c_1d_ ( plan_forward, n, in, out, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_forward )

    ! paso c) página 190 del Allen
    out = REAL(out) ** 2 + AIMAG(out) ** 2

    ! paso d) aplicar una fft inversa  a c
    call dfftw_plan_dft_c2r_1d_ ( plan_backward, n, out, c, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_backward )

    ! paso e) aplicar la normalizacion
    c = c / n

  endsubroutine calcula_corr_vel


  ! calcula la correlacion de velocidad de un vector temporal 3D
  subroutine calcula_corr_vel_3D (v, c)
    real(dp), dimension(:,:) :: v
    real(dp), dimension(:) :: c
    real(dp), allocatable :: corr(:,:)
    real(dp), allocatable :: v2(:,:)

    ! viene el resultado de la operacion
    allocate (corr(1:3, 1: 2 * size(v)))
    allocate (v2(1:3, 1: 2 * size(v)))
    v2 = 0
    corr = 0
    v2 (:,1:size(v)) = v

    call calcula_corr_vel (v2(1,:), corr(1,:))
    call calcula_corr_vel (v2(2,:), corr(2,:))
    call calcula_corr_vel (v2(3,:), corr(3,:))
    
    ! acumula las 3 dimensiones del vector velocidad
    c = (corr(1,:) + corr(2,:) + corr(3,:)) / 3

  endsubroutine calcula_corr_vel_3D 

  subroutine calcula_corr_vel_3D_b (v, c)
    ! Se puede pensar de otra manera
    ! acumular en un vector 1D cada una de las 3dimensiones
    ! del array de velocidades
    real(dp), dimension(:,:) :: v
    real(dp), dimension(:) :: c
    real(dp), allocatable :: corr(:)
    real(dp), allocatable :: v2(:)

    ! viene el resultado de la operacion
    allocate (corr(2 * size(v)))
    allocate (v2(2 * size(v)))
    v2 = 0
    corr = 0

    v2(1:size(v)) = (v(1,:) + v(2,:) + v(3,:)) / 3
    call calcula_corr_vel (v2(:), corr(:))

  endsubroutine calcula_corr_vel_3D_b 

#endif

end module mediciones
