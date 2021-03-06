module mediciones

#include "control.h"

  use types,      only: dp
  use globales,   only: gT, gNpart, gV, gM, gR, gF, gL, gSigma, gEpsil, gRc2, &
                        gPot_cut, gRho, gVol, gPot, gKin, gVir, gDt
  use ziggurat,   only: rnor

! Si se calcula la función g(r)
#ifdef CORR_PAR
  use globales,   only: gCorr_par, gNhist, gNgr, gDbin
#endif

! Si se usa el termostato de Langevin
#if THERM == 1
  use globales,   only: gGamma
#endif

! Si se utiliza openmp
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  private

  public :: calcula_pres, calcula_kin, calcula_fuerza, calcula_temp

contains
  
  !===============================================================================
  ! CALCULA LA FUERZA  
  !===============================================================================
  ! Calcula la fuerza ya la energía potencial del potencial de L-J

  subroutine calcula_fuerza()

    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: Fij       ! Módulo fuerza entre partículas i y j
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    integer                 :: i,j       
    real(dp)                :: cut4      ! Cuarta parte del potencial en r_c
#ifdef CORR_PAR
    real(dp)                :: r         ! Distancia entre partículas
    integer                 :: ind_bin   ! Indice de cada bin de la g(r)
    
    gNgr = gNgr + 1    ! Cuenta las veces que es llamada 
#endif

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
! Esto implicaría peor performance en caso de no utilizar openmp.
! Por ejemplo, para Npart = 100, Ntime = 10000
! los tiempos serían:  sin openmp  - 2.8
!                      1 thread    - 5.0
!                      2 threads   - 2.6
!                      4 threads   - 1.5
#ifdef _OPENMP
! Si se compila con OPENMP
! Se usa un loop corriendo sobre todas los ij. Se calcula la fuerza sólo una vez y se
! divide el potencial por 1/2 para cada partícula

#ifdef CORR_PAR
!$omp parallel &
!$omp shared (gNpart, gR, gL, gRc2, gF, cut4, gDbin, gCorr_par ) &
!$omp private (i, j, rij_vec, r2ij, r2in, r6in, Fij, r, ind_bin)
#else /* Si no se calcula la g(r) */
!$omp parallel &
!$omp shared (gNpart, gR, gL, gRc2, gF, cut4 ) &
!$omp private (i, j, rij_vec, r2ij, r2in, r6in, Fij)
#endif /* Fin CORR_PAR */
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
#ifdef CORR_PAR
          ! Calcula la función g(r) -  Ver pg 86 Frenkel
          r = sqrt(r2ij)
          if (r < gL/2.0_dp) then                           ! Sólo particulas a menos de gL/2
            ind_bin            = int(r/gDbin) + 1           ! En dónde cae la partícula
                                                            ! +1 vva porque definí indices 1:Nh
            gCorr_par(ind_bin) = gCorr_par(ind_bin) + 1     ! Actualizo contador del bin
          end if
#endif /* Fin CORR_PAR */
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
#ifdef CORR_PAR /* Si se calcula g(4) */
          ! Calcula la función g(r) - Ver pg 86 Frenkel
          r = sqrt(r2ij)
          if (r < gL/2.0_dp) then                           ! Sólo particulas a menos de gL/2
            ind_bin            = int(r/gDbin) + 1           ! En dónde cae la partícula
                                                            ! +1 vva porque definí indices 1:Nh
            gCorr_par(ind_bin) = gCorr_par(ind_bin) + 2     ! Actualizo contador del bin
          end if
#endif /* Fin CORR_PAR */
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

#if THERM == 1
  !===============================================================================
  ! TERMOSTATO DE LANGEVIN
  !===============================================================================
  ! Agrega la parte estocástica al a fuerza de cada partícula

  subroutine fuerza_langevin()

  integer  :: i, j

  do i = 1, gNpart
    do j = 1, 3
      gF(j,i) = gF(j,i) - gGamma * gV(j,i) + sqrt(2.0_dp*gT*gGamma/gDt) * rnor()
    end do
  end do

  end subroutine fuerza_langevin
#endif
 
  !===============================================================================
  ! CALCULA  PRESION 
  !===============================================================================
  ! Calcula la presion en base al teorema del virial (Ver 2.4 del Allen)

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
    
    real(dp), dimension(gNpart)   :: v2     ! Vector con la velocidad cuadratica
    integer                       :: i

!$omp parallel &
!$omp shared( gNpart, v2, gV) &
!$omp private ( i )
!$omp do
    do i = 1, gNpart
      v2(i) =  dot_product(gV(:,i),gV(:,i))  
    end do
!$omp end do
!$omp end parallel

!$omp parallel &
!$omp shared (gKin, gM, v2)
!$omp workshare
    gKin = 0.5_dp * gM * sum( v2 )
!$omp end workshare
!$omp end parallel

  end subroutine calcula_kin

  !===============================================================================
  ! CALCULA TEMPERATURA INSTANTANEA 
  !===============================================================================
  ! Calcula la anergia cinetica total del sistema

  subroutine calcula_temp(Temp_ins)
    
    real(dp), intent(out)         :: Temp_ins ! Temperatura instantánea

!$omp parallel &
!$omp shared (gKin, Temp_ins, gNpart)
!$omp workshare
    Temp_ins = 2.0_dp * gKin / (3*gNpart-3)
!$omp end workshare
!$omp end parallel

  end subroutine calcula_temp


end module mediciones
