module mediciones

#include "control.h"

  use types,      only: dp
  use globales,   only: gKmed, gV, gF, gR, gNpart, gKin, gNespecies,&
                        gCorrVfac_1, gCorrVfac_2, gCorrVfac_3, gCorrVver_1,&
                        gNCorrVfac_1, gNCorrVfac_2, gNCorrVfac_3, gNCorrVver_1,&
                        gNmodosVibra, gRho, gVir, gVol, gGamma, gDt,&
                        gPot, gTemperatura, gLado_caja, gCombEpsilon, &
                        gIndice_elemento, gMasa, gPot_cut, gRc2, gCombSigma
                         
  !use globales
  use ziggurat

! Si se utiliza openmp
#ifdef _OPENMP
  use omp_lib
#endif

! Si se calcula la función g(r)
#ifdef CORR_PAR
  use globales,      only: gCorr_par, gNhist, gNgr, gDbin 
#endif

  !#include  "mpif.h"
  !use mpi

  implicit none

  private

  public   :: calcula_fuerza, calcula_pres, calcula_kin, calcula_temp
#ifdef MODOS_VIB
  public   ::   acumula_velocidades_posicion, acumula_velocidades_equivalentes
#endif
#ifdef CORR_PAR
  public   :: normaliza_gr
#endif

contains

#ifdef LUIS 
  ! TODO: Esta se podría poner como subrutina aparte para que la calcule leyendo
  ! las configuraciones guardas en algún archivo. Ahí creo que se justifica volver
  ! a poner el doble loop de interacciones por separado.
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
#endif

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
    real(dp)                :: rc2       ! Cuadrado del radio de corte (ij)
    real(dp)                :: cut4      ! Cuarta parte del potencial en r_c (ij)
    integer                 :: ei, ej    ! Tipo de elemento (especie i ó j)
    real(dp)                :: epsil     ! epsilon de la interacción ei-ej
    real(dp)                :: sigma     ! sigma de la interacción ei-ej
#ifdef CORR_PAR
    real(dp)                :: r         ! Distancia entre partículas
    integer                 :: ind_bin   ! Indice de cada bin de la g(r)
    integer, dimension(2,2) :: igr       ! Matriz para ubicar la columna en g(r)
    integer                 :: ingr      ! Índice de la matriz igr
       
    gNgr = gNgr + 1    ! Cuenta las veces que es llamada 
    ! Creo una matriz para ubicar cada interacción ij en la clumna correcta
    ! de g(r) = [ g_11  , g12, g22]
    igr  = reshape( (/ 1, 2, 2, 3/), (/2,2/) ) 
#endif

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

! En OpenMP se deben definir las variables que se utilizan para calcular la g(r)
! Debo distinguir los dos casos:
#ifdef CORR_PAR /* Si se calcula la g(r) */
!$omp parallel &
!$omp shared (gNpart, gR, gLado_caja, gRc2, gF, gIndice_elemento, gCombEpsilon, gCombSigma,gCorr_par) &
!$omp private (i, j, rij_vec, r2ij, r2in, r6in, Fij,ei,ej,cut4,rc2,epsil,sigma,r,ind_bin)
#else /* Si no se calcula la g(r) */
!$omp parallel &
!$omp shared (gNpart, gR, gLado_caja, gRc2, gF, gIndice_elemento, gCombEpsilon, gCombSigma) &
!$omp private (i, j, rij_vec, r2ij, r2in, r6in, Fij,ei,ej,cut4,rc2,epsil,sigma)
#endif /* Fin CORR_PAR */
!$omp do reduction( + : gPot, gVir)
! do schedule(static,5) reduction( + : gPot, gVir)
! El static es casi irrelevante en este loop, porque no hay un desbalance de carga
! significativo. Serviría si se recorre el loop con i<j

     do i = 1, gNpart
      do j = 1, gNpart
        if (i /= j) then
          ei = gIndice_elemento(i)
          ej = gIndice_elemento(j)
          ! Cuarta parte del potencial evaluado en el radio de corte
          ! Para no hacer una multiplicación el el doble loop
          cut4 = gPot_cut(ei, ej) / 4.0_dp
          ! El radio de corte es siempore 2.5 por el sigma de cada interacción
          ! Lo comparo luego con la distancia dividida por el sigma de cada interacción
          rc2   = gRc2( ei, ej) 
          epsil = gCombEpsilon(ei,ej)
          sigma = gCombSigma(ei,ej)
          ! Distancia vectorial
          rij_vec = gR(:,i) - gR(:,j)
          ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
          ! Siempre en distancias relativas de sigma
          rij_vec = rij_vec - gLado_caja*nint(rij_vec/gLado_caja)
          ! Divido por el sigma de cada interacción
          ! Después sólo lo utilizo para calcular fuerza y potencial, siempre
          ! aparece la distancia dividida por el sigma.
          rij_vec = rij_vec / sigma
          r2ij   = dot_product( rij_vec , rij_vec )         ! Cuadrado de la distancia
          if ( r2ij < rc2 ) then
            r2in = 1.0_dp/r2ij                              ! Inversa al cuadrado
            r6in = r2in**3                                  ! Inversa a la sexta
            Fij  = r2in * r6in * epsil * (r6in - 0.5_dp)/sigma   ! Fuerza entre partículas
            gF(:,i) = gF(:,i) + Fij * rij_vec               ! Contribución a la partícula i
            gPot  = gPot + r6in * epsil * ( r6in - 1.0_dp) - cut4 ! Energía potencial
            gVir  = gVir + Fij * r2ij                       ! Término del virial para la presión
                                                            ! pg 48 de Allen W=-1/3 sum(r dv/dr)
          end if  ! Termina if del radio de corte
#ifdef CORR_PAR
          ! Calcula la función g(r) -  Ver pg 86 Frenkel
          ! Se debe a unidades absolutas (antes estaba r2ij dividiendo por el sigma)
          ingr = igr(ei,ej)                           ! Averigua en cuál g(r) se va a acumular 
          r = sqrt(r2ij) * sigma
          if (r < gLado_caja/2.0_dp) then             ! Sólo particulas a menos de gL/2
            ind_bin            = int(r/gDbin) + 1     ! En dónde cae la partícula
                                                      ! Va +1 porque definí indices 1:Nh
            gCorr_par(ingr, ind_bin) = gCorr_par(ingr, ind_bin) + 1  ! Actualizo contador del bin
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

! Se usa un loop corriendo sólo sobre los j<i. Se asigna la fuerza a dos partículas
! con signo contrario. Se calcula el potencial por cada interacción (sin repetir)

    do i = 1, gNpart - 1
      do j = i+1, gNpart
        ei = gIndice_elemento(i)
        ej = gIndice_elemento(j)
        ! Cuarta parte del potencial evaluado en el radio de corte
        ! Para no hacer una multiplicación el el doble loop
        cut4 = gPot_cut(ei, ej) / 4.0_dp
        ! El radio de corte es siempore 2.5 por el sigma de cada interacción
        ! Lo comparo luego con la distancia dividida por el sigma de cada interacción
        rc2   = gRc2( ei, ej) 
        epsil = gCombEpsilon(ei,ej)
        sigma = gCombSigma(ei,ej)
        ! Distancia vectorial
        rij_vec = gR(:,i) - gR(:,j)
        ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
        ! Siempre en distancias relativas de sigma
        rij_vec = rij_vec - gLado_caja*nint(rij_vec/gLado_caja)
        ! Divido por el sigma de cada interacción
        ! Después sólo lo utilizo para calcular fuerza y potencial, siempre
        ! aparece la distancia dividida por el sigma.
        rij_vec = rij_vec / sigma
        r2ij   = dot_product( rij_vec , rij_vec )         ! Cuadrado de la distancia
        if ( r2ij < rc2 ) then
          r2in = 1.0_dp/r2ij                               ! Inversa al cuadrado
          r6in = r2in**3                                   ! Inversa a la sexta
          Fij     = r2in * r6in * epsil *(r6in - 0.5_dp)/sigma    ! Fuerza entre partículas
          gF(:,i) = gF(:,i) + Fij * rij_vec                ! Contribución a la partícula i
          gF(:,j) = gF(:,j) - Fij * rij_vec                ! Contribucion a la partícula j
          gPot    = gPot + r6in * epsil * ( r6in - 1.0_dp) - cut4 ! Energía potencial
          gVir    = gVir + Fij * r2ij                      ! Término del virial para la presión
                                                           ! pg 48 de Allen W=-1/3 sum(r dv/dr)
        end if  ! Termina if del radio de corte
#ifdef CORR_PAR
          ! Calcula la función g(r) -  Ver pg 86 Frenkel
          ! Se debe a unidades absolutas (antes estaba r2ij dividiendo por el sigma)
          r = sqrt(r2ij) * sigma
          if (r < gLado_caja/2.0_dp) then                     ! Sólo particulas a menos de gL/2
            ind_bin            = int(r/gDbin) + 1     ! En dónde cae la partícula
                                                      ! Va +1 porque definí indices 1:Nh
            gCorr_par(igr, ind_bin) = gCorr_par(igr, ind_bin) + 2  ! Actualizo contador del bin
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

#ifdef CORR_PAR
  !===============================================================================
  ! NORMALIZACIÓN Y ESCRITURA DE LA g(r)
  !===============================================================================
  ! Esta subrutina normaliza e imprime la función g(r) ya calculada
  ! Agrega una cuarta columna con la g(r)_total, que es la correlación de pares
  ! sin importar la especie de partícula.

  subroutine normaliza_gr(grnor)
  ! Subrutina para normalizar y calcular la función g(r)     
  ! ver pg. 86 del Frenkel

    real(dp), dimension(:,:),intent(out) :: grnor  ! Función g(r)
    real(dp)                             :: dvol   ! Diferencie de volumen entre bines
    real(dp)                             :: nid    ! Parte de gas ideal en dvol
    integer                              :: i,j
    real(dp), parameter                  :: PI = 3.1415926535898_dp  ! pi
             

    do i = 1, gNhist
      grnor(1,i) = gDbin * (i + 0.5)                ! Posición centrada de cada bin
      dvol  = ( (i+1)**3 - i**3 ) * gDbin**3        ! Diferencia de vol entre los bin (i+1) e i
      nid   = (4.0_dp/3.0_dp) * PI * dvol * gRho    ! Parte del gas ideal en dvol 
      do j = 2, 4
        grnor(j,i) = real(gCorr_par(j,i),kind=dp) / (gNgr*gNpart*nid)      ! Función g(r) normalizada
      end do
        grnor(5,i) = real( sum(gCorr_par(:,i) ),kind=dp)/ (gNgr*gNpart*nid)  ! Función total
    end do

  end subroutine normaliza_gr 

#endif

  
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
        gF(j,i) = gF(j,i) - ma*gGamma * gV(j,i) + &
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

#ifdef  MODOS_VIB

  subroutine acumula_velocidades_posicion (VV)
    integer :: i, j
    real(dp), dimension(:,:) :: VV

    ! si todavia no se llegó al fin
    if (gNmodosVibra .lt. gKmed * gNpart) then
      !inicializo la velocidad de esa posicion a cero      
      gCorrVver_1(:, gNmodosVibra) = 0
      gCorrVfac_1(:, gNmodosVibra) = 0
      gCorrVfac_2(:, gNmodosVibra) = 0
      gCorrVfac_3(:, gNmodosVibra) = 0
      do i=1,gNpart / 4
        j = 4 * (i - 1)

        gCorrVver_1(:, gNmodosVibra) = gCorrVver_1(:, gNmodosVibra) + VV(:, j + 1)
        gCorrVfac_1(:, gNmodosVibra) = gCorrVfac_1(:, gNmodosVibra) + VV(:, j + 2)
        gCorrVfac_2(:, gNmodosVibra) = gCorrVfac_2(:, gNmodosVibra) + VV(:, j + 3)
        gCorrVfac_3(:, gNmodosVibra) = gCorrVfac_3(:, gNmodosVibra) + VV(:, j + 4)
  
      enddo
      gCorrVver_1(:, gNmodosVibra) = 4 * gCorrVver_1(:, gNmodosVibra) / gNpart
      gCorrVfac_1(:, gNmodosVibra) = 4 * gCorrVfac_1(:, gNmodosVibra) / gNpart
      gCorrVfac_2(:, gNmodosVibra) = 4 * gCorrVfac_2(:, gNmodosVibra) / gNpart
      gCorrVfac_3(:, gNmodosVibra) = 4 * gCorrVfac_3(:, gNmodosVibra) / gNpart
      gNmodosVibra = gNmodosVibra + 1
    endif  
  endsubroutine acumula_velocidades_posicion
  
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
         gCorrVfac_3(:, gNCorrVfac_3) = gV(:, j + 1)
         gNCorrVfac_3 = gNCorrVfac_3 + 1 
       endif  
             
       ! cara 1
       if (gIndice_elemento(j + 2) .eq. 1) then 
         gCorrVfac_1(:, gNCorrVfac_1) = gV(:, j + 2)
         gNCorrVfac_1 = gNCorrVfac_1  + 1 
       else 
         gCorrVfac_2(:, gNCorrVfac_2) = gV(:, j + 2)
         gNCorrVfac_2 = gNCorrVfac_2 + 1 
       endif 

       ! cara 2
       if (gIndice_elemento(j + 3) .eq. 1) then 
         gCorrVfac_1(:, gNCorrVfac_1) = gV(:, j + 3)
         gNCorrVfac_1 = gNCorrVfac_1 + 1 
       else 
         gCorrVfac_2(:, gNCorrVfac_2) =  gV(:, j + 3)
         gNCorrVfac_2 = gNCorrVfac_2 + 1 
       endif 

       ! cara 3
       if (gIndice_elemento(j + 4) .eq. 1) then 
         gCorrVfac_1(:, gNCorrVfac_1) = gV(:, j + 4)
         gNCorrVfac_1 = gNCorrVfac_1 + 1 
       else 
         gCorrVfac_2(:, gNCorrVfac_2) = gV(:, j + 4)
         gNCorrVfac_2 = gNCorrVfac_2 + 1 
       endif 
    enddo
  endsubroutine acumula_velocidades_equivalentes
#endif


end module mediciones
