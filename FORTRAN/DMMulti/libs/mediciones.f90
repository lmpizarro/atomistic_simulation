module mediciones
#include "control.h"
  use types,      only: dp
  use globales
  
  use ziggurat
!#include  "mpif.h"
  !use mpi
  implicit none


contains

  subroutine radial_distribution (switch, nhis, elemento)
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    integer :: i, j, switch, ig, i_elemento, j_elemento
    integer :: nhis ! numero total de bines
    integer :: elemento ! 0: total, 1: elemento tipo 1, 2: elemento tipo 2

    if (switch.eq.0) then !inicio
      gNgr = 0
      gDbin = gLado_caja / (2*nhis)
      allocate (gCorr_par(1:gNespecies * (gNespecies - 1) / 2 + gNespecies , 1:nhis))
      do i=1, nhis
        gCorr_par = 0
      enddo  
    else if (switch.eq.1) then ! muestreo       
      do i = 1, gNpart -1
        do j = 1, gNpart
          
          i_elemento = gIndice_elemento(i)
          j_elemento = gIndice_elemento(j)

          rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
          rij_vec = rij_vec - gLado_caja*nint(rij_vec/gLado_caja)
          r2ij   = sqrt(dot_product( rij_vec , rij_vec ))          ! Cuadrado de la distancia
          if (r2ij .lt. gLado_caja / 2) then
            ig = int(r2ij/gDbin)
            gCorr_par(i_elemento, ig) = gCorr_par(i_elemento, ig) + 2
          endif 
        enddo
      enddo
    else if (switch.eq.2) then ! resultados
     ! normalizacion de la gdr       
     gCorr_par = gCorr_par / maxval(gCorr_par)
    endif        
  endsubroutine radial_distribution


  subroutine kernel_fuerza(i, j, i_inter, j_inter)
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: Fij       ! Módulo fuerza entre partículas i y j
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    real(dp) :: cut4      ! Cuarta parte del potencial en r_c
    real(dp) :: rc2      ! Cuarta parte del potencial en r_c
    integer, intent(in) :: i, j, i_inter, j_inter

    cut4 = gPot_cut(i_inter, j_inter) / 4.0_dp
    rc2 = gRc2(i_inter, j_inter)

    ! aca hay que traer los parametros de sigma 
    ! y epsilon para calcular el potencial

    rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
    ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
    ! Siempre en distancias relativas de sigma
    rij_vec = rij_vec - gLado_caja*nint(rij_vec/gLado_caja)
    r2ij   = dot_product( rij_vec , rij_vec )          ! Cuadrado de la distancia
    if ( r2ij < rc2) then                
      r2in = (gCombSigma(i_inter,j_inter) **2)/r2ij                               ! Inversa al cuadrado
      r6in = r2in**3                                   ! Inversa a la sexta
      Fij     = r2in * r6in * (r6in - 0.5_dp)          ! Fuerza entre partículas
      gF(:,i) = gF(:,i) + Fij * rij_vec                ! Contribución a la partícula i
      gF(:,j) = gF(:,j) - Fij * rij_vec                ! Contribucion a la partícula j
      gPot    = gPot + gCombEpsilon(i_inter, j_inter) * r6in * &
                 ( r6in - 1.0_dp) - cut4  ! Energía potencial
      gVir    = gVir + Fij * r2ij                      ! Término del virial para la presión
                                                           ! pg 48 de Allen W=-1/3 sum(r dv/dr)
    end if  ! Termina if del radio de corte

  endsubroutine kernel_fuerza

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
       enddo
    enddo

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


#if THERM == 1
  !===============================================================================
  ! TERMOSTATO DE LANGEVIN
  !===============================================================================
  ! Agrega la parte estocástica al a fuerza de cada partícula

  subroutine fuerza_langevin()

    integer  :: i, j

    do i = 1, gNpart
      do j = 1, 3
        gF(j,i) = gF(j,i) - gGamma * gV(j,i) + sqrt(2.0_dp*gTemperatura*gGamma/gDt) * rnor()
      end do
    end do

  end subroutine fuerza_langevin
#endif
 
  !===============================================================================
  ! CALCULA  PRESION 
  !===============================================================================
  ! Calcula la presion en base al teorema del virial (Ver 2.4 del Allen)

  subroutine calcula_pres(presion)

    real(dp), intent(inout)    :: presion
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
    
    real(dp), allocatable   :: v2(:)    ! Vector con la velocidad cuadratica
    real(dp) :: masa
    integer :: inic, fin, i


    allocate(v2(1:gNpart))

    do i = 1, gNpart
      v2(i) =  dot_product(gV(:,i),gV(:,i))  
    end do

    gKin = 0.0_dp
    masa = 0.0_dp 
    inic = 1
    do i=1, gNespecies
      fin = gNp(i) + inic
      masa = gLj_param(i, 3)
      gKin = gKin + masa * sum(v2(inic:fin))
      inic = fin + 1
    enddo

    gKin = 0.5_dp * gKin

  endsubroutine calcula_kin
  !===============================================================================
  ! CALCULA TEMPERATURA INSTANTANEA 
  !===============================================================================
  ! Calcula la anergia cinetica total del sistema

  subroutine calcula_temp(Temp_ins)
    
    real(dp), intent(inout)         :: Temp_ins ! Temperatura instantánea

    Temp_ins = 2.0_dp * gKin / (3*gNpart-3)

  end subroutine calcula_temp


end module mediciones
