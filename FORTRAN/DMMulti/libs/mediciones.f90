module mediciones
#include "control.h"
  use types,      only: dp
  use globales
  
  use ziggurat
!#include  "mpif.h"
  !use mpi
  implicit none


contains

  subroutine kernel_fuerza(i, j, cut4, gRc2)
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: Fij       ! Módulo fuerza entre partículas i y j
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    real(dp), intent(in)    :: cut4      ! Cuarta parte del potencial en r_c
    real(dp), intent(in)    :: gRc2      ! Cuarta parte del potencial en r_c
    integer, intent(in) :: i,j

    rij_vec = gR(:,i) - gR(:,j)               ! Distancia vectorial
    ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
    ! Siempre en distancias relativas de sigma
    rij_vec = rij_vec - gLado_caja*nint(rij_vec/gLado_caja)
    r2ij   = dot_product( rij_vec , rij_vec )          ! Cuadrado de la distancia
    if ( r2ij < gRc2) then                
      r2in = 1.0_dp/r2ij                               ! Inversa al cuadrado
      r6in = r2in**3                                   ! Inversa a la sexta
      Fij     = r2in * r6in * (r6in - 0.5_dp)          ! Fuerza entre partículas
      gF(:,i) = gF(:,i) + Fij * rij_vec                ! Contribución a la partícula i
      gF(:,j) = gF(:,j) - Fij * rij_vec                ! Contribucion a la partícula j
      gPot    = gPot + r6in * ( r6in - 1.0_dp) - cut4  ! Energía potencial
      gVir    = gVir + Fij * r2ij                      ! Término del virial para la presión
                                                           ! pg 48 de Allen W=-1/3 sum(r dv/dr)
    end if  ! Termina if del radio de corte

  endsubroutine kernel_fuerza

  subroutine calcula_fuerza()

    !      
    ! elemento = especie, cada especie i tiene Npi partículas
    !
    real(dp):: cut4      ! Cuarta parte del potencial en r_c
    integer :: i, j, ke, je
    integer :: inic, fin
         
    integer :: inic_i, fin_i      
    ! calcula las fuerzas de interacción 
    ! de elementos iguales

    gF    = 0.0_dp
    gPot  = 0.0_dp
    gVir  = 0.0_dp

    ! calcula la fuerza para particulas de iguales elementos
    !
    inic = 1
    do ke=1, gNespecies  ! loop sobre los distintas elementos
      fin = gNp(ke) + inic - 1

#if DEBUG == 1
      print *, "calculo sobre los puros", inic, fin
#endif

      ! Por cuestiones de eficiencia. Se evita hacer una multiplicación dentro
      ! del loop anidado
      cut4 = gPot_cut(ke, ke) / 4.0_dp 

      do i= inic, fin - 1  ! para cada una de las partículas
        do j= i + 1, fin   ! se calcula la interaccion con la otra partícula
        !! calcula la fuerza

#if DEBUG == 1
          print *, "calculo puro ", ke, "i j ", i, j
#endif

          call kernel_fuerza (i, j, cut4, gRc2(ke,ke))

        enddo
      enddo
      inic = fin + 1
    enddo


    ! calculo de las fuerzas de interacción de  
    ! particulas distintas
    inic = 1
    ! ke apunta a una especie
    do ke=1, gNespecies - 1  ! para todos los elementos menos la última
      fin = gNp(ke) + inic - 1

#if DEBUG == 1
      print *, inic, fin
#endif

      inic_i = fin + 1

      ! para cada uno de los elementos
      ! calcular la interacción con los
      ! otros elementos
      do i=inic, fin   ! para todas las partículas de esa especie

        !je apunta a otra especie
        do je = ke + 1, gNespecies 
          cut4 = gPot_cut(ke, je) / 4.0_dp 

          fin_i = inic_i + gNp(je) - 1
            do j=inic_i, fin_i ! para todas las partículas de la otra especie
#if DEBUG == 1
              print *, "loop interno:", i, j
#endif              
              call kernel_fuerza (i, j, cut4, gRc2(ke, je))
            enddo 
          !inic_i = fin_i + 1
        enddo
      enddo

      inic = fin + 1
    enddo

    ! Constantes que faltaban en la energía
    gF = 48.0_dp * gF
    ! Constantes que faltaban en el potencial
    gPot =  4.0_dp * gPot
    ! Se agregan las constantes que faltan para el término del virial
    gVir = 48.0_dp*gVir / 3.0_dp

! Si se utiliza el termostato de Langevin
#if THERM == 1
    call fuerza_langevin()
#endif

#if DEBUG == 1
    do j=1,gNpart 
      print *, gF(1,j) , gF(2,j), gF(3,j)
    enddo  
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
