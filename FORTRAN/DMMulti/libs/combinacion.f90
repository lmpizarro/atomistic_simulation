!===============================================================================
! CALCULO DE PARAMETROS DE INTERACCIÓN PARA LA MEZCLA DE PARTÍCULAS 
!===============================================================================
!
module combinacion

  use types,       only: dp
  use globales

  implicit none

  private
  
  public    :: comb_epsilon, comb_sigma, corta_desplaza_pote

  integer :: i, j

contains

  !=============================================================================
  ! DEFINE EL EPSILON (DE L-J) ENTRE DISTINTAS ESPECIES 
  !=============================================================================

  subroutine comb_epsilon()
    
    ! TODO convendría que vaya en la parte de inicialización
    allocate (gCombEpsilon(1:gNespecies,1:gNespecies))

    ! asigna la diagonal
    do i=1, gNespecies
      gCombEpsilon(i,i) = gLj_param(i,1)
    enddo

    ! calcula los elementos fuera de la diagonal
    do i =1, gNespecies
      do j=i+1, gNespecies
         gCombEpsilon(i,j) = sqrt(gLj_param(i,1) *  gLj_param(j,1))
         gCombEpsilon(j,i) = gCombEpsilon(i,j) 
      enddo
    enddo
   
    gCombEpsilon = gCombEpsilon / gLj_param (1,1)

    ! imprime la matriz por pantalla
    print *, "gEpsilon"
    do i = 1, gNespecies
      write (*, 700) (gCombEpsilon (i,j), j = 1,gNespecies)
    enddo

    700 format (100g15.5)

  endsubroutine comb_epsilon
 
  !=============================================================================
  ! DEFINE EL SIGMA (DE L-J) ENTRE DISTINTAS ESPECIES 
  !=============================================================================

  subroutine comb_sigma()

    allocate (gCombSigma(1:gNespecies,1:gNespecies))

    ! asigna la diagonal
    do i=1, gNespecies
      gCombSigma(i,i) = gLj_param(i,2)
    enddo

    ! calcula los elementos fuera de la diagonal
    do i =1, gNespecies
      do j=i+1, gNespecies
         gCombSigma(i,j) =  (gLj_param(i,2) + gLj_param(j,2)) / 2
         gCombSigma(j,i) = gCombSigma(i,j)
      enddo
    enddo
 
    gCombSigma = gCombSigma / gLj_param(1,2)

    ! imprime la matriz por pantalla
    print *, "gSigma"
    do i = 1, gNespecies
      write (*, 700) (gCombSigma (i,j), j = 1,gNespecies)
    enddo

    700 format (100g15.5)

  endsubroutine comb_sigma

  !=============================================================================
  ! DEFINE EL RADIO DE CORTE Y SE EVALUA EL POTENCIA 
  !=============================================================================

  subroutine corta_desplaza_pote ()

    real(dp), parameter :: rc = 2.5_dp   ! Radio de corte en unidades de sigma

    allocate (gRc2(1:gNespecies, 1:gNespecies))
    allocate (gPot_Cut(1:gNespecies, 1:gNespecies))

    ! calcula para cada una de las especies 
    ! el radio de corte y el potencial de corte
    ! de acuerdo al cálculo en dinmod
    ! elementos de la diagonal
    do i = 1, gNespecies
       gRc2(i,i) = rc**2
       !gRc2(i,i) = (rc * gCombSigma(i,i) ) ** 2
       gPot_Cut(i,i) = 4.0_dp * gCombEpsilon(i,i) * & 
                     ( 1.0_dp / (gRc2(i,i)**6) - 1.0_dp / (gRc2(i,i)**3) ) 
    enddo

    ! calcula los elementos fuera de la diagonal
    do i =1, gNespecies
      do j=i, gNespecies
         gRc2(i,j) =  rc**2
         !gRc2(i,j) =  (rc * gCombSigma(i,j) ) ** 2
         gRc2(j,i) = gRc2(i,j)

         gPot_Cut(i,j) = 4.0_dp * gCombEpsilon(i,j) * &
                       ( 1.0_dp / (gRc2(i,j)**6) - 1.0_dp / (gRc2(i,j)**3) ) 
         gPot_Cut(j,i) = gPot_Cut(i,j) 
      enddo
    enddo

    ! Multiplico por el sgima en cada caso para tener el radio de corte correcto
    ! En lo anterior no hacía falta porque se simplificaba al evaluar en el potencial
   ! do i=1, gNespecies
   !  do j=i, gNespecies
   !    gRc2(i,j) = gRc2(i,j) * gCombSigma(i,j)**2
   !    gRc2(j,i) = gRc2(i,j)
   !  end do
   ! end do

    ! imprime la matriz de gRc2 por pantalla
    print *, "gRc2"
    do i = 1, gNespecies
      write (*, 700) (gRc2(i,j), j = 1,gNespecies)
    enddo

    ! imprime la matriz de gPot_Cut por pantalla
    print *, "gPot_Cut"
    do i = 1, gNespecies
      write (*, 700) (gPot_Cut (i,j), j = 1,gNespecies)
    enddo

    700 format (100g15.5)

  endsubroutine corta_desplaza_pote

endmodule combinacion
