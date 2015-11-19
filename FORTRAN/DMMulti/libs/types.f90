!
! Definición de tipos 
!
!
module types


     implicit none
     integer, parameter :: sp = selected_real_kind(6,37) ! single precision
     !integer, parameter :: dp = selected_real_kind(15,307) ! double precision
     
     ! Viene de las corridas de ising 
     !integer, parameter :: dp=kind(0.d0)
     integer, parameter :: DP=SELECTED_REAL_KIND( 12, 60 )

     !! Ejemplo de uso
     real(sp) :: r_sp = 1.0
     real(dp) :: r_dp = 1.0_dp

end module

